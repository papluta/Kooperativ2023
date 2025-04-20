#!/bin/sh

############# KOOPERATIV Bombus pascuorum samples from 2023 NGS pipeline ################ 


####################################################################################
##################### 	QUALITY FILTERING AND ADAPTER TRIMMING #####################
####################################################################################

#### using bbmap 38.92

# batch_number is the folder for each batch (1st_batch, 2nd_batch, 3rd_batch, 4th_batch, 5th_batch)

batch_number=2nd_batch

mkdir $batch_number/Clean

for i in `ls -1 $batch_number/*_R1_001.fastq.gz`; \
    do dname=$(dirname ${i}); name=$(basename ${i} _R1_001.fastq.gz); \
    bbmap/bbduk.sh in1=${dname}/${name}_R1_001.fastq.gz in2=${dname}/${name}_R2_001.fastq.gz \
    out1=$batch_number/Clean/${name}_clean_1.fastq.gz out2=$batch_number/Clean/${name}_clean_2.fastq.gz \
    ref=bbmap/resources/nextera.fa.gz ktrim=r k=17 mink=8 hdist=1 tpe tbo qtrim=rl trimq=30 ordered=t threads=112 stats=trimstats.txt ;
done

#ktrim=r right-trimming (3' adapters) k specifies the kmer size to use (must be at most the length of the adapters) 
#k=17 k specifies the kmer size to use (must be at most the length of the adapters) 
#mink=8  "mink" allows it to use shorter kmers at the ends of the read (for example, k=11 for the last 11 bases)
#hdist=1 "hdist" means "hamming distance"; this allows one mismatch. 
#tpe "tpe", which specifies to trim both reads to the same length
#tbo the flags "tbo", which specifies to also trim adapters based on pair overlap detection (which does not require known adapter sequences)
#qtrim=rl "qtrim=rl" means it will trim the left and right sides
#trimq=30 quality-trim to Q10 using the Phred algorithm
#ordered=t multiple threads are used, reads will not come out in the same order the went in, unless the “ordered” flag is used
#ref=adapters used by company (usually Nextera)


###########################################################################
##################### MAPPING TO THE REFERENCE GENOME #####################
###########################################################################

#### using bwa-mem2 2.1, samtools 1.11, Picard (GATK 4.2.2)

## make directories for aligned sequences
mkdir $batch_number/Aligned_GATK

## call the genome (using the curated GCF version), to download here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905332965.1/
cd reference_genomes/Bombus_pascuorum/
genome=GCF_905332965.1_iyBomPasc1.1_genomic.fna

## index genome
samtools faidx $genome
bwa-mem2 index $genome


### legend
# bwa-mem2 mem = aligning function (mem - maximal exact match; local alignment!), 
# -t 110 or -@ 110 = use 110 threads (PC efficiency)
# in a local alignment, you try to match your query with a substring (a portion) of your subject (reference).
# in a global alignment you perform an end to end alignment with the subject (and you may end up with a lot of gaps in global alignment if the sizes of query and subject are dissimilar)
# samtools view = print alignment in a specific format, -bSu = (bam format, S?, uncompressed); for flags -F: https://broadinstitute.github.io/picard/explain-flags.html
# samtools sort = sort the alignment by leftmost coordinates, -o = write into an output file
# samtools index = index the sorted alignment to make access to random parts easier
# MarkDuplicates = mark duplicated sequences (using Picard from GATK, could be also done by samtools markdup), gatk --java-options = gatk wrapper syntax for adding specific commands like memory allocation, "-Xmx100G" = use max 100G memory 
# rm -f = force remove intermediary files

cd $batch_number/Aligned_GATK

for i in `ls -1  $batch_number/Clean/*_clean_1.fastq.gz`
    do dname=$(dirname ${i}); name=$(basename ${i} _clean_1.fastq.gz)
    echo "name name is $name"
   
    in1=${dname}/${name}_clean_1.fastq.gz
    in2=${dname}/${name}_clean_2.fastq.gz
    bam=$batch_number/Aligned_GATK/${name}_aligned.bam
    sorted_bam=$batch_number/Aligned_GATK/${name}_aligned_sorted.bam
    rmd_bam=$batch_number/Aligned_GATK/${name}_rmd.bam
          
    bwa-mem2 mem -t 110 -R "@RG\tID:${name}\tSM:${name}\tPL:illumina\tLB:lib2\tPU:unit2" $genome $in1 $in2 | samtools view -@ 110 -bSu - > $bam
    samtools sort -@ 110 -o $sorted_bam $bam
    samtools index -@ 110 $sorted_bam    
    gatk-4.2.2.0/gatk --java-options "-Xmx80G" MarkDuplicates I=$sorted_bam O=$rmd_bam REMOVE_DUPLICATES=true M=${name}.duplicates.txt
    samtools index -@ 110 $rmd_bam

rm -f $bam
rm -f $sorted_bam
rm -f ${sorted_bam}.bai
done 


## get the total and the mapped only number of reads of a BAM file 
## legend
# samtools view -c = count reads, -F flag (filter out) reads: https://broadinstitute.github.io/picard/explain-flags.html
# -F 260 = filter out unmapped and secondary aligned reads
# -F 256 = filter out only secondary aligned reads


output_file="counts_summary.txt"

for i in $batch_number/Aligned_GATK/*_rmd.bam;
do
  base_name=$(basename "$i" _rmd.bam)
  count_total=$(samtools view -c "$i")
  count_mapped=$(samtools view -c -F 260 "$i")
  echo "$base_name, $count_total, $count_mapped" >> "$output_file"
 done


## get the genome coverage and error rate, run a script

# add exe permission
chmod +x bamGenomeCoverage2.sh

#run the script
bash bamGenomeCoverage2.sh 

## exclude samples with low coverage (< 90%) from further analysis (9 samples excluded)
## 705/714 samples retained
 
##########################################################################
##################### VARIANT CALLING USING BCFTOOLS #####################
##########################################################################

##### using bcftools 1.11

## move .bam files from all batches into one folder

mkdir all_batches
mv 1st_batch/Aligned_GATK/*_rmd.bam all_batches/
mv 2nd_batch/Aligned_GATK/*_rmd.bam all_batches/
mv 3rd_batch/Aligned_GATK/*_rmd.bam all_batches/
mv 4th_batch/Aligned_GATK/*_rmd.bam all_batches/
mv 5th_batch/Aligned_GATK/*_rmd.bam all_batches/

# Run mpileup per chromosome using parallel 
cat genome_chromosomes.txt | parallel -j 40 \
    "bcftools mpileup -Ou -f reference_genomes/Bombus_pascuorum/GCF_905332965.1_iyBomPasc1.1_genomic.fna \
    -b all_batches/all_bams.listCorrect \
    -r {} -q 20 -Q 30 -C 50 -d 250 -a FORMAT/DP,FORMAT/AD \
    -o all_batches/vcf/raw_mpileup_{}.bcf"
    
## and run variant calling separately for each chromosome
for i in $(cat genome_chromosomes2.txt)
do name=$(basename ${i})
    echo "processing $name"
bcftools call -mv -v -Oz --threads 80 -o all_batches/vcf/${name}_raw_snps.vcf.gz all_batches/vcf/raw_mpileup_${name}.bcf
done

## merge all vcfs from different chromosomes
ls -1 *_raw_snps.vcf.gz > vcf_list.txt
bcftools concat -Oz vcf_list.txt -o Goe_raw_snp.vcf.gz
tabix -p vcf Goe_raw_snp.vcf.gz

## check summary information stats from vcf file
bcftools stats  -s - Goe_raw_snp.vcf.gz > Goe_raw_snp.sumstats

#################################################################
##################### SNP QUALITY FILTERING #####################
#################################################################

## using vcftools v0.1.13, bcftools 1.11, plink v1.90b6.4

## prefilter the SNP dataset for missingness max. 0.5, min allele count 3, min phred 30
vcftools --gzvcf Goe_raw_snp.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --minDP 3 --recode --recode-INFO-all --out Goe_q5mac3q30dp3

## calculate relatedness, min. depth and missingness per individual
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --missing-inv --out Goe_q5mac3q30dp3
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --relatedness2 --out Goe_q5mac3q30dp3
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --depth-indv --out Goe_q5mac3q30dp3

## run R script to detect low quality/non-Bombus pascuorum individuals (phi < -0.5) and related individuals (phi > 0.2)
Rscript Scripts/relatedness.R

## remove samples from the vcf that are in the out_related.txt (127; H81 removed because it was related to H89 from another population, and more related to H89's population than to its own) and out_other_species.txt files (14)

vcftools --vcf --remove out_combined.txt Goe_raw_snp.recode.vcf --recode --recode-INFO-all --out Goe_BP_NR

## 564/714 samples retained

## zip and index
bgzip -c vcf Goe_BP_NR.recode.vcf > Goe_BP_NR.vcf.gz
tabix -p vcf Goe_BP_NR.vcf.gz

#### refilter with more stringent filters ####

### dataset for genetic diversity, differentiation and structure
vcftools --gzvcf Goe_BP_NR.vcf.gz --max-missing 0.95 --maf 0.05 --minQ 30 --minDP 4 --min-meanDP 5 \
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out Goe_BP_NR_g095maf005q30meandp5to50mindp4

### dataset for selection
vcftools --gzvcf Goe_BP_NR.vcf.gz --max-missing 0.85 --maf 0.05 --minQ 30 --minDP 4 --min-meanDP 5 \
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out Goe_BP_NR_g085maf005q30meandp5to50mindp4
	
#Get HWE by snp to exclude highly heterozygotsity sites (errors)
vcftools --vcf Goe_BP_NR_g095maf005q30meandp5to50mindp4.recode.vcf --hardy --out Goe_hardy_095
vcftools --vcf Goe_BP_NR_g085maf005q30meandp5to50mindp4.recode.vcf --hardy --out Goe_hardy_085

#run python script on both datasets
python Scripts/SnpsHe.py 

#remove snps identified with high heterozygosity
grep -Fwvf highHE_het60_095.txt Goe_BP_NR_g095maf005q30meandp5to50mindp4.recode.vcf > Goe_filtered095_Het60.recode.vcf
grep -Fwvf highHE_het60_085.txt Goe_BP_NR_g085maf005q30meandp5to50mindp4.recode.vcf > Goe_filtered085_Het60.recode.vcf

## pruning vcf files
# --set-missing-var-ids = assigns chromosome and position based IDs to the variants with missing IDs
# --indep-pairwise = window size (in kb), variant count to shift the window at the end of each step, paiwise r^2 threshold. Prunes the variants so that it excludes ones in linkage disequilibrium

plink --vcf Goe_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out Goe_filtered095_Het60_p1
plink --vcf Goe_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract Goe_filtered095_Het60_p1.prune.in --make-bed --out Goe_filtered095_Het60_p2
plink --bfile Goe_filtered095_Het60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out Goe_filtered095_Het60_pruned

plink --vcf Goe_filtered085_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out Goe_filtered085_Het60_p1
plink --vcf Goe_filtered085_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract Goe_filtered085_Het60_p1.prune.in --make-bed --out Goe_filtered085_Het60_p2
plink --bfile Goe_filtered085_Het60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out Goe_filtered085_Het60_pruned

## zip and index
bgzip -c Goe_filtered095_Het60_pruned.vcf > Goe_filtered095_Het60_pruned.vcf.gz
tabix -p vcf Goe_filtered095_Het60_pruned.vcf.gz

bgzip -c Goe_filtered085_Het60_pruned.vcf > Goe_filtered085_Het60_pruned.vcf.gz
tabix -p vcf Goe_filtered085_Het60_pruned.vcf.gz

## checking relatedness on a filtered (same as for diversity) and pruned dataset including related individuals
vcftools --gzvcf Goe_RELATED_filtered095_Het60_pruned.vcf.gz --relatedness2 --out Goe_RELATED_filtered095_Het60_pruned

## one more related individual detected and removed from not related dataset - SNPs refiltered

############################################################################
##################### DIVERISTY AND STRUCTURE ANALYSIS #####################
############################################################################

## plink v1.90b6.21

## Population structure with PCA
plink --vcf Goe_filtered095_Het60_pruned.vcf.gz --double-id --allow-extra-chr --pca --out Goe_filtered095_Het60_pruned_PCA

## Runs of Homozygosity - to do
plink --vcf Goe_filtered095_Het60_pruned.vcf.gz --allow-extra-chr --double-id --set-missing-var-ids @:#:\$1:\$2 --homozyg --homozyg-snp 50 --homozyg-kb 100 --homozyg-density 50 --homozyg-window-het 0 --homozyg-window-missing 3 --homozyg-group --out diversity_dataset/Goe_filtered095_Het60_pruned_Ho



