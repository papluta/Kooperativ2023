#!/bin/sh

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