#!/bin/sh

## prefilter the SNP dataset for missingness max. 0.5, min allele count 3, min phred 30
vcftools --gzvcf Goe_raw_snp.vcf.gz --max-missing 0.5 --mac 3 --minQ 30 --minDP 3 --recode --recode-INFO-all --out Goe_q5mac3q30dp3

## calculate relatedness, min. depth and missingness per individual
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --missing-inv --out Goe_q5mac3q30dp3
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --relatedness2 --out Goe_q5mac3q30dp3
vcftools --vcf Goe_q5mac3q30dp3.recode.vcf --depth --out Goe_q5mac3q30dp3

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
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out g095_dp4_not_related/Goe_BP_NR_g095maf005q30meandp5to50mindp4

### dataset for selection
vcftools --gzvcf Goe_BP_NR.vcf.gz --max-missing 0.90 --maf 0.01 --minQ 30 --minDP 4 --min-meanDP 5 \
    --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out Goe_BP_NR_g090maf001q30meandp5to50mindp4
	
#Get HWE by snp to exclude highly heterozygotsity sites (errors)
vcftools --vcf Goe_BP_NR_g095maf005q30meandp5to50mindp4.recode.vcf --hardy --out Goe_hardy_095
vcftools --vcf Goe_BP_NR_g0805maf001q30meandp5to50mindp4.recode.vcf --hardy --out Goe_hardy_090maf001

#run python script on both datasets
python Scripts/SnpsHe.py 

#remove snps identified with high heterozygosity
grep -Fwvf highHE_het60_095.txt Goe_BP_NR_g095maf005q30meandp5to50mindp4.recode.vcf > Goe_filtered095_Het60.recode.vcf
grep -Fwvf highHE_het60_085.txt Goe_BP_NR_g0805maf001q30meandp5to50mindp4.recode.vcf > Goe_filtered090maf001_Het60.recode.vcf

## pruning vcf files
# --set-missing-var-ids = assigns chromosome and position based IDs to the variants with missing IDs
# --indep-pairwise = window size (in kb), variant count to shift the window at the end of each step, paiwise r^2 threshold. Prunes the variants so that it excludes ones in linkage disequilibrium

plink --vcf Goe_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out Goe_filtered095_Het60_p1
plink --vcf Goe_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract Goe_filtered095_Het60_p1.prune.in --make-bed --out Goe_filtered095_Het60_p2
plink --bfile Goe_filtered095_Het60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out Goe_filtered095_Het60_pruned

plink --vcf Goe_filtered090maf001_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out Goe_filtered090maf001_Het60_p1
plink --vcf Goe_filtered090maf001_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract Goe_filtered090maf001_Het60_p1.prune.in --make-bed --out Goe_filtered085_Het60_p2
plink --bfile Goe_filtered090maf001_Het60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out Goe_filtered090maf001_Het60_pruned

## zip and index
bgzip -c Goe_filtered095_Het60_pruned.vcf > Goe_filtered095_Het60_pruned.vcf.gz
tabix -p vcf Goe_filtered095_Het60_pruned.vcf.gz

bgzip -c Goe_filtered090maf001_Het60_pruned.vcf > Goe_filtered090maf001_Het60_pruned.vcf.gz
tabix -p vcf Goe_filtered090maf001_Het60_pruned.vcf.gz

## checking relatedness on a filtered (same as for diversity) and pruned dataset including related individuals
vcftools --gzvcf Goe_RELATED_filtered095_Het60_pruned.vcf.gz --relatedness2 --out Goe_RELATED_filtered095_Het60_pruned

## one more related individual detected and removed from not related dataset - SNPs refiltered

plink --vcf Goe_neutral_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --indep-pairwise 50 5 0.5 --out Goe_filtered095_Het60_p1
plink --vcf Goe_neutral_filtered095_Het60.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --extract Goe_filtered095_Het60_p1.prune.in --make-bed --out Goe_filtered095_Het60_p2
plink --bfile Goe_filtered095_Het60_p2 --recode vcf --double-id --allow-extra-chr --set-missing-var-ids @:#:\$1:\$2 --out Goe_neutral_filtered095_Het60_pruned