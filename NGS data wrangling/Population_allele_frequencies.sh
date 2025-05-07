#!/bin/bash

# Convert VCF to PLINK format
echo "Converting VCF to PLINK format..."
vcftools --gzvcf Goe_filtered09maf001_Het60_pruned.vcf.gz \
         --plink \
         --out ./GOEcorrectplink \
         --chrom-map chrom_map.txt

# Impute missing alleles and convert to binary PLINK format
echo "Imputing missing alleles and converting to binary PLINK format..."
plink --file GOEcorrectplink \
      --fill-missing-a2 \
      --make-bed \
      --double-id \
      --allow-extra-chr \
      --chr-set 65 \
      --out GOEplink_imputed

# Recode binary PLINK format back to VCF
echo "Recoding binary PLINK format back to VCF..."
plink --allow-extra-chr \
      --bfile GOEplink_imputed \
      --recode vcf \
      --no-fid \
      --chr-set 65 \
      --out GOEimputed_output

# Subset VCF based on specific sample groups
echo "Subsetting VCF based on sample groups..."
for group in pops/*.txt; do
    basename=$(basename "$group" .txt)
    vcftools --vcf GOEimputed_output.vcf \
             --keep $group \
             --recode \
             --out pops/$basename
    echo "Subset created for $basename."
done

# Analyze allele frequencies for each subset
echo "Analyzing allele frequencies for each subset..."
for group in pops/*.recode.vcf; do
    basename=$(basename "$group" .recode.vcf)
    plink2 --vcf "$group" \
           --double-id \
           --allow-extra-chr \
           --chr-set 65 \
           --freq \
           --out pops/$basename
    echo "Allele frequency analysis completed for $group."
done

# Merge Allele Frequency Files 
python merge_afreq.py /pops ./allele_freq.txt
