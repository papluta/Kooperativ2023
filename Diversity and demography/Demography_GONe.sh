#!/bin/sh

# remove unassembled contigs
bcftools view -r NC_083488.1,NC_083489.1,NC_083490.1,NC_083491.1,NC_083492.1,NC_083493.1,NC_083494.1,NC_083495.1,NC_083496.1,NC_083497.1,NC_083498.1,NC_083499.1,NC_083500.1,NC_083501.1,NC_083502.1,NC_083503.1,NC_083504.1 \ 
Goe_demography_g095maf005.recode.vcf -Oz -o Goe_demography_filtered_17.vcf.gz

vcftools --vcf Goe_demography_filtered_17.vcf.gz --plink --chrom-map chrom_map.txt --out Goe_demography_filtered_17GONE

# replace ped first line with 1
awk 'BEGIN { OFS = "\t" } $1="1"' Goe_demography_filtered_17GONE.ped > Goe_demography_GONE_temp.ped
# add "-9" before genomic data 
awk 'BEGIN { OFS = "\t" } $6="-9"' Goe_demography_GONE_temp.ped > Goe_demography_filtered_17GONE.ped


bash script_GONE.sh Goe_demography_filtered_17GONE

