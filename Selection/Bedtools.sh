#!/bin/sh
GFF_FILE="/scratch/patrycja/reference_genomes/Bombus_pascuorum/annotation/genomic.gff"
VCF_FILE="/scratch/patrycja/reference_genomes/Bombus_pascuorum/annotation/Crop_BP_Goe.ann.vcf"
OUT_FILE="bedtools_output.csv"

# get the SNPs from the VCF file
bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\n' $VCF_FILE > snps.bed 
# merge the SNPs and create a BED file
sort -k1,1 -k2,2n snps.bed | bedtools merge -i - > merged_snps.bed 
# create a BED file with SNPs in the 5kb window around the genes
bedtools window -a merged_snps.bed -b $GFF_FILE -w 5000 > $OUT_FILE