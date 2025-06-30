#!/bin/sh

java -jar snpEff.jar build -gff3 -v Bpasc -noCheckCds -noCheckProtein -ud 0

snpEff Bpasc LST_SNPS_mis.vcf.gz > LST_BP_Goe.ann.vcf

snpEff Bpasc crop_SNPS_mis.vcf.gz > Crop_BP_Goe.ann.vcf
