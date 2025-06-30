library(dplyr)
library(tidyverse)

RDA <- read.csv("RDA/RDA_cand_crop_temp_xy.csv")
LFMM <- read.csv("LFMM/lfmm_k13.csv")


## overlaping snps in LFMM
LFMM <- rbind(LFMM_c, LFMM_t)
LFMM[duplicated(LFMM$snp), ]


common_snps <- RDA %>%
    inner_join(LFMM, by = join_by("snp", "predictor")) 

write.csv(common_snps, file = "RDA_LFMM_SNPS.csv", row.names = FALSE)


# hypergeometric test
total = 2038219
RDA_snps = nrow(RDA) # 4871 + 2592
RDA_snps = 7463
LFMM_snps = nrow(LFMM) # 293 + 6
LFMM_snps = 299
RDA_LFMM_snps = nrow(common_snps_only) # 221
RDA_LFMM_snps = 221

ht = phyper(RDA_LFMM_snps - 1, RDA_snps, total - RDA_snps, LFMM_snps, lower.tail = FALSE) #very low p-value
ht = phyper(221 - 1, 7463, total - 7463, 299, lower.tail = FALSE) #very low p-value


sanity <- RDA_old %>%
    full_join(RDA, by = join_by("snp", "predictor")) 