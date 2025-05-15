library(dplyr)
library(tidyverse)

RDA <- read.csv("RDA/RDA_crop_temp_xy.csv")
LFMM_c <- read.csv("LFMM/AF_lfmm_crop_k13.csv")
LFMM_t <- read.csv("LFMM/AF_lfmm_temp_k13.csv")


## overlaping snps in LFMM
LFMM <- rbind(LFMM_c, LFMM_t)
LFMM[duplicated(LFMM$snp), ]


common_snps_only <- RDA %>%
    select(-c(crop, LST)) %>%
    full_join(LFMM, by = join_by("snp", "predictor")) %>%
    filter(!is.na(q_value) & !is.na(correlation))

write.csv(common_snps_only, file = "RDA_LFMM_SNPS.csv", row.names = FALSE)
