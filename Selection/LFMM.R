library(lfmm)
library(data.table)
library(qvalue)

setwd("/scratch/patrycja/Goettingen2024/all_batches/vsf/")
frq <- fread("merged_frequencies.txt", header = TRUE, sep = "\t")
dim(frq)
frq <- frq[-14, ] # remove NOM16

gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("land_env.csv")
env <- env[env$Landscape != "NOM16", ]


### determine number of k ###
# pc
pc <- prcomp(gen)
plot(pc$sdev[1:20]^2, xlab = "PC", ylab = "Variance explained")
# points(5,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")

## crop and LST separately
crop <- env$crop

lfmm_crop <- lfmm_ridge(Y = gen, X = crop, K = 13)
crop.pv <- lfmm_test(Y = gen, X = crop, lfmm = lfmm_crop, calibrate = "gif")
crop.pv$gif #### SHOULD BE AROUND 1
crop.qv <- qvalue(crop.pv$calibrated.pvalue)$qvalues
length(which(crop.qv < 0.05)) ## h.w many SNPs have an FDR < 5%?
FDR.crop <- colnames(gen)[which(crop.qv < 0.05)]

LFMM_c <- data.frame(snp = colnames(gen), q_value = crop.qv) %>%
    filter(q_value < 0.05) %>%
    mutate(predictor = "crop")
rownames(LFMM_c) <- NULL


LST <- env$LST

lfmm_LST <- lfmm_ridge(Y = gen, X = LST, K = 13)
LST.pv <- lfmm_test(Y = gen, X = LST, lfmm = lfmm_LST, calibrate = "gif")
LST.pv$gif #### SHOULD BE AROUND 1
LST.qv <- qvalue(LST.pv$calibrated.pvalue)$qvalues
length(which(LST.qv < 0.05)) ## h.w many SNPs have an FDR < 5%?
FDR.LST <- colnames(gen)[which(LST.qv < 0.05)]

LFMM_t <- data.frame(snp = colnames(gen), q_value = LST.qv) %>%
    filter(q_value < 0.05) %>%
    mutate(predictor = "LST")
rownames(LFMM_t) <- NULL

LFMM = rbind(LFMM_c, LFMM_t)

write.csv(LFMM, "LFMM/lfmm_k13.csv", row.names = FALSE)

### all SNPs for Manhattan

LFMM_c_full <- data.frame(snp = colnames(gen), crop_q_value = crop.qv)
rownames(LFMM_c_full) <- NULL


LFMM_t_full <- data.frame(snp = colnames(gen), LST_q_value = LST.qv) 
rownames(LFMM_t_full) <- NULL

LFMM_full = inner_join(LFMM_c_full, LFMM_t_full, by = 'snp')
write.csv(LFMM_full, "LFMM/lfmm_full_k13.csv", row.names = FALSE)