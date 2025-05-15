library(lfmm)
library(data.table)
library(qvalue)


frq <- fread("merged_frequencies.txt", header = TRUE, sep = "\t")
dim(frq)
frq <- frq[-14, ] # remove NOM16

gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("land_rda.csv")
env <- env[env$Landscape != "NOM16", ]

pc <- prcomp(gen)
plot(pc$sdev[1:20]^2, xlab = "PC", ylab = "Variance explained")
# points(5,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")

## crop and LST separately
pred1 <- env$crop

lfmm_crop1 <- lfmm_ridge(Y = gen, X = pred1, K = 13)
AF.pv <- lfmm_test(Y = gen, X = pred1, lfmm = lfmm_crop1, calibrate = "gif")
AF.pv$gif #### SHOULD BE AROUND 1
AF.qv <- qvalue(AF.pv$calibrated.pvalue)$qvalues
length(which(AF.qv < 0.05)) ## h.w many SNPs have an FDR < 5%?
FDR.13.crop <- colnames(gen)[which(AF.qv < 0.05)]

LFMM_c <- data.frame(snp = colnames(gen), q_value = AF.qv) %>%
    filter(q_value < 0.05) %>%
    mutate(predictor = "crop")
rownames(LFMM_c) <- NULL


write.csv(LFMM_c, "LFMM/AF_lfmm_crop_k13.csv", row.names = FALSE)

pred2 <- env$LST

lfmm_LST1 <- lfmm_ridge(Y = gen, X = pred2, K = 13)
AF.pv <- lfmm_test(Y = gen, X = pred2, lfmm = lfmm_LST1, calibrate = "gif")
AF.pv$gif #### SHOULD BE AROUND 1
AF.qv <- qvalue(AF.pv$calibrated.pvalue)$qvalues
length(which(AF.qv < 0.05)) ## h.w many SNPs have an FDR < 5%?
FDR.13.LST <- colnames(gen)[which(AF.qv < 0.05)]

LFMM_t <- data.frame(snp = colnames(gen), q_value = AF.qv) %>%
    filter(q_value < 0.05) %>%
    mutate(predictor = "LST")
rownames(LFMM_t) <- NULL

write.csv(LFMM_t, "LFMM/AF_lfmm_LST_k13.csv", row.names = FALSE)
