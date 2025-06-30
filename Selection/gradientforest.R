library(data.table)
library(gradientForest)

setwd("/scratch/patrycja/Goettingen2024/all_batches/vcf/g090_maf001_dp4_not_related/")

frq <- fread("merged_frequencies.txt", header = TRUE, sep = "\t")
dim(frq)
frq <- frq[-14, ] # remove NOM16

# Extract Genotype Data and Environmental Data
gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("land_rda.csv")
env <- env[env$Landscape != "NOM16", ]
env <- env[, -1]

X11()
maxLevel <- log2(0.368 * nrow(env) / 2)
gf <- gradientForest(cbind(env, gen), predictor.vars = colnames(env), response.vars = colnames(gen), ntree = 10, maxLevel = maxLevel, corr.threshold = 0.50)
imp_matrix <- gf$imp.rsq

importance_scores_crop <- as.numeric(imp_matrix[1, ])
importance_scores_LST <- as.numeric(imp_matrix[4, ])
names(importance_scores_crop) <- colnames(imp_matrix)
names(importance_scores_LST) <- colnames(imp_matrix)
rsq_crop <- as.numeric(gf$imp.rsq[1, ])
rsq_LST <- as.numeric(gf$imp.rsq[4, ])
names(rsq_crop) <- colnames(gf$imp.rsq)
names(rsq_LST) <- colnames(gf$imp.rsq)

# Sort in descending order
sorted_scores_crop <- sort(importance_scores_crop, decreasing = TRUE)
sorted_scores_LST <- sort(importance_scores_LST, decreasing = TRUE)

sorted_rsq_crop <- sort(rsq_crop, decreasing = TRUE)
sorted_rsq_LST <- sort(rsq_LST, decreasing = TRUE)


top_1_crop <- sorted_scores_crop[1:round(length(sorted_scores_crop) * 0.01)]
top_1_LST <- sorted_scores_LST[1:round(length(sorted_scores_LST) * 0.01)]
write.table(top_1_percent, "top_snps_to_extract.txt", quote = FALSE, row.names = TRUE, col.names = FALSE)
