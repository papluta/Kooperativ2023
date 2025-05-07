library(data.table)
library(gradientForest)

frq <- fread("g090_maf001_dp4_not_related/merged_frequencies.txt", header = TRUE, sep = "\t")  
dim(frq)                                               
frq = frq[-14,] #remove NOM16

# Extract Genotype Data and Environmental Data
gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("g090_maf001_dp4_not_related/land_rda.csv")  
env = env[env$Landscape != 'NOM16',]
env = env[,-1]
 
X11()
maxLevel <- log2(0.368*nrow(env)/2)
gf <- gradientForest(cbind(env, gen), predictor.vars=colnames(env), response.vars=colnames(gen), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
imp_matrix <- gf$imp.rsq

importance_scores_crop <- as.numeric(imp_matrix[1, ]) 
importance_scores_LST <- as.numeric(imp_matrix[2, ])
names(importance_scores_crop) <- colnames(imp_matrix)
names(importance_scores_LST) <- colnames(imp_matrix)
imp_scores_crop <- as.numeric(gf$imp.rsq[1, ])
imp_scores_LST <- as.numeric(gf$imp.rsq[2, ])
names(imp_scores_crop) <- colnames(gf$imp.rsq)
names(imp_scores_LST) <- colnames(gf$imp.rsq)

# Sort in descending order
sorted_scores <- sort(imp_scores, decreasing = TRUE)
top_1_percent <- sorted_scores[1:round(length(sorted_scores) * 0.01)]
write.table(top_1_percent, "top_snps_to_extract.txt", quote = FALSE, row.names = TRUE, col.names = FALSE)
