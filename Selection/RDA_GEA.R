library(data.table)   
library(vegan)   
library(tidyverse)
library(dplyr)
library(ggplot2)   
library(psych)

# Allele Frequency Table
frq <- fread("merged_frequencies.txt", header = TRUE, sep = "\t")  
dim(frq)
frq = frq[-14,] #remove NOM16

# Extract Genotype Data and Environmental Data
gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("land_env.csv")  
env = env[env$Landscape != 'NOM16',]

# X11()
# pairs.panels(env[,2:5], scale=T)

# Redundancy Analysis Using Environmental and Geographic Variables
RDA <- rda(gen ~ crop + LST_histday + Condition(x + y), data = env, scale = TRUE) # LSThistday = from 2000 to 2023

RDA

RsquareAdj(RDA)

# screeplot(RDA)

scores(RDA, display = "bp")

# Extract and Plot Loadings for the First RDA Axis
load.rda <- scores(RDA, choices = c(1,2), display = "species")  # Extract species (SNP) scores
#hist(load.rda[, 1], main = "Loadings on RDA1")

# Save Loadings to CSV


# Define Function to Identify Outliers Based on Standard Deviation Threshold
outliers <- function(x, z) {
    lims <- mean(x) + c(-1, 1) * z * sd(x)  
    x[x < lims[1] | x > lims[2]]            
}

# Identify Outliers on the First RDA Axis (3 SD from Mean)
cand1 <- outliers(load.rda[, "RDA1"], 3)        # SNPs with loadings > 3 SD from mean on RDA1
cand2 <- outliers(load.rda[, "RDA2"], 3) 
# Count the Number of Outliers
ncand <- length(cand1) + length(cand2)
ncand

# Save Outliers as a Data Frame

### COMBINE CANDIDATE GENES

cand1.1 <- cbind.data.frame(predictor = rep("LST", times = length(cand1)), snp = names(cand1), loading = unname(cand1))
cand2.2 <- cbind.data.frame(predictor = rep("crop", times = length(cand2)), snp = names(cand2), loading = unname(cand2))
cand = rbind(cand1.1, cand2.2)  # Combine candidate SNPs from both axes

cand$snp <- as.character(cand$snp)                 # Ensure SNP column is character type
head(cand)

# Save Outliers
write.csv(cand, "RDA/RDA_cand_crop_temp_xy.csv", row.names = FALSE)

# Save all SNPs
load.rda  = as.data.frame(load.rda)
load.rda$snp = rownames(load.rda)
rownames(load.rda) = NULL
colnames(load.rda) = c("LST","crop", "snp")
write.csv(load.rda, "RDA/RDAload_all_crop_temp_xy.csv", row.names = FALSE)

# limits for plotting
lims1 <- mean(load.rda[, "LST"]) + c(1) * 3 * sd(load.rda[, "LST"]) # 0.0377047138230125
lims2 <- mean(load.rda[, "crop"]) + c(1) * 3 * sd(load.rda[, "crop"]) # 0.0322803435154258
lims1
lims2

######
# hm = as.data.frame(load.rda[,1])
# hm$snp = rownames(hm)
# rownames(hm) = NULL
# hm$outlier = ifelse(hm$snp %in% cand$snp, TRUE, FALSE)
# hm$pos = 1:nrow(hm)
# colnames(hm) = c('load','snp', 'outlier','pos')

# png(filename="outliers_rda.png", width = 1200, height = 480, units = "px")

# ggplot(hm, aes(pos, load))+
# geom_point(aes(col = outlier))

# dev.off()
# frq2 = frq
# colnames(frq2) = c(pop, 1:ncol(frq2))
# frq.plot = frq2 %>% pivot_longer(cols = 2:ncol(frq2), names_to = 'snp', values_to = 'loading') %>% mutate(outlier = ifelse(snp %in% cand1$snp, TRUE, FALSE))

# png(filename="freq_rda.png", width = 1800, height = 480, units = "px")

# ggplot(frq.plot, aes(snp, loading))+
# geom_point(aes(col = outlier))

# dev.off()