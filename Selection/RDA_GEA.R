library(data.table)   
library(vegan)   
library(tidyverse)
library(dplyr)
library(ggplot2)   
library(psych)  

# Allele Frequency Table
frq <- fread("g090_maf001_dp4_not_related/merged_frequencies.txt", header = TRUE, sep = "\t")  
dim(frq)                                               
frq = frq[-14,] #remove NOM16

# Extract Genotype Data and Environmental Data
gen <- frq[, 2:ncol(frq)]

# Environmental Variables
env <- read.csv("g090_maf001_dp4_not_related/land_rda.csv")  
env = env[env$Landscape != 'NOM16',]

X11()
pairs.panels(env[,2:5], scale=T)

# Redundancy Analysis Using Environmental and Geographic Variables
RDA <- rda(gen ~ crop + LST + Condition(x + y), data = env, scale = TRUE)

RDA

RsquareAdj(RDA)

screeplot(RDA)

# Extract and Plot Loadings for the First RDA Axis
load.rda <- scores(RDA, choices = c(1), display = "species")  # Extract species (SNP) scores for the first RDA axis
#hist(load.rda[, 1], main = "Loadings on RDA1")                # Plot histogram of loadings on RDA1

# Save Loadings to CSV
colnames(load.rda) <- c("loading")
write.csv(load.rda, "RDAload_crop_temp_xy.csv", row.names = FALSE)

# Define Function to Identify Outliers Based on Standard Deviation Threshold
outliers <- function(x, z) {
    lims <- mean(x) + c(-1, 1) * z * sd(x)   # Set upper and lower bounds (z SD from mean)
    x[x < lims[1] | x > lims[2]]             # Return SNP names in tails of distribution
}

# Identify Outliers on the First RDA Axis (3 SD from Mean)
cand1 <- outliers(load.rda[, 1], 3)        # SNPs with loadings > 3 SD from mean on RDA1
cand2 <- outliers(load.rda[, 2], 3) 
# Count the Number of Outliers
ncand <- length(cand1) + length(cand2)
ncand

# Save Outliers as a Data Frame

### COMBINE CANDIDATE GENES
cand <- cbind.data.frame(rep(1, times = length(cand1)), names(cand1), unname(cand1),
                         rep(1, times = length(cand2)), names(cand2), unname(cand2))
colnames(cand) <- c("axis", "snp", "loading")       # Add column names
cand$snp <- as.character(cand$snp)                 # Ensure SNP column is character type
head(cand)

# Save Outliers to CSV
write.csv(cand, "RDA_crop_temp_xy.csv", row.names = FALSE)

hm = as.data.frame(load.rda[,1])
hm$snp = rownames(hm)
rownames(hm) = NULL
hm$outlier = ifelse(hm$snp %in% cand$snp, TRUE, FALSE)
hm$pos = 1:nrow(hm)
colnames(hm) = c('load','snp', 'outlier','pos')

png(filename="outliers_rda.png", width = 1200, height = 480, units = "px")

ggplot(hm, aes(pos, load))+
geom_point(aes(col = outlier))

dev.off()
frq2 = frq
colnames(frq2) = c(pop, 1:ncol(frq2))
frq.plot = frq2 %>% pivot_longer(cols = 2:ncol(frq2), names_to = 'snp', values_to = 'loading') %>% mutate(outlier = ifelse(snp %in% cand1$snp, TRUE, FALSE))

png(filename="freq_rda.png", width = 1800, height = 480, units = "px")

ggplot(frq.plot, aes(snp, loading))+
geom_point(aes(col = outlier))

dev.off()