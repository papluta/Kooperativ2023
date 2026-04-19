library(vcfR)
library(adegenet)
library(hierfstat)
library(dartRverse)
library(StAMPP)

vcf <- read.vcfR("Data/variants_filtered_diversity.vcf.gz")  
gl <- vcfR2genlight(vcf)  
gl <- gl[indNames(gl) != "C49_C49", ] # removing NOM16 singular sample


pop.data <- read.table("Data/populations.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) %>% 
  mutate(sample = paste0(sample, "_", sample))

length(unique(pop.data$pop))

sample.order = gl@ind.names

pop.data.ordered = pop.data[match(sample.order, pop.data$sample), ]
sample.order
head(pop.data.ordered)

pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor

length(unique(pop.data.ordered$pop))

###### Assign Population Information ######
ploidy(gl) <- 2  
pop(gl) <- pop.data.ordered$pop  
gl  
tets = data.frame(gl@ind.names, gl@pop)
tets[is.na(tets),]

###### Compliance Check ######
glx <- gl.compliance.check(gl)
glx_nomono <- gl.filter.monomorphs(glx)


###### Heterozygosity Analysis ######


he <- gl.report.heterozygosity(glx_nomono, method = "pop")  

write.csv(he, "Data/heterozygosity.csv")  

###### Population Divergence Analysis ######

x2 <- as.matrix(gl)  
sample <- row.names(x2)  

## Merge with Population Information ##
pop.names <- pop.data.ordered$pop  
ploidy <- ploidy(gl)  
x2 = x2 * (1/ploidy)  # Convert allele counts to frequency
x2[is.na(x2)] <- NaN  # Replace missing values with NaN

## Format Data for StAMPP ##
format <- rep("freq", length(sample)) 
x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2)) 
geno <- stamppConvert(x.stampp, "r")  

## Compute Fst with 100 Bootstraps ##
fst <- stamppFst(geno, nboots = 100, percent = 95, nclusters = 1)

## Save Results ##
write.csv(fst$Bootstraps, "Data/genetic_distances_bootstraps.csv")
write.csv(fst$Pvalues, "Data/genetic_distances_Pvalues.csv")
write.csv(fst$Fsts, "Data/genetic_distances_Fsts.csv")



