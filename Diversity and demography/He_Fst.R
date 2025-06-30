###### Load Required Libraries ######
library(vcfR)
library(adegenet)
library(hierfstat)
library(dartRverse)
library(StAMPP)

###### Read and Process VCF File ######
vcf <- read.vcfR("Data/Goe_demography_g095maf005.recode.vcf.gz")  # Change to your VCF filename
gl <- vcfR2genlight(vcf)  # Convert to genlight object
gl <- gl[indNames(gl) != "C49_C49", ]

nam = data.frame(Sample = indNames(gl))
###### Read Population Data ######
pop.data <- read.table("Data/populations.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) %>% # Ensure this file has a column named "pop"
  mutate(sample = paste0(sample, "_", sample))

length(unique(pop.data$pop))

sample.order = gl@ind.names

pop.data.ordered = pop.data[match(sample.order, pop.data$sample), ]
sample.order
head(pop.data.ordered)

pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor

length(unique(pop.data$pop))

###### Assign Population Information ######
ploidy(gl) <- 2  
pop(gl) <- pop.data.ordered$pop  
gl  
tets = data.frame(gl@ind.names, gl@pop)
tets[is.na(tets),]

###### Compliance Check ######
glx <- gl.compliance.check(gl)
glx_nomono <- gl.filter.monomorphs(glx)

#save(glx, file = "Data/glx.RData")

###### Heterozygosity Analysis Using dartR ######

#######################################
######## outlier in NOM15 #############
######################################

he <- gl.report.heterozygosity(glx_nomono, method = "pop", nboots = 100, CI.type = "basic")  

write.csv(he, "Data/Goe_filtered095_Het60_pruned_heterozygosity_neutral.csv")  

###### Population Divergence Analysis (StAMPP) ######
## Convert genlight object to StAMPP-compatible format ##

x2 <- as.matrix(gl)  
sample <- row.names(x2)  

## Merge with Population Information ##
pop.names <- pop.data.ordered$pop  
ploidy <- ploidy(gl)  
x2 = x2 * (1/ploidy)  # Convert allele counts to frequency
x2[is.na(x2)] <- NaN  # Replace missing values with NaN

## Format Data for StAMPP ##
format <- rep("freq", length(sample))  # Define genotype format
x.stampp <- as.data.frame(cbind(sample, pop.names, ploidy, format, x2))  # Create dataframe for StAMPP
geno <- stamppConvert(x.stampp, "r")  # Convert to StAMPP format

## Compute Fst with 100 Bootstraps ##
fst <- stamppFst(geno, nboots = 100, percent = 95, nclusters = 1)

## Save Results ##
write.csv(fst$Bootstraps, "Results/Goe_filtered095_Het60_pruned_distances_bootstraps.csv")
write.csv(fst$Pvalues, "Results/Goe_filtered095_Het60_prunedc_distances_Pvalues.csv")
write.csv(fst$Fsts, "Results/Goe_filtered095_Het60_pruned_distances_Fsts.csv")

## plotting
ggplot(df, aes(pop, uHe))+
  geom_point()+
  geom_errorbar(aes(ymin = uHeLCI, ymax = uHeHCI))


