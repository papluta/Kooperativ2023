###### Load Required Libraries ######
library(vcfR)
library(adegenet)
library(hierfstat)
library(dartRverse)
library(StAMPP)

###### Read and Process VCF File ######
vcf <- read.vcfR("Data/Goe_filtered095_Het60_pruned.vcf.gz")  # Change to your VCF filename
gl <- vcfR2genlight(vcf)  # Convert to genlight object
gl <- gl[indNames(gl) != "C49_C49", ]
###### Read Population Data ######
pop.data <- read.table("Data/populations_filtered.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) # Ensure this file has a column named "pop"
pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor

length(unique(pop.data$pop))

sample.order = data.frame(sample = gl@ind.names)

pop.data.ordered = sample.order %>% left_join(pop.data, by = 'sample')

###### Assign Population Information ######
ploidy(gl) <- 2  
pop(gl) <- pop.data.ordered$pop  
gl  

###### Compliance Check ######
glx <- gl.compliance.check(gl)
glx_nomono <- gl.filter.monomorphs(glx)

###### Heterozygosity Analysis Using dartR ######
df <- gl.report.heterozygosity(glx, method = "pop",  nboots = 100)  
write.csv(df, "Results/Goe_filtered095_Het60_pruned_heterozygosity_boot.csv")  

###### Population Divergence Analysis (StAMPP) ######
## Convert genlight object to StAMPP-compatible format ##

x2 <- as.matrix(x)  
sample <- row.names(x2)  

## Merge with Population Information ##
pop.names <- pop.data.ordered$pop  
ploidy <- ploidy(x)  
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

print("Analysis complete! Results saved.")
