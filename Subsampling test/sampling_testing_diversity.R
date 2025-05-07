###### Load Required Libraries ######
library(vcfR)
library(adegenet)
library(hierfstat)
library(StAMPP)
library(dplyr)
library(dartR)

###### Read and Process VCF File ######
size = 'size21'
file = list.files(path = paste0('Data/',size, '/'),pattern="*.vcf.gz", full.names = T)
vcf.list = lapply(file, read.vcfR)
gl.list = lapply(vcf.list, vcfR2genlight)
correct_size = lapply(gl.list, function(x) length(x@ind.names))

###### Read Population Data ######
pop.data <- read.table("Data/populations_filtered.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) %>% filter(pop %in% c('NOM08', 'NOM31','NOM21'))

sample.order = data.frame(sample = gl.list[[1]]@ind.names)
pop.data.ordered = sample.order %>% left_join(pop.data, by = 'sample')

###### dartR loop ######

He_loop_fun = function(gl.list, pop) {
  sample.order = data.frame(sample = gl.list[[1]]@ind.names)
  pop.data.ordered = sample.order %>% left_join(pop, by = 'sample')
  glx = list()
  df = list()
  for (i in 1:length(gl.list)) {
    ploidy(gl.list[[i]]) <- 2  # Set ploidy level (diploid)
    pop(gl.list[[i]]) <- pop.data.ordered$pop 
    glx[[i]] <- gl.compliance.check(gl.list[[i]])
    df[[i]] <- gl.report.heterozygosity(glx[[i]], method = "pop", plot.out = FALSE)
  }
  df.comb = bind_rows(df, .id = 'i')
  return(df.comb)
}

df.inter = He_loop_fun(gl.list, pop.data.ordered)
df.comb = bind_rows(df.inter, .id = 'i')
m.df = df.comb %>% group_by(pop) %>% summarise(m.uHe = mean(uHe), sd.uHe = sd(uHe), m.He = mean(He), sd.He = sd(He))


write.csv(df.comb, file = paste0('Data/',size,'_results.csv'), row.names = F)

###### Population Divergence Analysis (StAMPP) ######
## Convert genlight object to StAMPP-compatible format ##
#x <- vcfR2genlight(vcf)  # Convert VCF again to ensure correct format

gl.m = lapply(gl.list, as.matrix)
head(gl.m[[1]][1:20,1:20])
sample.list = lapply(gl.m, row.names)

## Merge with Population Information ##
pop.names <- pop.data.ordered$pop  # Extract population names from file
ploidy <- lapply(gl.list, ploidy) 
gl.m2 = list()
for (i in 1:100) {
  gl.m2[[i]] = gl.m[[i]] * (1/ploidy[[i]])
}
gl.m2[is.na(gl.m2)] <- NaN  # Replace missing values with NaN


## Format Data for StAMPP ##
format <- rep("freq", length(sample))  # Define genotype format

x.stampp = list()
geno = list()
fst = list()
for (i in 1:100) {
x.stampp[[i]] <- as.data.frame(cbind(sample.list[[i]], pop.names, ploidy[[i]], format, gl.m2[[i]])) 
geno[[i]] <- stamppConvert(x.stampp[[i]], "r")  
fst[[i]] <- stamppFst(geno[[i]], nboots = 2, percent = 95, nclusters = 1)
}


## Save Results ##
for (i in 1:100){
write.csv(fst[[i]]$Bootstraps, paste0(size,"/i",i, "_gen_dist_bootstraps.csv"))
write.csv(fst[[i]]$Pvalues, paste0(size,"/i",i, "_gen_dist_Pvalues.csv"))
write.csv(fst[[i]]$Fsts, paste0(size,"/i",i, "_gen_dist_Fsts.csv"))
}

print("Analysis complete! Results saved.")