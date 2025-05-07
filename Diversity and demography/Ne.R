library(dplyr)
library(tidyverse)
library(dartRverse)
library(dartR.popgen)
library(ggplot2)
library(vcfR)
library(dplyr)
library(graph4lg)

vcf <- read.vcfR("Data/Goe_filtered095_Het60_pruned.vcf.gz")  # Change to your VCF filename
gl <- vcfR2genlight(vcf)
gl@chromosome

pop.data <- read.table("Data/populations_filtered.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) # Ensure this file has a column named "pop"
pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor

length(unique(pop.data$pop))

sample.order = data.frame(sample = gl@ind.names)

pop.data.ordered = sample.order %>% left_join(pop.data, by = 'sample')

###### Assign Population Information ######
ploidy(gl) <- 2  # Set ploidy level (diploid)
pop(gl) <- pop.data.ordered$pop  # Assign populations from external file
gl  # Check if populations are correctly assigned
gl@chromosome = as.factor(sub('\\.','_', as.character(gl@chromosome)))
unique(gl@chromosome)
gl@loc.names = sub('\\.','_', as.character(gl@loc.names))
unique(gl@loc.names)

#gl.sub = gl[gl@pop %in% c('NOM08', 'NOM21', 'NOM31')]
gl.sub = gl[,1:37507]

gen = gl2genepop(gl.sub, outfile = 'Goe_full.gen', outpath = getwd())




path.binaries <- "D://Programs/NeEstimator/"
path.data <- "C:/Users/patry/OneDrive/PhD/tralala/kooperativ/Kooperativ2023/Data/"

first = Sys.time()
nessep <- gl.LDNe(gl.sub,
                  outfile = "GoeLD2.txt", pairing="separate",outpath=getwd(),
                  neest.path = path.binaries,
                  critical = c(0, 0.05), singleton.rm = TRUE, mating = "random")
second = Sys.time()
save(nessep, file = 'Data/nessep.RData')

genepop = gl2genepop(gl, outfile = 'test.gen', outpath = getwd())


gen1 = readRDS('Data/gen1.rds')



genind_to_genepop(gen1, output = "gen1.txt")

########## for multiple files
size = 'size3'
file = list.files(path = paste0('D://Repositories/Kooperativ/',size, '/'),pattern="*.vcf.gz", full.names = T)
vcf.list = lapply(file, read.vcfR)
gl.list = lapply(vcf.list, vcfR2genlight)
names(gl.list) = 1:100

###### Assign Population Information ######
for (i in 1:length(gl.list)) {
ploidy(gl.list[[i]]) <- 2  
sample.order = data.frame(sample = gl.list[[i]]@ind.names)
pop.data.ordered = sample.order %>% left_join(pop.data, by = 'sample')
pop(gl.list[[i]]) <- pop.data.ordered$pop  
gl.list[[i]]@chromosome = as.factor(gsub('\\.','_', as.character(gl.list[[i]]@chromosome)))
gl.list[[i]]@loc.names = gsub('\\.','_', as.character(gl.list[[i]]@loc.names))
gl.list[[i]]@loc.names = gsub(':','_', as.character(gl.list[[i]]@loc.names))
}

### or ###

gen.list <- lapply(vcf.list, function(x) vcfR2genind(x))
gen.list[[1]]


for (i in 1:length(gen.list)) {
  ploidy(gen.list[[i]]) <- 2  
  sample.order = data.frame(sample = gen.list[[i]]@ind.names)
  pop.data.ordered = sample.order %>% left_join(pop.data, by = 'sample')
  pop(gen.list[[i]]) <- pop.data.ordered$pop  
  gl.list[[i]]@chromosome = as.factor(gsub('\\.','_', as.character(gen.list[[i]]@chromosome)))
  gl.list[[i]]@loc.names = gsub('\\.','_', as.character(gen.list[[i]]@loc.names))
  gl.list[[i]]@loc.names = gsub(':','_', as.character(gen.list[[i]]@loc.names))
}

lapply(gl.list, function(x) genind_to_genepop(x, output = paste0(size,'i',names(x),'.txt')))

#or

# path.binaries <- "D://Programs/NeEstimator/"
# path.data <- "C:/Users/patry/OneDrive/PhD/tralala/kooperativ/Kooperativ2023/Data/"
# 
# first = Sys.time()
# nessep = list()
# for (i in 1:length(gl.list)) {
# nessep[[i]] <- gl.LDNe(gl.list[[i]],
#                   outfile = paste0(size, 'i', names(gl.list[i]), ".txt"), pairing="separate", outpath=getwd(),
#                   neest.path = path.binaries,
#                   critical = c(0, 0.05), singleton.rm = TRUE, mating = "random")
# }
# 
# second = Sys.time()
# save(nessep, file = 'Data/nessep.RData')
# 
# 
# nessep[[2]] <- gl.LDNe(gl.list[[2]],
#                        outfile = paste0(size, 'i', names(gl.list[2]), ".txt"), pairing="separate", outpath=getwd(),
#                        neest.path = path.binaries,
#                        critical = c(0, 0.05), singleton.rm = TRUE, mating = "random")

############ results from NeEstimator

ne = read.csv('Data/NE_results_clean.csv')
ne.t = ne %>% pivot_longer(cols = 2:6, names_to = 'Ne_type', values_to = 'value') %>% 
  pivot_wider(names_from = Type, values_from = value)
