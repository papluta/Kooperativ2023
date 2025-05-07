library(dplyr)
library(tidyverse)

## all samples
pop = read.table('Data/populations.txt') %>% mutate(V1 = paste0(V1,'_',V1))

## to get the names of the samples that were included in analysis
samples = read.csv('Data/samples_NR.txt', sep = '', header = F) %>% 
  left_join(pop, by = join_by('V1')) 

#NOM40, NOM32, NOM21

## filter for analysed samples and from the three chosen sites
sampling.dat = samples %>% dplyr::select(V1, V2) %>% filter(V2 %in% c('NOM08', 'NOM31','NOM21'))
sam.dat = split(sampling.dat, sampling.dat$V2)

sampling_pops = function(data, n) {
  size = list()
  for (s in 1:length(data)){
    size[[s]] = list()
    for (i in 1:100) {
      size[[s]][[i]] = sample(data[[s]]$V1, n, F)
    }
  }
  size.c = list()
  for (i in 1:100) {
    size.c[[i]] = c(size[[1]][[i]], size[[2]][[i]], size[[3]][[i]])
  }
  return(size.c)
}

size3 = sampling_pops(sam.dat, n = 3)
size6 = sampling_pops(sam.dat, n = 6)
size9 = sampling_pops(sam.dat, n = 9)
size12 = sampling_pops(sam.dat, n = 12)
size15 = sampling_pops(sam.dat, n = 15)
size18 = sampling_pops(sam.dat, n = 18)
size21 = sampling_pops(sam.dat, n = 21)


## save files, for size3 only for now
for (i in 1:100) {
  write.table(size3[[i]],file = paste0('Data/size3/s3iter',i,'.txt'), row.names = F, quote = F, col.names = F)
}

for (i in 1:100) {
  write.table(size9[[i]],file = paste0('Data/size9/s9iter',i,'.txt'), row.names = F, quote = F, col.names = F)
}

for (i in 1:100) {
  write.table(size15[[i]],file = paste0('Data/size15/s15iter',i,'.txt'), row.names = F, quote = F, col.names = F)
}

for (i in 1:100) {
  write.table(size21[[i]],file = paste0('Data/size21/s21iter',i,'.txt'), row.names = F, quote = F, col.names = F)
}
