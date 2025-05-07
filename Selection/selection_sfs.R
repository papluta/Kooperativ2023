library(dplyr)
library(tidyverse)
library(ggplot2)

pop.data <- read.table("Data/populations_filtered.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) # Ensure this file has a column named "pop"
#pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor


pop.s = split(pop.data, pop.data$pop) %>% lapply(function(x) x %>% select(-pop) %>%  mutate(sample = paste0(sample,'_', sample)))


for (i in 1:length(pop.s)) {
  write.table(pop.s[[i]], file = paste0('Data/pops/', names(pop.s[i]), '.txt'), quote = F, row.names = F, col.names = F)
}



sfs = read.table('Data/merged_frequencies.txt', sep = '', header = T)
loci = colnames(sfs)[-1]
colnames(sfs) = c('pop', 1:(ncol(sfs)-1))
sfs2 = sfs
sfs2$pop = sub('pops/','', sfs2$pop)

sfs3 = sfs2 %>% pivot_longer(cols = 2:ncol(sfs2), names_to = 'chr', values_to = 'freq') %>% mutate(chr  = as.numeric(chr)) %>% filter(pop != 'NOM16')
sfs3.m = sfs3 %>% group_by(chr) %>% summarise(freq.m = mean(freq), freq.sd = sd(freq))

sfs3.m %>% ggplot(aes(chr, freq.m))+
  geom_bar(stat = 'identity')
