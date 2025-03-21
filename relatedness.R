library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)

# read relatedness, missingness and depth of coverage
rel = read.table('Data/Goe_RELATED_filtered095_Het60_pruned.relatedness2', header = T)
miss = read.table('Data/Goe_BP_RELATED_g095maf005q30meandp5to50mindp4.imiss', header = T) %>% select(INDV, F_MISS)
depth = read.table('Data/Goe_BP_RELATED_g095maf005q30meandp5to50mindp4.idepth', header = T) %>% select(INDV, MEAN_DEPTH)

### create file with low quality/non-BP samples

rel2 = rel %>% filter(RELATEDNESS_PHI != 0.5)

pop.raw = read.csv('Data/Kooperativ2023_samples.csv')

pop.raw$INDV = paste0(pop.raw$Sample_ID,'_', pop.raw$Sample_ID)


## other spec
other.species = rel2 %>% filter(RELATEDNESS_PHI < -0.5) %>% group_by(INDV1) %>% summarise(phi = mean(RELATEDNESS_PHI), n_comp = n()) %>% arrange(desc(n_comp)) %>% 
  left_join(pop.raw %>% select(INDV, Landscape_ID), by = join_by('INDV1' == 'INDV'))

write.table(other.species, file = 'other_species.txt', row.names = F, col.names = T, quote = F, sep = '\t')


## related individuals
related = rel2 %>% left_join(pop.raw %>% select(INDV, Landscape_ID), by = join_by('INDV1' == 'INDV')) %>% 
  rename(Pop1 = Landscape_ID) %>%
  left_join(pop.raw %>% select(INDV, Landscape_ID), by = join_by('INDV2' == 'INDV')) %>%
  rename(Pop2 = Landscape_ID) %>% select(INDV1, INDV2, Pop1, Pop2, RELATEDNESS_PHI) %>% filter(RELATEDNESS_PHI > 0.2) %>%
  mutate(Same = ifelse(Pop1 == Pop2, TRUE, FALSE)) %>%
  mutate(INDV1 = sub('_.*','', INDV1), INDV2 = sub('_.*','', INDV2)) %>%
  left_join(miss, by = join_by('INDV1' == 'INDV')) %>% rename(miss1 = F_MISS) %>%
  left_join(miss, by = join_by('INDV2' == 'INDV')) %>% rename(miss2 = F_MISS) %>%
  left_join(depth, by = join_by('INDV1' == 'INDV')) %>% rename(depth1 = MEAN_DEPTH) %>%
  left_join(depth, by = join_by('INDV2' == 'INDV')) %>% rename(depth2 = MEAN_DEPTH) %>%
  ungroup() %>% rowwise() %>%
  mutate(out_depth = ifelse(depth1 < depth2, INDV1, INDV2), out_miss = ifelse(miss1 > miss2, INDV1, INDV2)) %>%
  mutate(concordant = ifelse(out_depth == out_miss, T, F)) 

out_rel_test = unique(related$out_depth)
#out_rel_test = sub('_.*','',out_rel_test)

write.table(out_rel_test, file = 'out_related.txt', row.names = F, col.names = F, quote = F, sep = '\t')
write.table(related, file = 'out_related.txt', row.names = F, col.names = T, quote = F, sep = '\t')

rel_pop_n = rel_pop %>% filter(Same == TRUE) %>% group_by(INDV1, Pop1) %>% 
  summarise(n = n(), sib = paste(as.character(INDV2), collapse = ', '))

write.table(rel_pop_n2, file = 'sib_groups_c.txt', row.names = F, col.names = T, quote = F, sep = '\t')

### why do H81 from NOM42 and H89 from NOM09 have high relatedness???

h81 = rel2 %>% filter(INDV1 == 'H81_H81' | INDV2 == 'H81_H81') %>% filter(INDV1 != 'H89_H89' & INDV2 != 'H89_H89') %>%
  mutate(INDV_h81 = 'H81_H81', INDV_2 = ifelse(INDV1 == 'H81_H81', INDV2, INDV1)) %>% 
  left_join(pop.raw %>% select(INDV, Landscape_ID), by = join_by('INDV_2' == 'INDV')) %>%
  rename(Pop2 = Landscape_ID) %>%
  filter(Pop2 %in% c('NOM42', 'NOM09')) %>% group_by(Pop2) %>% summarise(phi_m = mean(RELATEDNESS_PHI))

h89 = rel2 %>% filter(INDV1 == 'H89' | INDV2 == 'H89') %>% filter(INDV1 != 'H81' & INDV2 != 'H81') %>%
  mutate(INDV_h89 = 'H89', INDV_2 = ifelse(INDV1 == 'H89', INDV2, INDV1)) %>% 
  left_join(pop.raw %>% select(Sample_ID, Landscape_ID), by = join_by('INDV_2' == 'Sample_ID')) %>%
  rename(Pop2 = Landscape_ID) %>%
  filter(Pop2 %in% c('NOM42', 'NOM09')) %>% group_by(Pop2) %>% summarise(phi_m = mean(RELATEDNESS_PHI))
