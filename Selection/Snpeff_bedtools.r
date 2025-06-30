library(dplyr)
library(tidyverse)
library(ggplot2)

LST_snpeff = read.table('Data/selection/LST_clean_ann.vcf', sep = '\t')

LST_snpeff_wide = LST_snpeff %>% separate_wider_delim(V3, delim = "|", names_sep = "_part", too_few = 'align_start')

LST_snpeff2 = LST_snpeff_wide[,1:10]
colnames(LST_snpeff2) = c('chr', 'pos', 'nuc', 'type', 'modifier', 'gen_loc', 'gen_loc2', 'trans', 'trans_loc', 'effect')
head(LST_snpeff2)

write.csv(LST_snpeff2, file = 'Results/LST_genes.csv', row.names = F)

LST_snpeff_plot = LST_snpeff2 %>% group_by(type) %>% summarise(n = n()) %>%
    arrange(n) %>% mutate(order = 1:nrow(.)) %>% mutate(type = ifelse(grepl('premature', type), '5_primer_UTR_PSCG_variant', type))


lst.plot = ggplot(LST_snpeff_plot, aes(n, as.factor(order))) +
  geom_bar(stat  = 'identity', col = 'black', fill = '#7b84ce')+
  theme_classic(base_size = 16)+
  labs(y = NULL, x = 'Count', title = "A")+
  theme(axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'))+
  scale_y_discrete(labels = gsub('_', ' ', LST_snpeff_plot$type))

#use plot from TopGO script
plot_grid(lst.plot, plot.lst.eff, ncol = 2, align = 'hv')

ggsave(file = 'LST_snpeff_bedtools_bar.pdf', height = 5.7, width = 13)

LST_bed_wide = LST_bed %>% separate_wider_delim(V1, delim = ";", names_sep = "p", too_few = 'align_start') %>% 
separate_wider_delim(V2, delim = ";", names_sep = "r", too_few = 'align_start') %>% select(1,2,3,6,7,8)

ncol(LST_bed_wide)
head(LST_bed_wide)

## misssense variants
missen = LST_snpeff2 %>% filter(type == 'missense_variant')
missense_variants = paste0(missen$chr,"_",missen$pos)


######## CROP

Crop_snpeff = read.table('Data/selection/Crop_clean_ann.vcf', sep = '\t')[,c(1,2,8)]

Crop_snpeff_wide = Crop_snpeff %>% separate_wider_delim(V8, delim = "|", names_sep = "_part", too_few = 'align_start')

Crop_snpeff2 = Crop_snpeff_wide[,1:10]
colnames(Crop_snpeff2) = c('chr', 'pos', 'nuc', 'type', 'modifier', 'gen_loc', 'gen_loc2', 'trans', 'trans_loc', 'effect')

write.csv(Crop_snpeff2, file = 'Results/Crop_genes.csv', row.names = F)


Crop_snpeff_plot = Crop_snpeff2 %>% group_by(type) %>% summarise(n = n()) %>%
  mutate(type = ifelse(grepl('premature', type), '5_primer_UTR_PSCG_variant', type)) %>% 
  right_join(data.frame(type = LST_snpeff_plot$type, order = LST_snpeff_plot$order), by = 'type') %>%
  mutate(n = ifelse(is.na(n), 0, n))

crop.plot = ggplot(Crop_snpeff_plot, aes(n, as.factor(order))) +
  geom_bar(stat  = 'identity', col = 'black', fill = '#7b84ce')+
  theme_classic(base_size = 16)+
  labs(y = NULL, x = 'Count', title = "B")+
  theme(axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'))+
  scale_y_discrete(labels = gsub('_', ' ', LST_snpeff_plot$type))+
  xlim(0,6)

library(cowplot)

plot_grid(lst.plot, crop.plot, nrow = 2, align = 'hv')

ggsave(file = "Results/snpeff_bar.png", height = 8, width = 6)

