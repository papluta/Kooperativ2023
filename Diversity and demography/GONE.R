library(tidyverse)

gone = read.table("Data/Output_Ne_Goe_demography_filtered_GONE", sep = '\t', header = T)


ggplot(gone[11:200,], aes(Generation, Geometric_mean))+
  annotate("rect", xmin = 131, xmax = 146, ymin = 0, ymax = max(gone[11:200,2]), alpha = 0.4, fill = 'grey')+
  geom_line(size = 1.3, col = "#6f4bb3")+
  #geom_vline(aes(xintercept = 147), linetype  = 'dashed', size = 1.3)+
  #geom_vline(aes(xintercept = 131), linetype  = 'dashed', size = 1.3)+
  #geom_vline(aes(xintercept = 7.3), linetype  = 'dashed', size = 1.3)+
  theme_bw(base_size = 16)+
  labs(y = expression('Mean Ne'), x = "Generation (back in time)")+
  scale_y_continuous()+
  theme(axis.text = element_text(color = 'black'), panel.grid = element_blank())

ggsave(file = "Results/Gone_allpops.png", width = 8, height = 4)
