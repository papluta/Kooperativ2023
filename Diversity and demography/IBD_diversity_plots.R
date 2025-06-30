library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(vegan)
library(geosphere)
library(ggplot2)
library(effects)


###### ISOLATION BY DISTANCE (FST VS GEOGRAPHIC DISTANCE)

fst_stampp = read.csv('Data/Goe_filtered095_Het60_pruned_genetic_distances_bootstraps.csv', header = T)

fst_stampp2 = fst_stampp[, c(2,3,104:107)] %>% filter(Population1 != 'NOM16' & Population2 != 'NOM16')


mean(fst_stampp2$Fst)
mean(fst_stampp2$Lower.bound.CI.limit)
mean(fst_stampp2$Upper.bound.CI.limit)

land = read.csv('Data/Kooperativ2023_land.csv', header = T) %>% 
  filter(Landscape_ID != 'NOM16' &  Landscape_ID != 'NOM11' &  Landscape_ID != 'NOM24')

# Compute pairwise geographic distances (great-circle distances in km)
geo_dist_matrix <- distm(land[, c("x", "y")]) / 1000  # Convert meters to km
rownames(geo_dist_matrix) <- land$Landscape_ID
colnames(geo_dist_matrix) <- land$Landscape_ID
geo_dist_matrix

# Load genetic distance bootstrap data
bootstrap_data <- fst_stampp2[,-c(3:5)]

# Ensure Population1 and Population2 columns exist
colnames(bootstrap_data) <- c("Population1", "Population2", "Fst")

# Create an empty Fst matrix
populations <- sort(unique(c(bootstrap_data$Population1, bootstrap_data$Population2)))
fst_matrix <- matrix(NA, nrow = length(populations), ncol = length(populations),
                     dimnames = list(populations, populations))

# Fill in the matrix with Fst values
for (i in 1:nrow(bootstrap_data)) {
  pop1 <- bootstrap_data$Population1[i]
  pop2 <- bootstrap_data$Population2[i]
  fst_value <- bootstrap_data$Fst[i]
  fst_matrix[pop1, pop2] <- fst_value
  fst_matrix[pop2, pop1] <- fst_value  # Ensure symmetry
}

# Save the ordered Fst matrix
#write.csv(fst_matrix, "AV_Fst_matrix_ordered.csv", row.names = TRUE)
fst_matrix

# Transform negative Fst values to 0
fst_matrix[fst_matrix < 0] = 0

# Perform Mantel test
# Linearize Fst using Fst / (1 - Fst)
fst_matrix <- fst_matrix / (1 - fst_matrix)

# Log-transform geographic distance
geo_dist_matrix <- log(geo_dist_matrix)


mantel_result <- mantel(as.dist(geo_dist_matrix), as.dist(fst_matrix), method = "pearson", permutations = 9999)
a = colnames(geo_dist_matrix)
b = colnames(fst_matrix)

# Print Mantel test result

print(mantel_result)

# Extract pairwise values for plotting
geo_dist_values <- as.vector(as.dist(geo_dist_matrix))
fst_values <- as.vector(as.dist(fst_matrix))

# Create a data frame for plotting
plot_data <- data.frame(GeographicDistance = geo_dist_values, GeneticDistance = fst_values)

# Plot Isolation by Distance
ggplot(plot_data, aes(x = exp(GeographicDistance), y = GeneticDistance)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_smooth(method = "lm",  color = "black") +
  labs(x = "Geographic Distance (km)", y = "Genetic Distance (Fst)",
       title = NULL) +
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))+
  geom_text(aes(x = 21, y = 0.004), label = 'r = 0.35, p = < 0.01')

write.csv(plot_data, "Results/Goe_IBD_exp_Fst.csv", row.names = TRUE)

ggsave(file = 'Results/Goe_IBD_exp.pdf', height = 5, width = 5)

mite.correlog3 <-  mantel.correlog(as.dist(fst_matrix), as.dist(geo_dist_matrix), nperm=1000, n.class= 30, mult = "bonferroni")
mite.correlog3
plot(mite.correlog3)


save('Results/correlogram.png', height = 4, width = 8)
########## HETEROZYGOSITY

het = read.csv('Data/Goe_filtered095_Het60_pruned_heterozygosity.csv') %>% filter(pop != 'NOM16')

mean(het$uHe)

het$label = ifelse(het$uHe > mean(het$uHe)+sd(het$uHe) |
                     het$uHe < mean(het$uHe)-sd(het$uHe), het$pop, NA)
het %>% 
  ggplot(aes(1, uHe))+
  geom_boxplot()+
  geom_point(aes(1, uHe, col = pop), size = 2, alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none')+
  geom_text(aes(label = label, col = pop), position = position_jitterdodge(dodge.width = 0.05, jitter.width = 0.2))+
  labs(x = NULL, y = 'Unbiased Heterogyzosity\n ', title = 'Average uHe = 0.228')+
  scale_x_continuous(expand = c(0.15,0,0,0.15))

mean + Z * (SD/sqrt(n))

het.ci = het %>% mutate(CI_upper = uHe + 1.960 * (uHeSD/sqrt(nLoc)),
               CI_lower = uHe - 1.960 * (uHeSD/sqrt(nLoc)))

het %>% ggplot(aes(pop, uHe))+
  geom_pointrange(aes(ymin = uHe - uHeSD, ymax = uHe + uHeSD))

het.ci %>% ggplot(aes(pop, uHe))+
  geom_pointrange(aes(ymin = CI_lower, ymax = CI_upper))

ggsave(file = 'Results/Goe_He_boxplot.png', height = 5, width = 5)
