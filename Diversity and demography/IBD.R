library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(vegan)
library(geosphere)
library(ggplot2)
library(effects)


###### ISOLATION BY DISTANCE (FST VS GEOGRAPHIC DISTANCE)

fst_stampp = read.csv('Data/genetic_distances_bootstraps.csv', header = T)

fst_stampp2 = fst_stampp[, c(2,3,104:107)]


mean(fst_stampp2$Fst)
mean(fst_stampp2$Lower.bound.CI.limit)
mean(fst_stampp2$Upper.bound.CI.limit)

land = read.csv('Data/coordinates.csv', header = T) %>% 
  filter(pop != 'NOM16')

# Compute pairwise geographic distances 
geo_dist_matrix <- distm(land[, c("x", "y")]) / 1000  # Convert meters to km
rownames(geo_dist_matrix) <- land$pop
colnames(geo_dist_matrix) <- land$pop
geo_dist_matrix

bootstrap_data <- fst_stampp2[,-c(3:5)]

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
  fst_matrix[pop2, pop1] <- fst_value  
}


# Transform negative Fst values to 0
fst_matrix[fst_matrix < 0] = 0

# Linearize Fst using Fst / (1 - Fst)
fst_matrix <- fst_matrix / (1 - fst_matrix)

# Log-transform geographic distance
geo_res_matrix = as.matrix(res)
geo_dist_matrix <- log(geo_dist_matrix)

mantel_result_res <- mantel(as.dist(geo_res_matrix), as.dist(fst_matrix), method = "pearson", permutations = 9999)
mantel_result_land <- mantel(as.dist(geo_dist_matrix), as.dist(fst_matrix), method = "pearson", permutations = 9999)
mantel_result_land <- mantel(as.dist(geo_dist_matrix), as.dist(fst_matrix), method = "pearson", permutations = 9999)

geo_dist_values <- as.vector(as.dist(geo_dist_matrix))
fst_values <- as.vector(as.dist(fst_matrix))

plot_data <- data.frame(GeographicDistance = geo_dist_values, GeneticDistance = fst_values)

# Plot Isolation by Distance
ggplot(plot_data, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.7, shape = 21) +
  geom_smooth(method = "lm",  color = "black") +
  labs(x = "Resistance", y = "Genetic Distance (Fst)",
       title = NULL) +
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))+
  annotate("text", x = 1.7e+5, y = 0.004, label = 'r = 0.43, p = < 0.01')

write.csv(plot_data, "Results/Goe_IBD_exp_Fst.csv", row.names = TRUE)

ggsave(file = 'Results/Goe_IBR.png', height = 5, width = 5)

mite.correlog3 <-  mantel.correlog(as.dist(fst_matrix), as.dist(geo_dist_matrix), nperm=1000, n.class= 30, mult = "bonferroni")
mite.correlog3
plot(mite.correlog3)

