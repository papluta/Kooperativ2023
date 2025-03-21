library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(vegan)
library(geosphere)
library(ggplot2)
library(effects)


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
# Perform Mantel test
# Linearize Fst using Fst / (1 - Fst)
fst_matrix <- fst_matrix / (1 - fst_matrix)

# Log-transform geographic distance
geo_dist_matrix <- log(geo_dist_matrix)


mantel_result <- mantel(as.dist(geo_dist_matrix), as.dist(fst_matrix), method = "pearson", permutations = 9999)
a = colnames(geo_dist_matrix)
b = colnames(fst_matrix)

a
b# Print Mantel test result

print(mantel_result)

# Extract pairwise values for plotting
geo_dist_values <- as.vector(as.dist(geo_dist_matrix))
fst_values <- as.vector(as.dist(fst_matrix))

# Create a data frame for plotting
plot_data <- data.frame(GeographicDistance = geo_dist_values, GeneticDistance = fst_values)

# Plot Isolation by Distance
ggplot(plot_data, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Geographic Distance (km)", y = "Genetic Distance (Fst)",
       title = paste("Isolation by Distance (r =", round(mantel_result$statistic, 2), ", p =", round(mantel_result$signif, 3), ")")) +
  theme_minimal()

write.csv(plot_data, "AF_Fst_Geo.csv", row.names = TRUE)

ggplot(pi, aes(pop, avg_pi)) +
  geom_boxplot()
