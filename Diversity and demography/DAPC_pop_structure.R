######## Load Required Libraries ########
library("vcfR")         # For reading VCF files in R
library("adegenet")     # For DAPC analysis and genetic data structures

######## Set Working Directory ########
#setwd("xxxxxxx")  # Set your working directory (optional)

######## Read and Convert VCF File ########
vcf <- read.vcfR("Data/Goe_filtered095_Het60_pruned.vcf.gz")   # Load compressed VCF file
gl <- vcfR2genlight(vcf)                             # Convert VCF data to a genlight object for compatibility with `adegenet`

######## Load Population Data ########
pop.data <- read.table("populations_filtered.txt", sep = "\t", header = F)
colnames(pop.data) = c('sample', 'pop') 
#pop.data$pop <- as.factor(pop.data$pop)  # Convert population column to factor

length(unique(pop.data$pop))

sample.order = gl@ind.names

pop.data.ordered = pop.data[match(sample.order, pop.data$sample), ]
sample.order
pop.data.ordered

pop(gl) <- as.factor(pop.data.ordered$pop)

######## Perform Clustering and DAPC ########
grp <- find.clusters(gl, max.n.clust = 36)          # Set the maximum number of clusters (sites), looking for a plateau in BIC plot

######## Optimize the Number of Principal Components ########
#dapc <- dapc(gl, n.da = 100, n.pca = 100)           # Run DAPC with 100 discriminant functions and 100 PCs
dapc <- dapc(gl, grp$grp)           # Run DAPC with 100 discriminant functions and 100 PCs
ascore <- a.score(dapc)                               # Assess the quality of the DAPC (a-score)
optm.dapc <- optim.a.score(dapc)                         # Optimize the number of PCs based on the a-score

######## Save and Plot Final DAPC with Optimal PCs ########
pdf("DAPCascore35.pdf")                          # Save plot to PDF
pnw.dapc <- dapc(gl, n.pca = 35, n.da = 2)          # Final DAPC analysis with 35 PCs and 2 discriminant functions
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = FALSE, 
        posi.leg = "bottomleft", cleg = 0.75)       # Plot DAPC scatter with specified colors and layout
dev.off()                                           # Close the PDF device

# Optional: Plot DAPC scatter in R viewer with specific colors for populations
scatter(pnw.dapc, col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
                          "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
                          "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", 
                          "black"))