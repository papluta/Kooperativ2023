library(vcfR)       
library(adegenet)     


vcf <- read.vcfR("Data/variants_filtered_diversity.vcf.gz")   
gl <- vcfR2genlight(vcf) 

pop.data <- read.table("Data/populations.txt", sep = "\t", header = F)  %>%
  rename(sample = 1, pop = 2) %>% 
  mutate(sample = paste0(sample, "_", sample))

length(unique(pop.data$pop))

sample.order = gl@ind.names

pop.data.ordered = pop.data[match(sample.order, pop.data$sample), ]
sample.order
head(pop.data.ordered)

pop.data$pop <- as.factor(pop.data$pop)

length(unique(pop.data.ordered$pop))

pop(gl) <- pop.data.ordered$pop

######## Perform Clustering and DAPC ########
grp <- find.clusters(gl, max.n.clust = 35) 

######## Optimize the Number of Principal Components ########
dapc <- dapc(gl, grp$grp)        
ascore <- a.score(dapc)                         
optm.dapc <- optim.a.score(dapc)

######## Save and Plot Final DAPC with Optimal PCs ########
pdf("DAPCascore35.pdf") 
pnw.dapc <- dapc(gl, n.pca = 35, n.da = 2)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = FALSE, 
        posi.leg = "bottomleft", cleg = 0.75)
dev.off() 

scatter(pnw.dapc, col = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", 
                          "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", 
                          "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", 
                          "black"))