library(data.table)
library(gradientForest)
library(tidyverse)


frq <- fread("allele_freq.txt", header = TRUE, sep = "\t")  
dim(frq)                                               


# Extract Genotype Data and Environmental Data
gen <- frq[, 2:ncol(frq)]
# head(gen)
gen[,1]

# Environmental Variables
land <- read.csv("landuse_1000m.csv")  
xy <- read.csv("coordinates.csv")
lst <- read.csv("lst.csv")
xy = xy[xy$pop != 'NOM16',]
land = land[land$Landscape != 'NOM16',]

xy = lst %>% left_join(xy, by = "pop") %>% rename(Landscape = pop)

env <- land %>% left_join(xy, by = "Landscape") %>% mutate(man_made = road + urban)

head(env)
nrow(env)

head(imp_matrix)
ncol(gen)


env = env %>% select(crop_conven, LST_historical_daytime, x, y, man_made)
head(env)
cor(env)

batch_size <- 100000
num_batches <- ceiling(ncol(gen) / batch_size)
num_batches

maxLevel <- ceiling(log2(0.368*nrow(env)/2))
maxLevel


for (i in 3:num_batches) {
  cat("Running batch", i, "of", num_batches, "\n")
  
  batch_snps <- gen[, ((i-1)*batch_size + 1) : min(i*batch_size, ncol(gen))]
  
  gf_model <- gradientForest(
    data = data.frame(env, batch_snps),
    predictor.vars = colnames(env),
    response.vars = colnames(batch_snps),
    ntree = 500,
    maxLevel = maxLevel,
    trace = FALSE
  )
  
    imp_weighted <- gf_model$imp.rsq  
    
  saveRDS(imp_weighted, file = file.path(output_dir, paste0("batch_", i, ".rds")))
  
  rm(gf_model, batch_snps, imp_weighted); gc()
}



# Combine all batches at the end

imp_all <- do.call(cbind, lapply(1:num_batches, function(i) {
  readRDS(file.path(output_dir, paste0("batch_", i, ".rds")))
}))

b1 <- readRDS(file.path(output_dir, paste0("batch_", 1, ".rds")))

dim(imp_all)

imp_all_t = t(imp_all)
imp_all_df = as.data.frame(imp_all_t)

dim(imp_all_t)
head(imp_all_df)

top_1_n <- round(nrow(imp_all_df) * 0.01)
top_1_n

crop_top1 <- arrange(imp_all_df, crop_conven) %>% top_n(top_1_n, crop_conven) %>% select(crop_imp = crop_conven) %>%
    rownames_to_column(var = "snp")
lst_top1 <- arrange(imp_all_df, crop_conven) %>% top_n(top_1_n, LST_historical_daytime) %>% select(lst_imp = LST_historical_daytime) %>%
    rownames_to_column(var = "snp")
man_top1 <- arrange(imp_all_df, crop_conven) %>% top_n(top_1_n, man_made) %>% select(man_imp = man_made) %>%
    rownames_to_column(var = "snp")

head(crop_top1)

write.csv(crop_top1, file = "gf_batches/crop_top1_importance.csv", row.names = F)
write.csv(lst_top1, file = "gf_batches/lst_top1_importance.csv", row.names = F)
write.csv(man_top1, file = "gf_batches/man_top1_importance.csv", row.names = F)
