
library(lfmm)
library(data.table)
library(qvalue)
library(dplyr)


frq <- fread("allele_freq.txt", header = TRUE, sep = "\t")
dim(frq)


gen <- frq[, 2:ncol(frq)]
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

# how many Ks?
pc <- prcomp(gen)
plot(pc$sdev[1:20]^2, xlab = "PC", ylab = "Variance explained")

pred_interest <- c("man_made", "LST_historical_daytime", "crop_conven")
X_env <- env %>%
    transmute(across(all_of(c("crop_conven", "LST_historical_daytime", "man_made")), ~ scale(.x)))
X_env

lfmm_all <- lfmm_ridge(Y = gen, X = X_env, K = 2)
lfmm.pv <- lfmm_test(Y = gen, X = X_env, lfmm = lfmm_all, calibrate = "gif")
lfmm.pv$gif 

lfmm.qv <- qvalue(lfmm.pv$calibrated.pvalue)$qvalues
lfmm.qv

lfmm_full <- as.data.frame(lfmm.qv)
lfmm_full$snp = rownames(lfmm.qv)
rownames(lfmm_full) = NULL

write.csv(lfmm_full, "LFMM/260306_lfmm_full_k2.csv", row.names = FALSE)

perm_pvals <- lapply(1:n_perm, function(i) {
  if (i %% 10 == 0) cat("  perm", i, "\n")
  
  Xp <- X_env
  idx <- sample.int(nrow(Xp))
  
  # Jointly permute all 3 predictors together (same row shuffle)
  Xp[,] <- Xp[idx, , drop = FALSE]
  
  modp <- lfmm_ridge(Y = gen, X = Xp, K = K)
  testp <- lfmm_test(Y = gen, X = Xp, lfmm = modp, calibrate = "gif")
  
  pp <- testp$calibrated.pvalue
  if (is.vector(pp)) pp <- matrix(pp, ncol = ncol(Xp))
  colnames(pp) <- colnames(Xp)
  
  pp
})


## Empirical thresholds per predictor (multiple tails)

thr_probs <- c("0.2pct" = 0.002, "0.1pct" = 0.001, "0.05pct" = 0.0005)

thr_by_pred <- rbindlist(lapply(names(thr_probs), function(lbl) {
  pr <- thr_probs[[lbl]]
  data.table(
    Predictor = pred_interest,
    Tail = lbl,
    Prob = pr,
    Threshold = sapply(pred_interest, function(p) {
      vec <- unlist(lapply(perm_pvals, function(pp) pp[, p]))
      as.numeric(quantile(vec, probs = pr, na.rm = TRUE))
    })
  )
}))

thr_file <- file.path(out_dir, "LFMM_Empirical_Thresholds_ByPredictor_MULTI.txt")
fwrite(thr_by_pred, thr_file, sep = "\t")
print(thr_by_pred)

chosen_tail <- "0.05pct"  

thr_chosen <- thr_by_pred[Tail == chosen_tail]
thr_map <- setNames(thr_chosen$Threshold, thr_chosen$Predictor)

results <- rbindlist(lapply(pred_interest, function(p) {
  data.table(
    SNP = colnames(Y),
    Predictor = p,
    PValue = p0[, p],
    Threshold = thr_map[[p]],
    Significant = p0[, p] <= thr_map[[p]]
  )
}), fill = TRUE)

res_file <- file.path(out_dir, paste0("LFMM_results_", chosen_tail, "_ByPredictor.txt"))
fwrite(results, res_file, sep = "\t")

thr_chosen <- thr_by_pred[thr_by_pred$Tail == chosen_tail]
thr_map <- setNames(thr_chosen$Threshold, thr_chosen$Predictor)

results <- rbindlist(lapply(pred_interest, function(p) {
  data.table(
    SNP = colnames(Y),
    Predictor = p,
    PValue = p0[, p],
    Threshold = thr_map[[p]],
    Significant = p0[, p] <= thr_map[[p]]
  )
}), fill = TRUE)

sig_only <- results[Significant == TRUE]

fwrite(
  sig_only,
  file.path(out_dir, paste0("LFMM_Significant_SNPs_", chosen_tail, "_ByPredictor.txt")),
  sep = "\t"
)