library(data.table)
library(tidyverse)
library(ggExtra)
library(cowplot)


lfmm <- read.table("LFMM/LFMM_results_0.05pct_ByPredictor.txt", sep = "\t", header = T)
gf_crop <- read.csv("gf_batches/crop_top1_importance.csv")
gf_lst <- read.csv("gf_batches/lst_top1_importance.csv")
gf_man <- read.csv("gf_batches/man_top1_importance.csv")

lfmm_perm <- lfmm %>% select(snp = SNP, Predictor, pval = PValue) %>% 
    pivot_wider(names_from = Predictor, values_from = pval)

lfmm_perm %>% 
    select(snp, q_value = crop_conven) %>% 
    filter(q_value <= threshold_crop) %>%
    inner_join(gf_crop, by = "snp")  %>% nrow()
lfmm_perm %>% 
    select(snp, q_value = LST_historical_daytime) %>%
    filter(q_value <= threshold_lst) %>%
    inner_join(gf_lst, by = "snp") %>% nrow()
lfmm_perm %>% 
    select(snp, q_value = man_made) %>%   
    filter(q_value <= threshold_man) %>%
    inner_join(gf_man, by = "snp")  %>% nrow()

methods_combined_crop <- lfmm_perm %>% 
    select(snp, crop_q_value = crop_conven) %>%    
    left_join(gf_crop, by = "snp")
methods_combined_lst <- lfmm_perm %>% 
    select(snp, LST_q_value = LST_historical_daytime) %>%    
    left_join(gf_lst, by = "snp")
methods_combined_man <- lfmm_perm` %>% 
    select(snp, man_q_value = man_made) %>%    
    left_join(gf_man, by = "snp")

head(methods_combined_crop)

crop_consensus <- methods_combined_crop %>% mutate(
    lfmm_flag = ifelse(crop_q_value <= 0.05, 1, 0),
    glm_flag = ifelse(!is.na(crop_p_value), 1, 0),
    gf_flag = ifelse(!is.na(crop_imp), 1, 0)) %>%
    mutate(no_flags = rowSums(across(ends_with("flag")))) %>%
    mutate(outlier = case_when(
        no_flags > 1 ~ "two or more",
        no_flags == 1 ~ "one",
        no_flags == 0 ~ "none",
        TRUE ~ "error"))


lst_consensus <- methods_combined_lst %>% mutate(
    lfmm_flag = ifelse(LST_q_value <= 0.05, 1, 0),
    glm_flag = ifelse(!is.na(lst_p_value), 1, 0),
    gf_flag = ifelse(!is.na(lst_imp), 1, 0)) %>%
    mutate(no_flags = rowSums(across(ends_with("flag")))) %>%
    mutate(outlier = case_when(
        no_flags > 1 ~ "two or more",
        no_flags == 1 ~ "one",
        no_flags == 0 ~ "none",
        TRUE ~ "error"))


man_consensus <- methods_combined_man %>% mutate(
    lfmm_flag = ifelse(man_q_value <= 0.05, 1, 0),
    glm_flag = ifelse(!is.na(man_p_value), 1, 0),
    gf_flag = ifelse(!is.na(man_imp), 1, 0)) %>%
    mutate(no_flags = rowSums(across(ends_with("flag")))) %>%
    mutate(outlier = case_when(
        no_flags > 1 ~ "two or more",
        no_flags == 1 ~ "one",
        no_flags == 0 ~ "none",
        TRUE ~ "error"))


threshold_crop = 0.000149854760244951
crop_consensus <- methods_combined_crop %>% mutate(
    outlier = case_when(
        crop_q_value <= threshold_crop & !is.na(crop_imp) ~ "consensus",
        crop_q_value > threshold_crop & !is.na(crop_imp) ~ "gf",
        crop_q_value <= threshold_crop & is.na(crop_imp) ~ "lfmm",
        crop_q_value > threshold_crop & is.na(crop_imp) ~ "none",
        TRUE ~ "error"))

threshold_lst = 0.000152254544035297
lst_consensus <- methods_combined_lst %>% mutate(
    outlier = case_when(
        LST_q_value <= threshold_lst & !is.na(lst_imp) ~ "consensus",
        LST_q_value > threshold_lst & !is.na(lst_imp) ~ "gf",
        LST_q_value <= threshold_lst & is.na(lst_imp) ~ "lfmm",
        LST_q_value > threshold_lst & is.na(lst_imp) ~ "none",
        TRUE ~ "error"))


threshold_man = 9.97999238592806e-05
man_consensus <- methods_combined_man %>% mutate(
    outlier = case_when(
        man_q_value <= threshold_man & !is.na(man_imp) ~ "consensus",
        man_q_value > threshold_man & !is.na(man_imp) ~ "gf",
        man_q_value <= threshold_man & is.na(man_imp) ~ "lfmm",
        man_q_value > threshold_man & is.na(man_imp) ~ "none",
        TRUE ~ "error"))


crop_plot_df <- crop_consensus %>%
    mutate(log_q = -log10(crop_q_value),
           chr = sub("(\\.1_).*",".1", snp),
          pos = sub(".*(\\.1_)","", snp) %>% sub("_.*","", .)) %>% 
    arrange(chr, pos) %>%
    mutate(chr_no = as.numeric(as.factor(chr))) %>%
    mutate(chr_no = ifelse(grepl("NW", chr), 18, chr_no)) %>% 
    group_by(chr_no) %>%
    mutate(pos_chr = 1:n()) %>% 
    ungroup() %>%
    arrange(chr_no, pos_chr) %>%
    mutate(pos_plot = 1:nrow(crop_consensus))

lst_plot_df <- lst_consensus %>%
    mutate(log_q = -log10(LST_q_value),
           chr = sub("(\\.1_).*",".1", snp),
          pos = sub(".*(\\.1_)","", snp) %>% sub("_.*","", .)) %>% 
    arrange(chr, pos) %>%
    mutate(chr_no = as.numeric(as.factor(chr))) %>%
    mutate(chr_no = ifelse(grepl("NW", chr), 18, chr_no)) %>% 
    group_by(chr_no) %>%
    mutate(pos_chr = 1:n()) %>% 
    ungroup() %>%
    arrange(chr_no, pos_chr) %>%
    mutate(pos_plot = 1:nrow(lst_consensus))

man_plot_df <- man_consensus %>%
    mutate(log_q = -log10(man_q_value),
           chr = sub("(\\.1_).*",".1", snp),
          pos = sub(".*(\\.1_)","", snp) %>% sub("_.*","", .)) %>% 
    arrange(chr, pos) %>%
    mutate(chr_no = as.numeric(as.factor(chr))) %>%
    mutate(chr_no = ifelse(grepl("NW", chr), 18, chr_no)) %>% 
    group_by(chr_no) %>%
    mutate(pos_chr = 1:n()) %>% 
    ungroup() %>%
    arrange(chr_no, pos_chr) %>%
    mutate(pos_plot = 1:nrow(man_consensus))


cutoff = 0
consensus_outlier <- man_plot_df %>% filter(outlier == "consensus")
gf_outlier <- man_plot_df %>% filter(outlier == "gf")
lfmm_outlier <- man_plot_df %>% filter(outlier == "lfmm")
none_outlier <- man_plot_df %>% filter(outlier == "none")

man_manhattan <- man_plot_df %>%  
  filter(log_q > cutoff) %>%
  ggplot(aes(x = pos_plot, y = log_q)) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 1) +
    geom_point(aes(fill = as.factor(chr_no)), size = 2, shape = 21, show.legend = F, stroke = NA) +
    geom_point(data = gf_outlier, col = "black", size = 2, alpha = 0.4) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 3) +
    scale_fill_manual(values = rep(c('#d0d0d0', '#ababab'), length.out = 18)) +
    labs(y = expression(atop(-log[10](italic(q)), "of man-made structure area")), x = NULL, col = 'Detected by', title = "C") +
#     scale_x_continuous(
#       breaks = chr_labels$center,
#       labels = chr_labels$chr_no,
#       expand = expansion(mult = c(0.01, 0.01))
#     ) +
    theme_classic(base_size = 22) +
    geom_hline(yintercept = -log10(threshold_man), col = '#c21328', linewidth = 1.5 ) +
    theme(
    axis.text.x = element_blank(), axis.text.y = element_text(size = 14, color = 'black'), 
    panel.grid = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = 'none'
    )

man_manhattan_hist <- ggMarginal(man_manhattan, type = "histogram", margins = "x", xparams = list(bins=60), fill = "#306EA1")

# ggsave("260308_man_manhattan.png", plot = man_manhattan_hist, width = 18, height = 5, dpi = 240)

cutoff = 0
consensus_outlier <- crop_plot_df %>% filter(outlier == "consensus")
gf_outlier <- crop_plot_df %>% filter(outlier == "gf")
lfmm_outlier <- crop_plot_df %>% filter(outlier == "lfmm")
none_outlier <- crop_plot_df %>% filter(outlier == "none")

crop_manhattan <- crop_plot_df %>%  
  filter(log_q > cutoff) %>%
  ggplot(aes(x = pos_plot, y = log_q)) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 1) +
    geom_point(aes(fill = as.factor(chr_no)), size = 2, shape = 21, show.legend = F, stroke = NA) +
    geom_point(data = gf_outlier, col = "black", size = 2, alpha = 0.4) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 3) +
    scale_fill_manual(values = rep(c('#d0d0d0', '#ababab'), length.out = 18)) +
    labs(y = expression(atop(-log[10](italic(q)), "of conventional crop area")), x = NULL, col = 'Detected by', title = "A") +
#     scale_x_continuous(
#       breaks = chr_labels$center,
#       labels = chr_labels$chr_no,
#       expand = expansion(mult = c(0.01, 0.01))
#     ) +
    theme_classic(base_size = 22) +
    geom_hline(yintercept = -log10(threshold_crop), col = '#c21328', linewidth = 1.5 ) +
    theme(
    axis.text.x = element_blank(), axis.text.y = element_text(size = 14, color = 'black'), 
    panel.grid = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = 'none'
    )

crop_manhattan_hist <- ggMarginal(crop_manhattan, type = "histogram", margins = "x", xparams = list(bins=60), fill = "#306EA1")

# ggsave("260308_crop_manhattan.png", plot = crop_manhattan_hist, width = 18, height = 5, dpi = 240)

cutoff = 0
consensus_outlier <- lst_plot_df %>% filter(outlier == "consensus")
gf_outlier <- lst_plot_df %>% filter(outlier == "gf")
lfmm_outlier <- lst_plot_df %>% filter(outlier == "lfmm")
none_outlier <- lst_plot_df %>% filter(outlier == "none")

lst_manhattan <- lst_plot_df %>%  
  filter(log_q > cutoff) %>%
  ggplot(aes(x = pos_plot, y = log_q)) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 1) +
    geom_point(aes(fill = as.factor(chr_no)), size = 2, shape = 21, show.legend = F, stroke = NA) +
    geom_point(data = gf_outlier, col = "black", size = 2, alpha = 0.4) +
    geom_point(data = consensus_outlier, col = "#306EA1", size = 3) +
    scale_fill_manual(values = rep(c('#d0d0d0', '#ababab'), length.out = 18)) +
    labs(y = expression(atop(-log[10](italic(q)), "of historical LST")), x = NULL, col = 'Detected by', title = "B") +
#     scale_x_continuous(
#       breaks = chr_labels$center,
#       labels = chr_labels$chr_no,
#       expand = expansion(mult = c(0.01, 0.01))
#     ) +
    theme_classic(base_size = 22) +
    geom_hline(yintercept = -log10(threshold_lst), col = '#c21328', linewidth = 1.5 ) +
    theme(
    axis.text.x = element_blank(), axis.text.y = element_text(size = 14, color = 'black'), 
    panel.grid = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = 'none'
    )

lst_manhattan_hist <- ggMarginal(lst_manhattan, type = "histogram", margins = "x", xparams = list(bins=60), fill = "#306EA1")

# ggsave("260308_lst_manhattan.png", plot = lst_manhattan_hist, width = 18, height = 5, dpi = 240)

combined_manhattan <- plot_grid(crop_manhattan_hist, lst_manhattan_hist, man_manhattan_hist, ncol = 1)

ggsave("260308_combined_manhattan.png", plot = combined_manhattan, width = 18, height = 15)

crop_outliers <- crop_plot_df %>% filter(outlier == "consensus") %>%
    select(chr, pos)
lst_outliers <- lst_plot_df %>% filter(outlier == "consensus") %>%
    select(chr, pos)
man_outliers <- man_plot_df %>% filter(outlier == "consensus") %>%
    select(chr, pos)

nrow(crop_outliers)
nrow(lst_outliers)
nrow(man_outliers)
head(crop_outliers)

write.table(crop_outliers, file = "260308crop_consensus.txt", row.names = F, col.names = F, quote = F,
            eol = "\n", sep = "\t")
write.table(lst_outliers, file = "260308lst_consensus.txt", row.names = F, col.names = F, quote = F,
            eol = "\n", sep = "\t")
write.table(man_outliers, file = "260308man_consensus.txt", row.names = F, col.names = F, quote = F,
            eol = "\n", sep = "\t")

crop_outliers <- read.table("260301crop_consensus.txt", header = F, sep = "\t")
lst_outliers <- read.table("260301lst_consensus.txt", header = F, sep = "\t")
man_outliers <- read.table("260301man_consensus.txt", header = F, sep = "\t")

common = crop_outliers %>%
    mutate(crop = 1) %>% 
    full_join(lst_outliers%>%
    mutate(lst = 1), by = c("V1", "V2")) %>% 
    full_join(man_outliers %>%
    mutate(man = 1), by = c("V1", "V2")) %>%
    select(-c(V1, V2))

rowSums(common, na.rm = T)

man_summary <- c(gf = sum(man_consensus$outlier == "consensus") + sum(man_consensus$outlier == "gf"),   
  lfmm = sum(man_consensus$outlier == "consensus") + sum(man_consensus$outlier == "lfmm"),      
  consensus = sum(man_consensus$outlier == "consensus"))
crop_summary <- c(gf = sum(crop_consensus$outlier == "consensus") + sum(crop_consensus$outlier == "gf"),   
  lfmm = sum(crop_consensus$outlier == "consensus") + sum(crop_consensus$outlier == "lfmm"),      
  consensus = sum(crop_consensus$outlier == "consensus"))
lst_summary <- c(gf = sum(lst_consensus$outlier == "consensus") + sum(lst_consensus$outlier == "gf"),   
  lfmm = sum(lst_consensus$outlier == "consensus") + sum(lst_consensus$outlier == "lfmm"),      
  consensus = sum(lst_consensus$outlier == "consensus"))
rbind(crop = crop_summary, lst = lst_summary, man = man_summary)

