library(topGO)
library(cowplot)
library(tidyverse)
library(VennDiagram)
library(grid)

## file prep for TopGO

Crop_raw  = read.table("Data/bedtools_crop_consensus.txt", header = F, sep = "\t")
Crop = Crop_raw %>% 
  mutate(gbkey = str_extract(V12, "gbkey=.*;"), ID = str_extract(V12, "ID=.*;"), 
         Parent = str_extract(V12, "Parent=.*;"), Genbank = str_extract(V12, "GenBank:.*;"), 
         Product = str_extract(V12, "product=.*;")) %>%
  mutate(across(gbkey:Product, function(x) sub(";.*","", x))) %>%
  mutate(across(gbkey:Product, function(x) sub(".*=","", x))) %>%
  mutate(Genbank = sub("GenBank:", "", Genbank)) %>%
  dplyr::select(-c(V3, V4, V9:V12)) %>%
  dplyr::rename(Chromosome = 1, Position = 2, Database = 3, Type = 4, Start = 5, End = 6)


LST_raw  = read.table("Data/bedtools_lst_consensus.txt", sep = "\t")
LST = LST_raw %>% 
  mutate(gbkey = str_extract(V12, "gbkey=.*;"), ID = str_extract(V12, "ID=.*;"), 
         Parent = str_extract(V12, "Parent=.*;"), Genbank = str_extract(V12, "GenBank:.*;"), 
         Product = str_extract(V12, "product=.*;")) %>%
  mutate(across(gbkey:Product, function(x) sub(";.*","", x))) %>%
  mutate(across(gbkey:Product, function(x) sub(".*=","", x))) %>%
  mutate(Genbank = sub("GenBank:", "", Genbank)) %>%
  dplyr::select(-c(V3, V4, V9:V12)) %>%
  dplyr::rename(Chromosome = 1, Position = 2, Database = 3, Type = 4, Start = 5, End = 6)

urban_raw  = read.table("Data/bedtools_urban_consensus.txt", sep = "\t")
urban = urban_raw %>% 
  mutate(gbkey = str_extract(V12, "gbkey=.*;"), ID = str_extract(V12, "ID=.*;"), 
         Parent = str_extract(V12, "Parent=.*;"), Genbank = str_extract(V12, "GenBank:.*;"), 
         Product = str_extract(V12, "product=.*;")) %>%
  mutate(across(gbkey:Product, function(x) sub(";.*","", x))) %>%
  mutate(across(gbkey:Product, function(x) sub(".*=","", x))) %>%
  mutate(Genbank = sub("GenBank:", "", Genbank)) %>%
  dplyr::select(-c(V3, V4, V9:V12)) %>%
  dplyr::rename(Chromosome = 1, Position = 2, Database = 3, Type = 4, Start = 5, End = 6)


geneID2GO <- readMappings("Data/merged_GO_terms.map") 

geneNames <- names(geneID2GO) 
head(geneNames)


LST_out <- LST %>%
  filter(grepl("cds", ID))
LST_ids <- paste0("lcl|", LST_out$Chromosome,"_",sub("-","_", LST_out$ID))

Crop_out <- Crop %>%
  filter(grepl("cds", ID))
Crop_ids <- paste0("lcl|", Crop_out$Chromosome,"_",sub("-","_", Crop_out$ID))

Urban_out <- urban %>%
  filter(grepl("cds", ID))
Urban_ids <- paste0("lcl|", Urban_out$Chromosome,"_",sub("-","_", Urban_out$ID))

go_function <- function(gene_id) {
geneList <- factor(as.integer(geneNames %in% gene_id)) 
names(geneList) <- geneNames

missing_genes <- gene_id[!gene_id %in% geneNames]
print(missing_genes)

GOdata_bio_proc <- new("topGOdata", description = "GO analysis BP", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata_bio_proc, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(GOdata_bio_proc, algorithm = "elim", statistic = "fisher")


goEnrichmentBP <- GenTable(GOdata_bio_proc,
                           elimFisher = resultFisher.elim,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500,
                           numChar = 200)

sig_goEnrichment_BP = goEnrichmentBP %>% mutate(p = as.numeric(elimFisher)) %>% filter(p < 0.05)
sig_goEnrichment_BP$Term = reorder(sig_goEnrichment_BP$Term, -log10(sig_goEnrichment_BP$p), FUN = sum, descending = TRUE)
sig_goEnrichment_BP = sig_goEnrichment_BP %>% ungroup() %>% arrange(desc(-log10(p))) 



####
GOdata_mol_func <- new("topGOdata", description = "GO analysis BP MF", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata_mol_func, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(GOdata_mol_func, algorithm = "elim", statistic = "fisher")

goEnrichment_MF <- GenTable(GOdata_mol_func,
                            elimFisher = resultFisher.elim,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500,
                           numChar = 200)

sig_goEnrichment_MF = goEnrichment_MF %>% mutate(p = as.numeric(elimFisher)) %>% filter(p < 0.05) %>%
  mutate(Term = sub(",.*", "", Term))
sig_goEnrichment_MF$Term = reorder(sig_goEnrichment_MF$Term, -log10(sig_goEnrichment_MF$p), FUN = sum, descending = TRUE)
sig_goEnrichment_MF = sig_goEnrichment_MF %>% ungroup() %>% arrange(desc(-log10(p))) 



GOdata_cellular_comp <- new("topGOdata", description = "GO analysis BP CC", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata_cellular_comp, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(GOdata_cellular_comp, algorithm = "elim", statistic = "fisher")

goEnrichment_CC <- GenTable(GOdata_cellular_comp,
                            elimFisher = resultFisher.elim,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500,
                           numChar = 200)

sig_goEnrichment_CC = goEnrichment_CC %>% mutate(p = as.numeric(elimFisher)) %>% filter(p < 0.05)
sig_goEnrichment_CC$Term = reorder(sig_goEnrichment_CC$Term, -log10(sig_goEnrichment_CC$p), FUN = sum, descending = TRUE)
sig_goEnrichment_CC = sig_goEnrichment_CC %>% ungroup() %>% arrange(desc(-log10(p)))

go_list <- list(BP = sig_goEnrichment_BP,
                MF = sig_goEnrichment_MF,
                CC = sig_goEnrichment_CC)
return(go_list)
}

lst_go_terms <- go_function(gene_id = LST_ids)
lst_go_terms_long <- bind_rows(lst_go_terms, .id = "go_type")
nrow(lst_go_terms_long)

crop_go_terms <- go_function(gene_id = Crop_ids)
crop_go_terms_long <- bind_rows(crop_go_terms, .id = "go_type")
nrow(crop_go_terms_long)
urban_go_terms <- go_function(gene_id = Urban_ids)
urban_go_terms_long <- bind_rows(urban_go_terms, .id = "go_type")
nrow(urban_go_terms_long)

write.csv(lst_go_terms_long, file = "Results/260308_lst_go_terms.csv", row.names = F)
write.csv(crop_go_terms_long, file = "Results/260308_crop_go_terms.csv", row.names = F)
write.csv(urban_go_terms_long, file = "Results/260308_urban_go_terms.csv", row.names = F)

plot.lst.BP = ggplot(sig_goEnrichment_BP_classic[1:10,], aes(-log10(p), Term))+
  geom_bar(stat = "identity", col = 'black', fill = "#84ce7b")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black', size = 14))+
  labs(x = expression(paste("-log"[10], " p-value")), y = NULL, title = "A")+
  scale_x_continuous(expand = c(0,0.1,0,0))

plot.lst.MF = ggplot(sig_goEnrichment_MF[1:10,], aes(-log10(p), Term))+
  geom_bar(stat = "identity", col = 'black', fill = "#84ce7b")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black', size = 14))+
  labs(x = expression(paste("-log"[10], " p-value")), y = NULL, title = "B")+
  scale_x_continuous(expand = c(0,0.1,0,0))

plot.lst.CC = ggplot(sig_goEnrichment_CC[1:10,], aes(-log10(p), Term))+
  geom_bar(stat = "identity", col = 'black', fill = "#84ce7b")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black', size = 14))+
  labs(x = expression(paste("-log"[10], " p-value")), y = NULL, title = "C")+
  scale_x_continuous(expand = c(0,0.1,0,0))

plot_grid(plot.lst.BP, plot.lst.MF, plot.lst.CC, nrow = 3, align = 'hv')

ggsave(file = "Results/GO_bar.pdf", height = 10, width = 10)

rbind(sig_goEnrichment_BP_classic, sig_goEnrichment_MF[1:10,], sig_goEnrichment_CC) %>%
  filter(grepl("...", Term))


plot.lst.BP = ggplot( sig_goEnrichment_BP_classic, aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#84ce7b")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 11))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "C")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot.lst.MF = ggplot(sig_goEnrichment_MF[1:10,],aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#7b84ce")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 11))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "D")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot.lst.CC = ggplot(sig_goEnrichment_CC, aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#ce7b84")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 11))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "E")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot_grid(plot.lst.BP, plot.lst.MF, plot.lst.CC, nrow = 1, align = 'hv')
ggsave(file = "Results/260218GO_bar.pdf", height = 8, width = 11)

sig_goEnrichment_BP_classic$type = "Biological process"
sig_goEnrichment_MF$type = "Molecular function"
sig_goEnrichment_CC$type = "Cellular component"

sig_all = rbind(sig_goEnrichment_BP_classic[,-c(6,7,8)], sig_goEnrichment_MF[,-c(6,7,8)],sig_goEnrichment_CC[,-c(6,7,8)])

write.csv(sig_all, file = 'Results/260218GO_terms_all.csv', row.names = F)

## venn

A <- crop_go_terms_long$GO.ID
B <- lst_go_terms_long$GO.ID
C <- urban_go_terms_long$GO.ID

overlaps <- calculate.overlap(list(A, B, C))

png("venn_diagramm.png")

venn.plot <- venn.diagram(
  x = list(A, B, C),
  category.names = c("Conventional\nagriculture", 
                     "Historical LST", 
                     "Man-made\nstructures"),
  filename = NULL
)

grid.newpage()
grid.draw(venn.plot)

add_labels <- function(features, x, y) {
  if (length(features) > 0 && length(features) < 10) {
    grid.text(paste(features, collapse = "\n"),
              x = x, y = y,
              gp = gpar(fontsize = 8))
  }
}

dev.off()
