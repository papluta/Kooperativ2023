# Load necessary library
library(topGO)
library(dplyr)
library(ggplot2)
library(cowplot)

## file prep for TopGO

Crop_raw  = read.table("Data/selection/Crop_out.txt", sep = "\t")
Crop = Crop_raw %>% 
  mutate(gbkey = str_extract(V12, "gbkey=.*;"), ID = str_extract(V12, "ID=.*;"), 
         Parent = str_extract(V12, "Parent=.*;"), Genbank = str_extract(V12, "GenBank:.*;"), 
         Product = str_extract(V12, "product=.*;")) %>%
  mutate(across(gbkey:Product, function(x) sub(";.*","", x))) %>%
  mutate(across(gbkey:Product, function(x) sub(".*=","", x))) %>%
  mutate(Genbank = sub("GenBank:", "", Genbank)) %>%
  dplyr::select(-c(V3, V4, V9:V12)) %>%
  dplyr::rename(Chromosome = 1, Position = 2, Database = 3, Type = 4, Start = 5, End = 6)


LST_raw  = read.table("Data/selection/LST_out.txt", sep = "\t")
LST = LST_raw %>% 
  mutate(gbkey = str_extract(V12, "gbkey=.*;"), ID = str_extract(V12, "ID=.*;"), 
         Parent = str_extract(V12, "Parent=.*;"), Genbank = str_extract(V12, "GenBank:.*;"), 
         Product = str_extract(V12, "product=.*;")) %>%
  mutate(across(gbkey:Product, function(x) sub(";.*","", x))) %>%
  mutate(across(gbkey:Product, function(x) sub(".*=","", x))) %>%
  mutate(Genbank = sub("GenBank:", "", Genbank)) %>%
  dplyr::select(-c(V3, V4, V9:V12)) %>%
  dplyr::rename(Chromosome = 1, Position = 2, Database = 3, Type = 4, Start = 5, End = 6)

#write.csv(LST, file = "LST_bedtools.csv", row.names = F)
#write.csv(Crop, file = "Crop_bedtools.csv", row.names = F)

geneID2GO <- readMappings("Data/selection/merged_GO_terms.map") #### load map file
geneID2GO ##chect it
str(head(geneID2GO)) ##check

geneNames <- names(geneID2GO) ##get gene names
head(geneNames)
ids <- read.csv("Data/selection/Temp.csv", header = F)

ids = as.vector(ids$V1) ## check ids
ids

geneList <- factor(as.integer(geneNames %in% ids)) #### use the outlier ids to create gene list of interest
names(geneList) <- geneNames
str(geneList)
geneList
table(geneList)
missing_genes <- ids[!ids %in% geneNames]
print(missing_genes)

GOdata_bio_proc <- new("topGOdata", description = "GO analysis BP", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
GOdata_bio_proc
resultFisher <- runTest(GOdata_bio_proc, algorithm = "classic", statistic = "fisher")
resultFisher
resultFisher.elim <- runTest(GOdata_bio_proc, algorithm = "elim", statistic = "fisher")


goEnrichmentBP <- GenTable(GOdata_bio_proc,
                           elimFisher = resultFisher.elim,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500)

sig_goEnrichment_BP_classic = goEnrichmentBP %>% mutate(p = as.numeric(classicFisher)) %>% filter(p < 0.05) %>%
  mutate(Term = ifelse(Term == "synaptic vesicle fusion to presynaptic a...", "synaptic vesicle fusion to\npresynaptic active zone membrane", Term)) %>%
  mutate(Term = ifelse(Term == "positive regulation of synaptic transmis...", "positive regulation of synaptic transmission", Term))

# write.csv(sig_goEnrichment_BP_classic, file = "GO_terms_LST.csv", row.names = F)
#sig_goEnrichment_BP_elim = goEnrichmentBP %>% mutate(p = as.numeric(elimFisher)) %>% filter(p < 0.05)

sig_goEnrichment_BP_classic$Term = reorder(sig_goEnrichment_BP_classic$Term, -log10(sig_goEnrichment_BP_classic$p), FUN = sum, descending = TRUE)
sig_goEnrichment_BP_classic = sig_goEnrichment_BP_classic %>% ungroup() %>% arrange(desc(-log10(p))) 


#sig_goEnrichment_CC %>% filter(grepl("presynaptic active zone cytoplasmic comp...", Term))



####
GOdata_mol_func <- new("topGOdata", description = "GO analysis BP MF", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata_mol_func, algorithm = "classic", statistic = "fisher")
resultFisher
goEnrichment_MF <- GenTable(GOdata_mol_func,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500)

sig_goEnrichment_MF = goEnrichment_MF %>% mutate(p = as.numeric(classicFisher)) %>% filter(p < 0.05) %>%
  mutate(Term = ifelse(Term == "protein tyrosine kinase collagen recepto...", "protein tyrosine kinase\ncollagen receptor activity", Term))

GOdata_cellular_comp <- new("topGOdata", description = "GO analysis BP CC", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, nodeSize = 5, gene2GO = geneID2GO)
resultFisher <- runTest(GOdata_cellular_comp, algorithm = "classic", statistic = "fisher")
resultFisher
goEnrichment_CC <- GenTable(GOdata_cellular_comp,
                           classicFisher = resultFisher,
                           orderBy = "classicFisher", topNodes = 500)

sig_goEnrichment_CC = goEnrichment_CC %>% mutate(p = as.numeric(classicFisher)) %>% filter(p < 0.05) %>%
  mutate(Term = ifelse(Term == "presynaptic active zone cytoplasmic comp...", "presynaptic active zone\ncytoplasmic component", Term))

sig_goEnrichment_BP_classic$Term = reorder(sig_goEnrichment_BP_classic$Term, -log10(sig_goEnrichment_BP_classic$p), FUN = sum, descending = TRUE)
sig_goEnrichment_BP_classic = sig_goEnrichment_BP_classic %>% ungroup() %>% arrange(desc(-log10(p))) 

sig_goEnrichment_MF$Term = reorder(sig_goEnrichment_MF$Term, -log10(sig_goEnrichment_MF$p), FUN = sum, descending = TRUE)
sig_goEnrichment_MF = sig_goEnrichment_MF %>% ungroup() %>% arrange(desc(-log10(p))) 

sig_goEnrichment_CC$Term = reorder(sig_goEnrichment_CC$Term, -log10(sig_goEnrichment_CC$p), FUN = sum, descending = TRUE)
sig_goEnrichment_CC = sig_goEnrichment_CC %>% ungroup() %>% arrange(desc(-log10(p)))


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

plot.lst.BP = ggplot( sig_goEnrichment_BP_classic[1:10,], aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#84ce7b")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 14))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "C")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot.lst.MF = ggplot(sig_goEnrichment_MF[1:10,],aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#7b84ce")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 14))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "D")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot.lst.CC = ggplot(sig_goEnrichment_CC[1:10,], aes(Term, -log10(p)))+
  geom_bar(stat = "identity", col = 'black', fill = "#ce7b84")+
  theme_classic(base_size = 16)+
  theme(panel.grid = element_blank(), axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', angle = 90, hjust = 0.98, vjust = 0.2, size = 14))+
  labs(y = expression(paste("-log"[10], " p-value")), x = NULL, title = "E")+
  scale_y_continuous(expand = c(0,0.1,0,0))

plot_grid(plot.lst.BP, plot.lst.MF, plot.lst.CC, nrow = 1, align = 'hv')
ggsave(file = "Results/GO_bar.pdf", height = 6, width = 11)

sig_goEnrichment_BP_classic$type = "Biological process"
sig_goEnrichment_MF$type = "Molecular function"
sig_goEnrichment_CC$type = "Cellular component"

sig_all = rbind(sig_goEnrichment_BP_classic[,-c(6,7)], sig_goEnrichment_MF,sig_goEnrichment_CC)

write.csv(sig_all, file = 'Results/GO_terms_all.csv', row.names = F)
## enrichment on x = -log(p), top 10? of the teTerm## enrichment on x = -log(p), top 10? of the terms on y, horizontal bar plots