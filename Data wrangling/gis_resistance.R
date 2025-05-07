library(dplyr)
library(tidyverse)


gis = read.csv('Data/GF_LEA_combined_full2023.csv')
lea.codes = read.csv('Data/lea_codes.csv')

lea = unique(gis$NC_FESTG)
lea.gs = gis %>% filter(fclass_ful != '') %>% group_by(NC_FESTG, fclass_ful) %>% summarise(n = n())
gf = unique(gis$fclass_ful)

lea.codes = lea.codes %>% filter(LEA_Ncode %in% lea)

write.csv(lea.codes, 'Data/lea_codes2_gf.csv', row.names = F)
order(lea)


####

gis = read.csv('Data/lea_g_aes2.csv')
uniq = gis %>% filter(gf_lea_cod != '' & FOERDERART != '') %>% group_by(gf_lea_cod, FOERDERART, landuse) %>% summarise(n = n())
trans = data.frame(FOERDERART = unique(gis$FOERDERART), meaning = c(
  'annual flower field','grassland AKUM','grassland AKUM','grassland AKUM', 'perennial flower field','grassland AKUM','organic','grassland AKUM','grassland AKUM',
  'grassland AKUM','perennial flower field','?','otherAES','SNHopen','otherAES', 'SNHopen',
  'grassland AKUM','fallow','MilanAES','grassland AKUM','grassland AKUM', 'otherAES','otherAES', 'otherAES','SNHwoody', 'NA'
))

trans = data.frame(FOERDERART = unique(gis$FOERDERART), meaning = c(
  ' BS1',' AKUM',' AKUM',' AKUM', ' BS2',' AKUM',' organic',' AKUM',' AKUM',
  ' AKUM',' BS2','?','AES',' SNHopen','AES', ' SNHopen',
  ' AKUM','fallow','AES',' AKUM',' AKUM', 'AES','AES', 'AES',' SNHwoody', 'NA'
))


write.csv(uniq, "Data/aes_table.csv", row.names = F)


aa = unique(gis$NC_FESTG)
aa
Crop = data.frame(Habitat = 'crops', NC_FESTG = c(52, 53, 60, 63, 64, 112:187, 210:292, 311:393, 411:433, 510:519, 601:687, 
                               701:710, 720:799, 802:806, 851:854, 860:866, 910, 912, 914, 941, 980, 992, 999, 48, 50, 51, 98))
Grassland = data.frame(Habitat = 'grassland', NC_FESTG = c(451, 452, 453, 972))
SNHopen = data.frame(Habitat = 'SNHopen', NC_FESTG = c(454, 462:467, 492, 592, 925, 991))
SNHwoody = data.frame(Habitat = 'SNHwoody', NC_FESTG = c(70:74, 480, 586, 587, 822, 982, 983, 99))
Fallow = data.frame(Habitat = 'fallow', NC_FESTG = c(54, 55, 57, 58, 62, 78, 573, 576, 583, 591, 920, 928))
Annual_flower_field = data.frame(Habitat = 'annual flower field', NC_FESTG = c(65, 574, 590, 594, 915, 49))
Structural_flower_field = data.frame(Habitat = 'structural flower field', NC_FESTG = c(575))
Perennial_flower_field = data.frame(Habitat = 'perennial flower field', NC_FESTG = c(66, 595, 918))
Plantation = data.frame(Habitat = 'plantation', NC_FESTG = c(59, 821, 823:847 ))

habitats = rbind(Crop, Grassland, SNHopen, SNHwoody, Fallow, Annual_flower_field, Structural_flower_field, Perennial_flower_field, Plantation)

write.csv(habitats, 'Data/gis/habitats.csv', row.names = F)
write.csv(trans, 'Data/gis/trans.csv', row.names = F)

exp = gis %>% left_join(habitats, by = 'NC_FESTG') %>% select(-c(22:25)) %>% 
  left_join(trans, by = 'FOERDERART') %>% 
  mutate(habitat.aes = paste0(Habitat, meaning)) %>% mutate(habitat.aes = sub('NA ', '', habitat.aes)) %>% mutate(habitat.aes = sub('NA', '', habitat.aes))


###

gis.raw = read.csv('Data/gis/Kooperativ_1000_raw2.csv') %>% mutate(landuse3 = ifelse(landuse3 == '', gf_lea_cod, landuse3)) %>% filter(area > 0)
count = gis.raw %>% group_by(landuse3) %>% summarise(n = n(), area = sum(area))
ll = unique(gis.raw$landuse3)
ll
names = c('SNHopen', 'settlement', 'grassland', 'crop', 'SNHwoody', 'water', 'forest', 'unclear', 'grassland AKUM', 'grassland organic', 'Kooperativ flower field', 'AESother',
          'crop organic', 'annual flower field', 'other', 'perennial flower field', 'annual flower field', 'AESother', 'structural flower field', 'organic', 'perennial flower field',
          'crop', 'SNHwoody', 'settlement', 'water', 'grassland','settlement', 'grassland', 'settlement', 'SNHwoody', 'SNHwoody', 'grassland', 'settlement', 'structural flower field',
          'water', 'settlement', 'settlement', 'settlement', 'SNHwoody','settlement', 'SNHopen', 'SNHwoody', 'crop organic', 'grassland AKUM', 'grassland organic', 'annual flower field', 'SNHopen', 'SNHwoody',
          'annual flower field', 'annual flower field', 'AESother', 'AESother', 'perennial flower field','AESother', 'grassland', 'SNHopen', 'SNHopen', 'perennial flower field', 'SNHopen', 'AESother', 'NA')
jj = data.frame(landuse3 = ll, landuse4 = names)

write.csv(jj, 'Data/gis/jj.csv', row.names = F)

gis = gis.raw %>% left_join(jj, by = 'landuse3')
unique(gis$landuse4)

koop1000 = gis %>% group_by(Landscape, landuse4) %>% summarise(area = round(sum(area)/(pi*1000^2),5)) %>% mutate(landuse4 = gsub(' ', '_', landuse4)) %>%
  pivot_wider(names_from = landuse4, values_from = area, values_fill = 0) %>%
  mutate(SNH = perennial_flower_field + structural_flower_field + annual_flower_field + Kooperativ_flower_field + SNHopen + SNHwoody + AESother + grassland_AKUM,
         pSNH = perennial_flower_field + structural_flower_field/2 + SNHopen + SNHwoody + AESother + grassland_AKUM,
         aSNH = structural_flower_field + annual_flower_field/2 + Kooperativ_flower_field)

write.csv(koop1000, file = 'Data/gis/Kooperativ_landuse1000.csv', row.names = F)
