# combine TERRG stable isotope data
# Sean Kinard
# 1-25-2022

library(tidyverse)
library(bio.infer)

# set work directory

# mac
# setwd("/Users/seankinard/Dropbox/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined")
# twin <- read_csv('data/TERRG_Jan2020_CNS_1-21.csv')
# tsum <- read_csv('data/TERGG_May2020_CN_37-52.csv')
# taut <- read_csv('data/TERRG_Oct2020_CN_22-36.csv')

# pc
setwd('C:\\Users\\s2kin\\Dropbox\\Research\\Manuscipts\\Diss.5_Isotope\\RAPID-TERRG-Combined')
twin <- read_csv('data\\TERRG_Jan2020_CNS_1-21.csv')
tsum <- read_csv('data\\TERGG_May2020_CN_37-52.csv')
taut <- read_csv('data\\TERRG_Sep2020_CN_22-36.csv')
tsum_key <- read_csv('data\\TERRG-May2020-CN_key.csv')

# Visual inspections
# twin 3 letter codes for fish and no species for sources
# tsum:  no max_length data, no notes included
filter(taut, ! is.na(notes)) %>%
  select(notes) %>% as.data.frame() # taut: notes indicate strange notes about rerun samples, fork lengths, vegetative species
 
#The goal is to create one dataset with UID, project, year, month, SampleID, site, type, species, replicate, forklength, C, N, H, S, C_percent, N_percent, H_percent, S_percent, notes

twin # UID = site, type (Detritus, Aquatic, Terrestrial, Fish, Invertebrate), species (fish 3-letter abbreviation, ALL CAPS), size, my_rep
 
tsum # UID = site, type (2 letters mix case), species (5 lowercase), replicate
# make my_rep by extract the last 2 string of UID, then replace . with ''
tsum$my_rep <- substr(tsum$UID, nchar(tsum$UID) - 2 + 1, nchar(tsum$UID)) %>%
  str_replace('-', '') %>%
  as.numeric()
  
# Explore october
colnames(taut)
taut$UID
taut[426:444,]
# UID = not standardized
# forklengths are not binned
# No replicate column, but the last .## of UID contains replicate # per site-species
taut$count
taut %>%
filter(!is.na(count)) %>%
  select(UID, count,  carbon, nitrogen, notes) %>%
  as.data.frame()
# Count probably refers to the number of crushed individuals in the sample
# Investigate notes
taut %>%
filter(!is.na(notes)) %>%
  select(UID, notes) %>%
  as.data.frame()
# numbers for blue crab and Cambarid 'Faxonius' taxons
# 'Rerun TERRG-21-36' refers to the fact that tray 23 was rerun as tray 36. In the former compile script, tray 23 was omitted and tray 36 was included.
# make my_rep by extract the last 2 string of UID, then replace . with ''
taut$my_rep <- substr(taut$UID, nchar(taut$UID) - 2 + 1, nchar(taut$UID)) %>%
  str_replace('.', '') %>%
  as.numeric()

# standardize column headers winter
twin <- rename(twin, S_percent = 'S%', N_percent = 'wt %N', C_percent = 'wt %C', my_rep = id, max_forklength = max_length)

# standardize column headers summer
tsum <- tsum %>%
  rename(carbon = "d13C (permil, vs VPDB)", nitrogen = "d15N (permil, vs AIR)", C_percent = CP, N_percent = NP, max_forklength = max_length) %>%
  select(UID, carbon, nitrogen, C_percent, N_percent, tray, site, type, species, max_forklength, details, my_rep)

# standardize column headers autumn
taut <- taut %>%
  rename(C_percent = carbon_percent, N_percent = nitrogen_percent, max_forklength = length) %>%
  select(UID, tray, site, type, species, max_forklength, carbon, nitrogen, C_percent, N_percent, my_rep)
taut$max_forklength <- as.numeric(taut$max_forklength)

# Create new columns: project, year, month, UID_new
twin <- twin %>%
  rename(old_id = UID) %>%
  mutate(new_id = paste('TERRG_2020_01', 
                        site, 
                        type, 
                        species, 
                        max_forklength, 
                        my_rep, 
                        sep = '_'))

tsum <- tsum %>%
  rename(old_id = UID) %>%
  mutate(new_id = paste('TERRG_2020_05', 
                        site, 
                        type, 
                        species, 
                        max_forklength, 
                        my_rep, 
                        sep = '_'))
# Fix species abbreviations by accessing weight key
tsum <- tsum_key %>%
  filter(Type != 'standard') %>%
  filter(Sample_ID != "USGS40") %>%
  select(Sample_ID, G.species, key_notes) %>%
  distinct() %>%
  filter(!duplicated(Sample_ID)) %>%
  right_join(tsum, by = c('Sample_ID' = 'old_id')) %>%
  select(-species, -details) %>%
  rename(old_id = Sample_ID, 
         species = G.species, 
         details = key_notes)

taut <- taut %>%
  rename(old_id = UID) %>%
  mutate(new_id = paste('TERRG_2020_09', 
                        site, 
                        type, 
                        species, 
                        max_forklength, 
                        my_rep, 
                        sep = '_')) 

# merge 
mster <- twin %>%
  full_join(tsum) %>%
  full_join(taut)

# Extract project, Year and month
mster %>%
  mutate(project = substr(mster$new_id, 1,5),
         year = substr(mster$new_id, 7,10),
         month = substr(mster$new_id, 12,13) )

# standardize row entries for the various factors (such as species)
# Use ifelse() and str_replace()
# Site codes
levels(as.factor(mster$site))
mster$site <- substr(mster$site, 1,2) # convert site to 2 letter code

# Types
levels(as.factor(mster$type))

mster <- mster%>%
  mutate(type = ifelse(type == 'Aquatic' | 
           type == 'Aquatic Source', 
         'Aquatic',
         ifelse(type == 'invertebrate' | 
                type == 'Invertebrate' | 
                type == 'Inverterbrate', 
                'Invertebrate',
         ifelse(type == 'Terrestrial Source' | 
                type =='Terrestrial', 
                'Terrestrial',
         ifelse(type == 'fish' | 
                type == 'Fish', 
                'Fish',
         ifelse(type == 'Detritus', 
                'Detritus', 
                NA))))))

# Species
levels(as.factor(mster$species)) # 200 levels is unwieldy, and formatting is type-dependent

# break mster into groups based on type, standardize levels, then merge back into mster
# Fish
my_fish <- filter(mster, type == 'Fish')
levels(as.factor(my_fish$species))

my_fish$species <- str_to_lower(my_fish$species)
my_fish$species <- str_to_title(my_fish$species)
my_fish$species <- str_replace(my_fish$species, ' ', '')
my_fish <- my_fish %>%
  mutate(species = ifelse(species %in% c('A.Mexicanus', 'Ame'), 'A.mexicanus',
                              ifelse(species %in% c('A.Natalis'), 'A.natalis',
                              ifelse(species == 'A.Rostrata', 'A.rostrata',
                              ifelse(species == 'A.Sayanus', 'A.sayanus',
                              ifelse(species == 'Ami', 'A.mitchilli',
                              ifelse(species == 'Asp', 'A.spatula',
                              ifelse(species %in% c('C.Lutrensis', 'Clu'), 'C.lutrensis',
                              ifelse(species %in% c('C.Variegatus', 'Cva'), 'C.variegatus',
                              ifelse(species %in% c('C.Venusta'), 'C.venusta',
                              ifelse(species == 'D.Cepedianum', 'D.cepedianum',
                              ifelse(species %in% c('D.Maculatus', 'Dma'), 'D.maculatus',
                              ifelse(species == 'Egr', 'E.gracile',
                              ifelse(species %in% c('F.Chrysotus'),'F.chrysotus',
                              ifelse(species %in% c('F.Grandis', 'Fgr'), 'F.grandis',
                              ifelse(species %in% c('F.Notatus', 'Fno'), 'F.notatus',
                              ifelse(species %in% c('G.Affinis', 'Gaf'), 'G.affinis',
                              ifelse(species %in% c('H.Cyanoguttatus', 'Hcy'), 'H.cyanoguttatus',
                              ifelse(species == 'HybridSunfish', 'L.hybrid',
                              ifelse(species %in% c('I.Punctatus', 'Ipu'), 'I.punctatus',
                              ifelse(species %in% c('L.Auritus', 'Lau'), 'L.auritus',
                              ifelse(species %in% c('L.Cyanellus', 'Lcy'), 'L.cyanellus',
                              ifelse(species %in% c('L.Gulosus', 'Lgu', 'L.gullosus'), 'L.gulosus',
                              ifelse(species %in% c('L.Macrochirus', 'Lma', 'L.macrochrius'), 'L.macrochirus',
                              ifelse(species %in% c('L.megalotus', 'L.Megalotis', 'Lme'), 'L.megalotis',
                              ifelse(species %in% c('L.Microlophus'), 'L.microlophus',
                              ifelse(species %in% c('L.Miniatus'), 'L.miniatus',
                              ifelse(species %in% c('L.Oculatus'), 'L.oculatus',
                              ifelse(species == 'Lhu', 'L.humilis',
                              ifelse(species %in% c('Mce'), 'M.cephalus',
                              ifelse(species %in% c('M.Salmoides', 'Msa'), 'M.salmoides',
                              ifelse(species == 'Mbe', 'M.berylina',
                              ifelse(species %in% c('N.Texanus'), 'N.texanus',
                              ifelse(species == 'Ngy', 'N.gyrinus', 
                              ifelse(species %in% c('O.Aureus', 'Oau'), 'O.aureus',
                              ifelse(species == 'P.cla', 'P.clarkii',
                              ifelse(species %in% c('P.Formosa'), 'P.formosa',
                              ifelse(species %in% c('P.Latipinna', 'Pla'), 'P.latipinna',
                              ifelse(species == 'P.Vigilax', 'P.vigilax', 
                              ifelse(species %in% c('Gdo', 'G.Dormitor'), 'G.dormitor', species ))))))))))))))))))))))))))))))))))))))))

# crayfish were included with fish. Switch to inverterbate data
mster <- mster %>%
  mutate(type = ifelse(species == 'P.clarkii', 'Invertebrate', type)) 

my_inv <- filter(mster, type == 'Invertebrate')
levels(as.factor(my_inv$species))

my_inv$species <- str_to_lower(my_inv$species)
my_inv$species <- str_replace(my_inv$species, 'dae', 'd')
my_inv$species <- str_remove(my_inv$species, 'sp2')
my_inv <- my_inv %>%
  mutate(species = ifelse(species %in% c('baetid', 'baetis'), 'baetid',
                          ifelse(species %in% c('cambarid', 'crambid', 'P.clarkii'), 'cambarid',
                          ifelse(species %in% c('coenagrion', 'coenagrionid'), 'coenagrionid',
                          ifelse(species %in% c('dysticid', 'dytiscid'), 'dysticid',
                          ifelse(species %in% c('elmid', 'elmidbeetle', 'elmidlarva'), 'elmid',
                          ifelse(species %in% c('halipid', 'haliplid', 'haplipid'), 'haliplid',
                          ifelse(species %in% c('hyalella', 'hyalellid', 'hyallelid'), 'hyalellid',
                          ifelse(species %in% c('hydrophilid', 'hydrophillid'), 'hydrophilid',
                          ifelse(species %in% c('hydrophsychid', 'hydropsychid'), 'Hydropsychid',
                          ifelse(species %in% c('invertef'), 'unknown',
                          ifelse(species %in% c('libelimid', 'libellid', 'libellulid'), 'Libellulid',
                          ifelse(species %in% c('physella', 'physid'), 'physid',
                          ifelse(species %in% c('veliid', 'vellid'), 'velliid', species))))))))))))))
                          
                          


my_src <- filter(mster, type %in% c('Aquatic', 
                                    'Terrestrial', 
                                    'Detritus'))
levels(as.factor(my_src$species))
my_src %>%
  filter(species %in% c('Aquat', 'AV', 'CPOM', 'Diato', 'Diatoms', 'FA', 'Filam', 'Filamentous Algae', 'FPOM', 'grass', 'leave', 'PERIPHYTON', 'TV')) %>%
  filter(!is.na(details)) %>%
  select(details) %>%
  as.data.frame()
filter(my_src, species == 'tv')  %>% as.data.frame()

my_src$species <- str_to_lower(my_src$species)
my_src <- my_src %>%
  mutate(species = ifelse(species == 'av', 'macrophyte',
                   ifelse(species %in% c('diatoms', 'periphyton'), 'periphyton',
                   ifelse(species %in% c('filamentous_algae', 'fa', 'filamentous algae'), 'filamentous_algae',
                   ifelse(species %in% c('leaves', 'shrub'), 'tree_leaf',
                   ifelse(species %in% c('tv'), 'mix_tv',
                   ifelse(species == 'fpom', 'fine_detritus',
                   ifelse(species == 'cpom', 'coarse_detritus', species ))))))))
                   
mster <- full_join(my_fish, my_inv) %>%
  full_join(my_src)

# replicates
mster %>%
  filter(is.na(my_rep)) %>%
  select(new_id, my_rep, details) %>%
  group_by(new_id) %>%
  summarize(n = n()) %>%
  filter(n > 1)

mster <- mster %>%
  mutate(my_rep = ifelse(is.na(my_rep), 1, my_rep))

which(proto$new_id == 'TERRG_2020_10_AR_invertebrate_Chironomidae_NA_NA')
# 1455 1456
which(proto$new_id == 'TERRG_2020_10_AR_Terrestrial Source_Salix spp._NA_NA')
# 1815 1822
which(proto$new_id == 'TERRG_2020_10_GC_Terrestrial Source_Fraxinus berlandieriana_NA_NA')
# 1803 1807
which(proto$new_id == 'TERRG_2020_10_MR_invertebrate_Chironomidae_NA_NA')
# 1440 1441
which(proto$new_id == 'TERRG_2020_10_SFC_invertebrate_Chironomidae_NA_NA')
# 1462 1463
mster$my_rep[c(1456, 1822, 1807, 1463, 1441)] <- 2

# Create new UID
mster <- mster %>%
  mutate(UID = paste(substr(new_id, 1, 13), site, type, species, max_forklength, my_rep, sep='_')) %>%
  select(-new_id, -old_id)

# replace NA with 'unkown' in species column
mster <- mster %>%
  mutate(species = ifelse(is.na(species), 'unknown', species)) %>%
  mutate(project = substr(UID, 1, 5),
         collection = as_date(paste(substr(UID, 7,10), substr(UID, 12,13), '15', sep = '-' )))

# export merged and finalized dataset
write_csv(mster, 'data/TERRG_CNS_combined.csv')




