# Marine-Derived Nutrient Assay
# Sean Kinard
# 4-18-2022

################
# Objectives  ##
################

## MDN_clean_v1'
### Setup: Import and tidy data
### Visualize data distributions
### Outlier Analysis
### export master data file

###################################################################
# Setup 
###################################################################

## Load Packages

library(tidyverse)

## set work directory

setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

###################################################################
## load, tidy data, merge data
###################################################################

### stable isotope data: filter TERRG winter samples (carbon, nitrogen, and sulfur)

miso <- read_csv('data/TERRG_CNS_combined.csv') %>%
  filter(collection == '2020-01-15') %>%
  select(-collection, -project) %>%
  mutate(species = str_to_upper(species)) %>%
  rename(guild = type)
  
colnames(miso) <- str_to_upper(colnames(miso))

#### fix fish species labels

miso <- miso %>% 
  filter(GUILD == 'Fish') %>%
  separate(SPECIES,
                  into = c('GENUS', 'SPECIES'),
                  sep = '[.]') %>%
  full_join(filter(miso, GUILD != 'Fish')) %>%
  mutate(SPECIES = str_to_lower(SPECIES))

(fish_species <- miso %>% 
  filter(GUILD == 'Fish') %>%
  pull(SPECIES) %>%
  unique())

#### Assign inverterbate family names
  
miso <- miso %>%
  filter(GUILD == 'Invertebrate') %>%
  mutate(SPECIES = str_to_title(SPECIES)) %>%
  mutate(SPECIES = str_replace_all(SPECIES, 'id', 'idae'))%>%
  mutate(FAMILY = ifelse(grepl('idae', SPECIES),
                         SPECIES,
                         NA)) %>%
  mutate(FAMILY = ifelse(SPECIES == 'Rharrisii',
                         'Panopeidae',
                  ifelse(SPECIES == 'Macrobrachium',
                         'Palaemonidae',
                  ifelse(SPECIES == 'Corbicula',
                         'Cyrenidae',
                  ifelse(SPECIES == 'Melanoidaees',
                         'Thiaridae',
                  ifelse(SPECIES == 'Csapidaeus',
                         'Portunidae',
                         FAMILY)))))) %>%
  filter(SPECIES != 'Unknown') %>%
  select(-SPECIES) %>%
  full_join(filter(miso, GUILD != 'Invertebrate'))

(invert_families <- miso %>%
    filter(GUILD == 'Invertebrate') %>%
    pull(FAMILY) %>%
    unique())

#### Consolidate Source identifiers

miso <- miso %>%
  filter(GUILD != 'Fish') %>%
  filter(GUILD != 'Invertebrate') %>%
  mutate(SOURCE_TYPE = SPECIES) %>%
  select(-SPECIES) %>%
  full_join(filter(miso, GUILD %in% c('Fish', 'Invertebrate')))

#### create universal lowest taxon identifier
miso <- miso %>%
  mutate(TAXON = ifelse(GUILD == 'Fish',
                        paste(GENUS, SPECIES, sep = '.'),
                 ifelse(GUILD == 'Invertebrate',
                        FAMILY,
                 SOURCE_TYPE)))

### Assign latitude, longitude, and nearest bay for sample sites

miso <- read_csv('data/MDN_site_lat_lon.csv') %>%
  rename(SITE = site_code,
         STANAME = site_name,
         LATITUDE = site_lat,
         LONGITUDE = site_lon,
         NEAREST_BAY = nearest_bay) %>%
  right_join(miso)

### USGS gauges ii: see gagesII_sept30_2011_var_desc.xlsx for meta data

miso <- read_csv('data/MDN_conterm.csv') %>%
  select(-STANAME) %>%
  right_join(miso)

### import streamer data, calculate distance to nearest bay (km) and clean
miso <- read_csv('data/streamer_data.csv') %>%
  select(site_code, position_between, over_under, trace_river_mile) %>%
  pivot_wider(names_from = 'over_under', 
              values_from = 'trace_river_mile') %>%
  mutate(baydist_km = 1.609344*(under + position_between*(over-under))) %>%
  select(site_code, baydist_km) %>%
  rename(SITE = site_code, BAYDIST_KM = baydist_km) %>%
  right_join(miso)

###################################################################
# MERGE with Functional Traits
###################################################################

### fish traits

#### Identified several taxa miscategorized as carnivores
not_PISC <- c("A.mitchilli", "C.variegatus", "F.grandis", "H.cyanoguttatus", "L.miniatus", "O.aureus")

miso <- read_csv('data/fish_trait matrix_22.csv') %>%
  separate(col = genus.species,
           into = c('genus', 'species'),
           sep = "[.]" ) %>%
  mutate(TAXON = paste(substr(genus, 1, 1), species, sep = '.')) %>%
  mutate(TAXON = str_to_title(TAXON)) %>%
  select(TAXON, genus, herbivore, planktivore, detrivore, invertivore, carnivore) %>%
  rename(GENUS_2 = genus,
         HERB = herbivore,
         PLAN = planktivore,
         DETR = detrivore,
         INVE = invertivore,
         PISC = carnivore
         ) %>%
  mutate(PISC = ifelse(TAXON %in% not_PISC, 0, PISC)) %>%
  right_join(miso)

### invertebrate trophic traits: PR, CG, HB, CF, SH, PA
miso <- read_csv('data/FreshwaterBioTraits_Transposed_20100927.csv') %>%
  filter(Family %in% invert_families) %>%
  select(Family, Feed_prim_abbrev) %>%
  unique() %>%
  filter(!is.na(Feed_prim_abbrev)) %>%
  mutate(ct = 1) %>%
  pivot_wider(names_from = Feed_prim_abbrev,
              values_from = ct,
              values_fill = 0) %>%
  rename(FAMILY = Family) %>%
  right_join(miso)

### Assign Consumer levels

#### Remove NAs within feeding group variables
miso <- miso %>%
  mutate(HERB = ifelse(is.na(HERB), 0, HERB),
         PR = ifelse(is.na(PR), 0, PR),
         CG = ifelse(is.na(CG), 0, CG),
         HB = ifelse(is.na(HB), 0, HB),
         CF = ifelse(is.na(CF), 0, CF),
         SH = ifelse(is.na(SH), 0, SH),
         PA = ifelse(is.na(PA), 0, PA),
         PLAN = ifelse(is.na(PLAN), 0, PLAN),
         DETR = ifelse(is.na(DETR), 0 , DETR),
         INVE = ifelse(is.na(INVE), 0 , INVE),
         PISC = ifelse(is.na(PISC), 0, PISC))

#### Create TROPHIC_LEVEL variable

##### Frequency distribution shows 3 size groups (<80, 80-119, =>120 mm)
miso %>%
  filter(GUILD == 'Fish') %>%
  ggplot(aes(MAX_FORKLENGTH)) +
  geom_freqpoly()
##### Assigning predator to known PISC taxa with Forklengths > 100mm

miso <- miso %>%
  mutate(TROPHIC_LEVEL = 
           ifelse(GUILD %in% c("Aquatic", 
                               "Detritus", 
                               "Terrestrial"), 0, 
                  ifelse(PISC > 0 & MAX_FORKLENGTH > 100, 3,
                  ifelse(PR > 0 | 
                  INVE > 0 | 
                  PLAN > 0, 2, 1))))

# Assign SITE_TYPE as Bay or Stream
miso <- miso %>% 
  mutate(SITE_TYPE = ifelse(SITE %in% c('BB', 'CB', 'HB', 'LB', 'NB'), 'Bay', 'Stream'))

# sample distribution: Guild
miso %>%
  filter(SITE_TYPE == 'Stream') %>%
  group_by(SITE, GUILD) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = SITE, values_from = n, values_fill = 0) %>%
  as.data.frame()

# Sample distribution: Taxon
miso %>%
  filter(SITE_TYPE == 'Stream') %>%
  group_by(SITE, TAXON) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = SITE, values_from = n, values_fill = 0) %>%
  as.data.frame()

# Find widespread species
(widespread <- miso %>%
  filter(SITE_TYPE == 'Stream') %>%
  group_by(SITE, GUILD, TAXON) %>%
  summarize(n=n()) %>%
  mutate(n = ifelse(n > 3, 1, 0)) %>%
  pivot_wider(names_from = SITE, values_from = n, values_fill = 0) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 3) %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  pull(TAXON) )

# Assign migrants as Diadromous and residents as Potamodromous:
mmigrants <- c("A.mitchilli", "A.spatula", "D.maculatus" , "F.grandis", "G.dormitor", "M.berylina", "M.cephalus", "Palaemonidae")

miso <- miso %>%
  mutate(MIGRANT = ifelse(TAXON %in% mmigrants, 'Diadromous',
                   ifelse(GUILD %in% c('Aquatic', 'Terrestrial', 'Detritus'), 'Local Source', 'Potamodromous')))

## Outlier analysis
# examine sources for outliers
# carbon
### standardize carbon for each site
miso <- miso %>%
  group_by(SITE, GUILD) %>%
  summarize(c_site_mean = mean(CARBON),
            c_site_std = sd(CARBON)) %>%
  full_join(miso) %>%
  mutate(carbon_s = ifelse(is.na(c_site_mean), 0, (CARBON-c_site_mean)/c_site_std )) 
### standardize nitrogen for each site
miso <- miso %>%
  group_by(SITE, GUILD) %>%
  summarize(n_site_mean = mean(NITROGEN),
            n_site_std = sd(NITROGEN)) %>%
  full_join(miso) %>%
  mutate(nitrogen_s = ifelse(is.na(n_site_mean), 0, (NITROGEN-n_site_mean)/n_site_std )) 
### standardize sulfur for each site
miso <- miso %>%
  group_by(SITE, GUILD) %>%
  summarize(s_site_mean = mean(SULFUR),
            s_site_std = sd(SULFUR)) %>%
  full_join(miso) %>%
  mutate(sulfur_s = ifelse(is.na(s_site_mean), 0, (SULFUR-s_site_mean)/s_site_std ))

### carbon outliers:
miso %>%
  ggplot(aes(carbon_s)) +
  geom_freqpoly()

(carbon_outliers <- miso %>% filter(carbon_s > 3 | carbon_s < -3)) # 2 entries

### nitrogen outliers:
miso %>%
  ggplot(aes(NITROGEN)) +
  geom_freqpoly()

(nitrogen_outliers <- miso %>% filter(nitrogen_s > 3 | nitrogen_s < -3)) # 2 entries

### sulfur outliers:
miso %>%
  ggplot(aes(SULFUR)) +
  geom_freqpoly()

(sulfur_outliers <- miso %>% filter(sulfur_s > 3 | sulfur_s < -3)) # 0 entries

### entries with stable isotope values that exceed 3 standard devations from the mean for a given site and type of material:
(outliers <- carbon_outliers %>% full_join(nitrogen_outliers) %>%
  select(CARBON, c_site_mean, NITROGEN, n_site_mean, UID) )

### remove outliers
miso <- miso %>%
  filter(! UID %in% outliers$UID)

## reorder sites along geography from West to East:
miso <- miso %>%
  mutate(SITE = factor(SITE)) %>%
  mutate(SITE = fct_relevel(SITE,c("BB","NB","CB","HB","LB",
                                   "TR","SF", "LN","UN","AR",                                          "PD","MR","PL","GC","WM", 
                                   "EM"))) %>%
  arrange(SITE)

# Export Master dataframe:
write_csv2(miso, "data/marine_master.csv")

# End


