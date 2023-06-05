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
miso %>%
  filter(SITE_TYPE == 'Stream') %>%
  group_by(SITE, GUILD, TAXON) %>%
  summarize(n=n()) %>%
  mutate(n = ifelse(n > 3, 1, 0)) %>%
  pivot_wider(names_from = SITE, values_from = n, values_fill = 0) %>%
  mutate(sum = rowSums(across(where(is.numeric)))) %>%
  filter(sum > 3) %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
# few species were caught at >3 sites : L. cyanellus, L. macrochirus, L.megalotis, cambarid, coenagrionid, hyalellid, physid

# create marine-migrant list:
mmigrants <- c("A.mitchilli", "A.spatula", "D.maculatus" , "F.grandis", "G.dormitor", "M.berylina", "M.cephalus", "macrobrachium")

d <- d %>%
  mutate(migrant = ifelse(species %in% mmigrants, 'amphidromous', 'freshwater'))

my_mmigrants <- my_coms %>%
  filter(species %in% mmigrants)

## Outlier analysis
# examine sources for outliers
# carbon
### standardize carbon for each site
d <- d %>%
  group_by(site_code, type) %>%
  summarize(c_site_mean = mean(carbon),
            c_site_std = sd(carbon)) %>%
  full_join(d) %>%
  mutate(carbon_s = ifelse(is.na(c_site_mean), 0, (carbon-c_site_mean)/c_site_std )) 
### standardize nitrogen for each site
d <- d %>%
  group_by(site_code, type) %>%
  summarize(n_site_mean = mean(nitrogen),
            n_site_std = sd(nitrogen)) %>%
  full_join(d) %>%
  mutate(nitrogen_s = ifelse(is.na(n_site_mean), 0, (nitrogen-n_site_mean)/n_site_std )) 
### standardize sulfur for each site
d <- d %>%
  group_by(site_code, type) %>%
  summarize(s_site_mean = mean(sulfur),
            s_site_std = sd(sulfur)) %>%
  full_join(d) %>%
  mutate(sulfur_s = ifelse(is.na(s_site_mean), 0, (sulfur-s_site_mean)/s_site_std ))

### carbon outliers:
d %>%
  ggplot(aes(carbon_s)) +
  geom_histogram()

carbon_outliers <- d %>% filter(carbon_s > 3 | carbon_s < -3) # 2 entries

### nitrogen outliers:
d %>%
  ggplot(aes(nitrogen)) +
  geom_histogram()

nitrogen_outliers <- d %>% filter(nitrogen_s > 3 | nitrogen_s < -3) # 2 entries

### sulfur outliers:
d %>%
  ggplot(aes(sulfur)) +
  geom_histogram()

sulfur_outliers <- d %>% filter(sulfur_s > 3 | sulfur_s < -3) # 0 entries

### entries with stable isotope values that exceed 3 standard devations from the mean for a given site and type of material:
outliers <- carbon_outliers %>% full_join(nitrogen_outliers) %>%
  select(carbon, c_site_mean, nitrogen, n_site_mean, UID)

### remove outliers
d <- d %>%
  filter(! UID %in% outliers$UID)

## reorder sites along geography from West to East:
d <- d %>%
  mutate(site_code = factor(site_code)) %>%
  mutate(site_code = fct_relevel(site_code,c("BB","NB","CB","HB","LB",
                                             "TR","SF", "LN","UN","AR",
                                             "PD","MR","PL","GC","WM",
                                             "EM"))) %>%
  arrange(site_code)


## add trophic level information from RAPID study
my_species <- d %>% ungroup %>% select(species, type) %>% unique()
my_species$species <- str_to_title(my_species$species)

RAPID_species <- RAPID_com %>%
  select(species, trophic_category, trophic_level_fbase) %>%
  unique()
RAPID_species$species <- str_to_title(RAPID_species$species)

RAPID_species %>%
  mutate(species = ifelse(species == 'Filamentous_algae', "Filamentous Algae" ,
                          ifelse(species == 'Macrophyte')))

## Data Visualization using Biplots
# Literature reference S&N MDN values

### MacAvoy, S. E., S. A. Macko, S. P. McIninch, and G. C. Garman. “Marine Nutrient Contributions to Freshwater Apex Predators.” Oecologia 122, no. 4 (March 1, 2000): 568–73. https://doi.org/10.1007/s004420050980.
## Resident freshwater C,S,N -25.6, 4.8, 12.2
## Anadromos Alosa     C,S,N -18.6, 18.4, 12.9
## Ictalurus furcatus  C,S,N -21.7, 10.6, 16.3

## Hicks, Brendan J., Mark S. Wipfli, Dirk W. Lang, and Maria E. Lang. “Marine-Derived Nitrogen and Carbon in Freshwater-Riparian Food Webs of the Copper River Delta, Southcentral Alaska.” Oecologia 144, no. 4 (August 1, 2005): 558–69. https://doi.org/10.1007/s00442-005-0035-2.
##             No spawners (C, N)        Spawners present
### Fish       -35:-29 , 4:8             -35:-20 , 6:15
### Aq inv     -40:-30 , 1:7             -34:-29,  2:9
### Aq veg     -31:-28 , 2:4             -30:-28   1:2
### Te Veg     -28:-27 , -2:-1           -28:-26   1:2

# region-wide analyses:


# overview boxplot
# carbon
d %>% ggplot(aes(x = site_code, y = carbon, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse')
# sulfur
d %>% ggplot(aes(x = site_code, y = sulfur, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse')
# nitrogen
d %>% ggplot(aes(x = site_code, y = nitrogen, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse')

# comparing amphidromous and freshwater taxa

# carbon 
d %>% ggplot(aes(x = migrant, y = carbon, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse') +
  facet_wrap(facets = d$site_code)
# sulfur 
d %>% ggplot(aes(x = migrant, y = sulfur, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse') +
  facet_wrap(facets = d$site_code)
# nitrogen
d %>% ggplot(aes(x = migrant, y = nitrogen, fill = PPTAVG_SITE)) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse') +
  facet_wrap(facets = d$site_code)

d %>%
  filter(type %in% c('Aquatic', 'Terrestrial')) %>%
  ggplot(aes(x = carbon, y = sulfur)) +
  geom_point(aes(color = site_type, fill = site_type), shape = 21, alpha = .7, color = 'black') +
  stat_ellipse(aes(color = site_type, fill = site_type), linetype=2) +
  ggtitle('Flora 34S and 13C Values @ Bays (pink) and Sampled Streams (blue)')
# Marine fauna contain greater proportions of 34S and 13C (increased values) compared to stream flora. This is consistent with trends exhibited in other marine-derived nutrient surveys.

d %>%
  ggplot(aes(x = carbon, y = sulfur)) +
  geom_point(data = filter(d, type %in% c('Fish','Invertebrate')), 
             aes(shape = type, fill = type), 
             alpha = .7, 
             color = 'black') +
  scale_shape_manual(values = 21:27) +
  stat_ellipse(data = filter(d, type %in% c('Aquatic', 'Terrestrial')), 
               aes(color = site_type, fill = site_type), 
               linetype=2) +
  geom_label(data = my_mmigrants, aes(x = carbon, y = sulfur), label = my_mmigrants$species) +
  ggtitle('Stream Fauna 34S and 13C Values @ Bays and Sampled Streams (amphidromous species are labeled)')
# All but one amphidromous samples have 34S and 13C values consistent with a mixed diet of marine and terrestrial sources. Stream fauna contain values of 34S and 13C that are within the 95% ellipses for 'sources' in this biplot.

# run mixing models to obtain % Terrestrial consumption. Then compare:
# MDN vs precip
# MDN vs species

# bay-specific analyses: 
# check that communities are within source value boundaries
# MDN vs precip
# MDN vs species

## Dam-specific analyses:
# check that communities are within source value boundaries
# community
# trophic
# species


# export data for simmr:
write_csv(d, 'data/marine_v1.csv')


