
### fish taxonomic information: match miso species format '<capitalized Genus 1st character><.><species lower case>'
### likely groupings: size class or taxonomic FAMILY

my_fish_sp <- read_csv('data/my_fish_species.csv')  %>%
  mutate(lowest_taxon = paste(str_sub(genus, start = 1, end = 1),
                              species,
                              sep = '.')) %>%
  select(lowest_taxon, family)
colnames(my_fish_sp) <- str_to_upper(colnames(my_fish_sp))

### load RAPID pilo stable isotope dataset: match levels for GUILD, FAMILY, SPECIES

RAPID_com <- read_csv('data/RAPID_com_data.csv') %>%
  select(trophic_category, species, trophic_level_fbase, guild) %>%
  mutate(species = str_to_upper(species)) %>%
  unique()
colnames(RAPID_com) <- str_to_upper(colnames(RAPID_com))

#### match guild levels

spp_terrg <- miso %>% 
  arrange(SPECIES) %>% 
  select(GUILD, SPECIES) %>%
  unique() 

spp_rapid <- RAPID_com %>%
  arrange(SPECIES) %>%
  select(GUILD, SPECIES) %>% 
  unique()

levels(as.factor(spp_terrg$GUILD))
levels(as.factor(spp_rapid$GUILD))

RAPID_com <- RAPID_com %>%
  mutate(GUILD = ifelse(GUILD == 'Aquatic primary producer',
                        'Aquatic',
                        ifelse(GUILD == 'fish',
                               'Fish',
                               ifelse(GUILD == 'invertebrate',
                                      'Invertebrate',
                                      ifelse(GUILD == 'organic material',
                                             'Detritus',
                                             ifelse(GUILD == 'Terrestrial primary producer',
                                                    'Terrestrial',
                                                    NA))))))

#### match species levels

##### Bio.infer TERRG benthic inverterbates

###### Obtain Taxonomic units
###### bio.infer formatting
icnt.TX <- miso %>% filter(GUILD == 'Invertebrate') %>%
  mutate(SPECIES = str_to_title(SPECIES)) %>%
  mutate(SPECIES = str_replace_all(SPECIES, 'id', 'idae')) %>%
  group_by(UID, SPECIES) %>%
  summarise(CountValue = n()) %>%
  rename(Taxon = SPECIES, SVN = UID) %>%
  ungroup() %>%
  as.data.frame()

###### match formatting provided in CRAN example data
icnt.TX$SVN <- as.factor(icnt.TX$SVN)
icnt.TX$Taxon <- as.factor(icnt.TX$Taxon)

str(icnt.TX)
str(bcnt.OR)

# i_spec <- get.taxonomic(icnt.TX)
# i_spec <- get.otu(i_spec)
# invert <- makess(i_spec)

###### get.taxonomic() returns data.frame with column headers and rows filled with <NA>



levels(as.factor(spp_terrg$SPECIES))
levels(as.factor(spp_rapid$SPECIES))

spp_rapid()



