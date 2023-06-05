# assessing sample distribution and variability within sites and across seasons for Terrg and rapid stable isotope data
# 2-2-2022
# Sean Kinard

# Load Packages
library(tidyverse)
library(gtable)
library(lubridate)

# Mac setup
setwd("/Users/seankinard/Dropbox/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined")
terrg <- read_csv('data/TERRG_CNS_combined.csv')
pilot <- read_csv('data/RAPID_June2019_master_isotope.csv')

# PC setup
setwd('C:\\Users\\s2kin\\Dropbox\\Research\\Manuscipts\\Diss.5_Isotope\\RAPID-TERRG-Combined')
terrg <- read_csv('data/TERRG_CNS_combined.csv')
pilot <- read_csv('data/RAPID_June2019_master_isotope.csv')

#

site_reps <- pilot %>%
  group_by(site, species, size_max) %>%
  summarize(n=n()) %>%
  filter(n>1) %>%
  as.data.frame() %>%
  mutate(temp_id = paste(site, species, size_max, sep = '_')) %>%
  select(temp_id, n)

levels(as.factor(pilot$guild))


todos <- pilot %>%
  rename(type = guild, max_forklength = size_max) %>%
  mutate(type = ifelse(type == 'Aquatic primary producer', 'Aquatic',
                       ifelse(type == 'fish', 'Fish', 
                              ifelse(type == 'invertebrate', 'Invertebrate', 
                                     ifelse(type == 'organic material', 'Detritus', 
                                            ifelse(type == 'Terrestrial primary producer', 'Terrestrial', NA)))))) %>%
  mutate(site = substr(site, 1, 2)) %>%
  mutate(my_rep = 1:length(pilot$trophic_category)) %>%
  mutate(UID = paste('RAPID_2019_06', site, type, species, max_forklength, my_rep)) %>% 
  select(UID, site, type, species, max_forklength, my_rep, carbon, nitrogen, hydrogen) %>%
  full_join(terrg)

# survey the number of samples at each site x type for each season
streams <- terrg %>%
  filter(type == 'Fish') %>%
  filter(site != 'LN' & site != 'UN') %>%
  select(site) %>%
  unique() %>%

# Creat Tables 
itemize <- function(x) {
  todos %>%
    filter(site %in% c("AR", "EM", "GC", "MR", "PD", "PL", "SF", "TR", "WM")) %>%
    filter(type == x) %>%
    mutate(event = substr(UID, 1, 13)) %>%
    group_by(event, site) %>%
    summarize(n=n()) %>%
    pivot_wider(names_from=site, values_from = n) }

item_Terrestrial <- itemize('Terrestrial')
item_Aquatic <- itemize('Aquatic')
item_Invertebrate <- itemize('Invertebrate')
item_Fish <- itemize('Fish')

# View Tables
item_Terrestrial # missing terrestrial in 2020_01 for 6 sites
item_Aquatic # missing aquatic in 2020_10 for 1 site
item_Invertebrate # invert samples < 10 at 7 sites for autumn
item_Fish # fish <10 at 1 site in 2020_01

# export tables
# View Tables
write_csv(item_Terrestrial, 'reports/terrestrial.csv')
write_csv(item_Aquatic, 'reports/aquatic.csv')
write_csv(item_Invertebrate, 'reports/Invertebrate.csv')
write_csv(item_Fish, 'reports/fish.csv')

# Is <1 sample per source a problem?

# What is the intra-site variability for a homogenized source?
my_subset<- todos %>%
  filter(site %in% c("AR", "EM", "GC", "MR", "PD", "PL", "SF", "TR", "WM")) %>%
  mutate(event = substr(UID, 1, 13))
my_subset$site <- factor( my_subset$site, levels = c('TR', 'SF', 'AR', 'PD', 'MR', 'GC', 'PL', 'WM', 'EM') )
my_subset$type <- factor( my_subset$type, levels = c('Aquatic', 'Terrestrial', 'Detritus', 'Invertebrate', 'Fish'))

filamentous <- my_subset %>%
  filter(species == 
  group_by(event,site,type) %>%
  summarize(n = n(),
         C_mean = mean(carbon),
         C_sd = sd(carbon),
         C_min = min(carbon),
         C_max = max(carbon),
         
         N_mean = mean(nitrogen),
         N_sd = sd(nitrogen),
         N_min = min(nitrogen),
         N_max = max(nitrogen) )


# Carbon variation
my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(event = substr(event, 9,13)) %>%
ggplot(aes(x=event, y = C_sd, fill = site, size = n)) +
  geom_point(shape = 23) +
  facet_wrap(facets = 'type') +
  scale_fill_viridis_d(direction = -1) +
  ylab('C13 standard deviation') + 
  ggtitle('Carbon variation within sites. 20-05 & 20-10 homogenized') +
  theme_classic(base_size = 12)

# Carbon over time table and plot
my_stats %>%
  mutate(event = paste('C_mean', substr(event, 9,13),  sep = '_')) %>%
  select(event, site, type, C_mean) %>%
  pivot_wider(names_from = event, values_from = C_mean) %>%
  mutate( may_jan = (C_mean_20_05 - C_mean_20_01),
          oct_may = (C_mean_20_10 - C_mean_20_05),
          oct_jan = (C_mean_20_10 - C_mean_20_01) ) %>%
  select(site, type, may_jan, oct_may, oct_jan) %>%
  filter(type %in% c('Aquatic', 'Terrestrial')) %>%
  arrange(type)
 
my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(year =  substr(event, 7,10),
         month = substr(event, 12, 13)) %>%
  mutate(my_date = as_date(paste(year, month, '01', sep = '-'))) %>%
ggplot(aes(x = my_date, y = C_mean, color = site)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = 'type') + 
  theme_classic(base_size = 12) + 
  scale_color_viridis_d()
  
# Nitrogen over time table and plots 
my_stats %>%
  mutate(event = paste('N_mean', substr(event, 9,13),  sep = '_')) %>%
  select(event, site, type, N_mean) %>%
  pivot_wider(names_from = event, values_from = N_mean) %>%
  mutate( may_jan = (N_mean_20_05 - N_mean_20_01),
          oct_may = (N_mean_20_10 - N_mean_20_05),
          oct_jan = (N_mean_20_10 - N_mean_20_01) ) %>%
  select(site, type, may_jan, oct_may, oct_jan) %>%
  filter(type %in% c('Aquatic', 'Terrestrial')) %>%
  arrange(type)
  
  my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(year =  substr(event, 7,10),
         month = substr(event, 12, 13)) %>%
  mutate(my_date = as_date(paste(year, month, '01', sep = '-'))) %>%
  ggplot(aes(x = my_date, y = N_mean, color = site)) +
  geom_point() +
  geom_line() +
  facet_wrap(facets = 'type')+ 
    theme_classic(base_size = 12) + 
    scale_color_viridis_d()





select(event, site, type, C_mean)

  ggplot(aes(x = event, y = C_mean, color = site)) +
  geom_point() +
  facet_wrap(facets = type)


 
my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(event = substr(event, 9,13)) %>%
  ggplot(aes(x=event, y = N_sd, fill = site, size = n)) +
  geom_point(shape = 23) +
  facet_wrap(facets = 'type') +
  scale_fill_viridis_d(direction = -1) +
  ylab('N15 standard deviation') + 
  ggtitle('Nitrogen variation within sites. 20-05 & 20-10 homogenized')

my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(event = substr(event, 9,13)) %>%
  ggplot(aes(x=event, y = C_mean, fill = site, size = n)) +
  geom_point(shape = 23) +
  facet_wrap(facets = 'type') +
  scale_fill_viridis_d(direction = -1) +
  ylab('C13 average') + 
  ggtitle('Carbon average across seasons')

my_stats %>%
  filter(type != 'Detritus') %>%
  mutate(event = substr(event, 9,13)) %>%
  ggplot(aes(x=event, y = N_mean, fill = site, size = n)) +
  geom_point(shape = 23) +
  facet_wrap(facets = 'type') +
  scale_fill_viridis_d(direction = -1) +
  ylab('N15 average') + 
  ggtitle('Nitrogen average across season')



        


# What is the intra-site variability for non-homogenized samples?



# What is the variability in sources across sample events within each site?

