Marine-Derived Nutrient Assay
# Sean Kinard
# 5-03-2022

################
# Objectives  ##
################

## 'MDN_assumption_checks_v1'
### Confirm absence of site effects on source stable isotope values prior to aggregating terrestrial and marine values.
### Confirm that mixtures are within source ranges

###################################################################
# Setup 
###################################################################

## Load Packages

library(tidyverse)

## set work directory

setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

d <- read_csv(
  "data/MDN_clean.csv")

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

# Data Visualization

## order site factor levels
d <- d %>%
  mutate(SITE = factor(SITE)) %>%
  mutate(SITE = fct_relevel(SITE,c("BB","NB","CB","HB","LB",
                                   "TR","SF", "LN","UN","AR",                                          "PD","MR","PL","GC","WM", 
                                   "EM"))) %>%
  arrange(SITE)


# carbon
d %>% 
  arrange(PPTAVG_SITE) %>%
  ggplot(aes(x = SITE, y = CARBON, fill = PPTAVG_SITE)) +
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