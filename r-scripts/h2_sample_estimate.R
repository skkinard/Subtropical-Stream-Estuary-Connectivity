# Hydrogen sampling outlook
# Sean Kinard
# 2-16-2022

# load packages
library(tidyverse)
library('sjmisc')
library(lubridate)
library(ggpubr)
library(rstatix)

# set work directory
setwd('C:\\Users\\s2kin\\Dropbox\\Research\\Manuscipts\\Diss.5_Isotope\\RAPID-TERRG-Combined\\')

# load data
d <- read_csv('data/TERRG_CNS_combined.csv') %>% 
  filter(site %in% c("AR","EM", "GC", "MR", "PD", "PL", "SF", "TR", "WM") )

d$site <- factor(d$site, levels = c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

# sample tables

all_samples <- d %>%
  group_by(site, collection, type) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples) %>%
    as.data.frame()

h2f_maysep <- d %>%
  filter(species %in% c('G.affinis', 'P.latipinna', 'H.cyanoguttatus', 'L.macrochirus', 'L.megalotis', 'L.gulosus', 'L.cyanellus', 'M.salmoides', 'L.oculatus') |
           type %in% c('Terrestrial', 'Aquatic')) %>%
  filter(collection != '2020-01-15') %>%
  group_by(collection, species, site) %>%
  summarize(samples = n() ) %>%
  mutate(samples = ifelse(samples > 3, 3, samples)) %>%
  ungroup() %>%
  group_by(species, site) %>%
  summarize(samples = sum(samples) ) %>%
  pivot_wider(names_from = site, values_from = samples) %>%
  filter(species != 'clover')

h2f_maysep[is.na(h2f_maysep)] <- 0

h2f_maysep_gtotal <- h2f_maysep %>% 
  rowwise() %>% 
  mutate(myTidySum = sum(c_across(PD:EM)),
         project = 'hydrogen') %>%
  group_by(project) %>%
  summarize(grandtotal = sum(myTidySum))

# The high-end estimate for the number of samples for hydrogen testing would be 347 plus any invertebrate powders that would be recovered and deemed suitable. This grand total includes a maximum of 3 replicates within each of nine species of fish and sources. Species were not distinguished by size. Patterns in mixing models from the 3-site pilot suggest that the basal resource patterns across the precipitation gradient are most evident at higher trophic levels; fish species from P.latipinna to L.gulosus had progressive trends with non-overlapping credible intervals.

# Aggregating seasons: The estimated 347 hydrogen samples includes powders from both summer and autumn sampling periods. A benefit for running both seasons is that not all fish were caught evenly through each sampling period, so an aggregated approach for hydrogen carbon mixing models could provide 9 estimates of % aquatic resources for a given species which could produce regressions. However, aggregating 2 seasons assumes temporal stability in resource isotopes. Given that hydrogen isotopes were more stable in sources spatially, and carbon does not consistently vary across seasons, I think aggregating hydrogen signals across the seasons would provide the most powerful dataset for terrestrial vs aquatic mixing models.

# Invertebrates: If we seek to connect the trophic interactions from upper to intermediate levels, we could run inverterbates. Invertebrate orders that are likely to contain leftover powder would include gastropoda, bivalva, odonata, decapoda, and hemiptera. Also, Hydrogen tins are filled with a much smaller amount of material, so recovering powders although tedious is feasible for more samples than Chris Groff and myself might initially think possible.

# size class: To assess whether basal resource dependency shifts with life-stage of common lepomis taxa, I would recommend running an additional (n = 9 sites, 2 additional Lepomis groups, 3 samples, 2 seasons) 112 samples. Thus, the top-end estimate This would make a total of 559 samples + inverts for hydrogen analysis.

## Available May samples
# missing periphyton AND filamentous algae at WM
# missing gambusia at TR AR MR
# missing gambusia AND poecilids at AR

# Running hydrogen on only May samples would involve about 177 samples + inverts (n~80) + sizeclass fish (n~100)

mafish_available <- d %>%
  filter(species %in% c('G.affinis', 'P.latipinna', 'H.cyanoguttatus', 'L.macrochirus', 'L.megalotis', 'L.gulosus', 'L.cyanellus', 'M.salmoides', 'L.oculatus') ) %>%
  filter(species != 'clover') %>%
  filter(collection == '2020-05-15') %>%
  group_by(species, site) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = site, values_from = samples) %>%
  as.data.frame()
mafish_available[is.na(mafish_available)] <- 0

masrc_available <- d %>%
  filter(type %in% c('Aquatic', 'Terrestrial') ) %>%
  filter(species != 'clover') %>%
  filter(collection == '2020-05-15') %>%
  group_by(type, species, site) %>%
  arrange(type) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = site, values_from = samples) %>%
  as_tibble() %>%
  relocate(TR, .before = SF)
masrc_available[is.na(masrc_available)] <- 0

mafish_total <- mafish_available %>%
  pivot_longer(cols = c(SF:MR), names_to = 'site', values_to = 'n') %>%
  mutate(n_cap = ifelse(n >= 3, 3, n),
         project = 'fish') %>%
  select(-n) %>%
  group_by(project) %>%
  summarize(n_tot = sum(n_cap))

masrc_total <- masrc_available %>%
  pivot_longer(cols = c(TR:EM), names_to = 'site', values_to = 'n') %>%
  mutate(project = 'source') %>%
  group_by(project) %>%
  summarise(n_tot = sum(n))

mafish_available
mafish_total
masrc_available  
masrc_total

## May small-mouthed sunfish (megalotis and macrochirus) size class sample distribution
# missing <50mm at TR MR GC
# missing <100mm at TR
# missing <200mm at SF

# recommend analyzing hydrogen in 100mm (n = 24) and 200mm (25)
May_lepomis_sm <- d %>%
  filter(species %in% c('L.macrochirus', 'L.megalotis') ) %>%
  filter(collection == '2020-05-15') %>%
  mutate(maxfl_bin = ifelse(max_forklength <= 50, '50mm',
                     ifelse(max_forklength > 50 & 
                            max_forklength <= 100 , '100mm',
                     ifelse(max_forklength > 100 &
                            max_forklength <= 200, '200mm',
                     ifelse(max_forklength > 200 &
                            max_forklength <= 300, '300mm',
                     ifelse(max_forklength > 300 &
                            max_forklength <= 500, '500mm', max_forklength)))))) %>%
  group_by(site, maxfl_bin) %>%
  summarize(available = n()) %>%
  pivot_wider(names_from = maxfl_bin, values_from = available) %>%
  relocate('50mm', .before='200mm') %>%
  relocate('100mm', .before='200mm')
May_lepomis_sm[is.na(May_lepomis_sm)] <- 0
May_lepomis_sm 

## May large-mouthed sunfish (megalotis and macrochirus) size class sample distribution
# missing <50mm at TR, SF, MR, PL, GC
# missing <100mm at TR
# missing <200mm at TR, PL, WM

# Recommend analyzing individuals in the <200 (n = 26) and <300 (n = 22)size class range
May_lepomis_lm <- d %>%
  filter(species %in% c('L.gulosus', 'L.cyanellus') ) %>%
  mutate(maxfl_bin = ifelse(max_forklength <= 50, '50mm',
                     ifelse(max_forklength > 50 & 
                            max_forklength <= 100 , '100mm',
                     ifelse(max_forklength > 100 &
                            max_forklength <= 200, '200mm',
                     ifelse(max_forklength > 200 &
                            max_forklength <= 300, '300mm',
                     ifelse(max_forklength > 300 &
                            max_forklength <= 500, '500mm', max_forklength)))))) %>%
  group_by(site, maxfl_bin) %>%
  summarize(available = n()) %>%
  pivot_wider(names_from = maxfl_bin, values_from = available) %>%
  relocate('50mm', .before='100mm')
May_lepomis_lm[is.na(May_lepomis_lm)] <- 0
May_lepomis_lm

## Available September samples
# September sources: missing periphyton AND filamentous algae at PL
# September gambusia missing at TR GC
# Gambusia AND Poecilid missing at GC

# Running only samples from September would entail 128 fish samples + 43 source samples for a total of 173 samples + inverts (3 spp -> n~80) + additional size class samples (n~100)

sefish_available <- d %>%
  filter(species %in% c('G.affinis', 'P.latipinna', 'H.cyanoguttatus', 'L.macrochirus', 'L.megalotis', 'L.gulosus', 'L.cyanellus', 'M.salmoides', 'L.oculatus') ) %>%
  filter(species != 'clover') %>%
  filter(collection == '2020-09-15') %>%
  group_by(species, site) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = site, values_from = samples) %>%
  as.data.frame()
sefish_available[is.na(sefish_available)] <- 0

sesrc_available <- d %>%
  filter(type %in% c('Aquatic', 'Terrestrial') ) %>%
  filter(species != 'clover') %>%
  filter(collection == '2020-09-15') %>%
  group_by(type, species, site) %>%
  arrange(type) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = site, values_from = samples) %>%
  as_tibble() %>%
  relocate(TR, .before = SF)
sesrc_available[is.na(sesrc_available)] <- 0

sefish_total <- sefish_available %>%
  pivot_longer(cols = c(SF:GC), names_to = 'site', values_to = 'n') %>%
  mutate(n_cap = ifelse(n >= 3, 3, n),
         project = 'fish') %>%
  select(-n) %>%
  group_by(project) %>%
  summarize(n_tot = sum(n_cap))

sesrc_total <- sesrc_available %>%
  pivot_longer(cols = c(PD:PL), names_to = 'site', values_to = 'n') %>%
  mutate(project = 'source') %>%
  group_by(project) %>%
  summarise(n_tot = sum(n))

sefish_available
sefish_total
sesrc_available  
sesrc_total


## September small-mouthed sunfish (megalotis and macrochirus) size class sample distribution
#  <50mm missing n at MR GC, low n at PL WM EM
#  <100mm low n at MR
#  <200mm low n at AR PD PL WM EM

# recommend analyzing hydrogen in 100mm (n = 19) and 200mm (13) in small-mouthed Lepomis September
September_lepomis_sm <- d %>%
  filter(species %in% c('L.macrochirus', 'L.megalotis') ) %>%
  filter(collection == '2020-09-15') %>%
  mutate(maxfl_bin = ifelse(max_forklength <= 50, '50mm',
                     ifelse(max_forklength > 50 & 
                            max_forklength <= 100 , '100mm',
                      ifelse(max_forklength > 100 &
                             max_forklength <= 200, '200mm',
                      ifelse(max_forklength > 200 &
                             max_forklength <= 300, '300mm',
                      ifelse(max_forklength > 300 &
                             max_forklength <= 500, '500mm', max_forklength)))))) %>%
  group_by(site, maxfl_bin) %>%
  summarize(available = n()) %>%
  pivot_wider(names_from = maxfl_bin, values_from = available) %>%
  relocate('50mm', .before='200mm') %>%
  relocate('100mm', .before='200mm')
September_lepomis_sm[is.na(September_lepomis_sm)] <- 0
September_lepomis_sm 

## September large-mouthed sunfish (megalotis and macrochirus) size class sample distribution
# <50mm missing n at GC, low n at SF, PD, MR, PL, EM
# <100mm missing n at SF, MR, low n at AR, PL, GC, EM
# <200mm missing n at SF, AR, WM, low n at MR, PL, GC, EM

# Recommend analyzing hydrogen in large-mouthed lepomis in the 100-300 size class range ( n = 18)
September_lepomis_lm <- d %>%
  filter(species %in% c('L.gulosus', 'L.cyanellus') ) %>%
  filter(collection == '2020-09-15') %>%
  mutate(maxfl_bin = ifelse(max_forklength <= 50, '50mm',
                     ifelse(max_forklength > 50 & 
                            max_forklength <= 100 , '100mm',
                      ifelse(max_forklength > 100 &
                             max_forklength <= 200, '200mm',
                      ifelse(max_forklength > 200 &
                             max_forklength <= 300, '300mm',
                      ifelse(max_forklength > 300 &
                             max_forklength <= 500, '500mm', max_forklength)))))) %>%
  group_by(site, maxfl_bin) %>%
  summarize(available = n()) %>%
  pivot_wider(names_from = maxfl_bin, values_from = available) %>%
  relocate('50mm', .before='100mm')
September_lepomis_lm[is.na(September_lepomis_lm)] <- 0
September_lepomis_lm




coarse_vs_leaf <- d %>%
  filter(species %in% c('periphyton', 'coarse_detritus', 'tree_leaf') ) %>%
  ggplot(aes(x=site, y = carbon, fill = species)) +
           geom_boxplot()






















