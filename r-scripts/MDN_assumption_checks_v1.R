# Marine-Derived Nutrient Assay
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
library(ggpubr)
library(gridExtra)
library(ggpmisc)

## set work directory

setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

d <- read_csv("data/MDN_clean.csv")

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
  arrange(SITE) %>%
  mutate(PPTAVG_SITE = ifelse(SITE == 'TR', # TR missing gauge data (import SF, shared watershed)
                              d %>% filter(SITE == 'SF') %>% pull(PPTAVG_SITE) %>% unique(),
                              PPTAVG_SITE) )

dlong <- d %>%
  pivot_longer(cols = all_of(c('CARBON', 'SULFUR', 'NITROGEN')),
               names_to = "ELEMENT",
               values_to = "PERMIL") %>%
  mutate(ELEMENT = factor(ELEMENT)) %>%
  mutate(ELEMENT = fct_relevel(ELEMENT,c("CARBON", "SULFUR", "NITROGEN"))) %>%
  arrange(PPTAVG_SITE)

dlong_scaled <- d %>% 
  mutate(CARBON = scale(CARBON),
         SULFUR = scale(SULFUR),
         NITROGEN = scale(NITROGEN) ) %>%
  pivot_longer(cols = all_of(c('CARBON', 'SULFUR', 'NITROGEN')),
               names_to = "ELEMENT",
               values_to = "PERMIL") %>%
  mutate(ELEMENT = factor(ELEMENT)) %>%
  mutate(ELEMENT = fct_relevel(ELEMENT,c("CARBON", "SULFUR", "NITROGEN"))) %>%
  arrange(PPTAVG_SITE)


# site comparisons:

## boxplots C N S
(box_elements <-  dlong_scaled %>%
    mutate(dlong_scaled, ELEMENT = str_to_title(ELEMENT)) %>%
    mutate(ELEMENT = fct_relevel(ELEMENT,c("Carbon", "Sulfur", "Nitrogen"))) %>%
  ggplot(aes(x = PERMIL, 
             y = SITE, 
             fill = PPTAVG_SITE)) +
  facet_wrap(facets = 'ELEMENT', nrow=1) +
  geom_boxplot() +
  scale_fill_viridis_c(trans = 'reverse') +
  theme_classic(base_size = 24) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    labs(fill = 'Rain (cm/y)') +
    ggtitle(expression(paste('Similar Patterns in ', sigma, 'C'^'13', ' and ', sigma, 'S'^'34', ' Across Sites'))))
      
      
      paste('Similar ', sigma, 'C'^'13',' and ', sigma, 'S'^'34', 'Patterns Across Sites')))

## Biplot C S
(biplot_site <- d %>%
  filter(SITE_TYPE == 'Stream') %>%
  ggplot(aes(x = CARBON, y = SULFUR, size = NITROGEN, fill = PPTAVG_SITE)) +
  stat_ellipse(aes(color = as.factor(SITE)), show.legend = FALSE) +
  geom_point(shape = 21, color = 'black') +
  scale_fill_viridis_c(trans = 'reverse') +
  scale_color_viridis_d(direction = -1) +
  theme_classic(base_size = 24) +
  xlab(expression(paste(sigma, 'C'^'13'))) +
  ylab(expression(paste(sigma, 'S'^'34'))) +
    labs(fill = 'Rain (cm/y)',
         size = expression(paste(sigma,'N'^'15'))) +
    ggtitle('Arid Stream Signatures Appear Differernt From Humid Streams'))

# Biplot C S for site Sources
(biplot_flora_site <- d %>%
    mutate(PPTAVG_SITE = round(PPTAVG_SITE, 0)) %>%
    filter(GUILD %in% c('Aquatic', 'Terrestrial')) %>%
    ggplot(aes(x = CARBON, 
               y = SULFUR)) +
    geom_point(aes(shape = TAXON,
                   fill = as.factor(PPTAVG_SITE)),
               size = 4) +
    scale_shape_manual(values = 21:25) +
    scale_fill_viridis_d(direction = -1,
                         labels=c('71', '82', '83', '93', '99', '106','108', '112', '113', 'Bay')) +
    guides(fill = guide_legend(
      override.aes=list(shape = 22, size = 8)),
      shape = guide_legend(
        override.aes = list(size = 7)) ) +
    theme_classic(base_size = 24) +
    labs(shape = 'Source Type',
         fill = 'Rainfall\n(cm/yr)') +
    xlab(expression(paste(sigma, 'C'^'13'))) +
    ylab(expression(paste(sigma, 'S'^'34'))) +
    scale_y_continuous(breaks=seq(-10,20,5)) +
    geom_text(data = filter(d, SULFUR < -4),
               label = filter(d, SULFUR < -4) %>% pull(UID),
               aes(x = CARBON,
                   y = SULFUR),
               nudge_x = 4) )

# Remove: Outlier macrophyte at Nueces Bay with sulfur < -4

d <- d %>%
  filter(SULFUR > -4)

(biplot_flora_site_2 <- d %>%
  mutate(PPTAVG_SITE = round(PPTAVG_SITE, 0)) %>%
  filter(GUILD %in% c('Aquatic', 'Terrestrial')) %>%
  ggplot(aes(x = CARBON, 
             y = SULFUR)) +
    geom_point(aes(color = SITE_TYPE),
               size = 3) +
  geom_point(aes(shape = TAXON,
                 fill = as.factor(PPTAVG_SITE)),
             size = 4) +
  stat_ellipse(aes(color = SITE_TYPE),
                 linetype=2, 
                 size = 2) +
  scale_shape_manual(values = 21:25) +
  scale_fill_viridis_d(direction = -1,
                       labels=c('71', '82', '83', '93', '99', '106','108', '112', '113', 'Bay')) +
  guides(fill = guide_legend(
    override.aes=list(shape = 22, size = 8)),
    shape = guide_legend(
      override.aes = list(size = 7)) ) +
  theme_classic(base_size = 24) +
  labs(shape = 'Source Type',
       fill = 'Rain (cm/y)',
       color = 'Site Type') +
  xlab(expression(paste(sigma, 'C'^'13'))) +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  scale_y_continuous(breaks=seq(-10,20,5)) +
  ggtitle('Estuarine Sources Have A Distinct Signature'))

# MANVOVA sources: site type 
my_sources <- c('Aquatic', 'Terrestrial', 'Detritus')

man_src_site_type <- manova(cbind(CARBON, SULFUR, NITROGEN) ~ SITE_TYPE , data = filter(d, GUILD %in% my_sources))

(summary.aov(man_src_site_type))

(box_flora_site_type <- dlong %>%
  filter(GUILD %in% my_sources) %>%
  ggplot(aes(x = SITE_TYPE, y = PERMIL, fill = SITE_TYPE)) +
  geom_boxplot(size = 1.5) +
  facet_wrap(facets = 'ELEMENT') +
  geom_signif(comparison = list(c("Bay",  "Stream")), 
              test = "t.test",
              textsize = 6) +
  ylim(c(-30,40)) +
  theme_classic(base_size = 24) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(fill = 'Locality') +
  ylab(expression(paste('Permil (', sigma, ')'))) +
  ggtitle('Estuarine Basal Resources Display Unique Signatures') +
  theme(legend.key.size = unit(1.5, 'cm'), 
        legend.key.height = unit(1.5, 'cm'),
         legend.key.width = unit(1.5, 'cm'),
         legend.title = element_text(size=24),
         legend.text = element_text(size=18)) )

# Do Sources in Stream vary by type?
d_src <- filter(d, GUILD %in% my_sources) %>%
  filter(SITE_TYPE == 'Stream') %>%
  mutate(TAXON = ifelse(TAXON =='Coarse Detritus',
                        'Detritus', TAXON)) %>%
  mutate(TAXON = factor(TAXON)) %>%
  mutate(TAXON = fct_relevel(TAXON,c("Periphyton", "Algae", "Macrophyte", "Detritus", "Riparian")))

man_flora_src_type <- manova(cbind(CARBON, SULFUR, NITROGEN) ~ TAXON , data = d_src)

(summary.aov(man_flora_src_type))

TukeyHSD(aov(CARBON ~ TAXON, 
             data = d_src), 
         conf.level=.95) 

# boxplot with specified mean comparisons
src_comparisons <- list( c("Periphyton", "Riparian"), 
                        c("Periphyton", "Algae"), 
                        c("Periphyton", "Detritus"),
                        c("Periphyton", "Macrophyte"))
sbox <- function(data, my_group, my_element) {
  ggboxplot(data, x = my_group, y = my_element) + 
     stat_compare_means(comparisons = src_comparisons,
                        label = "p.signif",
                        size = 5)+ # pairwise comparisons p-value
     # stat_compare_means() +    # global p-value
     theme_classic(base_size = 24) +
     theme(axis.title.x=element_blank(),
           axis.ticks.x=element_blank())
}

# Set C comparisons
src_comparisons <- list( c("Periphyton", "Algae"), 
                         c("Periphyton", "Macrophyte"), 
                         c("Periphyton", "Detritus"),
                         c("Periphyton", "Riparian"))
# C
box_src_c <- sbox(d_src, "TAXON", "CARBON") +
  ylab(expression(paste(sigma, 'C'^'13'))) +
  ylim(c(-40, 0)) +
  theme(axis.text.x=element_blank())

# Set S and N comparisons
src_comparisons <- list( c("Riparian", "Detritus"), 
                         c("Riparian", "Macrophyte"), 
                         c("Riparian", "Algae"),
                         c("Riparian", "Periphyton"))
# S
box_src_s <- sbox(d_src, "TAXON", "SULFUR") +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  ylim(c(-5,30)) +
  theme(axis.text.x=element_blank())

# N
box_src_n <- sbox(d_src, "TAXON", "NITROGEN") +
  ylab(expression(paste(sigma, 'N'^'15'))) +
  ylim(c(-5,40))

# Plot C, S, and N
(box_src <- grid.arrange(box_src_c, box_src_s, box_src_n, nrow = 3) )

## C: Periphyton differs from other sources
## S & N : Riparian and detritus differ from other sources
## numerous N outliers in riparian

###########################################################

# Inquiries:
## Is there differences among source types at bay sites

d_src_bay <- filter(d, GUILD %in% my_sources) %>%
  filter(SITE_TYPE == 'Bay') %>%
  mutate(TAXON = ifelse(TAXON =='Coarse Detritus',
                        'Detritus', TAXON)) %>%
  mutate(TAXON = factor(TAXON)) %>%
  mutate(TAXON = fct_relevel(TAXON,c("Algae", "Macrophyte", "Detritus")))

man_flora_src_type_b <- manova(cbind(CARBON, SULFUR, NITROGEN) ~ TAXON , data = d_src)

(summary.aov(man_flora_src_type_b))

TukeyHSD(aov(CARBON ~ TAXON, 
             data = d_src), 
         conf.level=.95) 

TukeyHSD(aov(SULFUR ~ TAXON, 
             data = d_src), 
         conf.level=.95) 

TukeyHSD(aov(NITROGEN ~ TAXON, 
             data = d_src), 
         conf.level=.95) 

dlong_src_bay <- filter(dlong, GUILD %in% my_sources) %>%
  filter(SITE_TYPE == 'Bay') %>%
  mutate(TAXON = ifelse(TAXON =='Coarse Detritus',
                        'Detritus', TAXON)) %>%
  mutate(TAXON = factor(TAXON)) %>%
  mutate(TAXON = fct_relevel(TAXON,c("Algae", "Macrophyte", "Detritus")))

bay_comparisons <- list( c("NB", "BB"), 
                         c("NB", "CB"),
                         c("NB", "HB"),
                         c("NB", "LB"))

dlong_src_bay %>%
  filter(TAXON %in% c('Algae', 'Macrophyte')) %>%
ggplot(aes(x = SITE, y = PERMIL)) +
  geom_boxplot() +
  geom_point(aes(shape = TAXON),
             size = 3) +
  scale_shape_manual(values = 21:22) +
  facet_wrap(~ ELEMENT) +
  geom_signif(comparison = bay_comparisons,
              textsize = 6,
              map_signif_level = TRUE,
              y_position = c(-10,-7,-4,-1),
              tip_length = 0) +
  ggtitle('Do Algae & Macrophyte values vary by bay site?') +
  ylab(expression(paste(sigma))) +
  theme_classic(base_size = 24)

# boxplot with specified mean comparisons
src_comparisons <- list( c("Algae", "Macrophyte"), 
                         c("Algae", "Detritus"), 
                         c("Macrophyte", "Detritus"))
# C
box_src_c_b <- sbox(d_src_bay, "TAXON", "CARBON") +
  ylab(expression(paste(sigma, 'C'^'13'))) +
  ylim(c(-40, 0)) +
  theme(axis.text.x=element_blank())

# S
box_src_s_b <- sbox(d_src_bay, "TAXON", "SULFUR") +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  ylim(c(-5,30)) +
  theme(axis.text.x=element_blank())

# N
box_src_n_b <- sbox(d_src_bay, "TAXON", "NITROGEN") +
  ylab(expression(paste(sigma, 'N'^'15'))) +
  ylim(c(-5,40))

# Plot C, S, and N
(box_src_bays <- grid.arrange(box_src_c_b, box_src_s_b, box_src_n_b, nrow = 3) )

# C Algae differs from macrophytes and detritus
# Other comparisons show no difference in means

d_src_bay %>%
  ggplot(aes(CARBON, SULFUR, shape = TAXON, color = SITE)) +
  geom_point(size = 6) +
  ggtitle('Do Sources Differ At Bay Sites')

dlong_src_bay %>%
  group_by(ELEMENT) %>%
  summarize(iso_mean = mean(PERMIL, na.rm=TRUE),
            iso_sd = sd(PERMIL, na.rm = TRUE),
            iso_n = n())

# remove detritus
dlong_src_bay %>%
  filter(TAXON != 'Detritus') %>%
  group_by(ELEMENT) %>%
  summarize(iso_mean = mean(PERMIL, na.rm=TRUE),
            iso_sd = sd(PERMIL, na.rm = TRUE),
            iso_n = n())

# Remove NB and Detritus
d_src_bay %>%
  filter(SITE != 'NB') %>%
  filter(TAXON != 'Detritus') %>%
  ggplot(aes(CARBON, SULFUR, shape = TAXON, color = SITE)) +
  geom_point(size = 6) +
  ggtitle('Do Sources Differ At Bay Sites')

dlong_src_bay %>%
  filter(SITE != 'NB') %>%
  filter(TAXON != 'Detritus') %>%
  group_by(ELEMENT) %>%
  summarize(iso_mean = mean(PERMIL, na.rm=TRUE),
            iso_sd = sd(PERMIL, na.rm = TRUE),
            iso_n = n())

# Bay end members (aggregating algae and macrophytes from 4/5 estuarine sites) (remove NB because values are different from other bays and the location might not have been as marine as we expected due to large freshwater discharge in the channel. Remove detritus since it cannot be confirmed where it originated)  
dlong_src_bay %>%
  filter(SITE != 'NB') %>%
  group_by(ELEMENT) %>%
  summarize(iso_mean = mean(PERMIL, na.rm=TRUE),
            iso_sd = sd(PERMIL, na.rm = TRUE),
            iso_n = n())

# Stream end member values (aggregating all sources and sites)
dlong %>%
  filter(SITE_TYPE == 'Stream') %>%
  filter(GUILD %in% c('Aquatic', 'Terrestrial')) %>%
  group_by(ELEMENT) %>%
  summarize(iso_mean = mean(PERMIL, na.rm=TRUE),
            iso_sd = sd(PERMIL, na.rm = TRUE),
            iso_n = n())
# end members will have to be corrected for uneven sampling distributions across source types and sites. Recommend generating normal distributions for each source type using mean and sd and pulling samples evenly across types with equivalent total n.
dlong %>%
  filter(GUILD %in% c('Aquatic', 'Terrestrial', 'Detritus')) %>%
  filter(SITE != 'NB') %>%
  group_by(SITE_TYPE, TAXON, ELEMENT) %>%
  summarize(n = n(),
            iso_mean = mean(PERMIL, na.rm=TRUE)) %>%
  pivot_wider(names_from = ELEMENT, values_from = iso_mean)

# Consider removing periphyton as an end member since it is between all other stream sources and estuarine source values. Could periphyton be influenced by marine signal?


#########################################################
# Is bay resource value different from terrestrial?

( box_src_site_type <- dlong %>%
    filter(GUILD %in% c('Aquatic', 
                        'Terrestrial', 
                        'Detritus')) %>%
    filter(SITE != 'NB') %>%
    filter(TAXON != 'Periphyton') %>%
  mutate(ELEMENT = str_to_title(ELEMENT)) %>%
  mutate(ELEMENT = factor(ELEMENT)) %>%
  mutate(ELEMENT = fct_relevel(ELEMENT,c("Carbon", "Sulfur", "Nitrogen"))) %>%
  ggplot(aes(x = SITE_TYPE, y = PERMIL)) +
  geom_boxplot() +
  facet_wrap(facets = 'ELEMENT') +
  geom_signif(comparison = list(c("Bay",  "Stream")), 
              test = "t.test",
              textsize = 6,
              y_position = -10,
              map_signif_level = TRUE) +
  theme_classic(base_size = 24) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab(expression(paste('Permil (', sigma, ')'))) +
  ggtitle('Do Sources differ between Bays and Streams?') )

## Aggregating all sources except periphyton. End members still need correction for uneven sampling.


# Correcting Source end-members for uneven sampling
options(digits = 4)
set.seed(10)
## summarize data to inform normal distributions
src_summary <- d %>%
  filter(GUILD %in% c('Aquatic','Terrestrial','Detritus')) %>%
  filter(TAXON != 'periphyton') %>%
  group_by(SITE_TYPE, TAXON) %>%
  summarize(n=n(),
            mean_c = mean(CARBON, na.rm=TRUE),
            sd_c =   sd(CARBON, na.rm=TRUE),
            mean_s = mean(SULFUR, na.rm = TRUE),
            sd_s =   sd(SULFUR, na.rm=TRUE),
            mean_n = mean(NITROGEN, na.rm=T),
            sd_n =   sd(NITROGEN, na.rm=T) )
## include total samples and types used for site types
src_summary <- src_summary %>%
  group_by(SITE_TYPE) %>%
  summarize(total = sum(n),
            srcs = n()) %>%
  right_join(src_summary)
## sample theoretical distributions

i_sources <- list()

for(i in 1:length(src_summary$n)) {
  i_sources$SITE_TYPE[[i]] <- src_summary$SITE_TYPE[i]
  
  i_sources$SRC[[i]] <- src_summary$TAXON[i]

  i_sources$C[[i]] <- rnorm(n = 200,
      mean = src_summary[i,] %>% pull(mean_c),
      sd = src_summary[i,] %>% pull(sd_c)) %>%
  sample(src_summary[i,] %>% 
           pull(total) / src_summary[i,] %>% 
           pull(srcs))

  i_sources$S[[i]] <- rnorm(n = 200,
           mean = src_summary[i,] %>% pull(mean_s),
           sd = src_summary[i,] %>% pull(sd_s)) %>%
  sample(src_summary[i,] %>% 
           pull(total) / src_summary[i,] %>% 
           pull(srcs))

  i_sources$N[[i]] <- rnorm(n = 200,
           mean = src_summary[i,] %>% pull(mean_n),
           sd = src_summary[i,] %>% pull(sd_n)) %>%
  sample(src_summary[i,] %>% 
           pull(total) / src_summary[i,] %>% 
           pull(srcs))
}

c_names <- src_summary %>%
  mutate(ID = paste('C',SITE_TYPE, TAXON, sep = "_")) %>%
  pull(ID)

s_names <- src_summary %>%
  mutate(ID = paste('S',SITE_TYPE, TAXON, sep = "_")) %>%
  pull(ID)

n_names <- src_summary %>%
  mutate(ID = paste('N',SITE_TYPE, TAXON, sep = "_")) %>%
  pull(ID)

i_bay <- tibble(i_sources$C[[1]],
       i_sources$C[[2]],
       i_sources$C[[3]],
       i_sources$S[[1]],
       i_sources$S[[2]],
       i_sources$S[[3]],
       i_sources$N[[1]],
       i_sources$N[[2]],
       i_sources$N[[3]])
colnames(i_bay) = c(c_names[1:3], s_names[1:3], n_names[1:3])
i_bay <- mutate(i_bay, SAMPLE = 1:length(i_bay$C_Bay_Algae))
i_bay <- i_bay %>% pivot_longer(cols = C_Bay_Algae:N_Bay_Macrophyte,
                                values_to = 'PERMIL',
                                names_to = 'ID') %>%
  separate(ID, c('ELEMENT', 'SITE_TYPE', 'TAXON'), sep = "_") %>%
  mutate(ELEMENT = ifelse(ELEMENT == 'C', 'CARBON',
                          ifelse(ELEMENT == 'S', 'SULFUR',
                                 'NITROGEN')))

i_stream <- tibble(i_sources$C[[4]],
                i_sources$C[[5]],
                i_sources$C[[6]],
                i_sources$C[[7]],
                i_sources$C[[8]],
                i_sources$S[[4]],
                i_sources$S[[5]],
                i_sources$S[[6]],
                i_sources$S[[7]],
                i_sources$S[[8]],
                i_sources$N[[4]],
                i_sources$N[[5]],
                i_sources$N[[6]],
                i_sources$N[[7]],
                i_sources$N[[8]])
colnames(i_stream) = c(c_names[4:8], s_names[4:8], n_names[4:8])
i_stream <- mutate(i_stream, SAMPLE = 1:length(i_stream$C_Stream_Algae))

i_stream <- i_stream %>% pivot_longer(cols = C_Stream_Algae:N_Stream_Riparian,
                                values_to = 'PERMIL',
                                names_to = 'ID') %>%
  separate(ID, c('ELEMENT', 'SITE_TYPE', 'TAXON'), sep = "_") %>%
  mutate(ELEMENT = ifelse(ELEMENT == 'C', 'CARBON',
                          ifelse(ELEMENT == 'S', 'SULFUR',
                                 'NITROGEN')))
## Theoretical source data
i_src <- full_join(i_bay, i_stream)

(biplot_i_src <- i_src %>%
  pivot_wider(names_from = ELEMENT,
              values_from = PERMIL) %>%
  ggplot(aes(x=CARBON, 
             y=SULFUR)) +
  stat_ellipse(aes(color = SITE_TYPE),
               linetype=2, 
               size = 2) +
  geom_point(aes(shape = TAXON,
                 fill = SITE_TYPE),
             size = 4) +
  scale_shape_manual(values = 21:26) +
  guides(fill = guide_legend(
    override.aes=list(shape = 22, size = 8)),
    shape = guide_legend(
      override.aes = list(size = 7)) ) +
  theme_classic(base_size = 24) +
  labs(shape = 'Source Type',
       color = 'Locality',
       fill = 'Locality') +
  xlab(expression(paste(sigma, 'C'^'13'))) +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  ggtitle('Theoretical Source Distribution') )

#################################################################
## Is periphyton different from bay resources
( box_bay_peri <- dlong %>%
    mutate(P_test = ifelse(TAXON == 'Periphyton', 'Periphyton',
                           ifelse(SITE_TYPE == 'Bay' & 
                                     TAXON %in% c('Algae', 'Macrophyte'), 
                                   'Bay Source', NA))) %>%
    filter(is.na(P_test)== FALSE) %>%
    mutate(ELEMENT = str_to_title(ELEMENT)) %>%
    mutate(ELEMENT = factor(ELEMENT)) %>%
    mutate(ELEMENT = fct_relevel(ELEMENT,c("Carbon", "Sulfur", "Nitrogen"))) %>%
    ggplot(aes(x = P_test, y = PERMIL)) +
    geom_boxplot() +
    facet_wrap(facets = 'ELEMENT') +
    geom_signif(comparison = list(c('Bay Source',  "Periphyton")), 
                test = "t.test",
                textsize = 6,
                y_position = -10,
                map_signif_level = TRUE) +
    theme_classic(base_size = 24) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab(expression(paste('Permil (', sigma, ')'))) +
    ggtitle('Is Stream Periphyton different from bay resources?') )

## Bay sources differ from stream periphyton for Sulfur but not Carbon or Nitrogen

## Is periphyton different from Stream resources
( box_bay_peri <- dlong %>%
    filter(SITE_TYPE == 'Stream') %>%
    filter(GUILD %in% c('Aquatic', 
                        'Terrestrial', 
                        'Detritus')) %>%
    mutate(P_test = ifelse(TAXON == 'Periphyton', 
                           'Periphyton',
                           'Stream Src')) %>%
    filter(is.na(P_test)== FALSE) %>%
    mutate(ELEMENT = str_to_title(ELEMENT)) %>%
    mutate(ELEMENT = factor(ELEMENT)) %>%
    mutate(ELEMENT = fct_relevel(ELEMENT,c("Carbon", "Sulfur", "Nitrogen"))) %>%
    ggplot(aes(x = P_test, y = PERMIL)) +
    geom_boxplot() +
    facet_wrap(facets = 'ELEMENT') +
    geom_signif(comparison = list(c('Stream Src',  "Periphyton")), 
                test = "t.test",
                textsize = 6,
                y_position = -10,
                map_signif_level = TRUE) +
    theme_classic(base_size = 24) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab(expression(paste('Permil (', sigma, ')'))) +
    ggtitle('Is Stream Periphyton different from stream resources?') )

## Bay sources differ from stream periphyton for Sulfur but not Carbon or Nitrogen

## Stream sources differ from periphyton for Carbon and Sulfur

# Periphyton should not be used when calculating end members for terrestrial vs estuarine signals since they are in between both. It is an interesting discussion point as periphyton may represent linkages in nutrient cycling between the 2 systems.

#########################################################################
## Are there site effects on source values?
d_str_src <- dlong %>%
  filter(SITE_TYPE == 'Stream') %>%
  filter(GUILD %in% c('Aquatic', 
                      'Terrestrial', 
                      'Detritus')) %>%
  filter(TAXON != 'Periphyton') %>%
    mutate(ELEMENT = str_to_title(ELEMENT)) %>%
    mutate(ELEMENT = factor(ELEMENT)) %>%
    mutate(ELEMENT = fct_relevel(ELEMENT,c("Carbon", "Sulfur", "Nitrogen")))

d_str_src %>%
  ggplot(aes(x = SITE, y = PERMIL)) +
    geom_boxplot() +
    facet_wrap(~ ELEMENT) +
    theme_classic(base_size = 24) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab(expression(paste('Permil (', sigma, ')'))) +
    ggtitle('Do Sources Vary By Site?') 




str(TukeyHSD(aov(PERMIL ~ SITE, 
             data = d_str_src %>% filter(ELEMENT == 'Carbon')), 
         conf.level=.95) )

t1 <- TukeyHSD(aov(PERMIL ~ SITE, 
             data = d_str_src %>% filter(ELEMENT == 'Carbon')), 
         conf.level=.95) 
as_tibble(t1$SITE) %>%
  mutate(site = rownames(t1$SITE)) %>%
  filter(`p adj` < 0.05) # TR is different from 4 sites Carbon


TukeyHSD(aov(PERMIL ~ SITE, 
             data = d_str_src %>% filter(ELEMENT == 'Sulfur')), 
         conf.level=.95)

t2 <- TukeyHSD(aov(PERMIL ~ SITE, 
                   data = d_str_src %>% filter(ELEMENT == 'Sulfur')), 
               conf.level=.95) 
as_tibble(t2$SITE) %>%
  mutate(site = rownames(t2$SITE)) %>%
  filter(`p adj` < 0.05) # 18 site pairs differ in Sulfur


t3 <- TukeyHSD(aov(PERMIL ~ SITE, 
                   data = d_str_src %>% filter(ELEMENT == 'Nitrogen')), 
               conf.level=.95) 
as_tibble(t2$SITE) %>%
  mutate(site = rownames(t3$SITE)) %>%
  filter(`p adj` < 0.05) # 18 site pairs differ in Nitrogen

## What sources should be used as the 'local' source signal
# all sources except periphyton
# all sites except TR


# Do migrants differ from local residents?

# manova C and S for migrant status
options(digits = 2)

man_migrant <- manova(cbind(CARBON, SULFUR, NITROGEN) ~ MIGRANT , data = filter(d, GUILD %in% c('Fish', 'Invertebrate')))

str(summary.aov(man_migrant))

summary.aov(man_migrant)[[2]][[5]][1]


(man_migrant_p_values <- tibble('Carbon' = summary.aov(man_migrant)[[1]][[5]][1] , 
                                'Sulfur' = summary.aov(man_migrant)[[2]][[5]][1] , 'Nitrogen' = summary.aov(man_migrant)[[3]][[5]][1] ) %>%
    pivot_longer(cols = Carbon:Nitrogen,
                 names_to = 'MANOVA',
                 values_to = 'p') %>%
    mutate(Sig = ifelse(p < .01 , '**',
                        ifelse(p < 0.05, "*",
                               'NS'))) %>%
    rename('p value' = p))

# Setup plot environment
biplot_migrant <- d %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  ggplot(aes(x = CARBON, 
             y = SULFUR,
             shape = MIGRANT,
             fill = MIGRANT)) +
  scale_shape_manual(values = 21:22) +
  theme_classic(base_size = 24) +
  guides(fill = guide_legend(
    override.aes=list(shape = 22, size = 8)),
    shape = guide_legend(
      override.aes = list(size = 7)) ) +
  theme_bw(base_size = 24) +
  xlab(expression(paste(sigma, 'C'^'13'))) +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  theme(legend.position = c(.9, .1)) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(colour = 'grey', 
                                         fill = 'white',
                                         linetype='solid'))

# Biplot:
biplot_migrant +
  geom_point(size = 4) +
  stat_ellipse(aes(color = MIGRANT),
               linetype = 5,
               size = 1,
               show.legend = FALSE) +
  annotate(geom = "table",
           x = -37,
           y = 15,
           label = list(man_migrant_p_values),
           size = 7) +
  ggtitle('Migrants Have Distinct Estuarine Signal') 

# plot estuarine boundary
biplot_migrant +
  geom_segment(aes(x=-35,xend=-20,y=10,yend=10),
               color = 'grey50') +
  geom_segment(aes(x=-20, xend=-20, y=-5, yend=10),
               color = 'grey50') +
  geom_point(size = 4) +
  stat_ellipse(aes(color = MIGRANT),
               linetype = 5,
               size = 1,
               show.legend = FALSE) +
  ggtitle('Some Residents Share Estuarine Signature') 

# label diadromous
diadromous <- d %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  filter(MIGRANT == 'Diadromous') %>%
  pull(TAXON) %>%
  unique()

d %>%
  filter(TAXON %in% diadromous) %>%
  ggplot(aes(x = CARBON, 
             y = SULFUR,
             fill = TAXON)) +
  geom_segment(aes(x=-35,xend=-20,y=10,yend=10),
               color = 'grey50') +
  geom_segment(aes(x=-20, xend=-20, y=-5, yend=10),
               color = 'grey50') +
  theme_classic(base_size = 24) +
  guides(fill = guide_legend(
    override.aes=list(size = 8)),
    shape = guide_legend(
      override.aes = list(size = 7)) ) +
  theme_bw(base_size = 24) +
  xlab(expression(paste(sigma, 'C'^'13'))) +
  ylab(expression(paste(sigma, 'S'^'34'))) +
  geom_label(label = d %>%
               filter(TAXON %in% diadromous) %>% 
               pull(TAXON),
             show.legend = FALSE) +
  ggtitle('Consistent Estuarine Signal Within Migrant Taxa')

# label estuarine-influenced Residents
influenced <- filter(d, CARBON > -20 | SULFUR > 10) %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  filter(MIGRANT == 'Potamodromous') %>%
  pull(TAXON) %>%
  unique()

bp <- function(my_sp, highlight) {
  d %>%
    filter(TAXON %in% influenced) %>%
    ggplot(aes(x = CARBON, 
               y = SULFUR,
               fill = TAXON)) +
    geom_segment(aes(x=-35,xend=-20,y=10,yend=10),
                 color = 'grey50') +
    geom_segment(aes(x=-20, xend=-20, y=-5, yend=10),
                 color = 'grey50') +
    theme_bw(base_size = 24) +
    xlab(expression(paste(sigma, 'C'^'13'))) +
    ylab(expression(paste(sigma, 'S'^'34'))) +
    ggtitle('Estuary Signature Varies Within Resident Taxa') +
    annotate(geom = 'text',
             x = -34,
             y = 14,
             label = highlight,
             size = 7) +
    stat_chull(data = filter(d, TAXON == highlight),
               alpha = .4,
               geom = 'polygon',
               show.legend = FALSE) +
    geom_label(data = filter(d, TAXON %in% my_sp),
               label = filter(d, TAXON %in% my_sp) %>%
                 pull(TAXON),
               show.legend = FALSE) }

bp(influenced[1],influenced[1])

p <- list()
for(i in 1:length(influenced)) {
  p[[i]] <- bp(influenced[1:i], influenced[i])
}

## labeling individual, estuarine-influenced taxa
p[1]
p[2]
p[3]
p[5]
p[6]
p[9]
p[10]


# Marine fauna contain greater proportions of 34S and 13C (increased values) compared to stream flora. This is consistent with trends exhibited in other marine-derived nutrient surveys.

# All but one amphidromous samples have 34S and 13C values consistent with a mixed diet of marine and terrestrial sources. One Mugil cephalus resembles potamodromous signatures implying that this estuarine individual feeds predominately in freshwater habitat. Stream fauna contain values of 34S and 13C that are within the 95% ellipses for 'sources' in this biplot.

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