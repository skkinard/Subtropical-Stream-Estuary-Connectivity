
library(tidyverse)
library('sjmisc')
library(lubridate)
library(ggpubr)
library(rstatix)


d <- read_csv('data/TERRG_CNS_combined.csv') %>% 
  filter(site %in% c("AR","EM", "GC", "MR", "PD", "PL", "SF", "TR", "WM") )

d$site <- factor(d$site, levels = c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

# sample tables

all_samples <- d %>%
  group_by(site, collection, type) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples) %>%
    as.data.frame()

fish <- d %>%
  filter(type == 'Fish') %>%
  group_by(site, collection) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples)

invert <- d %>%
  filter(type == 'Invertebrate') %>%
  group_by(site, collection) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples)

aquatic <- d %>%
  filter(type == 'Aquatic') %>%
  group_by(site, collection) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples)
 
terrestrial <-  d %>%
  filter(type == 'Terrestrial') %>%
  group_by(site, collection) %>%
  summarize(samples = n() ) %>%
  pivot_wider(names_from = collection, values_from = samples)

all_samples
fish
invert
aquatic
terrestrial

# Analytical reps
(arep_table <- d %>% 
  filter(details %in% c('"1/2"', '"2/2"')) %>%
  mutate(powder = paste(site, species, sep = '_'),
         analytical_rep = substr(details, 2, 2)) %>%
  select(carbon, nitrogen, powder, analytical_rep, species) %>%
  pivot_wider(names_from = analytical_rep, values_from = carbon:nitrogen) %>%
  mutate(c_diff = carbon_1 - carbon_2,
         n_diff = nitrogen_1 - nitrogen_2))

(arep_detritus_mu <- arep_table %>%
  filter(species == 'coarse_detritus') %>%
  select(c_diff, n_diff) %>%
  abs() %>%
  colMeans() %>%
  as_tibble() %>%
  t() 
colnames(arep_detritus_mu) <- c('c_mu', 'n_mu')
rownames(arep_detritus_mu) <- '')
  

(arep_detritus_variance_c <- arep_table %>%
  filter(species == 'coarse_detritus') %>%
  select(c_diff) %>%
  var()  %>%
  unlist() %>%
  unname() )
  
(arep_detritus_variance_n <- arep_table %>%
  filter(species == 'coarse_detritus') %>%
  select(n_diff) %>%
  var()  %>%
  unlist() %>%
  unname() )

arep_detritus_table <- cbind(rep('detritus',1),
                             round(arep_detritus_mu, 4), 
                             round(arep_detritus_variance_c, 4),
                                   round(arep_detritus_variance_n, 4) )
colnames(arep_detritus_table) <- c('type', 'c_mu', 'n_mu', 'c_var', 'n_var')
arep_detritus_table <- as_tibble(arep_detritus_table) %>%
  pivot_longer(cols = c(n_mu, n_var, c_mu,c_var), names_to = 'isotope', values_to = 'value') %>%
  mutate(statistic = substr(isotope, 3,4),
         isotope = substr(isotope,1,1) ) %>%
  mutate(statistic = ifelse(statistic == 'mu', 'mean_abs',
                            ifelse(statistic == 'va', 'variance', NA)),
         isotope = ifelse(isotope == 'n', 'nitrogen',
                          ifelse(isotope == 'c', 'carbon', NA))) %>%
  pivot_wider(names_from = statistic, values_from = value)

arep_detritus_table

# specific source comparisons
# filamentous_algae
var_table <- function(x) {
  d %>%
    filter(species == x) %>%
    group_by(site) %>%
    summarize(n=n(),
              c_mean = mean(carbon),
              c_var = var(carbon),
              n_mean = mean(nitrogen),
              n_var = var(nitrogen)) }
carbon_plot <- function(x) {
  d %>%
    filter(species == x) %>%
    ggplot(aes(x = collection, y = carbon)) +
    geom_point(size = 5, shape = 3) +
    geom_line(linetype=2) +
    facet_wrap(facets = 'site') +
    ggtitle(x) +
    xlab('Collection Date') +
    ylab('Carbon') +
    theme_grey(base_size = 16) +
    scale_x_date(date_labels="%b-%d",date_breaks  ="2 month") }

nitrogen_plot <- function(x) {
  d %>%
    filter(species == x) %>%
    ggplot(aes(x = collection, y = nitrogen)) +
    geom_point(size = 5, shape = 3) +
    geom_line(linetype=2) +
    facet_wrap(facets = 'site') +
    ggtitle(x) +
    xlab('Collection Date') +
    ylab('Carbon') +
    theme_grey(base_size = 16) +
    scale_x_date(date_labels="%b-%d",date_breaks  ="2 month") }


(filamentous_algae_table <- var_table(x = 'filamentous_algae'))
(filamentous_algae_cp <- carbon_plot(x = 'filamentous_algae'))
(filamentous_algae_np <- nitrogen_plot(x = 'filamentous_algae'))
  
(macrophyte_table <- var_table(x = 'macrophyte'))
(macrophyte_cp <- carbon_plot(x = 'macrophyte'))
(macrophyte_np <- nitrogen_plot(x = 'macrophyte'))

(periphyton_table <- var_table(x = 'periphyton'))
(periphyton_cp <- carbon_plot(x = 'periphyton'))
(periphyton_np <- nitrogen_plot(x = 'periphyton'))

(tree_leaf_table <- var_table(x = 'tree_leaf'))
(tree_leaf_cp <- carbon_plot(x = 'tree_leaf'))
(tree_leaf_np <- nitrogen_plot(x = 'tree_leaf'))

(grass_table <- var_table(x = 'grass'))
(grass_cp <- carbon_plot(x = 'grass'))
(grass_np <- nitrogen_plot(x = 'grass'))

(coarse_detritus_table <- var_table(x = 'coarse_detritus'))
(coarse_detritus_cp <- carbon_plot(x = 'coarse_detritus'))
(coarse_detritus_np <- nitrogen_plot(x = 'coarse_detritus'))

(mix_tv_table <- var_table(x = 'mix_tv'))
(mix_tv_cp <- carbon_plot(x = 'mix_tv'))
(mix_tv_np <- nitrogen_plot(x = 'mix_tv'))

# Biplot sources
src <- d %>%
  filter(type %in% c('Aquatic', 'Terrestrial')) %>%
  filter(species != 'clover')

src$species <- factor(src$species, levels = c('periphyton', 'filamentous_algae', 'macrophyte','coarse_detritus', 'mix_tv', 'tree_leaf', 'grass'))

biplot_allsources <- src %>%
  ggplot(aes(x = carbon, y = nitrogen, fill = species, shape = type)) +
  geom_point(size = 4, alpha = .5) +
  scale_shape_manual(values = c(21,24)) +
  scale_fill_viridis_d() +
  theme_classic(base_size = 16) +
  facet_wrap(facets=c('site')) +
  guides(fill = guide_legend(override.aes = list(shape = 21) ),
         shape = guide_legend(override.aes = list(fill = "black") ) )

biplot_limitedsources <- filter(src, species %in% c('periphyton', 'tree_leaf')) %>%
  ggplot(aes(x = carbon, y = nitrogen, fill = species, shape = type)) +
  geom_point(size = 4, alpha = .5) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = 23:24) +
  theme_classic(base_size = 16) +
  facet_wrap(facets=c('site')) +
  guides(fill = guide_legend(override.aes = list(shape = 21) ),
         shape = guide_legend(override.aes = list(fill = "black") ) )

# ANOVA section
# aquatic summary statistics
saq <- src %>%
  filter(species %in% c('filamentous_algae', 'periphyton'))      
saq$collection <- factor(saq$collection, levels = c('2020-01-15', '2020-05-15','2020-09-15')) 

(site_carbon_stats <- saq %>%
  group_by(site, collection) %>%
  get_summary_stats(carbon, type = "mean_sd"))

(site_nitrogen_stats <- saq %>%
    group_by(site, collection) %>%
    get_summary_stats(nitrogen, type = "mean_sd"))

# aquatic visualization
saq_cbox <- ggboxplot(
  saq, x = "site", y = "carbon",
  color = "collection", pallete = 'jco'
)
saq_cbox

saq_nbox <- ggboxplot(
  saq, x = "site", y = "nitrogen",
  color = "collection", pallete = 'jco'
)
saq_nbox

# aquatic outliers
(saq_2way_carbon_outlier <- saq %>% 
  group_by(site, collection) %>%
  identify_outliers(carbon) )

(saq_2way_nitrogen_outlier <- saq %>% 
    group_by(site, collection) %>%
    identify_outliers(carbon) )

# aquatic Build the linear model
saq_c_model  <- lm(carbon ~ site*collection,
             data = saq)
# aquatic Create a QQ plot of residuals
ggqqplot(residuals(saq_c_model))

# check normality assumption by groups
saq %>%
  group_by(site, collection) %>%
  shapiro_test(carbon)
# aquatic insufficient sample size groups

ggqqplot(saq, "carbon", ggtheme = theme_bw()) +
  facet_grid(carbon ~ site)

# Aquatic ANOVA
saqc_res.aov <- saq %>% anova_test(carbon ~ site * collection)
saqc_res.aov

saqn_res.aov <- saq %>% anova_test(nitrogen ~ site * collection)
saqn_res.aov

# terrestrial summary statistics
ste <- src %>%
  filter(species %in% c('tree_leaf', 'mix_tv' ))      
ste$collection <- factor(ste$collection, levels = c('2020-01-15', '2020-05-15','2020-09-15')) 

(site_carbon_stats <- ste %>%
    group_by(site, collection) %>%
    get_summary_stats(carbon, type = "mean_sd"))

(site_nitrogen_stats <- ste %>%
    group_by(site, collection) %>%
    get_summary_stats(nitrogen, type = "mean_sd"))

# terrestrial visualization
ste_cbox <- ggboxplot(
  ste, x = "site", y = "carbon",
  color = "collection", pallete = 'jco'
)
ste_cbox

ste_nbox <- ggboxplot(
  ste, x = "site", y = "nitrogen",
  color = "collection", pallete = 'jco'
)
ste_nbox

# terrestrial outliers
(ste_2way_carbon_outlier <- ste %>% 
    group_by(site, collection) %>%
    identify_outliers(carbon) )

(ste_2way_nitrogen_outlier <- ste %>% 
    group_by(site, collection) %>%
    identify_outliers(carbon) )

# remove extreme outliers
bad1 <- which(ste$site == 'GC' & 
      ste$collection == '2020-01-15' & 
      ste$species == 'mix_tv' & 
      near(ste$my_rep, 2)) # row = 2
bad2 <- which(ste$site == 'GC' & 
        ste$collection == '2020-01-15' & 
        ste$species == 'mix_tv' & 
        near(ste$my_rep, 5)) # row = 5
ste <- ste[-c(bad1, bad2), ]

# Aquatic ANOVA
stec_res.aov <- ste %>% anova_test(carbon ~ site * collection)
stec_res.aov

sten_res.aov <- ste %>% anova_test(nitrogen ~ site * collection)
sten_res.aov
























