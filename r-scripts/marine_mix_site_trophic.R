# mixing models: marine derived nutrients in trophic levels at each site
# Residents by trophic level
# Sean Kinard

# load packages
library(tidyverse)
library(simmr)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggpmisc)

# Set Work Directory

# linux
# setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

# mac
setwd('/Users/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

d <- read_csv("data/MDN_clean_final.csv")
ds <- read_csv("data/theoretical_source.csv") %>%
  pivot_wider(names_from = ELEMENT,
              values_from = PERMIL)
# Theoretical Source (E-Algae vs [F-riparian + F-detritus])
ds2 <- ds %>%
  filter(SITE_TYPE == 'Estuary' & TAXON == 'Algae') %>%
  full_join(filter(ds, SITE_TYPE == 'Fresh' & TAXON %in% c('Coarse Detritus', 'Riparian')))

# remove rows with na for either carbon or sulfur
d <- d %>%
  filter(! is.na(CARBON)) %>%
  filter(! is.na(SULFUR))

# order sites 
d$SITE = factor(d$SITE, levels = c('BB', 'CB', 'HB', 'LB', 'NB', 'LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

###########################
#  ad hoc SIMMR functions #
###########################
# prepares data for input to simmr
# Using: Sulfur
model_prep <- function(mix_data, team, src_data) {
  
  mix <- mix_data %>%
    select(SULFUR) %>%
    as.matrix()
  
  mix_grp <- mix_data %>%
    pull(team) %>%
    as.matrix()
  
  s_names <- c( "Estuary", "Fresh")
  
  s_means <- aggregate(x= select(src_data,
                                 SULFUR),
                       by= list(pull(src_data, SITE_TYPE)),
                       FUN=mean) %>%
    select(SULFUR) %>%
    as.matrix()
  
  s_sds <- aggregate(x= select(src_data,
                               SULFUR),
                     by= list(pull(src_data, SITE_TYPE)),
                     FUN=sd) %>%
    select(SULFUR) %>%
    as.matrix()
  
  my_list <- list(mix, mix_grp, s_names, s_means, s_sds, src_data)
  
  print(my_list) }

# extracts and aggregates mixing model statistical outputs.
mix_stat <- function(d) { # d=simmr_out
  
  # extract statistics and quantile statsitics for each group
  sumstat <- summary(d, 
                     type = c('statistics', 'quantiles'),
                     group = 1:length(levels(as.factor(prepared[[2]]))))
  # statistics
  my_stats <- as.data.frame(do.call(cbind, sumstat$statistics))
  # create vector of group names
  my_stats_groups <- sort(rep(levels(d$input$group),2))
  # rename columns to include groups
  colnames(my_stats) <- paste(my_stats_groups, 
                              colnames(my_stats), 
                              sep="_")
  # transpose dataframe
  my_stats <- t(my_stats)
  # add group and statistic columns
  my_stats <- mutate(as_tibble(my_stats), 
                     m_group = my_stats_groups,
                     statistic = rep(c('mean', 'sd'), 
                                     length(levels(d$input$group))) )
  # quantiles
  my_q <- as.data.frame(do.call(cbind, sumstat$quantiles))
  # create vector of group names
  my_q_groups <- sort(rep(levels(d$input$group),5))
  # rename columns to include groups
  colnames(my_q) <- paste(my_q_groups, 
                          colnames(my_q), 
                          sep="_")
  # transpose dataframe
  my_q <- t(my_q)
  # add group, site, and statistic columns
  my_q <- mutate(as_tibble(my_q), 
                 m_group = as.factor(my_q_groups),
                 statistic = rep(c('2.5%', 
                                   '25%', 
                                   '50%', 
                                   '75%', 
                                   '97.5%'), 
                                 length(levels(d$input$group))) )
  
  mix_output <- full_join(my_stats, my_q)
  return(mix_output) }

###########################
# Run SIMMR For Each Site #
###########################

# Assign sites
my_sites <- d %>% filter(SITE_TYPE == 'Fresh') %>% pull(SITE) %>% unique()

# Assign data, mixture groups, sources, and origin of interest
df <- d
grp <- 'TROPHIC_LEVEL'
src <- ds2
orgn <- 'Estuary'

# Create start tibble
my_table <- tibble( grp = df %>%
                      pull(grp)%>%
                      unique()) %>%
  mutate_at(1, as.character)
colnames(my_table) <- grp

# Loop: simmr grps at each site & merge stat outputs

for(var in my_sites) {
  
  # Filter data
  f_df <- df %>%
    filter(SITE == var) %>%
    filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
    filter(MIGRANT == 'Potamodromous')
  
  # Extract necessary data for simmr
  prepared <- model_prep(f_df, team = grp, src)
  # format Extracted data for simmr
  f_in =  simmr_load(mixtures = prepared[[1]],
                     group = prepared[[2]],
                     source_names = prepared[[3]],
                     source_means = prepared[[4]],
                     source_sds = prepared[[5]] )
  # Run mixing model (simmr)
  f_out = simmr_mcmc(f_in)
  # Extract stats from simmr output
  f_stat <- mix_stat(f_out)
  # extract group, mean, sd, 2.5%. 97.5%
  f_slim <- f_stat %>%
    select(all_of(orgn), m_group, statistic) %>%
    pivot_wider(values_from = all_of(orgn), 
                names_from = statistic) %>%
    mutate(SITE = rep(var))
  # Rename columns for easier recall        
  colnames(f_slim) <- str_to_upper(c(grp, 
                                     orgn,
                                     'sd',
                                     'CI_2.5', 
                                     'CI_25', 
                                     'CI_50', 
                                     'CI_75', 
                                     'CI_97.5',
                                     'SITE') )
  # Merge to previous
  my_table <- full_join(f_slim, my_table)
}


# Merge table to Environmental covariates
mix_site_trophic <- left_join(my_table, df%>%select(SITE, BAYDIST_KM:ELEV_SITE_M)%>%unique()) %>%
  filter(TROPHIC_LEVEL != '0')

# Generate Plots
plotter <- function(d, var) {
  
  d %>%
    ggplot(aes(x = {{var}},
               y = ESTUARY, 
               fill = TROPHIC_LEVEL)) +
    geom_errorbar(aes(ymin=CI_25,ymax=CI_75)) +
    geom_point(shape = 21,
               size = 4,
               show.legend = F) +
    geom_smooth(method = 'lm',
                show.legend = F) +
    stat_ma_eq(size = 6,
               label.x=.5,
               label.y = .99,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 6,
                     label.x=.5,
                     label.y = .93,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    facet_wrap(~TROPHIC_LEVEL) +
    ggtitle(expression(paste('% Estuary Signature In Potamodromous Fauna Using ', sigma,'S'^'34'))) +
    ylab('% Estuary') +
    theme_bw(base_size = 20) }


# Plots
( troph_PPTAVG_SITE <- plotter(mix_site_trophic, PPTAVG_SITE) +
    xlab('Rainfall (cm/yr)') )

( troph_BAYDIST_KM <- plotter(mix_site_trophic, BAYDIST_KM) +
    xlab('Estuary Distance (km)') )

( troph_ELEV_SITE_M <- plotter(mix_site_trophic, ELEV_SITE_M) +
    xlab('Elevation (m)') )

( troph_HIRES_LENTIC_DENS <- plotter(mix_site_trophic, HIRES_LENTIC_DENS) +
    xlab('Lentic Density') )

( troph_NPDES_MAJ_DENS <- plotter(mix_site_trophic, NPDES_MAJ_DENS) +
    xlab('Wastewater Density') )

( troph_RAW_DIS_NEAREST_MAJ_NPDES <- plotter(mix_site_trophic, RAW_DIS_NEAREST_MAJ_NPDES) +
    xlab('Wastewater Distance') )

( troph_RAW_DIS_NEAREST_DAM <- plotter(mix_site_trophic, RAW_DIS_NEAREST_DAM) +
    xlab('Dam Distance') )

# Trophic level 2 displays the strongest correlations: meaning that Invertivores display the greatest variation in Estuary influence among sites.

# Export plots
ggsave('Figures/troph_PPTAVG_SITE.png',
       troph_PPTAVG_SITE,
       'png')
ggsave('Figures/troph_BAYDIST_KM.png',
       troph_BAYDIST_KM,
       'png')
ggsave('Figures/troph_ELEV_SITE_M.png',
       troph_ELEV_SITE_M,
       'png')
ggsave('Figures/troph_HIRES_LENTIC_DENS.png',
       troph_HIRES_LENTIC_DENS,
       'png')
ggsave('Figures/troph_NPDES_MAJ_DENS.png',
       troph_NPDES_MAJ_DENS,
       'png')
ggsave('Figures/troph_RAW_DIS_NEAREST_MAJ_NPDES.png',
       troph_RAW_DIS_NEAREST_MAJ_NPDES,
       'png')
ggsave('Figures/troph_RAW_DIS_NEAREST_DAM.png',
       troph_RAW_DIS_NEAREST_DAM,
       'png')
