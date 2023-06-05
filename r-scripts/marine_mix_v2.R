# mixing models: marine derived nutrients
# guilds (resident vs migrants: fish and invertebrates)
# Sean Kinard

# load packages
library(tidyverse)
library(simmr)

# Set Work Directory

# linux
# setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

# mac
setwd('/Users/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

d <- read_csv("data/MDN_clean_final.csv")
ds <- read_csv("data/theoretical_source.csv") %>%
  pivot_wider(names_from = ELEMENT,
              values_from = PERMIL)

# carbon and sulfur are suitable for bayesian mixing model (see marine_v1 figures)

# remove rows with na for either carbon or sulfur
d <- d %>%
  filter(! is.na(CARBON)) %>%
  filter(! is.na(SULFUR))

# order sites 
d$SITE = factor(d$SITE, levels = c('BB', 'CB', 'HB', 'LB', 'NB', 'LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

# plots
ggplot(d, aes(x = CARBON, color = GUILD, y = SULFUR)) +
  facet_wrap(facets = 'SITE') +
  geom_point() +
  ggtitle('34 S vs 13 C for flora and fauna')

d %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  filter(SITE_TYPE == 'Fresh') %>%
  ggplot( aes(x = CARBON, color = MIGRANT, y = SULFUR)) +
  facet_wrap(facets = 'SITE') +
  geom_point() +
  ggtitle('34 S vs 13 C for stream fauna')

d %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  filter(SITE %in% c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM')) %>%
  ggplot(aes(x = CARBON, fill = PPTAVG_SITE, y = SULFUR)) +
  stat_ellipse(aes(color = as.factor(PPTAVG_SITE)), show.legend=FALSE) +
  scale_color_viridis_d(direction=-1) +
  scale_fill_viridis_c(trans = 'reverse') +
  geom_point(shape = 21, size = 2, alpha = .8) +
  ggtitle('34 S vs 13 C for stream fauna')

###########################
#  ad hoc SIMMR functions #
###########################
# prepares data for input to simmr
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
# Run SIMMR resident fauna #
###########################

# Prep data for mixing model
s_resident <- d %>%
  filter(GUILD %in% c('Fish', 'Invertebrate')) %>%
  filter(MIGRANT == 'Potamodromous') %>%
  filter(SITE_TYPE == 'Fresh')

# prep data
prepared <- model_prep(s_resident, team = 'SITE', ds)
# final prep and check for mixing model input        
simmr_in =  simmr_load(mixtures = prepared[[1]],
                       group = prepared[[2]],
                       source_names = prepared[[3]],
                       source_means = prepared[[4]],
                       source_sds = prepared[[5]] )
# Run mixing model
s_out = simmr_mcmc(simmr_in)
# AR stats
s_stat <- mix_stat(s_out)
# extract group, mean, sd, 2.5%. 97.5%
s_slim <- s_stat %>%
  select(Estuary, m_group, statistic) %>%
  pivot_wider(values_from = Estuary, names_from = statistic) %>%
  mutate(lower_CI = mean - sd,
         upper_CI = mean + sd) %>%
  rename(SITE = m_group)
s_slim$SITE <- factor(s_slim$SITE, levels = c('LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))
# plot  
( mix_residents_site <- s_slim %>%
    ggplot(aes(x = SITE, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `25%`, ymax= `75%`)) +
  ggtitle(expression(paste('% Estuary Signature In Potamodromous Fauna Using ', sigma,'S'^'34'))) +
    ylab('% Estuary') +
    xlab('Location'))

# Merge with site info
mster <- d %>%
  filter(SITE_TYPE == 'Fresh') %>%
  select(SITE, BAYDIST_KM, ELEV_SITE_M, NPDES_MAJ_DENS, RAW_DIS_NEAREST_MAJ_NPDES, RAW_DIS_NEAREST_DAM, PPTAVG_SITE, HIRES_LENTIC_DENS) %>%
  unique() %>%
  mutate(ELEV_SITE_M = ifelse(SITE == 'TR', 18,
                              ifelse(SITE == 'UN', 3,
                                     ifelse(SITE == 'LN', 2, ELEV_SITE_M)))) %>%
  full_join(s_slim)

# plot elevation, colored by site
( mix_residents_Elevation_1 <- mster %>%
    ggplot(aes(x = ELEV_SITE_M, y = mean, fill=PPTAVG_SITE)) +
    geom_errorbar(aes(ymin = `25%`, ymax= `75%`), size = 1) +
    geom_point(shape = 21,
               size = 4) +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(expression(paste('% Estuary Signature In Potamodromous Fauna Using ', sigma,'S'^'34'))) +
    ylab('% Estuary') +
    xlab('Elevation (m)') +
    theme_bw(base_size = 20) +
    theme(legend.position = c(.82, .9),
          legend.direction = 'horizontal',
          legend.key.width= unit(1, 'cm'),
          legend.background = element_rect(colour = 'grey', 
                                           fill = 'white',
                                           linetype='solid')) +
    labs(fill = 'Rainfall\n(cm/yr)') )

# elevation, colored site, regression
( mix_residents_Elevation_2 <- mix_residents_Elevation_1 +
    geom_smooth(method = 'lm', color = 'red', linetype = 2, se = F) +
    stat_ma_eq(size = 8,
               label.x=.6,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 8,
                     label.x=.6,
                     label.y = .85,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    geom_errorbar(aes(ymin = `25%`, ymax= `75%`), size = 1) +
    geom_point(shape = 21,
               size = 4) )

( mix_residents_PPT_1 <- mster %>%
    ggplot(aes(x = PPTAVG_SITE, y = mean)) +
    geom_errorbar(aes(ymin = `25%`, ymax= `75%`), size = 1) +
    geom_point(shape = 21,
               size = 4,
               fill = 'grey') +
    ggtitle(expression(paste('% Estuary Signature In Potamodromous Fauna Using ', sigma,'S'^'34'))) +
    ylab('% Estuary') +
    xlab('Annual Rainfall (cm)') +
    theme_bw(base_size = 20) )

( mix_residents_PPT_2 <- mix_residents_PPT_1 +
    geom_smooth(method = 'lm', color = 'red', linetype = 2, se = FALSE) +
    stat_ma_eq(size = 8,
               label.x=.6,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 8,
                     label.x=.6,
                     label.y = .85,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    geom_errorbar(aes(ymin = `25%`, ymax= `75%`), size = 1) +
    geom_point(shape = 21,
               size = 4, fill = 'grey') )

( Sulfur_Elevation <- d %>%
  filter(SITE_TYPE == 'Fresh') %>%
  filter(MIGRANT == 'Potamodromous') %>%
  ggplot(aes(x=ELEV_SITE_M, y = SULFUR, fill = PPTAVG_SITE)) +
  geom_smooth(method = 'lm', se = F, color = 'red') +
  geom_jitter(size = 4, shape = 21) +
  scale_fill_viridis_c(direction = -1) +
  ggtitle(expression(paste(sigma,'S'^'34',' In Potamodromous Fauna'))) +
  ylab('% Estuary') +
  xlab('Elevation (m)') +
  theme_bw(base_size = 20) +
  theme(legend.position = c(.5, .1),
        legend.direction = 'horizontal',
        legend.key.width= unit(1, 'cm'),
        legend.background = element_rect(colour = 'grey', 
                                         fill = 'white',
                                         linetype='solid')) +
  labs(fill = 'Rainfall\n(cm/yr)') +
    stat_ma_eq(size = 8,
               label.x=.75,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 8,
                     label.x=.75,
                     label.y=.85,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*')))  )
    
( Sulfur_Rainfall <- d %>%
  filter(SITE_TYPE == 'Fresh') %>%
  filter(MIGRANT == 'Potamodromous') %>%
  ggplot(aes(x=PPTAVG_SITE, y = SULFUR)) +
    stat_ma_eq(size = 8,
               label.x=.12,
               label.y=.1,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 8,
                     label.y = .1,
                     label.x=.42,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
  geom_smooth(method = 'lm',color = 'red', linetype = 1, se = F) +
  geom_jitter(size = 4, shape = 21, fill = 'grey') +
  scale_fill_viridis_c(direction = -1) +
  ggtitle(expression(paste('Estuary Influence (',sigma,'S'^'34', ') In Potamodromous Fauna'))) +
  ylab(expression(paste(sigma,'S'^'34'))) +
  xlab('Annual Rainfall (cm)') +
  theme_bw(base_size = 20) +
  theme(legend.position = c(.5, .1),
        legend.direction = 'horizontal',
        legend.key.width= unit(1, 'cm'),
        legend.background = element_rect(colour = 'grey', 
                                         fill = 'white',
                                         linetype='solid')) )
  
( Sulfur_baydist <- d %>%
    filter(SITE_TYPE == 'Fresh') %>%
    filter(MIGRANT == 'Potamodromous') %>%
    ggplot(aes(x=BAYDIST_KM, y = SULFUR)) +
    stat_ma_eq(size = 8,
               label.x=.15,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(size = 8,
                     label.x=.4,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    geom_smooth(method = 'lm',color = 'red', linetype = 3, se = F) +
    geom_jitter(size = 4, shape = 21, fill = 'grey') +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(expression(paste('Estuary Influence (',sigma,'S'^'34', ') In Potamodromous Fauna'))) +
    ylab(expression(paste(sigma,'S'^'34'))) +
    xlab('Distance To Estuary (km)') +
    theme_bw(base_size = 20) +
    theme(legend.position = c(.5, .1),
          legend.direction = 'horizontal',
          legend.key.width= unit(1, 'cm'),
          legend.background = element_rect(colour = 'grey', 
                                           fill = 'white',
                                           linetype='solid')) )  

ggsave('Figures/mix_residents_Elevation_1.png',
       mix_residents_Elevation_1,
       'png')

ggsave('Figures/mix_residents_Elevation_2.png',
       mix_residents_Elevation_2,
       'png')
ggsave('Figures/Sulfur_Elevation.png',
       Sulfur_Elevation,
       'png')
ggsave('Figures/Sulfur_Rainfall.png',
       Sulfur_Rainfall,
       'png')
ggsave('Figures/Sulfur_baydist.png',
       Sulfur_baydist,
       'png')






# check for algorithm convergence
summary(s_out, type = 'diagnostics') # values should not exceed 1.1
# Check model fit
# points should lie in the fitted value intervals (default 50%)
# post_pred = posterior_predictive(simmr_out, group = 1) 
# print(post_pred)
# unsolvable error "arguments imply differing number of rows". didn't work for groups 1:7
# plot prior and posteriors
prior_viz(s_out)
# simmr plots
# plot(simmr_out,type = c('density','matrix'),group = 2,title = 'F-2')
# plot(simmr_out, type = 'boxplot',group = 2,title = 'F-2')
# compare_groups(s_out, source = 'Aquatic Source', groups = 1:length(levels(as.factor(prepared[[2]]))))
# cleanup universals
remove(prepared)
remove(simmr_in)











