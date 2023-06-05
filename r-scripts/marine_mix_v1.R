# mixing models: marine derived nutrients
# guilds (resident vs migrants: fish and invertebrates)
# Sean Kinard

# load packages
library(tidyverse)
library(simmr)

# set work directory
setwd('/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

# load data
d <- read_csv('data/Marine_v1.csv')

# carbon and sulfur are suitable for bayesian mixing model (see marine_v1 figures)

# remove rows with na for either carbon or sulfur
d <- d %>%
  filter(! is.na(carbon)) %>%
           filter(! is.na(sulfur))

# order sites 
d$site_code = factor(d$site_code, levels = c('BB', 'CB', 'HB', 'LB', 'NB', 'LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

# plots
ggplot(d, aes(x = carbon, color = type, y = sulfur)) +
  facet_wrap(facets = 'site_code') +
  geom_point() +
  ggtitle('34 S vs 13 C for flora and fauna')

d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(site_type == 'stream') %>%
  ggplot( aes(x = carbon, color = migrant, y = sulfur)) +
  facet_wrap(facets = 'site_code') +
  geom_point() +
  ggtitle('34 S vs 13 C for stream fauna')

d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(site_code %in% c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM')) %>%
  ggplot(aes(x = carbon, fill = PPTAVG_SITE, y = sulfur)) +
  stat_ellipse(aes(color = as.factor(PPTAVG_SITE)), show.legend=FALSE) +
  scale_color_viridis_d(direction=-1) +
  scale_fill_viridis_c(trans = 'reverse') +
  geom_point(shape = 21, size = 2, alpha = .8) +
  ggtitle('34 S vs 13 C for stream fauna')

d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(site_code %in% c('LN', 'UN')) %>%
  mutate(site_code = ifelse(site_code == 'LN' , 'Below', 'Above')) %>%
  ggplot(aes(x = carbon, fill = site_code, y = sulfur)) +
  geom_point(shape = 21, size = 3, alpha = .8) +
  stat_ellipse(linetype = 2, aes(color = site_code), show.legend=FALSE) +
  ggtitle('34 S vs 13 C for fauna above and below the Nueces River dam')

d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(site_code %in% c('LN', 'UN')) %>%
  mutate(site_code = ifelse(site_code == 'LN' , 'Below', 'Above')) %>%
  ggplot(aes(x = carbon, fill = site_code, y = sulfur)) +
  geom_point(shape = 21, size = 3, alpha = .8) +
  geom_line(aes(group = species)) +
  stat_ellipse(linetype = 2, aes(color = site_code), show.legend=FALSE) +
  ggtitle('34 S vs 13 C for fauna above and below the Nueces River dam: line-species')

d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(site_code %in% c('LN', 'UN')) %>%
  mutate(site_code = ifelse(site_code == 'LN' , 'Below', 'Above')) %>%
  ggplot(aes(x = carbon, fill = site_code, y = sulfur)) +
  geom_point(shape = 21, size = 3, alpha = .8) +
  geom_line(aes(group = species)) +
  geom_label(label = d %>%
               filter(type %in% c('Fish', 'Invertebrate')) %>%
               filter(site_code %in% c('LN', 'UN')) %>%
               pull(species)) +
  stat_ellipse(linetype = 2, aes(color = site_code), show.legend=FALSE) +
  ggtitle('34 S vs 13 C for fauna above and below the Nueces River dam: line-species')

ggplot(d, aes(x = site_code, fill=PPTAVG_SITE, y = sulfur)) +
  geom_boxplot() +
  scale_fill_viridis_c(direction = -1)+
  ggtitle('Sulfur ratios for flora and fauna')

ggplot(d, aes(x = site_code, fill=PPTAVG_SITE, y = carbon)) +
  geom_boxplot() +
  scale_fill_viridis_c(direction = -1)+
  ggtitle('Carbon ratios for flora and fauna')


###########################
#  ad hoc SIMMR functions #
###########################
# prepares data for input to simmr
model_prep <- function(data, team) {
  
  mix <- data %>%
    filter(type == 'Fish' |
             type == 'Invertebrate') %>%
    select(carbon, sulfur) %>%
    as.matrix()
  
  mix_grp <- data %>%
    filter(type == 'Fish' |
             type == 'Invertebrate') %>%
    pull(team) %>%
    as.matrix()
  
  marine <- data %>% filter(site_type == 'bay' & type == 'Aquatic')
  ter <- data %>% filter(site_code %in% c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM') & type %in% c('Detritus', 'Terrestrial'))
  s_data <- full_join(marine, ter) %>%
    mutate(source_type = ifelse(site_type == 'bay' & type == 'Aquatic',
                                'Marine', 'Terrestrial'))
  
  s_names <- c( "Marine", "Terrestrial")
  
  s_means <- aggregate(x= select(s_data, 
                                 carbon, 
                                 sulfur),
                       by= list(pull(s_data, source_type)),
                       FUN=mean) %>%
    select(carbon, sulfur) %>%
    as.matrix()
  
  s_sds <- aggregate(x= select(s_data, carbon, sulfur),
                     by= list(pull(s_data, source_type)),
                     FUN=sd) %>%
    select(carbon,sulfur) %>%
    as.matrix()
  
  my_list <- list(mix, mix_grp,  s_names, s_means, s_sds, s_data)
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

# Prep data for mixing model of fauna vs site
s_resident <- d %>%
  filter(type %in% c('Fish', 'Invertebrate')) %>%
  filter(migrant == 'freshwater') %>%
  filter(site_type == 'stream')

marine <- d %>% filter(site_type == 'bay' & type == 'Aquatic')
ter <- d %>% filter(site_code %in% c('TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM') & type %in% c('Detritus', 'Terrestrial'))

combo <- full_join(s_resident, marine) %>% full_join(ter)
# prep data
prepared <- model_prep(data = combo, team = 'site_code')
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
  select(Marine, m_group, statistic) %>%
  pivot_wider(values_from = Marine, names_from = statistic) %>%
  mutate(lower_CI = mean - sd,
         upper_CI = mean + sd)
s_slim$m_group <- factor(s_slim$m_group, levels = c('LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))
# plot  
mix_residents <- ggplot(s_slim, aes(x = m_group, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_CI, ymax= upper_CI)) +
  ggtitle('% Marine-derived nutrients in resident fauna')
mix_residents 
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

# integrate site data from d
mega <- d %>% select(site_code, site_lat:ELEV_SITE_M) %>% 
  unique() %>%
  right_join(s_slim, by = c('site_code' = 'm_group'))

pairs(mega %>% select(mean, PPTAVG_SITE:ELEV_SITE_M))

ggplot(mega, aes(PPTAVG_SITE, mean)) +
  geom_point(size = 4) +
  geom_errorbar(width = .5, size = .8, aes(ymin = lower_CI, ymax = upper_CI)) +
  geom_smooth(method = 'lm', se = FALSE)

summary(lm(mean ~ PPTAVG_SITE, mega))

ggplot(mega, aes(ELEV_SITE_M, mean)) +
  geom_point(size = 4) +
  geom_errorbar(width = .5, size = .8, aes(ymin = lower_CI, ymax = upper_CI)) +
  geom_smooth(method = 'lm', se = FALSE)

summary(lm(mean ~ ELEV_SITE_M, mega))


# import streamer data
streamer <- read_csv('data/streamer_data.csv')

mega <- streamer %>%
  select(site_code, position_between, over_under, trace_river_mile) %>%
  pivot_wider(names_from = 'over_under', 
              values_from = 'trace_river_mile') %>%
  mutate(river_mile = under + position_between*(over-under)) %>%
  select(site_code, river_mile) %>%
  right_join(mega)

ggplot(mega, aes(river_mile, mean)) +
  geom_point(size = 4) +
  geom_errorbar(width = .5, size = .8, aes(ymin = lower_CI, ymax = upper_CI)) +
  geom_smooth(method = 'lm', se = FALSE)

summary(lm(mean ~ river_mile, mega))

mega <- mega[,c(1,2,9:15)]
mega <- scale(mega[,2:8])


full.model <- lm(mean ~ river_mile + PPTAVG_SITE + HIRES_LENTIC_DENS + NPDES_MAJ_DENS + RAW_DIS_NEAREST_MAJ_NPDES + RAW_DIS_NEAREST_DAM , ELEV_SITE_M, data=as.data.frame(mega) )
summary(full.model)










