setwd('/Users/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined')

d <- read_csv("data/MDN_clean_final.csv")
ds <- read_csv("data/theoretical_source.csv") %>%
  pivot_wider(names_from = ELEMENT,
              values_from = PERMIL)
# Theoretical Source (E-Algae vs [F-riparian + F-detritus])
ds2 <- ds %>%
  filter(SITE_TYPE == 'Estuary' & TAXON == 'Algae') %>%
  full_join(filter(ds, SITE_TYPE == 'Fresh' & TAXON %in% c('Coarse Detritus', 'Riparian')))

# order sites 
d$SITE = factor(d$SITE, levels = c('BB', 'CB', 'HB', 'LB', 'NB', 'LN', 'UN', 'TR', 'SF', 'AR', 'PD', 'MR', 'PL', 'GC', 'WM', 'EM'))

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

# Few species have sufficient samples across sites
my_sp2 <- c('Palaemonidae', 'G.affinis', 'L.megalotis', 'H.cyanoguttatus', 'L.macrochirus', 'P.latipinna')

# Generate Plots
plotter <- function(d, var) {
  
  d %>%
    ggplot(aes(x = {{var}},
               y = PERMIL, 
               fill = ELEMENT)) +
    geom_point(shape = 21,
               size = 4,
               show.legend = T) +
    geom_smooth(method = 'lm',
                show.legend = T) +
    stat_ma_eq(data = filter(d, ELEMENT == 'SULFUR'),
               size = 4,
               label.x=.1,
               label.y = .06,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    stat_correlation(data = filter(d, ELEMENT == 'SULFUR'),
                     size = 4,
                     label.x=.1,
                     label.y = .01,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    stat_ma_eq(data = filter(d, ELEMENT == 'CARBON'),
               size = 4,
               label.x=.1,
               label.y = .99,
               aes(label =  paste(after_stat(eq.label),
                                  sep = "~~italic(\"with\")~~"))) +
    
    stat_correlation(data = filter(d, ELEMENT == 'CARBON'),
                     size = 4,
                     label.x=.1,
                     label.y = .95,
                     mapping = aes(label = paste(after_stat(cor.label),
                                                 after_stat(p.value.label),
                                                 sep = '*"; "*'))) +
    ggtitle(expression(paste('Scaled ',sigma,'S'^'34', 'And ',sigma,'C'^'13', ' In Resident Freshwater Fauna'))) +
    ylab(expression(paste('Scaled ',sigma))) +
    theme_classic2(base_size = 20) }


my_d <- filter(dlong_scaled, GUILD == 'Fish') %>%
  filter(ELEMENT != 'NITROGEN')

( sc_rain <- plotter(my_d, PPTAVG_SITE) +
  xlab('Annual Rainfall (cm)') )
# while sulfur signatures increase with rainfall, carbon signatures decrease. These disparate trends may result from decreasing periphyton influence within freshwater communities as rainfall increases. But this still seems illogical given that carbon and sulfur covary according to the plot below. Granted, the autocorrelation has an R of .4 and an R^2 of .16, so the observed differnces within taxa might point towards some other factor.

# Scaled Sulfur vs Carbon: Autocorrelated
( svc <- dlong_scaled %>%
  pivot_wider(names_from = ELEMENT,
              values_from = PERMIL) %>%
    filter(GUILD == 'Fish') %>%
  ggplot(aes(x=CARBON, y = SULFUR, fill = PPTAVG_SITE)) +
  geom_point(shape = 21,
             size = 4,
             show.legend = F) +
  geom_smooth(method = 'lm',
              show.legend = T) +
  stat_correlation(size = 6,
                   label.x=.1,
                   label.y = .94,
                   mapping = aes(label = paste(after_stat(cor.label),
                                               after_stat(p.value.label),
                                               sep = '*"; "*'))) +
  stat_ma_eq(size = 6,
             label.x=.1,
             label.y = .99,
             aes(label =  paste(after_stat(eq.label),
                                sep = "~~italic(\"with\")~~"))) +
  ggtitle('Scaled Sulfur vs Scaled Carbon In Fish') )


ggsave('Figures/svc.png',
       svc,
       'png')
ggsave('Figures/sc_rain/png',
       sc_rain,
       'png')
