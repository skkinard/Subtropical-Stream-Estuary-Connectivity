setwd('C:\\Users\\s2kin\\Dropbox\\Research\\Manuscipts\\Diss.5_Isotope')

cn_sum19 <- read_csv('RAPID_sum_2019\\Data\\master_isotope.csv')
cn_sum20 <- read_csv('TERRG_sum_CN_22-36\\TERRG_CN_Cyc_03_22-36.csv')
cn_win20 <- read_csv('TERRG_win_CN_1-21\\TERRG_CN_01-21.csv')

my_table <- cn_sum19 %>%
  filter(guild == "Terrestrial primary producer" |
         guild == "Aquatic primary producer") %>%
  select(site, guild, species, carbon, nitrogen) %>%
  mutate(season = 'summer-2019') %>%
  rename(type = guild)

my_table$type <- str_replace(my_table$type, 'primary producer', 'Source')

my_table <- cn_sum20 %>%
  filter(type == "Terrestrial Source" |
           type == "Aquatic Source") %>%
  filter(site == 'AR' | site == 'GC' | site == 'SFC') %>%
  select(site, type, species, carbon, nitrogen) %>%
  mutate(season = 'summer-2020') %>%
  full_join(my_table)

my_table <- cn_win20 %>%
  filter(type == "Terrestrial Source" |
           type == "Aquatic Source") %>%
  filter(site == 'AR' | site == 'GC' | site == 'SFC') %>%
  rename(carbon = "d13C (permil, vs VPDB)", nitrogen = "d15N (permil, vs AIR)") %>%
  select(site, type, species, carbon, nitrogen) %>%
  mutate(season = 'winter-2020') %>%
  full_join(my_table)

my_table <- my_table %>%
    mutate(Rainfall = ifelse(site == 'SFC', '55 cm/yr',
                           ifelse(site =='AR', '70 cm/yr',
                                  ifelse(site == 'GC', '80 cm/yr', NA))) )
my_table$season <- factor(my_table$season, levels = c('summer-2019', 'winter-2020', 'summer-2020'))
                           
                           
my_table %>%
ggplot(aes(x=as.factor(season), y = nitrogen, fill = type)) +
  facet_wrap('Rainfall') +
  geom_violin() +
  geom_point(alpha = .5, position = position_dodge(width = .9), size = 2) +
  theme_bw(base_size = 24) +
  theme(legend.position = "bottom") +
  theme(legend.direction="horizontal") +
  theme(axis.title.x=element_blank()) +
  theme(legend.title = element_blank()) +
  ylab(expression(paste(sigma, 'N'^'15'))) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) 

my_table %>%
  ggplot(aes(x=as.factor(season), y = nitrogen, fill = type)) +
  facet_wrap('Rainfall') +
  geom_violin() +
  geom_point(alpha = .5, position = position_dodge(width = .9)) +
  theme(legend.position = "bottom") +
  theme(legend.direction="horizontal") +
  theme(legend.title = element_blank()) +
  xlab('Rainfall') +
  ylab(expression(paste(sigma, 'N'^'15')))


