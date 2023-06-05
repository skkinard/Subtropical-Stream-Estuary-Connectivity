# Extract site information from gages ii for MDN analyses:

library(tidyverse)

getwd()
setwd("/home/sean/Documents/Research/Manuscipts/Diss.5_Isotope/RAPID-TERRG-Combined")

# import data
streamer <- read_csv('data/streamer_data.csv')
basinid <- read_csv('data/conterm_basinid.txt')
climate <- read_csv('data/conterm_climate.txt')
hydromod_dams <- read_csv('data/conterm_hydromod_dams.txt')
hydromod_other <- read_csv('data/conterm_hydromod_other.txt')
landscape_pat <- read_csv('data/conterm_landscape_pat.txt')
topo <- read_csv('data/conterm_topo.txt')

# merge conterm gauge data
mega <- full_join(basinid, climate) %>%
  full_join(hydromod_dams) %>%
  full_join(hydromod_other) %>%
  full_join(landscape_pat) %>%
  full_join(topo)


# extract STAIDS
my_STAIDS <- streamer %>%
  pull(STAID) %>%
  unique()

mega$STAID <- as.numeric(mega$STAID)

my_mega <- mega %>%
  filter(STAID %in% my_STAIDS)

d <- my_mega %>%
  select(STAID,
         STANAME,
         LAT_GAGE,
         LNG_GAGE,
         PPTAVG_SITE,
         HIRES_LENTIC_DENS,
         NPDES_MAJ_DENS,
         RAW_DIS_NEAREST_MAJ_NPDES,
         RAW_DIS_NEAREST_DAM,
         ELEV_SITE_M )

# note that values of -999 are given to sites with no discharge or dams in their watershed
# Replace -999 with NA
d <- d %>%
  mutate(RAW_DIS_NEAREST_MAJ_NPDES = ifelse(RAW_DIS_NEAREST_MAJ_NPDES < 0 , NA, RAW_DIS_NEAREST_MAJ_NPDES),
         RAW_DIS_NEAREST_DAM = ifelse(RAW_DIS_NEAREST_DAM < 0, NA, RAW_DIS_NEAREST_DAM))

# export environmental variables
write_csv(d, 'data/MDN_conterm.csv')





