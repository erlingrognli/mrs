library(tidyverse); library(synthpop)

d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         lh_vPFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
         rh_vPFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area,
         .keep = 'unused') %>%

  select(all_of(c('age', 'ocd', 'gender', 'Left.Caudate', 'Right.Caudate', 
              "lh_lateraloccipital_area", "rh_lateraloccipital_area", 
              'lh_vPFC_area', 'rh_vPFC_area'))) %>%
  
  mutate(ocd = as.factor(ocd), gender = as.factor(gender))

d_syn <- syn(d)

d <- d_syn$syn
