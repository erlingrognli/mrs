d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         lh_vPFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
         rh_vPFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area,
         .keep = 'unused') # %>%

d_s <- select(d, c('ocd', 'gender', caudate, occipital_area, 'lh_vPFC_area', 'rh_vPFC_area'))
