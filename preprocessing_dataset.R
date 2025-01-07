# preprocessing of the dataset for analyses

library(tidyverse)

# define function for changing outliers to NA 
outlier_removal <- function(x, k){
  
  ub <- median(x) + k * mad(x)
  lb <- median(x) - k * mad(x)
  
  x[which(x>ub|x<lb)] <- NA
  
  return(x)}

# load demographic info on ocd participants

ocd_demog <- read.csv2('~/mrs_data/mrs_pasienter_demog.csv')

# remove excluded participants and irrelevant variables
# and fix "ok" text in subject ids.

ocd_demog <- select(filter(ocd_demog, Akseptert.deltagelse == 'J'), c(MR4K.Nr, F.dt, K.nn))

ocd_demog$MR4K.Nr <- gsub(pattern = " ok", replacement = '', x = ocd_demog$MR4K.Nr)

# recode and rename variables, add indicator variable for ocd  

ocd_demog <- mutate(ocd_demog, 
                    birth_date = as.Date(F.dt, format = '%d.%m.%Y'),
                    female = ifelse(K.nn == 'F', 1, 0),
                    SubjID = MR4K.Nr,
                    ocd = 1,
                    eos = 0,
                    .keep = 'none')

# load demographic info on healthy controls

hc_demog <- read.csv2('~/mrs_data/mrs_kontroller_demog.csv')

# remove irrelevant variables and recode/rename
# combine with ocd participants

hc_demog <- select(hc_demog, c(MR4K.Nr, F.dt, K.nn)) %>%
  
  mutate(SubjID = MR4K.Nr, 
         birth_date = as.Date(F.dt, format = '%d.%m.%Y'),
         female = ifelse(K.nn == 'J', 1, 0),
         ocd = 0,
         eos = 0,
         .keep = 'none')

demog <- rbind(hc_demog, ocd_demog)

# load mri dates, reformat, calculate age in days, 
# and convert to years with two decimals to comply with format of ytop data

demog <- right_join(demog, read.csv('~/mrs_data/studyDate.csv'), 
                    by = join_by(SubjID == SubjID)) %>%
  
  mutate(study_date = as.Date(as.character(StudyDate), format = '%Y%m%d'),
         age = round(as.double((study_date - birth_date)/365), digits = 2),
         .keep = 'unused')

# read in processed mri data, adding 

ocd_ytop <- read.csv('~/mrs_data/OCD_YTOP_data.csv')

ocd_ytop <- 
  
# filter out the rows of mri data containing ocd patients and hcs, 
# while removing age and sex columns with missing observations 
    
  select(filter(ocd_ytop, is.na(Age) == T), - c(Age, Sex))  %>%

# join filtered mri data with age and gender variables from demographics dataframe     

  full_join(select(demog, !study_date), 
            join_by(SubjID == SubjID))  %>%

# bind rows of mri data from ocd patients and hc with the ytop data, 
#  renaming age and gender columns to match
    
  bind_rows(
    mutate(
      filter(ocd_ytop, is.na(Age) == F),
      age = Age,
      female = ifelse(Sex == 'Female', 1, 0),
      ocd = 0,
      eos = ifelse(FullDx == 'EOS', 1, 0),
      .keep = 'unused')) %>%
  
  # remove EOP cases that are not EOS
  
  filter(!Dx%in%c('OTP','AFP')) %>%

  # remove brain measures that will not be used in the model
  
  select(!any_of(c("Dx", "FullDx", "Unknown", "Right.choroid.plexus", "Left.choroid.plexus", 
                   "Soft_Nonbrain_Tissue", "Fluid_Inside_Eyes", "CSF", "Skull", 
                   "X5th.Ventricle", "X4th.Ventricle",  "X3rd.Ventricle",
                   "Right.Inf.Lat.Vent", "Left.Inf.Lat.Vent",
                   #"Left.Lateral.Ventricle", "Right.Lateral.Ventricle",
                   "Optic.Chiasm", "Right.vessel", "Left.vessel",
                   "WM.hypointensities", "non.WM.hypointensities",
                   "Left.Cerebral.Cortex", "Right.Cerebral.Cortex",
                   "Left.Cerebellum.Cortex", "Right.Cerebellum.Cortex",
                   "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter",
                   "Left.Cerebral.White.Matter", "Right.Cerebral.White.Matter",
                   "Brain.Stem", "Left.VentralDC", "Right.VentralDC")))  %>% 
  
  filter(is.na(age) == F) %>% # remove participant without age data before making id consecutive
  
  mutate(id = 1:length(SubjID),
         .keep = 'unused')

# write processed data to csv file for outlier analysis etc

write_csv(ocd_ytop, file = '~/mrs_data/ocd_ytop_processed.csv')
  
# mark outliers by changing observations deviating more than 4 MAD to NA
  
ocd_ytop <- mutate(ocd_ytop, across(!c(id, female, age, ocd, eos, sbTIV), 
                ~ outlier_removal(.x, k=4)))
  
icv <- as.vector(scale(ocd_ytop$sbTIV)) # take intracranial volume out before pivoting to long format
                                 # as icv is used per individual, not per measurement

# define datasets for each analysis and pivot to long format

# subcortical volumes

d_v <- pivot_longer(ocd_ytop, cols = !c(id, female, age, ocd, eos, sbTIV), 
                         names_to = 'str_name', 
                         values_to = 'mri') %>%
  
  mutate(str_name = str_replace_all(str_name, c('[lr]h_' = '', 'Left.' = '', 'Right.' = '')),
         age_std = (age - min(age))/(max(age) - min(age)),
         gender = ifelse(female == 1, 0.5, -0.5),
         sbTIV = NULL,
         .keep = 'all') %>%
  
  # outlier observations are removed
  
  filter(is.na(mri)==FALSE & str_name%in%c('Amygdala', 'Hippocampus', 'Putamen', 
                                           'Pallidum', 'Caudate', 'Thalamus', 
                                           'Accumbens.area', 'Lateral.Ventricle')) %>%
  
  mutate(str_name = as_factor(str_name),
         str = as.integer(str_name))

# make and save plot illustrating that age can be represented by a linear effect

png(file = '~/mrs/plots/volume/subc_age_spline.png')

ggplot(d_v, aes(x = age, y = log(mri))) +
  geom_point() +
  geom_smooth(method = 'loess') +
  facet_wrap(vars(str_name)) +
  labs(title = 'MRI measurements across age', subtitle = '(With LOESS smooths)')

dev.off()

# save datafile for stan model in json format

cmdstanr::write_stan_json(
  file = '~/mrs_data/volume.json',
  data = list(
    N = length(unique(d_v$id)),
    n_obs = nrow(d_v),
    n_str = max(d_v$str),
    ind_id = d_v$id,
    ind_str = d_v$str,
    mri = d_v$mri,
    icv = icv,
    age = d_v$age_std,
    gender = d_v$gender,
    ocd = d_v$ocd,
    eos = d_v$eos,
    alpha_params = c(8, .1)))

# save data for plotting functions

write_rds(file = '~/mrs_data/volume_plot.rds',
          select(d_v, id, mri, str_name))

# cortical areas
  
d_a <- pivot_longer(ocd_ytop, cols = !c(id, female, age, ocd, eos, sbTIV), 
                  names_to = 'str_name', 
                  values_to = 'mri') %>%
  
  mutate(str_name = str_replace_all(str_name, c('[lr]h_' = '', 'Left.' = '', 'Right.' = '')),
         age_std = (age - min(age))/(max(age) - min(age)),
         gender = ifelse(female == 1, 0.5, -0.5),
         sbTIV = NULL,
         .keep = 'all') %>%
  
  # outlier observations are removed
  
  filter(is.na(mri)==FALSE & str_name%in%c('frontal_area', 'parietal_area', 'temporal_area', 
                                           'occipital_area', 'cingulate_area')) %>%
  
  mutate(str_name = as_factor(str_name),
         str = as.integer(str_name))

# make and save plot illustrating that age can be represented by a linear effect

png(file = '~/mrs/plots/area/areas_age_spline.png')

ggplot(d_a, aes(x = age, y = log(mri))) +
  geom_point() +
  geom_smooth(method = 'loess') +
  facet_wrap(vars(str_name)) +
  labs(title = 'MRI measurements across age', subtitle = '(With LOESS smooths)')

dev.off()

# save datafile for stan model in json format

cmdstanr::write_stan_json(
  file = '~/mrs_data/area.json',
  data = list(
    N = length(unique(d_a$id)),
    n_obs = nrow(d_a),
    n_str = max(d_a$str),
    ind_id = d_a$id,
    ind_str = d_a$str,
    mri = d_a$mri,
    icv = icv,
    age = d_a$age_std,
    gender = d_a$gender,
    ocd = d_a$ocd,
    eos = d_a$eos,
    alpha_params = c(9, .5)))

# save data for plotting functions

write_rds(file = '~/mrs_data/area_plot.rds',
          select(d_a, id, mri, str_name))

# cortical thickness

d_t <- pivot_longer(ocd_ytop, cols = !c(id, female, age, ocd, eos, sbTIV), 
                            names_to = 'str_name', 
                            values_to = 'mri') %>%
  
  mutate(str_name = str_replace_all(str_name, c('[lr]h_' = '', 'Left.' = '', 'Right.' = '')),
         age_std = (age - min(age))/(max(age) - min(age)),
         gender = ifelse(female == 1, 0.5, -0.5),
         sbTIV = NULL,
         .keep = 'all') %>%
  
  # outlier observations are removed
  
  filter(is.na(mri)==FALSE & str_name%in%c('frontal_thickness', 'parietal_thickness', 
                                           'temporal_thickness', 'occipital_thickness', 
                                           'cingulate_thickness')) %>%
  
  mutate(str_name = as_factor(str_name),
         str = as.integer(str_name))

# make and save plot illustrating that age can be represented by a linear effect

png(file = '~/mrs/plots/thickness/thickness_age_spline.png')

ggplot(d_t, aes(x = age, y = log(mri))) +
  geom_point() +
  geom_smooth(method = 'loess') +
  facet_wrap(vars(str_name)) +
  labs(title = 'MRI measurements across age', subtitle = '(With LOESS smooths)')

dev.off()

# save datafile for stan model in json format

cmdstanr::write_stan_json(
  file = '~/mrs_data/thickness.json',
  data = list(
    N = length(unique(d_t$id)),
    n_obs = nrow(d_t),
    n_str = max(d_t$str),
    ind_id = d_t$id,
    ind_str = d_t$str,
    mri = d_t$mri,
    icv = icv,
    age = d_t$age_std,
    gender = d_t$gender,
    ocd = d_t$ocd,
    eos = d_t$eos,
    alpha_params = c(.9, .2)))

# save data for plotting functions

write_rds(file = '~/mrs_data/thickness_plot.rds',
          select(d_t, id, mri, str_name))
