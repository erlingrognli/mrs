# preprocessing of the dataset for analyses

library(tidyverse)

# load demographic info on ocd participants

ocd_demog <- read.csv2('~/mrs_data/mrs_pasienter_demog.csv')

# remove excluded participants and irrelevant variables
# and fix "ok" text in subject ids.

ocd_demog <- select(filter(ocd_demog, consent == 1), c(MR4K.Nr, F.dt, K.nn))

ocd_demog$MR4K.Nr <- gsub(pattern = " ok", replacement = '', x = ocd_demog$MR4K.Nr)

# recode and rename variables, add indicator variable for ocd  

ocd_demog <- mutate(ocd_demog, 
                    birth_date = as.Date(F.dt, format = '%d.%m.%Y'),
                    female = ifelse(K.nn == 'F', 1, 0),
                    SubjID = MR4K.Nr,
                    ocd = 1,
                    scz = 0,
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
         scz = 0,
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
      scz = 0,
      .keep = 'unused')) %>%
  
  full_join(bind_rows(read_csv('~/mrs_data/OCD_sbTIV.csv', show_col_types = FALSE),
                      read_csv('~/mrs_data/YTOP_sbTIV.csv', show_col_types = FALSE)),
                      join_by(SubjID==SubjID)) %>%
  
  # remove brain measures that will not be used in the model
  
  select(!any_of(c("Unknown", "Right.choroid.plexus", "Left.choroid.plexus", 
                   "Soft_Nonbrain_Tissue", "Fluid_Inside_Eyes", "CSF", "Skull", 
                   "X5th.Ventricle", "X4th.Ventricle",  "X3rd.Ventricle",
                   "Right.Inf.Lat.Vent", "Left.Inf.Lat.Vent",
                   "Left.Lateral.Ventricle", "Right.Lateral.Ventricle",
                   "Optic.Chiasm", "Right.vessel", "Left.vessel",
                   "WM.hypointensities", "non.WM.hypointensities",
                   "Left.Cerebral.Cortex", "Right.Cerebral.Cortex",
                   "Left.Cerebellum.Cortex", "Right.Cerebellum.Cortex",
                   "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter",
                   "Left.Cerebral.White.Matter", "Right.Cerebral.White.Matter",
                   "Brain.Stem", "Left.VentralDC", "Right.VentralDC")))
  

write_csv(ocd_ytop, file = '~/mrs_data/mrs_wf_data.csv')

# identify and record outliers and values

outlier_id <- function(x, k, id){
  
  ub <- median(x) + k * mad(x)
  lb <- median(x) - k * mad(x)
  
  return(id[which(x>ub|x<lb)])
  
}

outliers_all <- list('> 5 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, scz)), 
                                               outlier_id, k=5, id = ocd_ytop$SubjID)),
                     '> 4 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, scz)), 
                                               outlier_id, k=4, id = ocd_ytop$SubjID)),
                     '> 3.5 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, scz)), 
                                                 outlier_id, k=3.5, id = ocd_ytop$SubjID))) %>%
  
  lapply(enframe, name = "region", value = "id") %>%
  
  enframe(name = 'distance') %>%
  
  unnest('value') %>%
  
  mutate(region = str_remove(region, pattern = '[1234567]'),
         mri = NA,
         median = NA)

for(i in 1:nrow(outliers_all)){
  outliers_all$mri[i] <-
    select(
      filter(
        ocd_ytop, SubjID == outliers_all$id[i]),
      outliers_all$region[i])
  
  outliers_all$median[i] <- sapply(select(ocd_ytop, outliers_all$region[i]), median)
  outliers_all$mri <- as.double(outliers_all$mri)}

write.csv(outliers_all, file = '~/mrs_data/mrs_outliers_mad.csv')

# remove outliers > 4 mad

outlier_removal <- function(x, k){
  
  ub <- median(x) + k * mad(x)
  lb <- median(x) - k * mad(x)
  
  x[which(x>ub|x<lb)] <- NA
  
  return(x)
}

ocd_ytop <- mutate(ocd_ytop, across(!c(SubjID, female, age, ocd, scz), 
                                    ~ outlier_removal(.x, k=4)))

# pivot data to long format, remove outlier observations and write to file

write.csv(
  
 pivot_longer(ocd_ytop,
    cols = !c(SubjID, female, age, ocd, scz), 
    names_to = 'area_name', 
    values_to = 'mri') %>%
    
  filter(is.na(mri) == FALSE),
  
  file = '~/mrs_data/mrs_lf_data.csv')



