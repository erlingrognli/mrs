library(tidyverse)

d <- read.csv('~/mrs_data/Synthetic_neuroimaging_data.csv') %>%
  
    select(SubjID, Age, Sex, Diagnosis, Volume_1:Volume_5) %>%
    
    mutate(id = 1:length(SubjID), 
           age = Age, 
           female = ifelse(Sex == 'male', 0, 1), 
           schiz = ifelse(Diagnosis == 'CTR', 0, 1), 
           .keep = 'unused') %>%
  
    pivot_longer(cols = Volume_1:Volume_5, 
                 names_to = 'volume_id', 
                 names_prefix = "Volume_", 
                 values_to = 'volume', 
                 names_transform = list(volume_id = as.integer))

library(brms)

mod <- brmsformula(volume ~ s(age) + (1|id) + (1|volume_id) + female + schiz,
                   sigma ~ (1|volume), # letting errors vary by volume
                   family = 'gaussian')

get_prior(mod, data = d)

pr <- prior('normal(0, 5)', 'b')

fit <- brm(mod, data = d, prior = pr)

d <- readRDS('synth_data.rds')

d <- d$syn

d <- mutate(d, id = 1:nrow(d)) %>%
  
  pivot_longer(cols = !c(id, gender, age, ocd, scz), 
               names_to = 'area_name', 
               values_to = 'mri')
  
  


