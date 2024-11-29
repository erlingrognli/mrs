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

libary(tidyverse)

d <- readRDS('synth_data.rds')

d <- mutate(d, id = 1:nrow(d)) %>%
  
  pivot_longer(cols = !c(id, gender, age, ocd), 
               names_to = 'area_name', 
               values_to = 'mri')

vPFC_area <- c("rh_lateralorbitofrontal_area", 
               "lh_lateralorbitofrontal_area", 
               "rh_medialorbitofrontal_area", 
               "lh_medialorbitofrontal_area",
               'lh_vPFC_area',
               'rh_vPFC_area')

caudate <- c('Left.Caudate', 'Right.Caudate')

occipital_area <- c("lh_lateraloccipital_area", 
                    "rh_lateraloccipital_area")

d <- mutate(d, ind_area = 
              case_when(area_name %in% caudate ~ 1,
                        area_name %in% vPFC_area ~ 2,
                        area_name %in% occipital_area ~ 3),
            .keep = 'all')

ind_pred <- mutate(d, 
                   gender = ifelse(gender == 1, -.5, .5),
                   vpfc_area_ocd = ifelse(area_name %in% vPFC_area & ocd == 1, 1, 0),
                   caudate_ocd = ifelse(area_name %in% caudate & ocd == 1, 1, 0),
                   .keep = 'none')

dat <- list(N = length(unique(d$id)),
            n_obs = nrow(d),
            n_area = max(d$ind_area),
            n_beta = ncol(ind_pred),
            ind_id = d$id,
            ind_area = d$ind_area,
            ind_pred = as.matrix(ind_pred),
            mri = d$mri,
            age = d$age,
            n_knots = 5,
            age_knots = c(11.0, 13.0, 15.0, 17.0, 19.0))

library(cmdstanr); library(bayesplot);library(posterior)

options(mc.cores = 4)

m <- cmdstan_model('mri_mod.stan')

fit <- m$sample(data = dat)

saveRDS(fit, 'fit.rds')

png(file = 'pairs.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(variables = c('area_icpt', 'beta', 'knot_values')))

dev.off()

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

png(file = 'ppc_violin.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$ind_area,
                   y_draw = 'both')

dev.off()

ppc_pit_ecdf(y = log(d$mri), yrep = ppc_draws)

betas <- summarise_draws(fit$draws(variables = c('beta', 'sigma')))

areas <- summarise_draws(fit$draws(variables = c('area_icpt', 'var_area_icpt')))

ages <- summarise_draws(fit$draws(variables = c('knot_values')))


