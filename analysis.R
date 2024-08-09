library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

caudate <- c('Left.Caudate', 'Right.Caudate')

vPFC_thickness <- c("rh_lateralorbitofrontal_thickness", 
                    "lh_lateralorbitofrontal_thickness", 
                    "rh_medialorbitofrontal_thickness", 
                    "lh_medialorbitofrontal_thickness")

vPFC_area <- c("rh_lateralorbitofrontal_area", 
               "lh_lateralorbitofrontal_area", 
               "rh_medialorbitofrontal_area", 
               "lh_medialorbitofrontal_area")

vPFC <- c(vPFC_area, vPFC_thickness)

occipital_thickness <- c("lh_lateraloccipital_thickness",
                         "rh_lateraloccipital_thickness")

occipital_area <- c("lh_lateraloccipital_area", 
                    "rh_lateraloccipital_area")

occipital <- c(occipital_area, occipital_thickness)

all_areas <- c(caudate, vPFC, occipital)

d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         .keep = 'unused') %>%
  
  pivot_longer(cols = !c(id, gender, age, ocd, scz), 
               names_to = 'area_name', 
               values_to = 'mri') %>%
  
  filter(area_name %in% all_areas)

d <- mutate(d, 
            ind_beta = 
              case_when(ocd == 1 & area_name %in% vPFC ~ 1,
                        ocd == 1 & area_name %in% caudate ~ 2,
                        ocd == 1 & area_name %in% occipital ~ 3,
                        ocd == 0 ~ 3),
            ind_area = 
              case_when(area_name %in% vPFC_area ~ 1,
                        area_name %in% vPFC_thickness ~ 2,
                        area_name %in% caudate ~ 3,
                        area_name %in% occipital_area ~ 4,
                        area_name %in% occipital_thickness ~ 5),
            .keep = 'all')

ggplot(data = d,
       mapping = aes(age, log(mri))) +
  geom_point() + 
  geom_smooth()

dat <- list(N = length(unique(d$id)),
            n_obs = nrow(d),
            n_area = max(d$ind_area),
            n_beta = max(d$ind_beta),
            ind_id = d$id,
            ind_beta_ocd = d$ind_beta,
            ind_area = d$ind_area,
            ind_gender = d$gender,
            mri = d$mri,
            age = d$age,
            n_knots = 5,
            age_knots = c(11.0, 13.0, 15.0, 17.0, 19.0))
      
m <- cmdstan_model('mri_mod.stan')

fit <- m$sample(data = dat)

mcmc_pairs(fit$draws(variables = c('sigma', 'beta_ocd', 'beta_gender')))

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$ind_area,
                   y_draw = 'both')

ppc_pit_ecdf(y = log(d$mri), yrep = ppc_draws)

betas <- summarise_draws(fit$draws(variables = c('beta', 'beta_female', 'sigma')))

areas <- summarise_draws(fit$draws(variables = c('area_icpt', 'var_area_icpt')))

ages <- summarise_draws(fit$draws(variables = c('knot_values')))


