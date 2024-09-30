library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

caudate <- c('Left.Caudate', 'Right.Caudate')

vPFC_thickness <- c("rh_lateralorbitofrontal_thickness", 
                    "lh_lateralorbitofrontal_thickness", 
                    "rh_medialorbitofrontal_thickness", 
                    "lh_medialorbitofrontal_thickness")

vPFC_area <- c("rh_lateralorbitofrontal_area", 
               "lh_lateralorbitofrontal_area", 
               "rh_medialorbitofrontal_area", 
               "lh_medialorbitofrontal_area",
               'lh_vPFC_area',
               'rh_vPFC_area')

vPFC <- c(vPFC_area, vPFC_thickness)

occipital_area <- c("lh_lateraloccipital_area", 
                    "rh_lateraloccipital_area")

all_areas <- c(caudate, vPFC, occipital_area)

d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         lh_vPFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
         rh_vPFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area,
         .keep = 'unused') %>%
  
  pivot_longer(cols = !c(id, female, age, ocd, scz), 
               names_to = 'area_name', 
               values_to = 'mri') %>%
  
  filter(area_name %in% all_areas)

d <- mutate(d, ind_area = 
              case_when(area_name %in% vPFC_thickness ~ 1,
                        area_name %in% caudate ~ 2,
                        area_name %in% vPFC_area ~ 3,
                        area_name %in% occipital_area ~ 4),
            .keep = 'all')

ind_pred <- mutate(d, 
                   #female = female,
                   vpfc_area_ocd = ifelse(area_name %in% vPFC_area & ocd == 1, 1, 0),
                   vpfc_thickness_ocd = ifelse(area_name %in% vPFC_thickness & ocd == 1, 1, 0),
                   caudate_ocd = ifelse(area_name %in% caudate & ocd == 1, 1, 0),
                   .keep = 'none')

ggplot(data = filter(d, ind_area == 1),
       mapping = aes(age, log(mri))) +
  geom_point() + 
  geom_smooth()



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


      
m <- cmdstan_model('mri_mod.stan')

fit <- m$sample(data = dat, 
                       iter_warmup = 500,
                       iter_sampling = 500)

loo_output <- fit$loo(moment_match = TRUE)

pareto_k <- loo_output[['diagnostics']][['pareto_k']]

wide_d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
                     show_col_types = F) %>%
  
  filter(is.na(age) == F) 

bad_k <- filter(wide_d, pareto_k > .7) # not immediately clear why these are problematic/influential


png(file = 'pairs.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(variables = c('area_icpt', 'beta', 'knot_values')))

dev.off()

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$ind_area,
                   y_draw = 'both')

ppc_pit_ecdf(y = log(d$mri), yrep = ppc_draws)

# overall adequate fit

# PPC plot of vPFC area separately for ocd and controls

ppc_violin_grouped(y = log(d$mri[which(d$ind_area==3)]),
                   yrep = ppc_draws[,which(d$ind_area==3)],
                   group = d$ocd[which(d$ind_area==3)],
                   y_draw = 'both')


ppc_violin_grouped(y = log(d$mri[which(d$ind_area!=1)]),
                   yrep = ppc_draws[,which(d$ind_area!=1)],
                   group = d$female[which(d$ind_area!=1)],
                   y_draw = 'both')



# there is something going on here, in the lower tail for the ocd patients

fit$init_model_methods()

betas <- fit$draws(variables = 'beta')

betas <- summarise_draws(fit$draws(variables = c('beta', 'sigma')))

mcmc_areas(fit$draws(variables = 'beta'))

areas <- summarise_draws(fit$draws(variables = c('area_icpt', 'var_area_icpt')))

ages <- summarise_draws(fit$draws(variables = c('knot_values')))

exp(betas$mean)
