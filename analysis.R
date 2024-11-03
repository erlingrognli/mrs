library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

# define groupings of columns to fit effect of ocd or psychosis for

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

caudate <- c('rh_caudate',
             'lh_caudate')

# read in dataset, make id variable consecutive

d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%  
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         .keep = 'unused') %>%
  
# merge measures of individual areas/volumes

  mutate(lh_vPFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
         rh_vPFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area
         ) %>%

# pivot to long format
  
  pivot_longer(cols = !c(id, female, age, ocd, scz), 
               names_to = 'region', 
               values_to = 'mri') %>%
  
# define variable for thickness(1), area (2) or volume (3), for spline function
  
  mutate(ind_spl = 
           case_when(endsWith(region, 'thickness') ~ 1,
                     endsWith(region, 'area') ~ 2,
                     endsWith(region, 'Cortex') ~ 2,
                     .default = 3))
 
    
# d <- mutate(d, ind_re = 
#               case_when(area_name %in% vPFC_thickness ~ 1,
#                         area_name %in% caudate ~ 2,
#                         area_name %in% vPFC_area ~ 3,
#                         area_name %in% occipital_area ~ 4),
#             .keep = 'all')

ind_pred <- mutate(d, 
                   female = female,
                   vpfc_area_ocd = ifelse(region %in% vPFC_area & ocd == 1, 1, 0),
                   vpfc_thickness_ocd = ifelse(region %in% vPFC_thickness & ocd == 1, 1, 0),
                   caudate_ocd = ifelse(region %in% caudate & ocd == 1, 1, 0),
                   .keep = 'none')

d$region <- str_remove(d$region, '[lr]h_')
d$region <- str_remove(d$region, 'Left.')
d$region <- str_remove(d$region, 'Right.')

d$region <- as.integer(as.factor(d$region))



ggplot(data = filter(d, region == 3),
       mapping = aes(age, log(mri))) +
  geom_point() + 
  geom_smooth()



dat <- list(N = length(unique(d$id)),
            n_obs = nrow(d),
            n_region = max(d$region),
            n_beta = ncol(ind_pred),
            n_spl = as.vector(table(d$ind_spl)),
            ind_id = d$id,
            ind_region = d$region,
            ind_spl = d$ind_spl,
            ind_pred = as.matrix(ind_pred),
            mri = d$mri,
            age = d$age,
            n_knots = 3,
            age_knots = c(11.0, 15.0, 19.0))


      
m <- cmdstan_model('mri_mod.stan')


fit <- m$sample(data = dat, 
                       iter_warmup = 1000,
                       iter_sampling = 1000)

np <- nuts_params(fit)

png(file = 'pairs%d.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(variables = c('beta', 'alpha', 'sigma')),
           np = np,
           max_treedepth = 10)

mcmc_pairs(fit$draws(variables = c('thick_knot_values', 'area_knot_values', 'volume_knot_values', 'alpha')),
           np = np,
           max_treedepth = 10)

mcmc_pairs(fit$draws(variables = c('region_icpt_raw'[1:20])),
           np = np,
           max_treedepth = 10)

dev.off()

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

png(file = 'violin_grouped.png',
    width = 75,
    height = 45,
    units = 'cm',
    res = 100)

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$region,
                   y_draw = 'both')

dev.off()

png(file = 'violin_grouped_spl%d.png',
    width = 75,
    height = 45,
    units = 'cm',
    res = 100)

ppc_violin_grouped(y = log(d$mri)[which(d$ind_spl==1)], 
                   yrep = ppc_draws[,which(d$ind_spl==1)], 
                   group = d$region[which(d$ind_spl==1)],
                   y_draw = 'both')

ppc_violin_grouped(y = log(d$mri)[which(d$ind_spl==2)], 
                   yrep = ppc_draws[,which(d$ind_spl==2)], 
                   group = d$region[which(d$ind_spl==2)],
                   y_draw = 'both')

ppc_violin_grouped(y = log(d$mri)[which(d$ind_spl==3)], 
                   yrep = ppc_draws[,which(d$ind_spl==3)], 
                   group = d$region[which(d$ind_spl==3)],
                   y_draw = 'both')

dev.off()

ppc_pit_ecdf(y = log(d$mri), yrep = ppc_draws)

# overall adequate fit

loo_output <- fit$loo(moment_match = TRUE)

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

summarise_draws(fit$draws(variables = c('beta_exp', 'sigma_exp')))

mcmc_areas(fit$draws(variables = 'beta_exp'))

areas <- summarise_draws(fit$draws(variables = c('area_icpt', 'var_area_icpt')))

ages <- summarise_draws(fit$draws(variables = c('knot_values')))
