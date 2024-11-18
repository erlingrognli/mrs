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

caudate <- c('Right.Caudate',
             'Left.Caudate')

# define function for changing outliers to NA 
outlier_removal <- function(x, k){
  
  ub <- median(x) + k * mad(x)
  lb <- median(x) - k * mad(x)
  
  x[which(x>ub|x<lb)] <- NA
  
  return(x)}

# read in dataset, make id variable consecutive

d <- read_csv(file = '~/mrs_data/mrs_wf_data.csv',
              show_col_types = F) %>%  
  
  filter(is.na(age) == F) %>%
  
  mutate(id = 1:length(SubjID),
         .keep = 'unused') %>%
  
# mark outliers by changing observations deviating more than 4 MAD to NA
  
  mutate(across(!c(id, female, age, ocd, scz), 
                        ~ outlier_removal(.x, k=4))) %>% 
  
# merge measures of individual areas/volumes

  mutate(lh_PFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
         rh_PFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area
         ) %>%

# pivot to long format and remove outliers marked as NA
  
  pivot_longer(cols = !c(id, female, age, ocd, scz), 
               names_to = 'region', 
               values_to = 'mri') %>%
  
  filter(is.na(mri) == FALSE) %>%
  
# define measurement variable for thickness(1), area (2) or volume (3)
  
  mutate(measure = 
           case_when(endsWith(region, 'thickness') ~ 1,
                     endsWith(region, 'area') ~ 2,
                     endsWith(region, 'Cortex') ~ 2,
                     .default = 3))

# construct predictor matrix
 
pred <- mutate(d,
               vpfc_area_ocd = ifelse(region %in% vPFC_area & ocd == 1, 1, 0),
               vpfc_thickness_ocd = ifelse(region %in% vPFC_thickness & ocd == 1, 1, 0),
               caudate_ocd = ifelse(region %in% caudate & ocd == 1, 1, 0),
               .keep = 'none')


# remove markings of left or right hemisphere and change region variable to integer coding

d$region <- str_remove(d$region, '[lr]h_')
d$region <- str_remove(d$region, 'Left.')
d$region <- str_remove(d$region, 'Right.')

d$region_number <- as.integer(as_factor(d$region))

# make index of measure by region

ind_measure <- select(d, region, measure) %>%
  
  reframe(measure = unique(measure), .by=region) %>%
  
  arrange(region)

dat <- list(N = length(unique(d$id)),
            n_obs = nrow(d),
            n_region = max(d$region_number),
            n_beta = ncol(pred),
            ind_id = d$id,
            ind_region = d$region_number,
            pred = as.matrix(pred),
            ind_measure = ind_measure$measure,
            mri = d$mri)

lpr <- function(mu, sd){ round(exp(qnorm(c(.001, .05, .5, .95, .999), mu, sd)), digits = 2)}

      
m <- cmdstan_model('mri_mod_nospline.stan')


fit <- m$sample(data = dat, 
                       iter_warmup = 500,
                       iter_sampling = 500)

np <- nuts_params(fit)

png(file = 'pairs.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(starts_with('beta_exp'), sigma, starts_with('measure_icpt'), 'region_icpt[30]', 'region_icpt[60]', 'region_icpt[75]' ),
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
                   y_draw = 'violin')

dev.off()

png(file = 'violin_grouped_spl%d.png',
    width = 75,
    height = 45,
    units = 'cm',
    res = 100)

ppc_violin_grouped(y = log(d$mri)[which(d$measure==1)], 
                   yrep = ppc_draws[,which(d$measure==1)], 
                   group = d$region[which(d$measure==1)],
                   y_draw = 'violin')

ppc_violin_grouped(y = log(d$mri)[which(d$measure==2)], 
                   yrep = ppc_draws[,which(d$measure==2)], 
                   group = d$region[which(d$measure==2)],
                   y_draw = 'violin')

ppc_violin_grouped(y = log(d$mri)[which(d$measure==3)], 
                   yrep = ppc_draws[,which(d$measure==3)], 
                   group = d$region[which(d$measure==3)],
                   y_draw = 'violin')

dev.off()

ppc_pit_ecdf_grouped(log(d$mri), 
                     yrep = ppc_draws, 
                     group = d$measure,
                     plot_diff = TRUE)




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

summarise_draws(fit$draws(variables = c('beta_exp', 'measure_icpt', 'sigma')))

mcmc_areas(fit$draws(variables = c('beta_exp')))
