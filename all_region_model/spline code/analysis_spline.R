library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

options(posterior.num_args = c('digits' = 2),
        mc.cores = 4)


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
  
  # mutate(lh_vPFC_area = lh_lateralorbitofrontal_area + lh_medialorbitofrontal_area,
  #        rh_vPFC_area = rh_lateralorbitofrontal_area + rh_medialorbitofrontal_area
  #        ) %>%
  
  # pivot to long format and remove outliers marked as NA
  
  pivot_longer(cols = !c(id, female, age, ocd, scz), 
               names_to = 'region', 
               values_to = 'mri') %>%
  
  filter(is.na(mri) == FALSE) %>%
  
  # define measurement variable for thickness(1), area (2) or volume (3)
  # remove markings of left or right hemisphere
  
  mutate(measure = case_when(endsWith(region, 'thickness') ~ 1,
                             endsWith(region, 'area') ~ 2,
                             endsWith(region, 'Cortex') ~ 2,
                             .default = 3),
         region = str_replace_all(region, c('[lr]h_' = '', 'Left.' = '', 'Right.' = '')))


# make index of measure by region, selecting three reference areas as region 1-3

d$region_number <- as_factor(d$region) %>%
  
  fct_relevel('lateraloccipital_thickness', 'lateraloccipital_area', 'Brain.Stem') %>%
  
  as.integer()

ind_measure <- select(d, region, measure, region_number) %>%
  
  reframe(measure = unique(measure), region_number = unique(region_number), .by=region) %>%
  
  arrange(region_number)

# construct predictor matrix
# define groupings of columns to fit effect of ocd or psychosis for

vPFC_thickness <- c("lateralorbitofrontal_thickness", 
                    "medialorbitofrontal_thickness")

vPFC_area <- c('lateralorbitofrontal_area', 
               'medialorbitofrontal_area')

pred <- mutate(d,
               vpfc_area_ocd = ifelse(region %in% vPFC_area & ocd == 1, 1, 0),
               vpfc_thickness_ocd = ifelse(region %in% vPFC_thickness & ocd == 1, 1, 0),
               caudate_ocd = ifelse(region == 'Caudate' & ocd == 1, 1, 0),
               .keep = 'none')


dat <- list(N = length(unique(d$id)),
            n_obs = nrow(d),
            n_knots = 3,
            n_spl = as.vector(table(d$measure)),
            n_region = max(d$region_number),
            n_beta = ncol(pred),
            ind_id = d$id,
            ind_region = d$region_number,
            pred = as.matrix(pred),
            ind_measure = ind_measure$measure,
            mri = d$mri,
            age = d$age,
            age_knots = c(11, 15, 19))

lpr <- function(mu, sd){ round(exp(qnorm(c(.001, .05, .5, .95, .999), mu, sd)), digits = 2)}

m <- cmdstan_model('mri_mod_spline.stan')


pf <- m$pathfinder(data = dat,
                   save_cmdstan_config = TRUE)

write_rds(pf, file = 'pf_spline.rds')

pf <- read_rds(pf, file = 'pf_spline.rds')

fit <- m$sample(data = dat, 
                       iter_warmup = 300,
                       iter_sampling = 300,
                init = pf)

np <- nuts_params(fit)

png(file = 'pairs_spl.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(starts_with('knot_values_raw'), icpt, 'region_icpt[30]', 'region_icpt[60]', 'region_icpt[75]', 'id_icpt_raw[10]', 'id_icpt_raw[30]', 'measure_icpt[2]', 'measure_icpt[3]', starts_with('sigma')))

dev.off()

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

bayesplot_theme_update(axis.text.x = element_text(angle = 90, hjust = 1))

png(file = 'violin_thickness_spl.png',
    width = 45,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri)[which(d$measure==1)], 
                   yrep = ppc_draws[,which(d$measure==1)], 
                   group = d$region[which(d$measure==1)],
                   y_draw = 'violin')
dev.off()

png(file = 'violin_area_spl.png',
    width = 45,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri)[which(d$measure==2)], 
                   yrep = ppc_draws[,which(d$measure==2)], 
                   group = d$region[which(d$measure==2)],
                   y_draw = 'violin')
dev.off()

png(file = 'violin_volume_spl.png',
    width = 13,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri)[which(d$measure==3)], 
                   yrep = ppc_draws[,which(d$measure==3)], 
                   group = d$region[which(d$measure==3)],
                   y_draw = 'violin')
dev.off()

png(file = 'pit_ecdf_spl.png',
    width = 90,
    height = 30,
    units = 'cm',
    res = 100)

  ppc_pit_ecdf_grouped(log(d$mri), 
                     yrep = ppc_draws, 
                     group = d$measure,
                     plot_diff = TRUE)
dev.off()


loo_output <- fit$loo(moment_match = TRUE)

write_rds(loo_output, file = 'loo_spline.rds')

regpar <- summarise_draws(fit$draws(variables = c('beta_exp', 'sigma', 'alpha')))

icpts <- summarise_draws(fit$draws(variables = c('measure_icpt_raw', 'region_icpt_raw')))