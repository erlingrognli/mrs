library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

options(posterior.digits = 2,
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
  
  mutate(across(!c(id, female, age, ocd, scz, sbTIV), 
                        ~ outlier_removal(.x, k=4)))

  icv <- as.vector(scale(d$sbTIV))
  
# pivot to long format
  
d <- pivot_longer(d, cols = !c(id, female, age, ocd, scz, sbTIV), 
               names_to = 'str_name', 
               values_to = 'mri') %>%
    
  mutate(str_name = str_replace_all(str_name, c('[lr]h_' = '', 'Left.' = '', 'Right.' = '')),
         age_std = (age - min(age))/(max(age) - min(age)),
         sbTIV = NULL,
         .keep = 'all') %>%
    
  filter(is.na(mri)==FALSE & str_name%in%c('Amygdala', 'Hippocampus', 'Putamen', 'Pallidum', 'Caudate', 'Thalamus', 'Accumbens.area', 'Lateral.Ventricle')) %>%
  
  mutate(str_name = as_factor(str_name) %>% fct_relevel('Thalamus'), # make Thalamus reference region for intercepts
         str = as.integer(str_name))

# make and save plot illustrating that age can be represented by a linear effect

png(file = '~/mrs/plots/subc_age_spline.png')
  
  ggplot(d, aes(x = age, y = log(mri))) +
  geom_point() + 
  geom_smooth(method = 'loess') + 
  facet_wrap(vars(str_name)) + 
  labs(title = 'MRI measurements across age', subtitle = '(With LOESS smooths)')

  dev.off()

dat = dat <- list(
            N = length(unique(d$id)),
            n_obs = nrow(d),
            n_str = max(d$str),
            ind_id = d$id,
            ind_str = d$str,
            mri = d$mri,
            icv = icv,
            age = d$age_std,
            female = d$female,
            ocd = d$ocd,
            alpha_params = c(8, .1))


m <- cmdstan_model('mod_mrs.stan')

fit <- m$sample(data = dat, 
                iter_warmup = 1000,
                iter_sampling = 1500)

fit$cmdstan_diagnose()

setwd('~/mrs/plots')

png(file = 'pairs_hyper.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(beta_icv, sigma_alpha_id,
                                    mu_beta_age, mu_beta_female,
                                    mu_beta_ocd),
           np = np,
           max_treedepth = 10)

dev.off()

png(file = 'pairs_icpt.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars('beta_icv', 'sigma_alpha_id', starts_with('alpha_str_raw'), 
                                    'alpha_id[10]', 'alpha_id[15]', 'alpha_id[30]'),
           np = np,
           max_treedepth = 10)
dev.off()

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

bayesplot_theme_update(axis.text.x = element_text(angle = 90, hjust = 1))

png(file = 'violin_ppc.png',
    width = 45,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$str_name,
                   y_draw = 'both') +
  labs(title = 'Posterior predictive plot')

dev.off()

png(file = 'pit_ecdf.png',
    width = 36,
    height = 12,
    units = 'cm',
    res = 100)

  ppc_pit_ecdf_grouped(log(d$mri), 
                     yrep = ppc_draws, 
                     group = d$str_name,
                     plot_diff = TRUE)
dev.off()

png(file = 'ocd_betas.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'exp_beta_ocd') %>%
  
  rename_variables('Thalamus' = 'exp_beta_ocd[1]', 
                   'Amygdala' = 'exp_beta_ocd[2]', 
                   'Hippocampus' = 'exp_beta_ocd[3]',
                   'Accumbens(area)' = 'exp_beta_ocd[4]',
                   'Putamen' = 'exp_beta_ocd[5]',
                   'Lateral Ventricle' = 'exp_beta_ocd[6]',
                   'Pallidum' = 'exp_beta_ocd[7]',
                   'Caudate' = 'exp_beta_ocd[8]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD (posterior distributions)', 
       subtitle = 'Ajusted for age, gender and intracranial volume' )

  dev.off()


setwd('~/mrs')

loo_output <- fit$loo(moment_match = TRUE)

fit$summary(variables = c('mu_beta_female', 'mu_beta_age', 'beta_icv'))
