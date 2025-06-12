library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

options(posterior.num_args = list(digits = 2),
        mc.cores = 4)

# m <- cmdstan_model('~/mrs/mod_mrs.stan')
# 
# fit_thickness <- m$sample(data = 'R:/Prosjekter_VVHF/MRS_1000368/mrs_data/thickness.json',
#                 sig_figs = 9)

fit_thickness <- read_rds('fit_thickness.rds')

fit_thickness$cmdstan_diagnose() # no problems detected

setwd('~/mrs/plots/thickness')

# inspect sampling with pairs plots 

png(file = 'pairs_hyper.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit_thickness$draws(), pars = vars(beta_icv, sigma_alpha_id,
                                    sigma_alpha_str,
                                    mu_beta_age, mu_beta_female,
                                    mu_beta_ocd, mu_beta_eos))

dev.off()

png(file = 'pairs_icpt.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit_thickness$draws(), pars = vars('beta_icv', 'sigma_alpha_id', starts_with('alpha_str_raw'),
                                    'alpha_id[10]', 'alpha_id[15]', 'alpha_id[30]'))

dev.off()

# graphical posterior predictive checking

bayesplot_theme_update(axis.text.x = element_text(angle = 90, hjust = 1))

ppc_draws <- fit_thickness$draws(variables = 'ppc', format = 'draws_matrix')

d <- read_rds('R:/Prosjekter_VVHF/MRS_1000368/mrs_data/thickness_plot.rds')

png(file = 'violin_ppc_thickness.png',
    width = 30,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$str_name,
                   y_draw = 'violin') +
  labs(title = 'Posterior predictive plot')

dev.off()

png(file = 'violin_ppc_thickness_diagn.png',
    width = 30,
    height = 15,
    units = 'cm',
    res = 200)

ppc_violin_grouped(y = log(d$mri), 
                   yrep = ppc_draws, 
                   group = d$diagnosis,
                   y_draw = 'violin') +
  labs(title = 'Posterior predictive plot') 
  
  dev.off()

png(file = 'pit_ecdf_thickness.png',
    width = 36,
    height = 12,
    units = 'cm',
    res = 100)

  ppc_pit_ecdf_grouped(log(d$mri), 
                     yrep = ppc_draws,
                     group = d$str_name,
                     plot_diff = TRUE)
dev.off()

# loo-pit plot

loo_thickness <- fit_thickness$loo(moment_match = TRUE,
                                  save_psis = TRUE)

write_rds(loo_thickness, file = '~/mrs/loo_thickness.rds')

png(file = 'loo_pit_thickness.png',
    width = 12,
    height = 12,
    units = 'cm',
    res = 100)

ppc_loo_pit_overlay(log(d$mri), 
                     yrep = ppc_draws,
                     psis_object = loo_thickness$psis_object)
dev.off()

# plot estimates from model

png(file = 'ocd_betas_thickness.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

  fit_thickness$draws(variables = 'beta_ocd', format = 'draws_df') %>%
    
   mutate(across(starts_with('beta_ocd'), exp)) %>%
  
   rename_variables('Frontal cortex' = 'beta_ocd[1]', 
                   'Parietal cortex' = 'beta_ocd[2]',
                   'Temporal cortex' = 'beta_ocd[3]',
                   'Occipital cortex' = 'beta_ocd[4]',
                   'Cingulate cortex' = 'beta_ocd[5]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD on cortical thickness per lobe', 
       subtitle = 'Ajusted for age, gender and intracranial volume') + 
  vline_at(1)

  dev.off()

png(file = 'eos_betas_thickness.png',
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)

  fit_thickness$draws(variables = 'beta_eos', format = 'draws_df') %>%
    
    mutate(across(starts_with('beta_eos'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_eos[1]', 
                     'Parietal cortex' = 'beta_eos[2]',
                     'Temporal cortex' = 'beta_eos[3]',
                     'Occipital cortex' = 'beta_eos[4]',
                     'Cingulate cortex' = 'beta_eos[5]') %>%
    mcmc_areas() +
    labs(title = 'Multiplicative effect of EOS on cortical thickness per lobe', 
         subtitle = 'Ajusted for age, gender and intracranial volume' ) +
    vline_at(1)
  
  dev.off()
  
png(file = 'female_betas_thickness.png', 
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)
  
  fit_thickness$draws(variables = 'beta_female', format = 'draws_df') %>%
    
    mutate(across(starts_with('beta_female'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_female[1]', 
                     'Parietal cortex' = 'beta_female[2]',
                     'Temporal cortex' = 'beta_female[3]',
                     'Occipital cortex' = 'beta_female[4]',
                     'Cingulate cortex' = 'beta_female[5]') %>%
    
    mcmc_areas() +
    labs(title = 'Multiplicative effect of female gender on cortical thickness per lobe', 
         subtitle = 'Ajusted for age, diagnosis and intracranial volume' ) +
    vline_at(1)
  
  dev.off()

png(file = 'age_betas_thickness.png', 
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)
  
  fit_thickness$draws(variables = 'beta_age', format = 'draws_df') %>%
    
    mutate(across(starts_with('beta_age'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_age[1]', 
                     'Parietal cortex' = 'beta_age[2]',
                     'Temporal cortex' = 'beta_age[3]',
                     'Occipital cortex' = 'beta_age[4]',
                     'Cingulate cortex' = 'beta_age[5]') %>%
    
    mcmc_areas() +
    labs(title = 'Multiplicative effect of aging 7.3 years from 11.5 on cortical thickness per lobe',
         subtitle = 'Ajusted for gender, diagnosis and intracranial volume') +
    vline_at(1)
  
dev.off()
  
str_names <- levels(fct_recode(d$str_name, 
                               'Frontal cortex' = 'frontal_thickness', 
                               'Parietal cortex' = 'parietal_thickness',
                               'Temporal cortex' = 'temporal_thickness',
                               'Occipital cortex' = 'occipital_thickness',
                               'Cingulate cortex' = 'cingulate_thickness'))

ppd <- 
  
  bind_rows(
    fit_thickness$draws(variables = 'ppd_ocd', format = 'draws_df') %>% set_variables(str_names),
    fit_thickness$draws(variables = 'ppd_ctr', format = 'draws_df') %>% set_variables(str_names),
    fit_thickness$draws(variables = 'ppd_eos', format = 'draws_df') %>% set_variables(str_names),
    .id = 'Diagnosis') %>%
  
  select(!c(.chain, .iteration, .draw)) %>%
  
  mutate(Diagnosis = fct_recode(as_factor(Diagnosis), OCD = '1', Control = '2', EOS = '3'), .keep = 'unused') %>%
  
  pivot_longer(!Diagnosis, names_to = 'structure', values_to = 'thickness')

ppd_generated <- m$generate_quantities(fitted_params = fit_thickness, 
                                       data = 'R:/Prosjekter_VVHF/MRS_1000368/mrs_data/thickness.json')
ppd <- 
  
  bind_rows(
    ppd_generated$draws(variables = 'ppd_ocd', format = 'draws_df') %>% set_variables(str_names),
    ppd_generated$draws(variables = 'ppd_ctr', format = 'draws_df') %>% set_variables(str_names),
    ppd_generated$draws(variables = 'ppd_eos', format = 'draws_df') %>% set_variables(str_names),
    .id = 'dx') %>%
  
  select(!c(.chain, .iteration, .draw)) %>%
  
  mutate(Diagnosis = fct_recode(as_factor(dx), OCD = '1', Control = '2', EOS = '3'), .keep = 'unused') %>%
  
  pivot_longer(!Diagnosis, names_to = 'structure', values_to = 'thickness') %>%
  
  bind_rows(ppd)



png(file = 'ppd_thickness.png',
  width = 30,
  height = 20,
  units = 'cm',
  res = 400)

ggplot(data = ppd, aes(x = Diagnosis, y = thickness,
                       fill = Diagnosis)) + 
  geom_violin(alpha = .4,
              draw_quantiles = c(.05, .5, .95)) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(title = 'Posterior predictive distributions of cortical thickness across lobes', 
       y = 'Measurements in millimeters',
       x = NULL) + 
# annotate('text', x = 1:3, label = c(paste(mean(thickness), '(', quantile(thickness, c(.05, .95)), ')'))) +
  theme(legend.position = 'none') + 
  facet_wrap(vars(structure))

dev.off()

# writing posterior summary to file

prob_increase <- function(x) {length(which(x>1))/length(x)}
prob_decrease <- function(x) {length(which(x<1))/length(x)}

estimates <- fit_thickness$draws(format = 'draws_df') %>%
  
  select(starts_with('beta_'), starts_with('sigma'), starts_with('mu')) %>%
  
  mutate(across(everything(), exp)) %>%
  
  summarise_draws(mean, sd, quantile2, 
                  prob_increase, prob_decrease, 
                  ess_bulk, ess_tail, rhat) %>%
  
  mutate(across(!variable, \(x) round(x, digits = 2))) %>%
  
  mutate(variable = str_replace_all(variable, 
                                    c('\\[1\\]' = '_frontal',
                                      '\\[2\\]' = '_parietal',
                                      '\\[3\\]' = '_temporal',
                                      '\\[4\\]' = '_occipital',
                                      '\\[5\\]' = '_cingulate')))

write_csv(estimates, file = '~/mrs/thickness_estimates.csv')
