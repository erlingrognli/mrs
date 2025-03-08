library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

options(posterior.digits = 2,
        mc.cores = 4)

m <- cmdstan_model('mod_mrs.stan', compile_model_methods = TRUE)

fit <- m$sample(data = '~/mrs_data/area.json', iter_sampling = 1000, sig_figs = 9)

fit$cmdstan_diagnose()

setwd('~/mrs/plots/area')

np <- nuts_params(fit)

png(file = 'pairs_hyper.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(beta_icv, sigma_alpha_id,
                                    sigma_alpha_str, sigma_beta_age,
                                    sigma_beta_female, sigma_beta_ocd,
                                    sigma_beta_eos,
                                    mu_beta_age, mu_beta_female,
                                    mu_beta_ocd, mu_beta_eos), 
           np = np,
           max_treedepth = 12)

dev.off()

png(file = 'pairs_icpt.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars('alpha', 'beta_icv', 'sigma_alpha_id', 'sigma_alpha_str', starts_with('alpha_str'),
                                    'alpha_id[10]', 'alpha_id[15]', 'alpha_id[30]'))

dev.off()

# graphical posterior predictive checking

bayesplot_theme_update(axis.text.x = element_text(angle = 90, hjust = 1))

d <- read_rds('~/mrs_data/area_plot.rds')

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

png(file = 'violin_ppc_area.png',
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

png(file = 'pit_ecdf_area.png',
    width = 36,
    height = 12,
    units = 'cm',
    res = 100)

  ppc_pit_ecdf_grouped(log(d$mri), 
                     yrep = ppc_draws, 
                     group = d$str_name,
                     plot_diff = TRUE)
dev.off()


# diagnostics using loo and loo-pit-plot

loo_area <- fit$loo(moment_match = TRUE, 
                    save_psis = TRUE)

write_rds(loo_area, file = '~/mrs/loo_area.rds')

png(file = 'loo_pit_area.png',
    width = 12,
    height = 12,
    units = 'cm',
    res = 100)

ppc_loo_pit_overlay(log(d$mri), 
                    yrep = ppc_draws,
                    psis_object = loo_area$psis_object)
dev.off()


# plot estimates from model

png(file = 'ocd_betas_area.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'beta_ocd', format = 'draws_df') %>%
  
  mutate(across(starts_with('beta_ocd'), exp)) %>%
  
  rename_variables('Frontal cortex' = 'beta_ocd[1]', 
                   'Parietal cortex' = 'beta_ocd[2]',
                   'Temporal cortex' = 'beta_ocd[3]',
                   'Occipital cortex' = 'beta_ocd[4]',
                   'Cingulate cortex' = 'beta_ocd[5]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD on cortical lobe area',
       subtitle = 'Ajusted for age, gender and intracranial volume' ) + 
  vline_at(1)

  dev.off()

png(file = 'eos_betas_area.png',
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)

fit$draws(variables = 'beta_eos', format = 'draws_df') %>%
  
    mutate(across(starts_with('beta_eos'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_eos[1]', 
                     'Parietal cortex' = 'beta_eos[2]',
                     'Temporal cortex' = 'beta_eos[3]',
                     'Occipital cortex' = 'beta_eos[4]',
                     'Cingulate cortex' = 'beta_eos[5]') %>%
    mcmc_areas() +
    labs(title = 'Multiplicative effect of EOS on cortical lobe area', 
         subtitle = 'Ajusted for age, gender and intracranial volume' ) +
    vline_at(1)
  
  dev.off()
  
png(file = 'female_betas_area.png', 
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)
  
  fit$draws(variables = 'beta_female', format = 'draws_df') %>%
    
    mutate(across(starts_with('beta_female'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_female[1]', 
                     'Parietal cortex' = 'beta_female[2]',
                     'Temporal cortex' = 'beta_female[3]',
                     'Occipital cortex' = 'beta_female[4]',
                     'Cingulate cortex' = 'beta_female[5]') %>%
    
    mcmc_areas() +
    labs(title = 'Multiplicative effect of female gender on cortical lobe area',
         subtitle = 'Ajusted for age, diagnosis and intracranial volume' ) +
    vline_at(1)
  
  dev.off()
  
png(file = 'age_betas_area.png', 
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)
  
  fit$draws(variables = 'beta_age', format = 'draws_df') %>%
    
    mutate(across(starts_with('beta_age'), exp)) %>%
    
    rename_variables('Frontal cortex' = 'beta_age[1]', 
                     'Parietal cortex' = 'beta_age[2]',
                     'Temporal cortex' = 'beta_age[3]',
                     'Occipital cortex' = 'beta_age[4]',
                     'Cingulate cortex' = 'beta_age[5]') %>%
    
    mcmc_areas() +
    labs(title = 'Multiplicative effect of aging 7.3 years from 11.5 on cortical lobe area',
         subtitle = 'Ajusted for gender, diagnosis and intracranial volume') +
    vline_at(1)
  
  dev.off()
  
  str_names <- levels(fct_recode(d$str_name, 
                                 'Frontal cortex' = 'frontal_area', 
                                 'Parietal cortex' = 'parietal_area',
                                 'Temporal cortex' = 'temporal_area',
                                 'Occipital cortex' = 'occipital_area',
                                 'Cingulate cortex' = 'cingulate_area'))

  ppd <- 
  
  bind_rows(
    fit$draws(variables = 'ppd_ocd', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_ctr', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_eos', format = 'draws_df') %>% set_variables(str_names),
    .id = 'dx') %>%
    
  select(!c(.chain, .iteration, .draw)) %>%
  
  mutate(Diagnosis = fct_recode(as_factor(dx), OCD = '1', Control = '2', EOS = '3'), .keep = 'unused') %>%
  
  pivot_longer(!Diagnosis, names_to = 'structure', values_to = 'volume')

png(file = 'ppd_area.png',
  width = 30,
  height = 20,
  units = 'cm',
  res = 400)

ggplot(data = ppd, aes(x = Diagnosis, y = volume,
                       fill = Diagnosis)) + 
  geom_violin(alpha = .4,
              draw_quantiles = c(.05, .5, .95)) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(title = 'Posterior predictive distributions of cortical lobe areas', 
       x = NULL,
       y = 'Measurements in square millimeters') + 
  theme(legend.position = 'none') + 
  facet_wrap(vars(structure), scales = 'free')

dev.off()

setwd('~/mrs')
