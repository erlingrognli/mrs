library(cmdstanr); library(tidyverse); library(ggplot2); library(bayesplot); library(posterior)

options(posterior.digits = 2,
        mc.cores = 4)

m <- cmdstan_model('mod_mrs.stan')

fit <- m$sample(data = '~/mrs_data/thickness.json',
                sig_figs = 9)

fit$cmdstan_diagnose() # no problems detected

setwd('~/mrs/plots/thickness')

# inspect sampling with pairs plots 

png(file = 'pairs_hyper.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(beta_icv, sigma_alpha_id,
                                    sigma_alpha_str,
                                    mu_beta_age, mu_beta_gender,
                                    mu_beta_ocd, mu_beta_eos))

dev.off()

png(file = 'pairs_icpt.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars('beta_icv', 'sigma_alpha_id', starts_with('alpha_str_raw'),
                                    'alpha_id[10]', 'alpha_id[15]', 'alpha_id[30]'))

dev.off()

# graphical posterior predictive checking

bayesplot_theme_update(axis.text.x = element_text(angle = 90, hjust = 1))

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

d <- read_rds('~/mrs_data/thickness_plot.rds')

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

setwd('~/mrs')

# compare with model without predictors using loo

m_null <- cmdstan_model('null_mod_mrs.stan')

null_fit <- m_null$sample(data = '~/mrs_data/thickness.json',
                          sig_figs = 9)

null_fit$cmdstan_diagnose()

loo_thickness <- list(fit$loo(moment_match = TRUE),
                   null_fit$loo(moment_match = TRUE))

# plot estimates from model

setwd('~/mrs/plots/thickness')

png(file = 'ocd_betas_thickness.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'exp_beta_ocd') %>%
  
  rename_variables('Frontal cortex' = 'exp_beta_ocd[1]', 
                   'Parietal cortex' = 'exp_beta_ocd[2]',
                   'Temporal cortex' = 'exp_beta_ocd[3]',
                   'Occipital cortex' = 'exp_beta_ocd[4]',
                   'Cingulate cortex' = 'exp_beta_ocd[5]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD on cortical thickness (posterior distributions)', 
       subtitle = 'Ajusted for age, gender and intracranial volume' ) + 
  vline_at(1)

  dev.off()

png(file = 'eos_betas_thickness.png',
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)

fit$draws(variables = 'exp_beta_eos') %>%
    
    rename_variables('Frontal cortex' = 'exp_beta_eos[1]', 
                     'Parietal cortex' = 'exp_beta_eos[2]',
                     'Temporal cortex' = 'exp_beta_eos[3]',
                     'Occipital cortex' = 'exp_beta_eos[4]',
                     'Cingulate cortex' = 'exp_beta_eos[5]') %>%
    mcmc_areas() +
    labs(title = 'Multiplicative effect of EOS on cortical thickness (posterior distributions)', 
         subtitle = 'Ajusted for age, gender and intracranial volume' ) +
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
    fit$draws(variables = 'ppd_ctr', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_ocd', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_eos', format = 'draws_df') %>% set_variables(str_names),
    .id = 'Diagnosis') %>%
  
  select(!c(.chain, .iteration, .draw)) %>%
  
  mutate(Diagnosis = fct_recode(as_factor(Diagnosis), Control = '1', OCD = '2', EOS = '3'), .keep = 'unused') %>%
  
  pivot_longer(!Diagnosis, names_to = 'structure', values_to = 'volume')

png(file = 'ppd_thickness.png',
  width = 30,
  height = 30,
  units = 'cm',
  res = 200)

ggplot(data = ppd, aes(x = volume, 
                       colour = Diagnosis,
                       fill = Diagnosis)) + 
  geom_density(alpha = .4) + 
  scale_fill_discrete() + 
  labs(title = 'Posterior predictive distributions of thickness across cortical areas', 
       x = 'Measurements in millimeters',
       y = NULL) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  facet_wrap(vars(structure), scales = 'free')

dev.off()

