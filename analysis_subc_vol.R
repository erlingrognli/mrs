library(cmdstanr); library(bayesplot); library(posterior)

options(posterior.digits = 2,
        mc.cores = 4)

m <- cmdstan_model('~/mrs/mod_mrs.stan')

fit <- m$sample(data = '~/mrs_data/volume.json', 
                iter_sampling = 1000,
                sig_figs = 9)

fit$cmdstan_diagnose()

# inspect sampling with pairs plots

setwd('~/mrs/plots/volume')

png(file = 'pairs_hyper.png',
    width = 45,
    height = 45,
    units = 'cm',
    res = 100)

mcmc_pairs(fit$draws(), pars = vars(beta_icv, sigma_alpha_id,
                                    sigma_alpha_str, sigma_beta_age,
                                    sigma_beta_gender, sigma_beta_ocd,
                                    sigma_beta_eos,
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

d <- read_rds('~/mrs_data/volume_plot.rds')

ppc_draws <- fit$draws(variables = 'ppc', format = 'draws_matrix')

png(file = 'violin_ppc_volume.png',
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

# loo-pit plot

loo_volume <- fit$loo(moment_match = TRUE, 
                      save_psis = TRUE)

write_rds(loo_volume, file = '~/mrs/loo_volume.rds')

png(file = 'loo_pit_volume.png',
    width = 12,
    height = 12,
    units = 'cm',
    res = 100)

ppc_loo_pit_overlay(log(d$mri), 
                    yrep = ppc_draws,
                    psis_object = loo_volume$psis_object)
dev.off()

# plot model estimates

png(file = 'ocd_betas_volume.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'beta_ocd', format = 'draws_df') %>%
  
  mutate(across(starts_with('beta_ocd'), exp)) %>%
  
  rename_variables('Amygdala' = 'beta_ocd[1]', 
                   'Hippocampus' = 'beta_ocd[2]',
                   'Accumbens(area)' = 'beta_ocd[3]',
                   'Putamen' = 'beta_ocd[4]',
                   'Lateral Ventricle' = 'beta_ocd[5]',
                   'Pallidum' = 'beta_ocd[6]',
                   'Caudate' = 'beta_ocd[7]',
                   'Thalamus' = 'beta_ocd[8]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD on subcortical volumes', 
       subtitle = 'Ajusted for age, gender and intracranial volume' ) + 
  vline_at(1)

  dev.off()

png(file = 'eos_betas_volume.png',
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)

fit$draws(variables = 'beta_eos', format = 'draws_df') %>%
  
    mutate(across(starts_with('beta_eos'), exp)) %>%
    
    rename_variables('Amygdala' = 'beta_eos[1]',
                     'Hippocampus' = 'beta_eos[2]',
                     'Accumbens(area)' = 'beta_eos[3]',
                     'Putamen' = 'beta_eos[4]',
                     'Lateral Ventricle' = 'beta_eos[5]',
                     'Pallidum' = 'beta_eos[6]',
                     'Caudate' = 'beta_eos[7]',
                     'Thalamus' = 'beta_eos[8]') %>%
    mcmc_areas() +
    labs(title = 'Multiplicative effect of EOS on subcortical volumes', 
         subtitle = 'Ajusted for age, gender and intracranial volume' ) +
    vline_at(1)
  
  dev.off()
  
png(file = 'female_betas_volume.png', 
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'beta_female', format = 'draws_df') %>%
  
  mutate(across(starts_with('beta_female'), exp)) %>%
  
  rename_variables('Amygdala' = 'beta_female[1]',
                   'Hippocampus' = 'beta_female[2]',
                   'Accumbens(area)' = 'beta_female[3]',
                   'Putamen' = 'beta_female[4]',
                   'Lateral Ventricle' = 'beta_female[5]',
                   'Pallidum' = 'beta_female[6]',
                   'Caudate' = 'beta_female[7]',
                   'Thalamus' = 'beta_female[8]') %>%
  
  mcmc_areas() +
  labs(title = 'Multiplicative effect of female gender on subcortical volumes', 
       subtitle = 'Ajusted for age, diagnosis and intracranial volume' ) +
  vline_at(1)

dev.off()

png(file = 'age_betas_volume.png', 
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'beta_age', format = 'draws_df') %>%
  
  mutate(across(starts_with('beta_age'), exp)) %>%
  
  rename_variables('Amygdala' = 'beta_age[1]',
                   'Hippocampus' = 'beta_age[2]',
                   'Accumbens(area)' = 'beta_age[3]',
                   'Putamen' = 'beta_age[4]',
                   'Lateral Ventricle' = 'beta_age[5]',
                   'Pallidum' = 'beta_age[6]',
                   'Caudate' = 'beta_age[7]',
                   'Thalamus' = 'beta_age[8]') %>%
  
  mcmc_areas() +
  labs(title = 'Multiplicative effect of aging 7.3 years from 11.5 on subcortical volumes', 
       subtitle = 'Ajusted for age, diagnosis and intracranial volume' ) +
  vline_at(1)

dev.off()

str_names <- levels(fct_recode(d$str_name, 'Accumbens (area)' = 'Accumbens.area', 'Lateral Ventricle' = 'Lateral.Ventricle'))

ppd <- 
  
  bind_rows(
    fit$draws(variables = 'ppd_ocd', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_ctr', format = 'draws_df') %>% set_variables(str_names),
    fit$draws(variables = 'ppd_eos', format = 'draws_df') %>% set_variables(str_names),
    .id = 'dx') %>%
  
  select(!c(.chain, .iteration, .draw)) %>%
  
  mutate(Diagnosis = fct_recode(as_factor(dx),  OCD = '1', Control = '2', EOS = '3'), .keep = 'unused') %>%
  
  pivot_longer(!Diagnosis, names_to = 'structure', values_to = 'volume')

png(file = 'ppd_volume.png',
  width = 30,
  height = 20,
  units = 'cm',
  res = 400)

ggplot(data = ppd, aes(y = volume, x = Diagnosis)) + 
  geom_violin(alpha = .4, 
              draw_quantiles = c(.05, .5, .95),
              aes(fill = Diagnosis)) + 
  viridis::scale_fill_viridis(discrete = TRUE) +
  labs(title = 'Posterior predictive distributions across subcortical structures', 
       x = NULL,
       y = 'Measurements in cubic/square millimeters') + 
  theme(legend.position = 'none') + 
  facet_wrap(vars(structure), scales = 'free')

dev.off()

setwd('~/mrs')
