library(cmdstanr); library(bayesplot); library(posterior)

options(posterior.digits = 2,
        mc.cores = 4)

m <- cmdstan_model('mod_mrs.stan')

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

# comparing with null model using loo

m_null <- cmdstan_model('~/mrs/null_mod_mrs.stan')

null_fit <- m_null$sample(data = '~/mrs_data/volume.json', 
                          iter_sampling = 1000, 
                          sig_figs = 9)

null_fit$cmdstan_diagnose()

loo_volume <- list(fit$loo(moment_match = TRUE),
                 null_fit$loo(moment_match = TRUE))

write_rds(loo_volume, file = '~/mrs/loo_volume.rds')

# plot model estimates

png(file = 'ocd_betas.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'exp_beta_ocd') %>%
  
  rename_variables('Amygdala' = 'exp_beta_ocd[1]', 
                   'Hippocampus' = 'exp_beta_ocd[2]',
                   'Accumbens(area)' = 'exp_beta_ocd[3]',
                   'Putamen' = 'exp_beta_ocd[4]',
                   'Lateral Ventricle' = 'exp_beta_ocd[5]',
                   'Pallidum' = 'exp_beta_ocd[6]',
                   'Caudate' = 'exp_beta_ocd[7]',
                   'Thalamus' = 'exp_beta_ocd[8]') %>%
  mcmc_areas() +
  labs(title = 'Multiplicative effect of OCD (posterior distributions)', 
       subtitle = 'Ajusted for age, gender and intracranial volume' ) + 
  vline_at(1)

  dev.off()

png(file = 'eos_betas.png',
      width = 20,
      height = 20,
      units = 'cm',
      res = 200)

fit$draws(variables = 'exp_beta_eos') %>%
    
    rename_variables('Amygdala' = 'exp_beta_eos[1]',
                     'Hippocampus' = 'exp_beta_eos[2]',
                     'Accumbens(area)' = 'exp_beta_eos[3]',
                     'Putamen' = 'exp_beta_eos[4]',
                     'Lateral Ventricle' = 'exp_beta_eos[5]',
                     'Pallidum' = 'exp_beta_eos[6]',
                     'Caudate' = 'exp_beta_eos[7]',
                     'Thalamus' = 'exp_beta_eos[8]') %>%
    mcmc_areas() +
    labs(title = 'Multiplicative effect of EOS (posterior distributions)', 
         subtitle = 'Ajusted for age, gender and intracranial volume' ) +
    vline_at(1)
  
  dev.off()
  
png(file = 'gender_betas.png', 
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)

fit$draws(variables = 'beta_gender', format = 'draws_df') %>%
  
  mutate(across(starts_with('beta_gender'), exp)) %>%
  
  rename_variables('Amygdala' = 'beta_gender[1]',
                   'Hippocampus' = 'beta_gender[2]',
                   'Accumbens(area)' = 'beta_gender[3]',
                   'Putamen' = 'beta_gender[4]',
                   'Lateral Ventricle' = 'beta_gender[5]',
                   'Pallidum' = 'beta_gender[6]',
                   'Caudate' = 'beta_gender[7]',
                   'Thalamus' = 'beta_gender[8]') %>%
  
  mcmc_areas() +
  labs(title = 'Multiplicative differences by gender, female compared to male participants (posterior distributions)', 
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

fit$summary(variables = c('mu_beta_gender', 'sigma_beta_gender', 'mu_beta_age', 'sigma_beta_age', 'beta_icv'))
