library(cmdstanr); library(bayesplot); library(showtext); library(ggplot2); library(magrittr); library(posterior)

# fit_volume <- readRDS('fit_volume.rds')
# fit_area <- readRDS('fit_area.rds')
# 
# ocd_volume <- fit_volume$draws(variables = 'beta_ocd', format = 'draws_df')
# eos_volume <- fit_volume$draws(variables = 'beta_eos', format = 'draws_df')
# ocd_area <- fit_area$draws(variables = 'beta_ocd', format = 'draws_df')
# eos_area <- fit_area$draws(variables = 'beta_eos', format = 'draws_df')

# save('ocd_volume', 'ocd_area', 'eos_volume', 'eos_area', file = 'posterplot_data.Rdata')

load()

font_add("calibri", "C:/Windows/Fonts/calibri.ttf")
showtext_auto()

ahus_lysblaa <- "#D5EEF6"
ahus_mellomblaa <- "#6CADDF"
ahus_klarblaa <- "#00529B"
ahus_morkblaa <- "#002042"

bayesplot_theme_set(theme_default(base_size = 36, base_family = "calibri") +
                      theme(
                        plot.background = element_rect(fill = ahus_lysblaa),  # Background color
                        panel.background = element_rect(fill = 'white'),     # Panel background color
                        axis.text = element_text(size = 22, family = "calibri"),  # Axis text font & size
                        axis.title = element_text(size = 28, family = "calibri"), # Axis title font & size
                        axis.line = element_line(color = ahus_morkblaa, linewidth = 1),  # Axis line color
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_line(color = ahus_mellomblaa, linewidth = 0.1)    # Minor grid line color
                      ))


png(file = 'ocd_betas_volume_poster.png',
    width = 30,
    height = 60,
    units = 'cm',
    res = 200)


  mutate(across(starts_with('beta_ocd'), exp)) %>%
  
  rename_variables('Amygdala' = 'beta_ocd[1]', 
                   'Hippocampus' = 'beta_ocd[2]',
                   'Nucleus accumbens(area)' = 'beta_ocd[3]',
                   'Putamen' = 'beta_ocd[4]',
                   'Lateral Ventricle' = 'beta_ocd[5]',
                   'Pallidum' = 'beta_ocd[6]',
                   'Caudate' = 'beta_ocd[7]',
                   'Thalamus' = 'beta_ocd[8]') %>%
  select('Thalamus', 'Amygdala', 'Nucleus accumbens(area)', 'Putamen', 'Amygdala', 'Caudate') %>%
  mcmc_areas() +
  scale_x_continuous(breaks = seq(0, 1.1, by = 0.01)) +
  labs(title = 'Multiplicative effect of OCD on subcortical volumes', 
       subtitle = 'Ajusted for age, gender and intracranial volume' )

dev.off()

png(file = 'eos_betas_volume.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)


  
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



png(file = 'ocd_betas_area.png',
    width = 20,
    height = 20,
    units = 'cm',
    res = 200)


  
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

fit_area$draws(variables = 'beta_eos', format = 'draws_df') %>%
  
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