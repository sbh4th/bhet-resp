#  program:  bhet-resp-analysis.R
#  task:     estimate models for respiratory outcomes
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2023-09-30


##  0 Load needed packages ----
library(here)
library(tidyverse)
library(tidybayes)
library(haven)
library(kableExtra)
library(modeldb)
library(brms)
library(cmdstanr)
library(modelr)
library(osfr)
library(modelsummary)
library(bayesplot)
library(patchwork)
library(marginaleffects)

# Use the cmdstanr backend for Stan
# You need to install the cmdstanr package first
# (https://mc-stan.org/cmdstanr/) and then run cmdstanr::install_cmdstan() to
# install cmdstan on your computer.
options(mc.cores = 4,
        brms.backend = "cmdstanr")

## download data files from OSF (data-clean component)
# dir.create("data-clean")
bhet_project <- osf_retrieve_node("b4wze")
bhet_project %>%
  osf_ls_files("Master Dataset (Season 4)",
               pattern = "dta") %>%
  osf_download(path = here("data-clean"),
    conflicts = "overwrite")

# set path to upload model fits to OSF
# code component
u2s_fits <- osf_retrieve_node("cv9qg")



## 1 Read in dataset, limit to resp vars ----
d <- read_dta(here("data-clean", 
                   "BHET_master_data_6Sep2023.dta"), 
  col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, 
                ban_status_2019, ban_status_2020, 
                ban_status_2021, ban_status_no, 
                ban_status_composite,
                freq_cough, freq_phlegm,
                freq_wheezing, freq_breath,
                freq_no_chest)) 

# define main outcome
d1 <- d %>% drop_na() %>%
  mutate(resp = if_else(
    freq_cough < 3 |
    freq_phlegm < 3 |
    freq_wheezing < 3 |
    freq_breath < 3 |
    freq_no_chest < 3, 1, 0),
    year = if_else(wave==1, 2018, 
      if_else(wave==2, 2019,
        if_else(wave==4, 2021, 0))),
    cohort_year = if_else(
      ban_status_composite==1, 2019, 
      if_else(ban_status_composite==2, 2020, 
              if_else(ban_status_composite==3, 2021, 2022))),
    treat = ifelse(year >= cohort_year, 1, 0),
    cohort_year = ifelse(cohort_year == 2022,-Inf, 
                         cohort_year)) %>%
  # relabel last cohort year 
  # treatment cohort dummies
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2019,2020,2021), 
    remove_original = F) %>%
  # wave dummies
  add_dummy_variables(year, 
    values=c(2018,2019,2021), remove_original = F)

# check the coding
vars <- c("freq_cough", "freq_phlegm",
  "freq_wheezing", "freq_breath",
  "freq_no_chest")

map(vars, ~ 
    d1 %>% 
      select(resp, starts_with("freq")) %>%
      group_by(resp, across(all_of(.x))) %>%
      tally() %>%
      kbl() %>% kable_styling()
   )

map(vars, ~ 
    d1 %>% 
      select(resp, starts_with("freq")) %>%
      group_by(resp, across(all_of(.x))) %>%
      tally() %>%
      kbl() %>% kable_styling()
   )

w

d1 %>% ggplot(aes(y = as_factor(freq_cough))) + 
  geom_bar(aes(x = (..count..)/sum(..count..))) +
  labs(y = "", x = "Proportion")




# overall prevalence
d1 %>% group_by(resp) %>% tally()

# limit sample

d2 <- d1 %>% 
  select(starts_with(c("year","cohort")), 
  "ID_VILLAGE","resp","treat") %>%
  
  # limit to complete cases
  drop_na() %>% 
  
  # create unique continuous village Id
  group_by(ID_VILLAGE) %>%
  mutate(v_id = cur_group_id()) %>%
  ungroup()

b1 <-
  brm(data = d2, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id),
      prior = c(prior(normal(0, 1), class = Intercept), 
                prior(exponential(1), class = sd)),        
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 287,
      file = "code/fits/bhet-resp-b1")

## check the chains
mcmc_trace(b1, pars="b_Intercept") +
  theme_classic() + ylab("")

## visualize prior draws on probability scale
prior_draws(b1) %>%
  # rescale to absolute probabilities
  mutate(p1_s = inv_logit_scaled(Intercept)) %>%
  ggplot(aes(x = p1_s)) + 
    stat_halfeye(color= "black", fill =  '#1b9e77') +
    labs(x = "Probability of poor respiratory health", 
      y = NULL, subtitle = "") +
  theme_classic()

b2 <-
  brm(data = d2, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + 
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 3975,
      file = "code/fits/bhet-resp-b2")

## check the chains
mcmc_trace(b2, pars=c("b_Intercept", "b_year_2021")) +
  theme_classic() + ylab("")

## extract prior draws
pd <- prior_draws(b2) %>%
  # rescale to absolute probabilities
  mutate(pd0 = inv_logit_scaled(Intercept),
         pd1 = inv_logit_scaled(Intercept + b)) %>%
  # treatment effect
  mutate(diff = pd1 - pd0)

## plot for intercept
plot_pd1 <- pd %>%
  ggplot(aes(x = pd1)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Probability of poor lung health", 
    y = NULL, subtitle = "Prior for Intercept") +
  theme_classic()

## plot for treatment effect
plot_pdd <- pd %>%
  ggplot(aes(x = diff)) + 
  stat_halfeye(color= "black", fill =  '#d95f02') +
  labs(x = "Difference in treated vs. control", 
    y = NULL, subtitle = "Prior for treatment effect") +
  theme_classic()

# Combined plot for priors
priors_lung <- (plot_pd1 | plot_pdd) +
  plot_annotation(
    title = "Priors for poor lung symptoms",
                  theme = theme_classic())

# Marginal predictions for Bayesian ETWFE (simple)

## load brms model
b2 <- readRDS(here("code/fits", 
  "bhet-resp-b2.rds"))

bme_pred <- predictions(
  b2, 
  newdata   = subset(d2, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# plot of predicted probabilities by treatment
bme_pred_p <- bme_pred |>
  posterior_draws() |>
  ggplot(aes(x = draw, fill=factor(treat))) +
    stat_halfeye(slab_alpha = .5) + 
    annotate("text", x = 0.57, y = 0.7, 
           label="Control", color='#1b9e77') +
    annotate("text", x = 0.49, y = 0.95, 
           label="Treated", color='#d95f02') +
  scale_x_continuous(
    "Probability of poor respiratory symptoms", 
    limits=c(0.4,0.7)) +
  scale_y_continuous("Posterior Density") +
  scale_fill_manual(values = c('#1b9e77','#d95f02')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=12))

bme_avg <- slopes(
  b2, 
  newdata   = subset(d2, treat==1),
  variables = "treat", 
  by        = "treat"
  ) 

# plot of treatment effect
bme_avg_p <- bme_avg |>
  posterior_draws() |>
  ggplot(aes(x = draw)) +
    stat_halfeye(slab_alpha = .5, fill = "#7570b3") +
    annotate("text", x = -0.075, y = 0.95, 
           label="Difference", color = '#7570b3') +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60") +
    scale_x_continuous("Marginal Effect", limits=c(-0.25,0.1)) +
    scale_y_continuous("") +
    theme_classic() + 
    theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=12))


f2 <- bme_pred_p + bme_avg_p + 
  plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Posterior distributions of marginal predictions: poor respiratory symptoms")
f2


# wrangle aggregate ATTs for model summary table
bmp <- data.frame(
  term = bme_pred$treat,
  estimate = bme_pred$estimate,
  conf.low = bme_pred$conf.low,
  conf.high = bme_pred$conf.high,
  std.error = abs(bme_pred$conf.high - 
                    bme_pred$conf.low) / (2 * 1.96)
)

bti <- data.frame(
  term = 2,
  estimate = bme_avg$estimate,
  conf.low = bme_avg$conf.low,
  conf.high = bme_avg$conf.high,
  std.error = abs(bme_avg$conf.high - 
                    bme_avg$conf.low) / (2 * 1.96)
)

bta <- bind_rows(bmp,bti) %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

gl <- data.frame()

betwfe_me_avg <- list(tidy = bta, glance = gl)
class(betwfe_me_avg) <- "modelsummary_list"

# Freq analysis
library(etwfe)
library(fixest)

fm <- glm(resp ~  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
        family = "binomial", data = d2)

library(sandwich)
library(lmtest)
clustered_se <- vcovCL(fm, cluster = d2$cohort_year)

etwfe = fixest::feols(
  resp ~ treat:i(cohort_year, i.year, ref=-Inf, ref2 = 2018) | cohort_year + year,
  data = d2,
  vcov = ~v_id
)

me_pred <- predictions(
  fm, 
  newdata   = subset(d2, treat==1),
  variables = "treat", 
  by        = "treat",
  vcov = "HC3"
  )

me_avg <- slopes(
  fm, 
  newdata   = subset(d2, treat==1),
  variables = "treat", 
  by        = "treat",
  vcov = "HC3"
  )

# wrangle aggregate ATTs for model summary table
fmp <- data.frame(
  term = me_pred$treat,
  estimate = me_pred$estimate,
  conf.low = me_pred$conf.low,
  conf.high = me_pred$conf.high,
  std.error = abs(me_pred$conf.high - 
                    me_pred$conf.low) / (2 * 1.96)
)

fti <- data.frame(
  term = 2,
  estimate = me_avg$estimate,
  conf.low = me_avg$conf.low,
  conf.high = me_avg$conf.high,
  std.error = abs(me_avg$conf.high - 
                    me_avg$conf.low) / (2 * 1.96)
)

fta <- bind_rows(fmp,fti) %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

gl <- data.frame()

fetwfe_me_avg <- list(tidy = fta, glance = gl)
class(fetwfe_me_avg) <- "modelsummary_list"

# wrangle aggregate ATTs for model summary table
ti <- data.frame(
  term = paste("ATT(", me_avg$term, ")", sep = ""),
  estimate = me_avg$estimate,
  std.error = me_avg$std.error,
  conf.low = me_avg$estimate - 1.96*me_avg$std.error,
  conf.high = me_avg$estimate + 1.96*me_avg$std.error)

gl <- data.frame()

etwfe_me_avg <- list(tidy = ti, glance = gl)
class(etwfe_me_avg) <- "modelsummary_list"

modelsummary(list("Simple Average" = etwfe_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')




# by cohort

bme_c <- slopes(
  b2, 
  newdata   = subset(d2, treat & cohort_year),
  variables = "treat", 
  by        = "cohort_year"
  )

bti_c <- data.frame(
  term = paste("ATT(g", bme_c$cohort_year, ")", sep=""), 
  estimate = bme_c$estimate,
  conf.low = bme_c$conf.low,
  conf.high = bme_c$conf.high,
  std.error = abs(bme_c$conf.high - 
                    bme_c$conf.low) / (2 * 1.96)
)

betwfe_me_c <- list(tidy = bti_c, glance = gl)
class(betwfe_me_c) <- "modelsummary_list"

modelsummary(list("Bayesian Cohort Average" = betwfe_me_c),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')



## first stack datasets with sample all treated (treated = 1)
## and all untreated (treated = 0)
mf_0 <- d2 %>% 
  mutate(treat = 0)
mf_1 <- d2 %>% 
  # select(school_3, blockf, treated, school_1, h_id) %>%
  mutate(treat = 1)
mf_me <- bind_rows(mf_0, mf_1)

# marginal predicted risks
mpp <- add_epred_draws(b2, newdata=mf_me,
    allow_new_levels = TRUE, seed = 294) %>%
    mutate(Tx = recode_factor(treat, 
      `0` = "Control", `1` = "Treated")) %>%
  group_by(Tx, .draw) %>%
  summarise(`Pr(y)` = mean(`.epred`))

# marginal predictions by treatment status
mpr <- mpp %>% group_by(Tx) %>%
  median_hdi(`Pr(y)`, .width = 0.95)
mpr

# marginal treatment effect
me <- mpp %>% 
  pivot_wider(names_from = Tx, values_from = `Pr(y)`) %>%
  mutate(diff = `Treated` - `Control`) %>%
  median_hdi(diff, .width = 0.95) 
me