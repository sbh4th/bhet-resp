#  program:  bhet-resp-analysis.R
#  task:     mediation models for respiratory outcomes
#  input:    bhet-master; indoor temp files
#  output:   
#  project:  BHET
#  author:   sam harper \ 2024-01-18


##  0 Load needed packages ----
library(here)
library(tidyverse)
library(tidybayes)
library(haven)
library(readxl)
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


# load data

## 1 Read in dataset, limit to resp vars ----
dresp <- read_rds(here("data-clean", "bhet-resp-data.rds"))

# limit sample

d3 <- d2 %>% 
  select(starts_with(c("year","cohort")), 
  "ID_VILLAGE","resp","treat", "med_h",
  "min_h", "max_h", "ipm25", "district_2",
  "district_3", "district_4") %>%
  
  # limit to complete cases
  drop_na(-ipm25, -med_h, -min_h, -max_h) %>% 
  
  # create unique continuous village Id
  group_by(ID_VILLAGE) %>%
  mutate(v_id = cur_group_id()) %>%
  ungroup()

te_med <-
  brm(data = d3, 
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
      seed = 2847,
      file = "code/fits/bhet-resp-te")

# load if model already exists
te_med <- readRDS("code/fits/bhet-resp-te.rds")

## check the chains
mcmc_trace(te_med, pars=c("b_Intercept", "b_year_2021")) +
  theme_classic() + ylab("")

## extract prior draws
pd <- prior_draws(te_med) %>%
  # rescale to absolute probabilities
  mutate(pd0 = inv_logit_scaled(Intercept),
         pd1 = inv_logit_scaled(Intercept + b)) %>%
  # treatment effect
  mutate(diff = pd1 - pd0)

## plot for intercept prior
plot_pd1 <- pd %>%
  ggplot(aes(x = pd1)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Probability of poor lung health", 
    y = NULL, subtitle = "Prior for Intercept") +
  theme_classic()

## plot for treatment effect prior
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
#b2 <- readRDS(here("code/fits", 
 # "bhet-resp-b2.rds"))

bme_pred <- predictions(
  te_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# plot of predicted probabilities by treatment
bme_pred_p <- bme_pred |>
  posterior_draws() |>
  ggplot(aes(x = draw, fill=factor(treat))) +
    stat_halfeye(slab_alpha = .5) + 
    annotate("text", x = 0.57, y = 0.6, 
           label="Control", color='#1b9e77') +
    annotate("text", x = 0.50, y = 0.95, 
           label="Treated", color='#d95f02') +
  scale_x_continuous(
    "Probability of poor respiratory symptoms", 
    limits=c(0.4,0.8)) +
  scale_y_continuous("Posterior Density") +
  scale_fill_manual(values = c('#1b9e77','#d95f02')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=16))

bme_avg <- slopes(
  te_med, 
  newdata   = subset(d3, treat==1),
  variables = "treat", 
  by        = "treat"
  ) 

# plot of treatment effect
bme_avg_p <- bme_avg |>
  posterior_draws() |>
  ggplot(aes(x = draw)) +
    stat_halfeye(slab_alpha = .5, fill = "#7570b3") +
    annotate("text", x = -0.07, y = 0.95, 
           label="Difference", color = '#7570b3') +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60") +
    scale_x_continuous("Marginal Effect", limits=c(-0.3,0.1)) +
    scale_y_continuous("") +
    theme_classic() + 
    theme(legend.position = "none", 
        axis.text.y = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=16))


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

modelsummary(list("Bayesian Simple Average" = betwfe_me_avg),
  shape = term ~ model + statistic, 
  statistic = c("({std.error})", "conf.int"),
  gof_omit ='._*')


# by cohort
bme_c <- slopes(
  te_med, 
  newdata   = subset(d3, treat & cohort_year),
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


# mediation model
# model for outcome (including missing indoor PM2.5)
bf_resp <- bf(resp ~ 1 + (1 | v_id) + mi(ipm25) +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        treat:cohort_year_2019:year_2019:mi(ipm25) + 
        treat:cohort_year_2019:year_2021:mi(ipm25)+
        treat:cohort_year_2020:year_2021:mi(ipm25) +
        treat:cohort_year_2021:year_2021:mi(ipm25) +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
        family = bernoulli() )

# model for missing PM2.5
bf_pm <- bf(ipm25 | mi() ~ 1, family = lognormal())

# combined multivariate model
pmod <- bf_resp + bf_pm + set_rescor(FALSE)

# priors
# get_prior(data=d3, pmod)

# priors
priormed <- c(
  prior(normal(0, 1), class = Intercept, resp = resp),
  prior(normal(0, 1), class = b, resp = resp),
  prior(exponential(1), class = sd, resp = resp))

cde_pm_ln <- 
  brm(data = d3, 
      pmod,  # here we insert the model
      prior = priormed,
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed =257,
      file = "code/fits/cde_pm_ln")

# new data for estimated predictions
ndpm_med <- subset(d3, treat==1) %>%
  mutate(ipm25 = 55)

# counterfactual predictions (median PM)
bcde_pred_median <- predictions(
  cde_pm_ln, 
  newdata   = ndpm_med,
  variables = "treat", 
  by        = "treat"
  )

bcde_med <- bcde_pred_median %>% 
  filter(group == "resp") %>% 
  pivot_wider(names_from = treat, 
    values_from = c(estimate, conf.low, conf.high)) %>% 
  mutate(estimate_2 = estimate_1 - estimate_0, 
    se_0 = abs(conf.high_0 - conf.low_0) / (2 * 1.96), 
    se_1 = abs(conf.high_1 - conf.low_1) / (2 * 1.96)) %>% 
  mutate (se_2 = sqrt(se_0^2 + se_1^2),
          conf.low_2 = estimate_2 - 1.96 * se_2,
          conf.high_2 = estimate_2 + 1.96 * se_2) %>%
  pivot_longer(!group, names_to = c(".value", "treat"),
    names_pattern = "(.*)_(.)")

ndpm_q25 <- subset(d3, treat==1) %>%
  mutate(ipm25 = 32)

# counterfactual predictions (25%ile of PM)
bcde_pred_q25 <- predictions(
  cde_pm_ln, 
  newdata   = ndpm_q25,
  variables = "treat", 
  by        = "treat"
  )

bcde_q25 <- bcde_pred_q25 %>% 
  filter(group == "resp") %>% 
  pivot_wider(names_from = treat, 
    values_from = c(estimate, conf.low, conf.high)) %>% 
  mutate(estimate_2 = estimate_1 - estimate_0, 
    se_0 = abs(conf.high_0 - conf.low_0) / (2 * 1.96), 
    se_1 = abs(conf.high_1 - conf.low_1) / (2 * 1.96)) %>% 
  mutate (se_2 = sqrt(se_0^2 + se_1^2),
          conf.low_2 = estimate_2 - 1.96 * se_2,
          conf.high_2 = estimate_2 + 1.96 * se_2) %>%
  pivot_longer(!group, names_to = c(".value", "treat"),
    names_pattern = "(.*)_(.)")

  
# wrangle aggregate ATTs for model summary table
bmp_med <- data.frame(
  term = bcde_med$treat,
  estimate = bcde_med$estimate,
  conf.low = bcde_med$conf.low,
  conf.high = bcde_med$conf.high,
  std.error = bcde_med$se
)

bmp_med <- bmp_med %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

bmp_q25 <- data.frame(
  term = bcde_q25$treat,
  estimate = bcde_q25$estimate,
  conf.low = bcde_q25$conf.low,
  conf.high = bcde_q25$conf.high,
  std.error = bcde_q25$se
)

bmp_q25 <- bmp_q25 %>%
  mutate(term = recode_factor(term,
    `0` = "Untreated", `1` = "Treated",
    `2` = "Difference"))

gl <- data.frame()

bcde_me_med <- list(tidy = bmp_med, glance = gl)
class(bcde_me_med) <- "modelsummary_list"

bcde_me_q25 <- list(tidy = bmp_q25, glance = gl)
class(bcde_me_q25) <- "modelsummary_list"


# all estimates in one table
modelsummary(list("Total Effect" = betwfe_me_avg,
  "CDE (50th %ile)" = bcde_me_med,
  "CDE (25th %ile)" = bcde_me_q25),
  shape = term ~ model + statistic, 
  statistic = "conf.int",
  gof_omit ='._*')



# Supplementary analyses
# include product terms for cohort and time
bf_resp_full <- bf(resp ~ 1 + (1 | v_id) + 
        mi(ipm25) * (treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021), 
        family = bernoulli() )

# drop product terms for exposure-mediation interaction
bf_resp_np <- bf(resp ~ 1 + (1 | v_id) + 
        mi(ipm25) + treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021, 
        family = bernoulli() )

# combined multivariate model
# pmod_full <- bf_resp_full + bf_pm + set_rescor(FALSE)
pmod_np <- bf_resp_np + bf_pm + set_rescor(FALSE)

cde_pm_ln_full <- 
  brm(data = d3, 
      pmod_full,  # insert the model
      prior = priormed, # priors defined above
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed =3779,
      file = "code/fits/cde_pm_ln_full")

cde_pm_ln_np <- 
  brm(data = d3, 
      pmod_np,  # insert the model
      prior = priormed, # priors defined above
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed =3779,
      file = "code/fits/cde_pm_ln_np")



cde_temp_min <-
  brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + min_h +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        treat:cohort_year_2019:year_2019:min_h + 
        treat:cohort_year_2019:year_2021:min_h +
        treat:cohort_year_2020:year_2021:min_h +
        treat:cohort_year_2021:year_2021:min_h +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(0, 1), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 265,
      file = "code/fits/bhet-resp-cde_med_min")



# does the policy affect indoor PM?
bf_treat <- bf(ipm25_10 ~ 1 + (1 | v_id) +
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021)

bf_pm <- bf(ipm25_10 | mi() ~ 1)

bmodel <- bf_treat + bf_pm + set_rescor(FALSE)

te_pm <-
  brm(data = d3, 
      family = gaussian(),

      prior = c(prior(normal(0, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes", 
      seed = 265,
      file = "code/fits/te_temp_min")

# load if already run
te_temp_min <- readRDS("code/fits/te_temp_min.rds")

ndc <- subset(d3, treat==1) %>%
  mutate(mean_h=20)





# mediator model
m_model <- bf(min_h ~ 1 + (1 | v_id) +
                treat:cohort_year_2019:year_2019 + 
                treat:cohort_year_2019:year_2021 +
                treat:cohort_year_2020:year_2021 +
                treat:cohort_year_2021:year_2021 +
                cohort_year_2019 + cohort_year_2020 +
                cohort_year_2021 + year_2019 + year_2021) +
                gaussian()

# outcome model
y_model <- bf(resp ~ 1 + (1 | v_id) + min_h + 
                treat:cohort_year_2019:year_2019 + 
                treat:cohort_year_2019:year_2021 +
                treat:cohort_year_2020:year_2021 +
                treat:cohort_year_2021:year_2021 +
                cohort_year_2019 + cohort_year_2020 +
                cohort_year_2021 + year_2019 + year_2021) +
                bernoulli()

# priors
priormed <- c(
  prior(normal(10, 5), class = Intercept, resp = minh),
  prior(normal(0, 5), class = b, resp = minh),
  prior(exponential(1), class = sd, resp = minh),
  prior(normal(0, 1), class = Intercept, resp = resp),
  prior(normal(0, 1), class = b, resp = resp),
  prior(exponential(1), class = sd, resp = resp))

medfit <- brm(
  m_model + y_model + set_rescor(FALSE),
  data = d3, iter = 2000, warmup = 1000,
  prior = priormed,
  chains = 4, cores = 4,
  sample_prior = "yes")


brm(data = d3, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + min_h +
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
      seed = 265,
      file = "code/fits/bhet-resp-cde_med_min_ni")


x <- seq(0.01, 10, length.out = 1000)
mu <- 0
sigma <- 0.5
y <- dlnorm(x, meanlog = mu, sdlog = sigma)



te_wt <-
  brm(data = dh, 
      family = gaussian(),
      weight ~ 1 + (1 | v_id) + height +  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(65, 10), class = Intercept),
        prior(normal(0, 5), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes")

