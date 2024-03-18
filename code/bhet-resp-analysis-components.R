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

# model for coughing
te_cough <- 
  brm(data = dresp, 
      family = bernoulli(),
      cough ~ 1 + (1 | v_id) + 
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
      seed = 3957,
      file = "code/fits/bhet-resp-cough")

# load if model already exists
# te_cough <- readRDS("code/fits/bhet-resp-cough.rds")

## check the chains
mcmc_trace(te_cough, pars=c("b_Intercept", "b_year_2021")) +
  theme_classic() + ylab("")

## extract prior draws
pd <- prior_draws(te_cough) %>%
  # rescale to absolute probabilities
  mutate(pd0 = inv_logit_scaled(Intercept),
         pd1 = inv_logit_scaled(Intercept + b)) %>%
  # treatment effect
  mutate(diff = pd1 - pd0)

## plot for intercept prior
plot_pd1 <- pd %>%
  ggplot(aes(x = pd1)) + 
  stat_halfeye(color= "black", fill =  '#1b9e77') +
  labs(x = "Probability of cough", 
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



# model for phlegm
te_phlegm <- 
  brm(data = dresp, 
      family = bernoulli(),
      phlegm ~ 1 + (1 | v_id) + 
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
      seed = 3857,
      file = "code/fits/bhet-resp-phlegm")


# model for wheezing
te_wheeze <- 
  brm(data = dresp, 
      family = bernoulli(),
      wheeze ~ 1 + (1 | v_id) + 
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
      seed = 3857,
      file = "code/fits/bhet-resp-wheeze")

# model for breath
te_breath <- 
  brm(data = dresp, 
      family = bernoulli(),
      breath ~ 1 + (1 | v_id) + 
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
      seed = 3857,
      file = "code/fits/bhet-resp-breath")

# model for chest
te_nochest <- 
  brm(data = dresp, 
      family = bernoulli(),
      nochest ~ 1 + (1 | v_id) + 
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
      seed = 3857,
      file = "code/fits/bhet-resp-nochest")

te_csmoke <- 
  brm(data = dresp, 
      family = bernoulli(),
      csmoke ~ 1 + (1 | v_id) + 
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
      seed = 75,
      file = "code/fits/bhet-resp-csmoke")


te_resp_a <- 
  brm(data = dresp, 
      family = bernoulli(),
      resp ~ 1 + (1 | v_id) + 
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021 + 
        csmoke + fsmoke + age_health + male + temp,
      prior = c(prior(normal(0, 1), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes")

