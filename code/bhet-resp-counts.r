#  program:  bhet-resp-analysis.R
#  task:     estimate models for respiratory outcomes
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2024-03-04


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

# bhet_project <- osf_retrieve_node("b4wze")
# set path to upload model fits to OSF
# code component
#
# u2s_fits <- osf_retrieve_node("cv9qg")



## 1 Read in dataset, limit to resp vars ----
dc <- read_rds(here("data-clean", "bhet-resp-data.rds"))

# empty model to evaluate priors

set.seed(2479)

prior <-
  tibble(i = 1:n,
         a = rnorm(n, mean = 3, sd = 0.5),
         b = rnorm(n, mean = 0, sd = 0.2)) %>% 
  expand_grid(x = seq(from = log(100), to = log(200000), length.out = 100))

# left
p1 <-
  prior %>% 
  ggplot(aes(x = x, y = exp(a + b * x), group = i)) +
  geom_line(linewidth = 1/4, alpha = 2/3,
            color = wes_palette("Moonrise2")[4]) +
  labs(subtitle = expression(beta%~%Normal(0*', '*0.2)),
       x = "log population",
       y = "total tools") +
  coord_cartesian(xlim = c(log(100), log(200000)),
                  ylim = c(0, 500))
# right
p2 <-
  prior %>% 
  ggplot(aes(x = exp(x), y = exp(a + b * x), group = i)) +
  geom_line(linewidth = 1/4, alpha = 2/3,
            color = wes_palette("Moonrise2")[4]) +
  labs(subtitle = expression(beta%~%Normal(0*', '*0.2)),
       x = "population",
       y = "total tools") +
  coord_cartesian(xlim = c(100, 200000),
                  ylim = c(0, 500))

# combine
p1 | p2




bc1 <-
  brm(data = dc,
      family = poisson(),
      cresp | trunc(lb = 5, ub = 20) ~ 1 + (1 | v_id),
      prior = c(prior(normal(2, 0.5), class = Intercept),
                prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 2037,
      file = "code/fits/bhet-cresp-bc1")


## check the chains
mcmc_trace(bc1, pars="b_Intercept") +
  theme_classic() + ylab("")

## visualize prior draws on probability scale
prior_draws(bc1_bounds) %>%
  # rescale to absolute probabilities
  mutate(count = exp(Intercept)) %>%
  ggplot(aes(x = count)) + 
    stat_halfeye(color= "black", fill =  '#1b9e77') +
    labs(x = "Probability of poor respiratory health", 
      y = NULL, subtitle = "") +
  theme_classic()

bc2 <-
  brm(data = dc,
      family = poisson(),
      cresp | trunc(lb = 5, ub = 20) ~ 1 + (1 | v_id) + 
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
      prior = c(prior(normal(2, 0.5), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      sample_prior = "yes",
      seed = 4867,
      file = "code/fits/bhet-cresp-bc2")

## check the chains
mcmc_trace(bc2, pars=c("b_Intercept", "b_year_2021")) +
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
bc2 <- readRDS(here("code/fits", 
  "bhet-resp-bc2.rds"))

bme_pred <- predictions(
  bc2, 
  newdata   = subset(dc, treat==1),
  variables = "treat", 
  by        = "treat"
  )

# plot of predicted probabilities by treatment
bme_pred_p <- bme_pred |>
  posterior_draws() |>
  ggplot(aes(x = draw, fill=factor(treat))) +
    stat_halfeye(slab_alpha = .5) + 
    annotate("text", x = 16.6, y = 0.6, 
           label="Control", color='#25681A') +
    annotate("text", x = 17, y = 1, 
           label="Treated", color='#830223') +
  scale_x_continuous(
    "Respiratory symptom scale (higher = better)", 
    limits=c(15.5, 18)) +
  scale_y_continuous("Posterior Density") +
  scale_fill_manual(values = c('#25681A','#830223')) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=12))

bme_avg <- slopes(
  bc2, 
  newdata   = subset(dc, treat==1),
  variables = "treat", 
  by        = "treat"
  ) 

# plot of treatment effect
bme_avg_p <- bme_avg |>
  posterior_draws() |>
  ggplot(aes(x = draw)) +
    stat_halfeye(slab_alpha = .5, fill = "#7570b3") +
    annotate("text", x = 0.4, y = 0.95, 
           label="Difference", color = '#7570b3') +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray60") +
    scale_x_continuous("Marginal Effect", limits=c(-1,2)) +
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


dv <- dc %>% group_by(v_id, wave) %>% summarize(y = sum(cough), ban_status = mean(ban_status_composite), pop = n()) %>%
  mutate(    
    year = if_else(wave==1, 2018, 
      if_else(wave==2, 2019,
        if_else(wave==4, 2021, 0))),
    cohort_year = if_else(
      ban_status==1, 2019, 
      if_else(ban_status==2, 2020, 
              if_else(ban_status==3, 2021, 2022))),
    treat = ifelse(year >= cohort_year, 1, 0),
    cohort_year = ifelse(cohort_year == 2022,-Inf, 
                         cohort_year),
    lnpop = log(pop)) %>%
  # relabel last cohort year 
  # treatment cohort dummies
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2019,2020,2021), 
    remove_original = F) %>%
  # wave dummies
  add_dummy_variables(year, 
    values=c(2018,2019,2021), remove_original = F) 

pm <- glm(y ~  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021, 
        offset = lnpop,
        family = "poisson", data=dv)

pm_pred <- predictions(
  pm,
  type = "response",
  newdata   = me_data,
  variables = "treat", 
  by        = "treat",
  vcov = "HC3"
  )

pmr <- glm(yrate ~  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021,
        family = "gaussian", data=dv)

pml <- glm(y ~  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021 + 
          lnpop,
        family = "poisson", data=dv)


df <- df %>% mutate(lnpop = log(population))
library(RStata)
options("RStata.StataVersion" = 16)
options("RStata.StataPath"= '/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp')

m1 <- glm(deaths ~ gender + age_group + year, offset=lnpop, data = df)

# Bootstrap to get standard error of RR
bsme <- function(splits) {
  x <- analysis(splits)
  model <- glm(y ~  
        treat:cohort_year_2019:year_2019 + 
        treat:cohort_year_2019:year_2021 +
        treat:cohort_year_2020:year_2021 +
        treat:cohort_year_2021:year_2021 +
        cohort_year_2019 + cohort_year_2020 +
        cohort_year_2021 + year_2019 + year_2021, 
        offset = lnpop,
        family = "poisson", data=x)
  me_0 <- x %>% filter(treat==1) %>%
    mutate(treat=0, lnpop=0)
  me_1 <- x %>% filter(treat==1) %>%
    mutate(treat=1, lnpop=0)
  amp_0 <- mean(predict(model, newdata = me_0, 
                   type = "response"))
  amp_1 <- mean(predict(model, newdata = me_1, 
                   type = "response"))
  amp_1 - amp_0
}

set.seed(3846)
bs_samples <- rsample::bootstraps(dv, times = 1000)
# iterate over each bootstrap sample and compute statistic
bs_samples$ame <- map_dbl(bs_samples$splits, bsme)
quantile(bs_samples$ame, probs = c(.05, .5, .95))










# set treated pop to untreated, offset to 0
mf_0 <- dv %>% filter(treat==1) %>%
    mutate(treat=0, lnpop=0)
# set treated pop to treated, offset to 0
mf_1 <- dv %>% filter(treat==1) %>%
    mutate(treat=1, lnpop=0)
# combine data
me_data <- bind_rows(mf_0, mf_1)
me_data %>% 
  add_predictions(pm, type = 'response') %>% 
  group_by(treat) %>% 
  pivot_wider(names_from = "treat", 
    names_prefix = "t", values_from = pred) %>% 
  mutate(me = t1 - t0) %>% 
  summarize(ame = mean(me), amp0 = mean(t0),
            amp1 = mean(t1))


me_data %>%
  add_predictions(pm_fe, type = 'response') %>%
  group_by(treat) %>%
  pivot_wider(names_from = "treat",
    names_prefix = "t", values_from = pred) %>%
  mutate(me = t1 - t0) %>%
  summarize(ame = mean(me), amp0 = mean(t0),
            amp1 = mean(t1))


set.seed(13)

x <- rnorm(100, sd = 0.1)
y <- rpois(100, exp(5 * x))

e <- rpois(100, 5) + 1       # this is your offset/exposure
y_weighted <- y / e          # weighting by your offset/exposure

### --- Using offset(.)

mod_1 <- glm(y ~ x + offset(log(e)), family = 'poisson')

### --- Using the weighted outcome

mod_2 <- glm(y_weighted ~ x, family = 'poisson', weights = e)

mod_3 <- glm(y_weighted ~ x, family = 'poisson')

round(mod_1$coefficients, 3)
round(mod_2$coefficients, 3)
round(mod_3$coefficients, 3)
