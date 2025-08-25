#  program:  bhet-resp-event-study-plot.R
#  task:     alternative models for self-reported outcomes
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2025-08-21

## 0  packages to load
pkgs <- c('here', 'tidyverse', 'modelsummary', 
          'fixest', 'marginaleffects',
          'patchwork', 'estimatr',
          'tinytable')

#load all packages at once
lapply(pkgs, library, character.only=TRUE)


## 1 create analytic dataset

# data
esd <- read_rds(here("data-clean", 
                     "bhet-resp-data.rds")) %>%
  drop_na(cresp:farm_4, -bmi) %>%
  mutate(
    # add event indicators
    event = year - cohort_year,
    # add indicator for ever treated
    evtreat = if_else(cohort_year!=-Inf,1,0))

## 2 logit model for any symptom

# GLM model
esm <- glm(
  resp ~ treat:cohort_year_2019:year_2019 +
    treat:cohort_year_2019:year_2021 + 
    treat:cohort_year_2020:year_2019 +
    treat:cohort_year_2020:year_2021 +
    treat:cohort_year_2021:year_2019 +
    treat:cohort_year_2021:year_2021 + 
    cohort_year_2019 + cohort_year_2020 + cohort_year_2021 + 
    year_2019 + year_2021,
  data = esd, family = "binomial")

# Marginal effects for ETWFE
# aggregates ATTs by time since treatment
me_esm <- slopes(
  esm, 
  newdata   = subset(esd, treat & event>=0),
  variables = "treat", 
  by        = "event",
  # make sure to use cluster-robust SEs
  vcov      = ~v_id)


## 3 Estimate the model for pre-periods

# GLM for any symptom
esm_pre <- glm(
  resp ~ evtreat:cohort_year_2020:year_2019 +
    evtreat:cohort_year_2021:year_2019 + 
    cohort_year_2019 + cohort_year_2020 + 
    cohort_year_2021 + year_2019 + year_2021,
  # restrict to only untreated periods
  data = subset(esd, treat==0), family = "binomial")

# marginal effects for pre-periods
me_esm_pre <- 
  slopes(esm_pre, 
         variables = "evtreat", 
         newdata = subset(esd, treat==0 & 
                            event > -3 & event < 0), 
         by = "event",
         # make sure to use cluster-robust SEs
         vcov = ~v_id)

# can test here whether all pre-treatment 
# effects are zero
hypotheses(me_esm_pre, joint = TRUE)

## 4 make the plot

# put pre-treatment and dynamic effects together
ete <- bind_rows(me_esm_pre, me_esm)

# make the plot
ete %>% 
  mutate(`Intervention time` = 
           if_else(event < 0, "Pre", "Post")) %>% 
  ggplot(aes(x = event, y = estimate, ymin = conf.low, 
             ymax = conf.high, color = `Intervention time`)) +
  geom_hline(yintercept = 0) +
  geom_pointrange() + theme_classic() +
  scale_color_brewer(palette="Set1") +
  labs(x = "Time period relative to treatment",
       y = "Average ATT")
