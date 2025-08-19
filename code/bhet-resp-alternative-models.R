#  program:  bhet-resp-alternative-models.R
#  task:     alternative models for self-reported outcomes
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2025-08-19

## 0  packages to load
pkgs <- c('here', 'tidyverse', 'modelsummary', 
          'fixest', 'marginaleffects',
          'patchwork', 'estimatr',
          'MASS', 'nnet', 'modeldb')

#load all packages at once
lapply(pkgs, library, character.only=TRUE)


## 1 create analytic dataset
d <- read_csv(here("data-clean", 
                   "BHET_master_data_04Oct2024.csv"), 
              col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, ID_COUNTY, 
                            ban_status_2019, ban_status_2020, 
                            ban_status_2021, ban_status_no, 
                            ban_status_composite,
                            freq_cough, freq_phlegm,
                            freq_wheezing, freq_breath,
                            freq_no_chest, height, weight, 
                            age_health, gender_health,
                            smoking, temp, p_usable_pm, 
                            PM25conc_exposureugm3,
                            lived_with_smoker, years_with_smoker,
                            height, weight, occupation,
                            freq_farming, freq_exercising,
                            freq_drink))

d2 <- d %>% mutate(
  # total symptoms
  cresp = rowSums(across(c(freq_breath, freq_cough, 
                           freq_phlegm, freq_wheezing, freq_no_chest))),
  
  # separate symptoms
  cough = if_else(freq_cough < 3, 1, 0),
  phlegm = if_else(freq_phlegm < 3, 1, 0),
  wheeze = if_else(freq_wheezing < 3, 1, 0),
  breath = if_else(freq_breath < 3, 1, 0),
  nochest = if_else(freq_no_chest < 3, 1, 0),
  resp = if_else(
    freq_cough < 3 |
      freq_phlegm < 3 |
      freq_wheezing < 3 |
      freq_breath < 3 |
      freq_no_chest < 3, 1, 0),
  nosym = if_else(freq_cough==4 & freq_phlegm==4 & 
    freq_no_chest==4 & freq_breath==3 & 
    freq_wheezing==5, 1, 0),
  
  # change symptoms to factors
#  across(c("freq_breath", "freq_cough", "freq_phlegm",
#    "freq_wheezing", "freq_no_chest"), ~ factor(.x, ordered = TRUE)),
  
  # covariates
  male = if_else(gender_health == 1, 1, 0),
  # csmoke = if_else(smoking == 1, 1, 0),
  # fsmoke = if_else(smoking == 2, 1, 0),
  bmi = weight / (height/100)^2,
  
  occ = case_when(
    occupation == 1 ~ "Agriculture",
    occupation %in% c(2,11,12) ~ "Manual",
    occupation %in% c(3:5) ~ "Non-manual",
    occupation %in% c(6:10) ~ "Other"),
  # add dummies for occupation
  occ_2 = if_else(
    occ=="Manual", 1, 0),
  occ_3 = if_else(
    occ=="Non-manual", 1, 0),
  occ_4 = if_else(
    occ=="Other", 1, 0),
  
  ets = case_when(
    smoking == 1 ~ "Current smoker",
    smoking == 2 ~ "Former smoker",
    smoking == 3 & lived_with_smoker %in% c(2,3) ~ 
      "Never smoker lived with smoker",
    smoking == 3 & lived_with_smoker == 1 ~ 
      "No smoking exposure"),
  # add dummies for smoking
  ets_2 = if_else(
    ets=="Former smoker", 1, 0),
  ets_3 = if_else(
    ets=="Never smoker lived with smoker", 1, 0),
  ets_4 = if_else(
    ets=="No smoking exposure", 1, 0),
  
  drink = case_when(
    freq_drink == 1 ~ "Never",
    freq_drink %in% c(2:5) ~ "Occasional",
    freq_drink %in% c(6:8) ~ "Regular",
    freq_drink == 9 ~ "Everyday"),
  # add dummies for drinking
  drink_2 = if_else(
    drink=="Occasional", 1, 0),
  drink_3 = if_else(
    drink=="Regular", 1, 0),
  drink_4 = if_else(
    drink=="Everyday", 1, 0),
  
  farm = case_when(
    freq_farming == 1 ~ "Never",
    freq_farming %in% c(2:3) ~ "Occasional",
    freq_farming == 4 ~ "Regular",
    freq_farming == 5 ~ "Everyday"),
  # add dummies for farming
  farm_2 = if_else(
    farm=="Occasional", 1, 0),
  farm_3 = if_else(
    farm=="Regular", 1, 0),
  farm_4 = if_else(
    farm=="Everyday", 1, 0)) %>%
  
  # treatment related
  mutate(year = if_else(wave==1, 2018, 
                        if_else(wave==2, 2019,
                                if_else(wave==4, 2021, 0))),
         cohort_year = if_else(
           ban_status_composite==1, 2019, 
           if_else(ban_status_composite==2, 2020, 
                   if_else(ban_status_composite==3, 2021, 2022))),
         treat = ifelse(year >= cohort_year, 1, 0),
         cohort_year = ifelse(cohort_year == 2022,-Inf, 
                              cohort_year),
         ppm25 = case_when(
           PM25conc_exposureugm3 <= 0 ~ NA,
           p_usable_pm == 1 ~ PM25conc_exposureugm3,
           p_usable_pm== 0 ~ NA)) %>%
  rename(district = ID_COUNTY) %>%
  
  # add dummies for district
  mutate(
    district_2 = if_else(
      district==24576, 1, 0),
    district_3 = if_else(
      district==38376, 1, 0),
    district_4 = if_else(
      district==23494, 1, 0)) %>%
  # relabel last cohort year 
  # treatment cohort dummies
  add_dummy_variables(cohort_year, 
                      values=c(-Inf,2019,2020,2021), 
                      remove_original = F) %>%
  # wave dummies
  add_dummy_variables(year, 
                      values=c(2018,2019,2021), remove_original = F) %>%
  # create unique continuous village Id
  group_by(ID_VILLAGE) %>%
  mutate(v_id = cur_group_id()) %>%
  ungroup()

## Ordinal

# Define labels
labels_list <- list(
  freq_cough    = c("Most days a week",
                    "Several days a week",
                    "Only with a chest infection",
                    "Not at all"),
  
  freq_phlegm   = c("Most days a week",
                    "Several days a week",
                    "Only with a chest infection",
                    "Not at all"),
  
  freq_wheezing = c("Most days a week",
                    "Several days a week",
                    "A few days a month",
                    "Only with a chest infection",
                    "Not at all"),
  
  freq_breath   = c("Most days a week",
                    "Several days a week",
                    "Not at all"),
  
  freq_no_chest = c("No good days",
                    "A few good days",
                    "Most days are good",
                    "Every day is good")
)

d2 <- d2 %>%
  mutate(across(
    all_of(names(labels_list)),
    ~ {
      var <- cur_column()
      f <- factor(.x,
                  levels = seq_along(labels_list[[var]]),
                  labels = labels_list[[var]])
    }
  ))

# set reference levels
ref_levels <- list(
  freq_cough    = 4,
  freq_phlegm   = 4,
  freq_no_chest = 4,
  freq_breath   = 3,
  freq_wheezing = 5
)

d2 <- d2 %>%
  mutate(
    across(
      all_of(names(ref_levels)),
      ~ relevel(as.factor(.x, ordered = TRUE), 
                ref = ref_levels[[cur_column()]])
    )
  )

o_cough <- polr(
  freq_cough ~ treat:cohort_year_2019:year_2019 + 
    treat:cohort_year_2019:year_2021 + 
    treat:cohort_year_2020:year_2021 + 
    treat:cohort_year_2021:year_2021 + 
    cohort_year_2019 + cohort_year_2020 + 
    cohort_year_2021 + year_2019 + year_2021, 
  data = d2, Hess = TRUE)

o_cough_mp <- marginaleffects::avg_predictions(
  o_cough, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

o_cough_me <- marginaleffects::avg_slopes(
  o_cough, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)



o_phlegm <- polr(
  freq_phlegm ~ treat:cohort_year_2019:year_2019 + 
    treat:cohort_year_2019:year_2021 + 
    treat:cohort_year_2020:year_2021 + 
    treat:cohort_year_2021:year_2021 + 
    cohort_year_2019 + cohort_year_2020 + 
    cohort_year_2021 + year_2019 + year_2021, 
  data = d2, Hess = TRUE)

o_phlegm_mp <- marginaleffects::avg_predictions(
  o_phlegm, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

o_phlegm_me <- marginaleffects::avg_slopes(
  o_phlegm, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

o_wheeze <- polr(
  freq_wheezing ~ treat:cohort_year_2019:year_2019 + 
    treat:cohort_year_2019:year_2021 + 
    treat:cohort_year_2020:year_2021 + 
    treat:cohort_year_2021:year_2021 + 
    cohort_year_2019 + cohort_year_2020 + 
    cohort_year_2021 + year_2019 + year_2021, 
  data = d2, Hess = TRUE)

o_wheeze_mp <- marginaleffects::avg_predictions(
  o_wheeze, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

o_wheeze_me <- marginaleffects::avg_slopes(
  o_wheeze, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

mn_wheeze <- multinom(
  freq_wheezing ~ treat:cohort_year_2019:year_2019 + treat:cohort_year_2019:year_2021 + treat:cohort_year_2020:year_2021 + treat:cohort_year_2021:year_2021 + cohort_year_2019 + cohort_year_2020 + cohort_year_2021 + year_2019 + year_2021, data = d2, Hess = TRUE)


o_nochest <- polr(
  freq_no_chest ~ treat:cohort_year_2019:year_2019 + 
    treat:cohort_year_2019:year_2021 + 
    treat:cohort_year_2020:year_2021 + 
    treat:cohort_year_2021:year_2021 + 
    cohort_year_2019 + cohort_year_2020 + 
    cohort_year_2021 + year_2019 + year_2021, 
  data = d2, Hess = TRUE)

o_nochest_mp <- marginaleffects::avg_predictions(
  o_nochest, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)

o_nochest_me <- marginaleffects::avg_slopes(
  o_nochest, 
  newdata = subset(d2, treat==1), 
  var = "treat", vcov = ~v_id)


mn_breath <- multinom(
  freq_breath ~ treat:cohort_year_2019:year_2019 + treat:cohort_year_2019:year_2021 + treat:cohort_year_2020:year_2021 + treat:cohort_year_2021:year_2021 + cohort_year_2019 + cohort_year_2020 + cohort_year_2021 + year_2019 + year_2021, data = d2, Hess = TRUE)

o_breath <- polr(
  freq_breath ~ treat:cohort_year_2019:year_2019 + treat:cohort_year_2019:year_2021 + treat:cohort_year_2020:year_2021 + treat:cohort_year_2021:year_2021 + cohort_year_2019 + cohort_year_2020 + cohort_year_2021 + year_2019 + year_2021, data = d2, Hess = TRUE)
