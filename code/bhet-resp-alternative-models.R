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
          'MASS', 'nnet', 'modeldb',
          'tinytable')

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
      f <- factor(.x, ordered = TRUE,
                  levels = seq_along(labels_list[[var]]),
                  labels = labels_list[[var]])
    }
  ))


## 2 Function to run models, 
# estimate marginal effects ----

# Function to estimate ordered logit for multiple outcomes
estimate_ologit_did <- function(outcome_vars, predictor_vars) {
  models <- outcome_vars %>%
    set_names() %>%
    map(~ {
      formula <- as.formula(paste(.x, "~", 
        paste(predictor_vars, collapse = " + ")))
      polr(formula, data = dresp_cc, Hess = TRUE)
    })
  
  return(models)
}


## 3 Run models ----

# binary outcomes
b_ord_out <- c("freq_cough", "freq_phlegm", 
  "freq_wheezing", "freq_breath", "freq_no_chest")

# basic DiD specification with fixed effects
# for cohort and year
rhs_did <- c("treat:cohort_year_2019:year_2019",
             "treat:cohort_year_2019:year_2021", 
             "treat:cohort_year_2020:year_2021",
             "treat:cohort_year_2021:year_2021", "cohort_year_2019",
             "cohort_year_2020", "cohort_year_2021",
             "year_2019", "year_2021")

dresp_cc <- d2 %>%
  # limit to complete cases (ignoring bmi)
  drop_na(cresp:farm_4, -bmi)



## gather estimates across models ----
## and write results to outputs folder

# basic DiD
ologit_did <- estimate_ologit_did(b_ord_out, rhs_did)
write_rds(ologit_did, file = here(
  "outputs/ologit-did.rds"))

# estimate marginal predictions (simple average)
ologit_mp <- lapply(ologit_did, 
  marginaleffects::avg_predictions, 
  newdata = subset(dresp_cc, treat==1), 
  variables = "treat", by = "treat",
  # make sure to use cluster-robust SEs
  vcov = ~v_id)
write_rds(ologit_mp, file = here(
  "outputs/ologit-mp.rds"))

# estimate marginal effects (simple average)
# from the basic DiD
ologit_me <- lapply(ologit_did, 
  marginaleffects::slopes, 
  newdata = subset(dresp_cc, treat==1), 
  variables = "treat", by = "treat",
  # make sure to use cluster-robust SEs
  vcov = ~v_id)
write_rds(ologit_me, file = here(
  "outputs/logit-me.rds"))


# grab estimates and SEs from DiD results
ologit_tab1 <- bind_rows(ologit_mp,
  .id = "outcome") %>%
  mutate(
    category = group,
    treat = treat,
    est = estimate * 100,
    se = std.error * 100,
    ll = est - 1.96 * se,
    ul = est + 1.96 * se,
    ci = paste("(", sprintf('%.1f', ll), ", ",
    sprintf('%.1f', ul), ")", sep="")) %>%
  
  # keep relevant columns
  dplyr::select(outcome, category, 
    treat, est, ci) %>%
  
  # reshape exposure to wide
  pivot_wider(names_from = "treat",
    values_from = c("est", "ci"),
    names_vary = "slowest")


ologit_tab2 <- bind_rows(ologit_me,
  .id = "outcome") %>%
  mutate(
    category = group,
    est = estimate * 100,
    se = std.error * 100,
    ll = est - 1.96 * se,
    ul = est + 1.96 * se,
    ci = paste("(", sprintf('%.1f', ll), ", ",
      sprintf('%.1f', ul), ")", sep="")) %>%
  
    # keep relevant columns
  dplyr::select(outcome, category, 
    est, ci) 

ologit_table <- ologit_tab1 %>%
  left_join(ologit_tab2, join_by(outcome, category)) %>%
  mutate(across(c(est_0, est_1, est), ~ round(.x, 1))) %>%
  dplyr::select(-outcome)

colnames(ologit_table) <- c("Category", "%", "(95% CI)", 
  "%", "95% CI", "ATT (%)", "95% CI")

tt(ologit_table,
   notes = "Note: Average marginal predictions and ATT from ETWFE ordered logit models without covariates.") %>%
  group_tt(
    j = list("Treated" = 2:3, 
             "Untreated" = 4:5,
             "Difference" = 6:7),
    i = list("Frequency of coughing" = 1,
             "Frequency of phlegm" = 5,
             "Frequency of wheezing" = 9,
             "Frequency of trouble breathing" = 14,
             "Frequency of days with little chest trouble" = 17)) |>
  style_tt(i = c(1, 6, 11, 17, 21), align = "l", italic=T) 

  
