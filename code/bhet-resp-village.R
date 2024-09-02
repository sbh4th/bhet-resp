#  program:  bhet-resp-village.R
#  task:     village-level analysis
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2024-08-28


##  0 Load needed packages ----
library(here)
library(tidyverse)
library(haven)
library(modeldb)
library(modelr)
library(osfr)
library(modelsummary)
library(tinytable)

## 1 Read in dataset, limit to resp vars ----
dv <- read_dta(here("data-clean", 
                   "BHET_master_data_22Jul2024.dta"), 
  col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, ID_COUNTY, 
                ban_status_2019, ban_status_2020, 
                ban_status_2021, ban_status_no, 
                ban_status_composite,
                freq_cough, freq_phlegm,
                freq_wheezing, freq_breath,
                freq_no_chest)) 

# define main outcomes and covariates
dv1 <- dv %>%
  mutate(
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
    freq_no_chest < 3, 1, 0)) %>%
  
    # limit to complete cases across outcomes
    drop_na(resp, cough, phlegm, 
      wheeze, breath, nochest) %>%

    # create unique continuous village Id
    group_by(ID_VILLAGE) %>%
    mutate(v_id = cur_group_id()) %>%
    ungroup()
    
# aggregate to village level
dv2 <- dv1 %>% 
  group_by(v_id, wave) %>% 
    summarise(
      across(c("resp", "cough", 
        "phlegm", "wheeze", "breath", "nochest"), 
        ~ sum(.x, na.rm = TRUE)), 
      pop = n(), 
      ban = max(ban_status_composite)) %>%
  
    # treatment related variables
    mutate(
      year = if_else(wave==1, 2018, 
      if_else(wave==2, 2019,
        if_else(wave==4, 2021, 0))),
    cohort_year = if_else(
      ban==1, 2019, 
      if_else(ban==2, 2020, 
              if_else(ban==3, 2021, 2022))),
    treat = ifelse(year >= cohort_year, 1, 0),
    cohort_year = ifelse(cohort_year == 2022,-Inf, 
                         cohort_year)) %>% 
  
  # treatment cohort dummies
  add_dummy_variables(cohort_year, 
    values=c(-Inf,2019,2020,2021), 
    remove_original = F) %>%
  
  
  # wave dummies
  add_dummy_variables(year, 
    values=c(2018,2019,2021), remove_original = F) 

## 3 village level analysis ----

# Function to estimate GLM for multiple outcomes
estimate_v_did <- function(outcome_vars, predictor_vars) {
  models <- outcome_vars %>%
    set_names() %>%
    map(~ {
      formula <- as.formula(
        paste0("cbind(", .x, ", pop - ", .x, ") ~", 
        paste(predictor_vars, collapse = " + ")))
      glm(formula, data = dv2, weights = pop,
          family = "binomial")
    })
  
  return(models)
}

# binary outcomes
b_out <- c("resp", "cough", "phlegm", "wheeze",
              "breath", "nochest")

# basic DiD specification covariates
rhs_did <- c("treat:cohort_year_2019:year_2019",
  "treat:cohort_year_2019:year_2021", 
  "treat:cohort_year_2020:year_2021",
  "treat:cohort_year_2021:year_2021", "cohort_year_2019",
  "cohort_year_2020", "cohort_year_2021",
  "year_2019", "year_2021")

# gather estimates across models
v_did <- estimate_v_did(b_out, rhs_did)

# estimate marginal effects (simple average)
v_did_me <- lapply(v_did, 
  marginaleffects::slopes, 
  newdata = subset(dv2, treat==1), 
  variables = "treat", 
  by = "treat",
  # heteroskedasticity robust standard errors
  vcov = "HC0")

write_rds(v_did_me, file = here(
  "outputs/v-did-me.rds"))

# grab estimates and SEs from DiD results
did_t1 <- v_did_me %>% {
  tibble(
    # table = 1,
    outcome = names(v_did_me),
    est = map_dbl(., "estimate") * 100,
    stderror = map_dbl(., "std.error") * 100,
    ll = est - 1.96 * stderror,
    ul = est + 1.96 * stderror,
    ci = paste("(", sprintf("%.1f", ll), ", ",
      sprintf("%.1f", ul), ")", sep="")
  )
}

# grab number of observations
did_obs <- v_did %>% {
  tibble(
    outcome = names(v_did),
    nobs = map_int(., "df.null") + 1
  )
}

# add observations to DiD results
did_t1 <- did_t1 %>%
  inner_join(did_obs) %>%
  select(-outcome) %>%
  relocate(nobs)

# grab estimates and SEs from 
# individual DiD results
me_i <- read_rds(logit_me, file = here(
  "outputs/logit-me.rds"))

did_t2 <- me_i %>% {
  tibble(
    # table = 2,
    outcome = names(logit_me),
    esti = map_dbl(., "estimate") * 100,
    stderrori = map_dbl(., "std.error") * 100,
    lli = esti - 1.96 * stderrori,
    uli = esti + 1.96 * stderrori,
    cii = paste("(", sprintf("%.1f", lli), ", ",
      sprintf("%.1f", uli), ")", sep="")
  )
}

logit_i <- read_rds(logit_did, file = here(
  "outputs/logit-did.rds"))

# grab number of observations
did_obs2 <- logit_i %>% {
  tibble(
    outcome = names(logit_i),
    nobsi = map_int(., "df.null") + 1
  )
}

did_t2 <- did_t2 %>%
  inner_join(did_obs2) %>%
  relocate(outcome, nobsi)

didt <- cbind(did_t2, did_t1) %>%
  mutate(dest = esti - est,
    dse = sqrt(stderrori^2 + stderror^2),
    dci = paste("(", sprintf("%.1f", dest - 1.96 * dse), ", ",
      sprintf("%.1f", dest + 1.96 * dse), ")", sep="")) %>%
  dplyr::select(nobsi, esti, cii, 
    nobs, est, ci, dest, dci) %>%
  mutate(outcome = c("Any symptom", "Coughing",
    "Phlegm", "Wheezing attacks", "Trouble breathing",
    "Chest trouble")) %>%
  # rename for combining with other tables
  # rename(`obs_1` = "nobsi", `estimate_1` = "esti", 
  #  `ci_1` = "cii", `obs_2` = "nobs", 
  #  `estimate_2` = "est", `ci_2` = "ci") %>%
  relocate(outcome)

write_rds(didt, file = here("outputs",
  "did-compare-v-i.rds"))

# code for table
colnames(didt) <- c(" ", "Obs", "ATT", 
  "(95% CI)", "Obs", "ATT", "(95% CI)",
  "ATT", "(95% CI)")

tt(didt,
   digits = 2,
   #width = c(3.5, 3, 1, 0.5, 2, 0.5, 2, 0.5, 2, 0.5, 2),
   notes = list("Note: ATT = Average Treatment Effect on the Treated, CI = confidence interval, Obs = observations", a = list(i=0, j=3, text = "Standard errors clustered by village"), b = list(i=0, j=6, text = "Village-level models weighted by cluster size with robust standard errors"))) %>%
  group_tt(
    j = list("Individual-level" = 2:4, 
             "Village-level" = 5:7,
             "Difference" = 8:9)) %>%
  style_tt(j = 1:9, align = "lcccccccc") %>%
  format_tt(escape = TRUE) %>%
  format_tt(j=c(3,6,8), sprintf = "%.1f") 

