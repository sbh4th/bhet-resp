#  program:  bhet-resp-analysis.R
#  task:     estimate models for respiratory outcomes
#  input:    bhet-master
#  output:   
#  project:  BHET
#  author:   sam harper \ 2024-01-17


##  0 Load needed packages ----
library(here)
library(tidyverse)
library(haven)
library(kableExtra)
library(modeldb)
library(modelr)
library(osfr)
library(modelsummary)


## download data files from OSF (data-clean component)
# dir.create("data-clean")
bhet_project <- osf_retrieve_node("vxur5")
bhet_project %>%
  osf_ls_files("Master Dataset (Seasons 1-4)",
               pattern = "dta") %>%
  osf_download(path = here("data-clean"),
    conflicts = "overwrite")
#
## set path to upload model fits to OSF
## code component
#u2s_fits <- osf_retrieve_node("cv9qg")



## 1 Read in dataset, limit to resp vars ----
d <- read_dta(here("data-clean", 
                   "BHET_master_data_22Jul2024.dta"), 
  col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, ID_COUNTY, 
                ban_status_2019, ban_status_2020, 
                ban_status_2021, ban_status_no, 
                ban_status_composite,
                freq_cough, freq_phlegm,
                freq_wheezing, freq_breath,
                freq_no_chest, height, weight, 
                age_health, gender_health,
                smoking, temp, PM25_exp_remove, 
                PM25conc_exposureugm3)) 

# define main outcome
d1 <- d %>%
  mutate(
    # total symptoms
    cresp = rowSums(across(c(freq_breath, freq_cough, 
      freq_phlegm, freq_wheezing, freq_no_chest))),
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
    male = if_else(gender_health == 1, 1, 0),
    csmoke = if_else(smoking == 1, 1, 0),
    fsmoke = if_else(smoking == 2, 1, 0),
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
    rename(ppm25 = PM25conc_exposureugm3,
           district = ID_COUNTY) %>%
  # add dummies for district
  add_dummy_variables(district,
  values = c(1,2,3,4),
  remove_original = T) %>%
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


write_rds(d1, here("data-clean", "bhet-resp-data.rds"))

