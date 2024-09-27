db <- read_dta(here("data-clean", 
                   "BHET_master_data_15Jan2024.dta"), 
  col_select= c(hh_id, ptc_id, wave, ID_VILLAGE, ID_COUNTY, 
                ban_status_2019, ban_status_2020, 
                ban_status_2021, ban_status_no, 
                ban_status_composite,
                freq_cough, freq_phlegm,
                freq_wheezing, freq_breath,
                freq_no_chest, PM25_indoor_seasonal_hs,
                height, weight, temp, sys_brachial,
                dia_brachial)) 


db1 <- db %>% drop_na(-PM25_indoor_seasonal_hs) %>%
  mutate(
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
    rename(ipm25 = PM25_indoor_seasonal_hs,
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
    values=c(2018,2019,2021), remove_original = F)

dt <- readxl::read_excel(here("data-clean", 
  "HH_indoor_temperature_HEATING_SEASON.xlsx")) %>%
  select(ptc_id, wave, min_h, med_h, max_h,
         sd_h, iqr_h)

# merge indoor temp data with resp data
db2 <- db1 %>%
  left_join(dt, by = join_by(ptc_id, wave))

library(fixest)

mm <- fixest::feols(
  temp ~ treat | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_te <- fixest::feols(
  dia_brachial ~ treat | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_cde <- fixest::feols(
  dia_brachial ~ treat + temp | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_cde_i <- fixest::feols(
  dia_brachial ~ treat * temp | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_cde_ic <- fixest::feols(
  dia_brachial ~ treat * c_temp | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_te_h <- fixest::feols(
  dia_brachial ~ treat:i(cohort_year, i.year, ref=-Inf, ref2 = 2018)  | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)

om_cde_h <- fixest::feols(
  dia_brachial ~ treat:i(cohort_year, i.year, ref=-Inf, ref2 = 2018) + temp  | 
    cohort_year + year,
  data = db2,
  vcov = ~ID_VILLAGE)




regmedint_obj <- regmedint(data = db2,
                           ## Variables
                           yvar = "dia_brachial",
                           avar = "treat",
                           mvar = "c_temp",
                           cvar = c(c("year_2019", "year_2021", "cohort_year_2019", "cohort_year_2020", "cohort_year_2021")),
                           # eventvar = "event",
                           ## Values at which effects are evaluated
                           a0 = 0,
                           a1 = 1,
                           m_cde = 0,
                           c_cond = c(0.34,0.38,0.21,0.14,0.05),
                           ## Model types
                           mreg = "linear",
                           yreg = "linear",
                           ## Additional specification
                           interaction = TRUE,
                           casecontrol = FALSE)

    
    
    etwfe_c = fixest::feols(
  y ~ treat:i(cohort_year, i.year, ref=Inf, ref2 = 2011) | 
    cohort_year + year,
  data = d,
  vcov = ~state
)
)
