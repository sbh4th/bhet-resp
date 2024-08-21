# packages to load
pkgs <- c('here', 'tidyverse', 'modelsummary', 
          'fixest', 'marginaleffects',
          'kableExtra', 'patchwork', 'estimatr')

#load all packages at once
lapply(pkgs, library, character.only=TRUE)

## 1 Read in dataset, limit to resp vars ----
dresp <- read_rds(here("data-clean", 
  "bhet-resp-data.rds")) 

## 2 Estimate DD across several binary outcomes

# Function to estimate GLM for multiple outcomes
estimate_logit_did <- function(outcome_vars, predictor_vars) {
  models <- outcome_vars %>%
    set_names() %>%
    map(~ {
      formula <- as.formula(paste(.x, "~", 
        paste(predictor_vars, collapse = " + ")))
      glm(formula, data = dresp_cc, family = "binomial")
    })
  
  return(models)
}

# binary outcomes
b_out <- c("resp", "cough", "phlegm", "wheeze",
              "breath", "nochest")

# basic DiD specification
rhs_did <- c("treat:cohort_year_2019:year_2019",
  "treat:cohort_year_2019:year_2021", 
  "treat:cohort_year_2020:year_2021",
  "treat:cohort_year_2021:year_2021", "cohort_year_2019",
  "cohort_year_2020", "cohort_year_2021",
  "year_2019", "year_2021")

rhs_dida <- c(rhs_did, "age_health", "male", 
  "csmoke", "fsmoke")

dresp_cc <- dresp %>%
  # limit to complete cases
  drop_na(b_out, age_health, male, csmoke, fsmoke)

# gather estimates across models
# basic DiD
logit_did <- estimate_logit_did(b_out, rhs_did)

# with covariates
logit_dida <- estimate_logit_did(b_out, rhs_dida)

# estimate marginal effects (simple average)
# basic DiD
logit_me <- lapply(logit_did, marginaleffects::slopes, 
  newdata = subset(dresp, treat==1), 
  variables = "treat", by = "treat",
  # make sure to use cluster-robust SEs
  vcov = ~v_id)

# heterogenous treatment effects
logit_me_het <- lapply(logit_did, 
  marginaleffects::slopes, 
  newdata = subset(dresp, treat==1), 
  variables = "treat", 
  by = c("cohort_year", "year"),
  # make sure to use cluster-robust SEs
  vcov = ~v_id)

# adjusted DiD
logit_mea <- lapply(logit_dida, marginaleffects::slopes, 
  newdata = subset(dresp, treat==1), 
  variables = "treat", by = "treat",
  # make sure to use cluster-robust SEs
  vcov = ~v_id)

# heterogenous treatment effects
logit_mea_het <- lapply(logit_dida, 
  marginaleffects::slopes, 
  newdata = subset(dresp, treat==1), 
  variables = "treat", 
  by = c("cohort_year", "year"),
  # make sure to use cluster-robust SEs
  vcov = ~v_id)

# contrasts of ATTs
# heterogenous treatment effects
logit_mea_het_c <- lapply(logit_dida, 
  marginaleffects::slopes, 
  newdata = subset(dresp, treat==1), 
  variables = "treat", 
  by = c("cohort_year", "year"),
  hypothesis = c("b1 - b2 =0", 
    "b1 - b3 = 0", "b1 - b4 = 0"),
  # make sure to use cluster-robust SEs
  vcov = ~v_id)

# tests of heterogeneity
# Define a function to apply hypotheses function
apply_hypotheses <- function(variable_name) {
  hypotheses(logit_mea_het_c[[variable_name]], 
             joint = TRUE)
}

het_tests <- setNames(lapply(b_out, 
  apply_hypotheses), b_out)

# write heterogeneity test results to file
write_rds(het_tests, file = here("outputs/models",
  "r_het_tests.rds"))

# grab estimates and SEs from DiD results
did_t1 <- logit_me %>% {
  tibble(
    # table = 1,
    outcome = names(logit_me),
    est = map_dbl(., "estimate") * 100,
    stderror = map_dbl(., "std.error") * 100,
    ll = est - 1.96 * stderror,
    ul = est + 1.96 * stderror,
    ci = paste("(", sprintf("%.2f", ll), ", ",
      sprintf("%.2f", ul), ")", sep="")
  )
}

# grab number of observations
did_obs <- logit_did %>% {
  tibble(
    outcome = names(logit_did),
    nobs = map_int(., "df.null") + 1
  )
}

did_t1 <- did_t1 %>%
  inner_join(did_obs) %>%
  relocate(outcome, nobs)

# grab estimates and SEs from adjusted DiD results
did_t2 <- logit_mea %>% {
  tibble(
    # table = 2,
    # outcome = names(logit_me),
    esta = map_dbl(., "estimate") * 100,
    stderrora = map_dbl(., "std.error") * 100,
    lla = esta - 1.96 * stderrora,
    ula = esta + 1.96 * stderrora,
    cia = paste("(", sprintf("%.2f", lla), ", ",
      sprintf("%.2f", ula), ")", sep="")
  )
}

# grab number of observations
did_obs2 <- logit_dida %>% {
  tibble(
    outcome = names(logit_dida),
    nobs = map_int(., "df.null") + 1
  )
}

did_t2 <- did_t2 %>%
  inner_join(did_obs2) %>%
  relocate(outcome, nobs) %>%
  select

didt <- cbind(did_t1, did_t2) %>%
  dplyr::select(est, ci, esta, cia) %>%
  mutate(outcome = c("Any symptom", "Coughing",
    "Phlegm", "Wheezing attacks", "Trouble breathing",
    "Chest trouble"),
    category = "Self-reported (pp)") %>%
  # rename for combining with other tables
  rename(`estimate_1` = "est", `ci_1` = "ci",
         `estimate_2` = "esta", `ci_2` = "cia") %>%
  relocate(category, outcome)
  

# write results table to dataset
write_rds(didt, file = here("data-clean", 
  "resp-did.rds"))

# table of marginal effects
modelsummary(list("Any symptom" = logit_me$resp,
  "Cough" = logit_me$cough, "Phlegm" = logit_me$phlegm,
  "Wheezing" = logit_me$wheeze,
  "Shortness of breath" = logit_me$breath,
  "Chest trouble" = logit_me$nochest,
  "Any symptom" = logit_mea$resp,
  "Cough" = logit_mea$cough, "Phlegm" = logit_mea$phlegm,
  "Wheezing" = logit_mea$wheeze,
  "Shortness of breath" = logit_mea$breath,
  "Chest trouble" = logit_mea$nochest),
  shape = model ~ term + statistic,
  statistic = 'conf.int')

didtable <- modelsummary(list("Any" = logit_me$resp, "Cough" = logit_me$cough), shape = model ~ term + statistic,
  statistic = 'conf.int')

modelsummary(list("Any" = logit_me$resp, "Cough" = logit_me$cough), shape =  ~ "" + statistic,
  statistic = 'conf.int')

didatable <- modelsummary(list("Any" = logit_mea$resp, "Cough" = logit_mea$cough), shape = model ~ term + statistic,
  statistic = 'conf.int')

data_tables <- data.frame(good_table = didtable, 
                          bad_table = didatable)

# any respiratory outcome
resp_het_any <- bind_rows(logit_mea$resp, 
                          logit_mea_het$resp) %>% 
  mutate(outcome = "resp") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_any, file = here("outputs", 
  "resp-het-resp.rds"))

# heterogeneity table for any respiratory outcome

# cough
resp_het_cough <- bind_rows(logit_mea$cough, 
                          logit_mea_het$cough) %>% 
  mutate(outcome = "cough") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_cough, file = here("outputs", 
  "resp-het-cough.rds"))

# phlegm
resp_het_phlegm <- bind_rows(logit_mea$phlegm, 
                          logit_mea_het$phlegm) %>% 
  mutate(outcome = "phlegm") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_phlegm, file = here("outputs", 
  "resp-het-phlegm.rds"))


# wheeze
resp_het_wheeze <- bind_rows(logit_mea$wheeze, 
                          logit_mea_het$wheeze) %>% 
  mutate(outcome = "wheeze") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_wheeze, file = here("outputs", 
  "resp-het-wheeze.rds"))

# shortness of breah
resp_het_breath <- bind_rows(logit_mea$breath, 
                          logit_mea_het$breath) %>% 
  mutate(outcome = "breath") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_breath, file = here("outputs", 
  "resp-het-breath.rds"))

# no chest trouble
resp_het_nochest <- bind_rows(logit_mea$nochest, 
                          logit_mea_het$nochest) %>% 
  mutate(outcome = "chest") %>%
  select(outcome, estimate, conf.low, conf.high, 
         cohort_year, year) %>%
  mutate_at(vars(c(cohort_year,year)), ~ recode(., 
         `2019` = "2019",
         `2020` = "2020",
         `2021` = "2021",
         .missing = "All")) %>%
  relocate(outcome, cohort_year, year) %>%
  mutate(ci = paste("(", sprintf('%.2f', conf.low), ", ",
    sprintf('%.2f', conf.high), ")", sep="")) %>%
  select(-conf.low, -conf.high) 

# write table to data
write_rds(resp_het_nochest, file = here("outputs", 
  "resp-het-nochest.rds"))

# pre-trends plots
# set theme for pre-trends
theme_pt <- function() {
  theme_classic() + 
    theme(axis.title = element_text(size=14),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.subtitle = element_text(size = 12))
}


# Define the function to estimate the model, gather statistics, and create a plot
run_analysis <- function(outcome, data, 
  year_var = "year", treatment_var = "et", 
  cohort_var = "cohort_year_2021", 
  vcov_cluster = "v_id") {
  
  # Estimate the model for the given outcome
  formula <- as.formula(paste(outcome, "~", 
    year_var, "*", treatment_var))
  pt_model <- glm(formula, data = data)
  
  # Compute average comparisons for the time trends by treatment status
  pt_comparisons <- avg_comparisons(pt_model, 
    var = year_var, by = cohort_var, 
    vcov = as.formula(paste("~", vcov_cluster)))
  
  # Compute the difference in time trends by treatment
  pt_trend_diff <- avg_comparisons(pt_model, 
    var = year_var, by = cohort_var, 
    hypothesis = "b2 - b1 = 0", 
    vcov = as.formula(paste("~", vcov_cluster)))
  
  # Gather statistics for the difference in pre-trends
  pt_text <- paste("Difference in trend (SE) for treated",
  "vs. untreated villages:", sep = "\n") 
  pt_stats <- paste(sprintf("%.1f", 
    pt_trend_diff$estimate),
    " (", sprintf("%.1f", pt_trend_diff$std.error), 
    ")", ", 95% CI: ", sprintf("%.1f", pt_trend_diff$conf.low),
    ", ", sprintf("%.1f", pt_trend_diff$conf.high), sep = "")
  pt_test <- paste(pt_text, pt_stats, sep = "\n")
  
  # Plot the trends and estimates
  pt_plot <- plot_predictions(pt_model, 
    condition = c(year_var, treatment_var)) + 
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(breaks = c(2018, 2019)) +
    labs(subtitle = pt_test, 
      y = paste("Probability of", outcome),, x = "") +
    scale_color_manual(name = "Treated in 2021?",
                       labels = c("No", "Yes"),
                       values = c("#e41a1c", "#377eb8")) +
    scale_fill_manual(name = "Treated in 2021?",
                      labels = c("No", "Yes"),
                      values = c("#e41a1c", "#377eb8")) + 
    theme_pt() # Using a basic theme, adjust as needed

  return(list(model = pt_model, 
    comparisons = pt_comparisons, plot = pt_plot, 
    test_text = pt_test))
}

# List of outcomes
b_out <- c("resp", "cough", "phlegm", "wheeze", 
           "breath", "nochest")

# Corresponding descriptive y-axis labels
y_labels <- c("Probability of any symptom", 
              "Probability of coughing", 
              "Probability of phlegm", 
              "Probability of wheezing", 
              "Probability of breathlessness", 
              "Probability of no chest symptoms")

# List of outcomes
b_out <- c("resp", "cough", "phlegm", 
           "wheeze", "breath", "nochest")

# Apply the function across all outcomes and store results
results <- lapply(b_out, 
  function(outcome) run_analysis(outcome, data = drpt))

# Access plots, models, and comparison results
plots <- lapply(results, function(res) res$plot)
models <- lapply(results, function(res) res$model)
comparisons <- lapply(results, function(res) res$comparisons)
tests <- lapply(results, function(res) res$test_text)

combined_plot <- wrap_plots(plots) + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("images", "resp-pretrends.png"), 
       plot=combined_plot, width=8.5, height=11)

for (outcome in seq_along(b_out)) {
  nobs(logit_did$outcome)
}

