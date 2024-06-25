# Load basic packages -----------------------------------------------------
library( tidyverse)
library( lubridate)
library( adjustedCurves)
library( survival)
library( survminer)
library(patchwork)
library( ggpmisc)

# Universal functions -----------------------------------------------------
temp_func <- function( raw_model){
  
  temp <- summary( raw_model )
  
  p <- temp$coefficients %>% as.data.frame()%>% janitor::clean_names() %>% pull( pr_z)
  
  Cox_results <- 
    temp$conf.int %>% as_tibble( rownames = "Vars") %>% janitor::clean_names() %>%   
    mutate( coef = temp$coefficients[,"coef"],
            se = temp$coefficients[,4],
            n = temp$n,
            nevent = temp$nevent,
            HR_CI = paste(format(round(exp_coef, 2), nsmall = 2, scientific = FALSE), " (",
                          format(round(lower_95, 2), nsmall = 2), " to ",
                          format(round(upper_95, 2), nsmall = 2), ")",sep = ""),
            p_value = p) 
  return( Cox_results)
  
}
missingness_func <- function( df){
  
  df %>% summarise( across( everything(), ~ sum( is.na(.))/ n()))
  
  
}

# Import Covid-19 testing data --------------------------------------------
#England Covid-19 test results
covid19_result_england <- 
  data.table::fread( input = "D:/DPhil/UK_Biobank_opioid_application/COVID_test_results/Version_20230405/covid19_result_england.txt") %>% 
  mutate( specdate = lubridate::dmy( specdate)) %>% 
  select( eid, specdate, result, origin) 

covid19_result_scotland <- 
  data.table::fread( input = "D:/DPhil/UK_Biobank_opioid_application/COVID_test_results/Version_20230405/covid19_result_scotland.txt") %>% 
  mutate( specdate = lubridate::dmy( specdate)) %>%  
  select( eid, specdate, result)

covid19_result_wales <- 
  data.table::fread( input = "D:/DPhil/UK_Biobank_opioid_application/COVID_test_results/Version_20230405/covid19_result_wales.txt") %>% 
  mutate( specdate = lubridate::dmy( specdate))%>% 
  select( eid, specdate, result)


# Import baseline raw datasets ---------------------------------------------------------
# UKBB England participants
Population_characteristics___Baseline_characteristics___Indices_of_Multiple_Deprivation <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/Main_dataset/Population_characteristics___Baseline_characteristics___Indices_of_Multiple_Deprivation.rds")
England_participants_ID <- filter( Population_characteristics___Baseline_characteristics___Indices_of_Multiple_Deprivation, !is.na(f.26410.0.0)) %>% rename( eid =  f.eid)

UK_Biobank_Assessment_Centre___Touchscreen___Sociodemographics___Education <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/Main_dataset/UK_Biobank_Assessment_Centre___Touchscreen___Sociodemographics___Education.rds")

baseline_df <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/R_datasets/COVID-19_MAIN_df/baseline_df.rds")

# Import lifestyle --------------------------------------------------------
lancet_lifestyle_tidy_df <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Covid_PRS/lancet_lifestyle_tidy_df.rds")
# Import HES --------------------------------------------------------------

# (attention: different HES version)
hesin <- read_delim("D:/DPhil/UK_Biobank_opioid_application/HES/Version_20230405/hesin.txt",
                    delim = "\t",
                    escape_double = FALSE,
                    col_types = cols(epistart = col_date(format = "%d/%m/%Y"),
                                     epiend = col_date(format = "%d/%m/%Y"),
                                     elecdate = col_date(format = "%d/%m/%Y"),
                                     admidate = col_date(format = "%d/%m/%Y"),
                                     disdate = col_date(format = "%d/%m/%Y"),
                                     admisorc = col_character(),
                                     disdest = col_character()),
                    guess_max = 5000,
                    trim_ws = TRUE)

hesin_diag <- read_delim("D:/DPhil/UK_Biobank_opioid_application/HES/Version_20230405/hesin_diag.txt",
                         delim = "\t",
                         escape_double = FALSE,
                         col_types = cols(diag_icd9 = col_character(),
                                          diag_icd9_nb = col_character(),
                                          diag_icd10_nb = col_character()),
                         trim_ws = TRUE)

# Import death registry ---------------------------------------------------
death <- read_delim("D:\\DPhil\\UK_Biobank_opioid_application\\Death\\Version_20230405\\death.txt", 
                    delim = "\t",
                    col_types = cols(date_of_death = col_date(format = "%d/%m/%Y")))

death_cause <- read_delim("D:\\DPhil\\UK_Biobank_opioid_application\\Death\\Version_20230405\\death_cause.txt", delim = "\t")


# get rid of duplication
death <- 
  death %>% 
  arrange(date_of_death) %>% 
  group_by( eid) %>% 
  filter( row_number() == 1) %>% 
  ungroup() 

# Import GP vaccination (diagnosis source) -----------------------------------------
vaccine_diagnosis_source <- readRDS("D:/DPhil/UK_Biobank_opioid_application/UKB_based_projects/Sub_project/Breakthrough_infection/Derived_dataset/vaccine_diagnosis_source.rds")

full_vaccination_df <- 
  vaccine_diagnosis_source %>% 
  filter( code == "Dose_2") %>% 
  arrange( event_dt) %>% 
  group_by( eid) %>% 
  filter( row_number() == 1) %>% 
  ungroup() %>% 
  filter( event_dt >= as.Date("2020/12/01"), event_dt <= as.Date("2021/10/01")) %>% 
  rename( vaccination_dt = event_dt)

# Curate all population ---------------------------------------------------
All <- 
  baseline_df %>% 
  left_join( death, by = "eid") %>% 
  mutate( index_date = as.Date("2020/03/15")) %>% 
  filter( is.na(date_of_death) | date_of_death > index_date) %>% 
  left_join( select(lancet_lifestyle_tidy_df, 1:10) , by = c("eid" = "f.eid")) %>% # lifestyle factor
  left_join( baseline_education, by = "eid") %>%  # addtional baseline information
  mutate( Tidy_BMI = case_when( bmi >= 30 ~ 1, bmi < 30 ~ 0, TRUE ~ NA_integer_)) 


# Curate target population cohort (COVID-19)-----------------------------------------
target_population <- 
  covid19_result_combined %>% 
  filter( result == 1) %>% 
  arrange( specdate) %>% 
  group_by( eid) %>% 
  filter( row_number() == 1) %>%  
  ungroup() %>% 
  filter( specdate <= as.Date("2022/03/01")) %>% 
  left_join( death, by = "eid") %>% 
  left_join( full_vaccination_df, by = "eid") %>% 
  rename( index_date = specdate) %>% 
  mutate( index_date = case_when( origin == 1 ~ index_date - sample(1:7, n(), replace = TRUE), TRUE ~ index_date)) %>% 
  mutate( index_date_update = index_date + 30) %>% 
  filter( index_date <= date_of_death | is.na(date_of_death), eid != "6024915") %>% # important QC 
  select( eid, index_date, vaccination_dt, index_date_update, date_of_death, origin) 

# baseline charateristics
linked_baseline <- 
  target_population %>% 
  left_join( select(lancet_lifestyle_tidy_df, 1:10) , by = c("eid" = "f.eid")) %>% # lifestyle factor
  left_join( baseline_df, by = "eid") %>%  # baseline information
  left_join( baseline_education, by = "eid")   # addtional baseline information


# Modelling with COX ---------------------------------------------------------------
Model_df_overall <- 
  linked_baseline %>% 
  left_join( history_cohort, by = "eid") %>% 
  left_join( outcome_cohort, by = "eid") %>% 
  mutate( across(.cols = starts_with("history_"), ~ case_when( is.na(.x) ~ 0, TRUE ~ .))) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + 210, as.Date( "2022-10-01"), na.rm = TRUE), 
                                          TRUE ~ pmin( date_of_death, index_date + 210, as.Date( "2022-10-01"), na.rm = TRUE)),
          across( .cols = starts_with("occur_date"), ~ case_when( !is.na(.x) & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
          across( .cols = starts_with("occur_date"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}"),
          across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x ))) 

Model_df_acute <- 
  linked_baseline %>% 
  left_join( history_cohort, by = "eid") %>% 
  left_join( outcome_cohort, by = "eid") %>% 
  mutate( across(.cols = starts_with("history_"), ~ case_when( is.na(.x) ~ 0, TRUE ~ .))) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + 30, as.Date( "2022-10-01")), 
                                          TRUE ~ pmin( date_of_death, index_date + 30, as.Date( "2022-10-01"))),
          across( .cols = starts_with("occur_date"), ~ case_when( !is.na(.x) & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
          across( .cols = starts_with("occur_date"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}"),
          across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x ))) 

Model_df_postacute <- 
  linked_baseline_sub %>% 
  left_join( history_cohort_sub, by = "eid") %>% 
  left_join( outcome_cohort_sub, by = "eid") %>% 
  mutate( across(.cols = starts_with("history_"), ~ case_when( is.na(.x) ~ 0, TRUE ~ .))) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date_update + 180, as.Date( "2022-10-01")), 
                                          TRUE ~ pmin( date_of_death, index_date_update + 180, as.Date( "2022-10-01"))),
          across( .cols = starts_with("occur_date"), ~ case_when( !is.na(.x) & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
          across( .cols = starts_with("occur_date"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date_update), .names = "follow_up_days_{.col}"),
          across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x ))) 



# estimate HR (Unfavourable vs Intermediate vs favourable)
Model_lifestyle_category_func <- function( df = Model_df_overall, 
                                           outcomes = "Pulmonary", 
                                           exposure = "Tidy_total_score_CAT"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0)
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( exposure, "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=2)
  
  output_statistics <- 
    complete_df %>% 
    group_by( Tidy_total_score_CAT) %>% 
    summarise( IR = sum(.data[[paste("incident_outcome_occur_date_", outcomes, sep = "")]])/sum(.data[[paste("follow_up_days_occur_date_", outcomes, sep = "")]]),
               group_n = n(),
               group_event = sum(.data[[paste("incident_outcome_occur_date_", outcomes, sep = "")]])) %>% 
    ungroup() %>% 
    mutate( IR_year_1000 = IR*365*1000) %>% 
    mutate( vars  = paste( "Tidy_total_score_CAT", Tidy_total_score_CAT, sep = "")) %>% 
    select( -Tidy_total_score_CAT)
  
  output_all <- 
    output_statistics %>% 
    left_join( output_model, by = "vars") %>% 
    mutate( IRD_year_1000 = (exp_coef - 1) * .data$IR_year_1000[1],
            IRD_year_1000_lower = (lower_95 - 1) * .data$IR_year_1000[1],
            IRD_year_1000_upper = (upper_95 - 1) * .data$IR_year_1000[1])
  
  return(output_all)
  
}

# all cohort for sequelae 
lifestyle_category_overall_HR <- 
  map_df( c(unique(Final_code_list$organ_system), "any", "hospital", "death") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = Model_df_overall,
          .id = "target_outcome") 

lifestyle_category_acute_HR <- 
  map_df( c(unique(Final_code_list$organ_system), "any", "hospital", "death") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = Model_df_acute,
          .id = "target_outcome")

lifestyle_category_postacute_HR <- 
  map_df( c(unique(Final_code_list$organ_system), "any", "hospital", "death") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_category_func,
          df = Model_df_postacute,
          .id = "target_outcome") 



# estimate HR (Unfavourable vs Intermediate vs favourable)
Cali_lifestyle_category_func <- function( df = Model_df_overall, 
                                           outcomes = "any", 
                                           exposure = "Tidy_total_score_CAT"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0) 
  
  Control <- 
    complete_df %>% 
    filter( .data[[exposure]] == 1) %>% 
    summarise( outcome = sum( incident_outcome_occur_date_any), no_outcome = n() - outcome)
  
  Exposure_moderate <- 
    complete_df %>% 
    filter( .data[[exposure]] == 2) %>% 
    summarise( outcome = sum( incident_outcome_occur_date_any), no_outcome = n() - outcome)
  
  Exposure_healthy <- 
    complete_df %>% 
    filter( .data[[exposure]] == 3) %>% 
    summarise( outcome = sum( incident_outcome_occur_date_any), no_outcome = n() - outcome)
  
  Calibrate_moderate<- episensr::misclassification( matrix(c(unname(unlist(Exposure_moderate[1,])), unname(unlist(Control[1,]))),
                                  dimnames = list(c("Case", "No_case"),
                                                  c("Exposure", "Control")),
                                  nrow = 2, byrow = FALSE),
                                  type = "exposure",
                                  bias_parms = c(0.95, 0.95, 0.90, 0.90))
  
  Moderate_lifestyle <- 
  tibble( Analysis = c("Original", "Cali"), 
          OR = c(Calibrate_moderate$obs.measures[2,1], Calibrate_moderate$adj.measures[2,1]), 
          OR_lower = c(Calibrate_moderate$obs.measures[2,2], Calibrate_moderate$adj.measures[2,2]), 
          OR_upper = c(Calibrate_moderate$obs.measures[2,3], OR_upper = Calibrate_moderate$adj.measures[2,3])) 

  Calibrate_healthy <- episensr::misclassification(matrix(c(unname(unlist(Exposure_healthy[1,])), unname(unlist(Control[1,]))),
                                     dimnames = list(c("Case", "No_case"),
                                                     c("Exposure", "Control")),
                                     nrow = 2, byrow = FALSE),
                                     type = "exposure",
                                     bias_parms = c(0.95, 0.95, 0.90, 0.90))
  
  Healthy_lifestyle <- 
  tibble( Analysis = c("Original", "Cali"), 
          OR = c(Calibrate_healthy$obs.measures[2,1], Calibrate_healthy$adj.measures[2,1]), 
          OR_lower = c(Calibrate_healthy$obs.measures[2,2], Calibrate_healthy$adj.measures[2,2]), 
          OR_upper = c(Calibrate_healthy$obs.measures[2,3], OR_upper = Calibrate_healthy$adj.measures[2,3])) 

  output_df <- bind_rows( list( Moderate_lifestyle = Moderate_lifestyle, Healthy_lifestyle = Healthy_lifestyle), .id = "Group")
  
  return( output_df)
  
}

# all cohort for sequelae 
Cali_lifestyle_overall_OR <- 
  map_df( c("any", "hospital", "death") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Cali_lifestyle_category_func,
          df = Model_df_overall,
          .id = "target_outcome") 


# Modelling with mediator -------------------------------------------------
# generate mediator variables
function(){
  
  # Outcome records before index date
  M_history_sequelae_cohort <- 
    target_population %>% 
    left_join( select( baseline_df, eid, recruitment_date), by = "eid") %>% 
    left_join( all_outcome_record, by = "eid") %>% 
    filter( admidate < index_date, admidate >= recruitment_date) %>% 
    arrange( desc(admidate)) %>% 
    group_by( eid, organ_system) %>% 
    filter( row_number() == 1) %>% 
    mutate( index = 1) %>% 
    pivot_wider( id_cols = eid, names_from = organ_system, values_from = index, names_prefix = "history_") %>% 
    ungroup() %>% 
    rowwise() %>%
    mutate(history_any = min(c_across(-1), na.rm = TRUE)) %>%
    ungroup()
  
  M_history_hospital_cohort <- 
    target_population %>% 
    left_join( select( baseline_df, eid, recruitment_date), by = "eid") %>% 
    left_join( hesin, by = "eid") %>% 
    filter( admidate < index_date, admidate >= recruitment_date) %>% 
    arrange( desc(admidate)) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    mutate( history_hospital = 1) %>% 
    select( eid, history_hospital)
  
  M_history_cohort <- 
    M_history_sequelae_cohort %>% 
    full_join( M_history_hospital_cohort, by = "eid") %>% 
    mutate( history_death = NA) 
  
  # Outcome records after and on index date
  M_outcome_sequelae_cohort <- 
    target_population %>% 
    left_join( all_outcome_record, by = "eid") %>% 
    filter( admidate >= index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid, organ_system) %>% 
    filter( row_number() == 1) %>% 
    pivot_wider( id_cols = eid, names_from = organ_system, values_from = admidate, names_prefix = "occur_date_") %>% 
    ungroup() %>% 
    rowwise() %>%
    mutate(occur_date_any = min(c_across(-1), na.rm = TRUE)) %>%
    ungroup()
  
  M_outcome_hospital_cohort <- 
    target_population %>% 
    left_join( hesin, by = "eid") %>% 
    filter( admidate >= index_date) %>% 
    arrange( admidate) %>% 
    group_by( eid) %>% 
    filter( row_number() == 1) %>% 
    mutate( occur_date_hospital = admidate) %>% 
    select( eid, occur_date_hospital) %>% 
    ungroup()
  
  M_outcome_death_cohort <- 
    target_population %>% 
    mutate( occur_date_death = date_of_death) %>% 
    select( eid, occur_date_death)
  
  M_outcome_cohort <- 
    M_outcome_sequelae_cohort %>% 
    full_join(M_outcome_hospital_cohort, by = "eid") %>%
    full_join(M_outcome_death_cohort, by = "eid") 
  
  
}

M_Model_df_overall <- 
  linked_baseline %>% 
  left_join( M_history_cohort, by = "eid") %>% 
  left_join( M_outcome_cohort, by = "eid") %>% 
  mutate( across(.cols = starts_with("history_"), ~ case_when( is.na(.x) ~ 0, TRUE ~ .))) %>% 
  mutate( follow_up_end_date = case_when( is.na( date_of_death) ~ pmin( index_date + 210, as.Date( "2022-10-01"), na.rm = TRUE), 
                                          TRUE ~ pmin( date_of_death, index_date + 210, as.Date( "2022-10-01"), na.rm = TRUE)),
          across( .cols = starts_with("occur_date"), ~ case_when( !is.na(.x) & .x <= follow_up_end_date ~ 1, TRUE ~ 0), .names = "incident_outcome_{.col}"),
          across( .cols = starts_with("occur_date"), ~ as.numeric( pmin( follow_up_end_date, .x, na.rm = TRUE) - index_date), .names = "follow_up_days_{.col}"),
          across( .cols = starts_with("follow_up_days_"), ~ case_when( .x == 0 ~ .x + 1, TRUE ~ .x ))) 



# estimate HR (Unfavourable vs Intermediate vs favourable)
M_Model_lifestyle_category_func <- function( df = M_Model_df_overall, 
                                             outcomes = "any", 
                                             exposure = "Tidy_total_score_CAT"){
  
  complete_df <- df 

  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( exposure, paste("history_", outcomes, sep = ""),"age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=3)
  
  return(output_model)
  
}

# all cohort for sequelae 
Direct_effect_lifestyle_category_overall_HR <- 
  map_df( c(unique(Final_code_list$organ_system), "any") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          M_Model_lifestyle_category_func,
          df = M_Model_df_overall,
          .id = "target_outcome") 

# causal mediation estimation
Causal_medication_analysis_func <- function( df = M_Model_df_overall, 
                                             outcomes = "Diabetes", 
                                             exposure = "Tidy_total_score_BI"){
  
  complete_df <- df 
  
  # run model
  med_fit_formula_input <- 
    as.formula( paste( 
      paste( paste("history_", outcomes, sep = ""), "~", sep = ""), 
      paste( c( exposure, "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  out_fit_formula_input <- 
    as.formula( paste( 
      paste( "incident_outcome_occur_date_", outcomes, "~", sep = ""), 
      paste( c( exposure, paste("history_", outcomes, sep = ""), "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  med_fit <- glm( med_fit_formula_input, data = complete_df, family = binomial())
  out_fit <-glm( out_fit_formula_input, data = complete_df, family = binomial())
  
  # Quasi-Bayesian Monte Carlo 
  med_out_12 <-mediation::mediate( med_fit, 
                                   out_fit, 
                                   treat = "Tidy_total_score_BI", 
                                   mediator = paste("history_", outcomes, sep = ""), 
                                   control.value = 1, 
                                   treat.value = 2,
                                   robustSE = TRUE, 
                                   sims = 100) 
  
  # med_out_13 <-mediation::mediate( med_fit, 
  #                                  out_fit, 
  #                                  treat = "Tidy_total_score_CAT", 
  #                                  mediator = paste("history_", outcomes, sep = ""), 
  #                                  control.value = 1, 
  #                                  treat.value = 3,
  #                                  robustSE = TRUE, 
  #                                  sims = 100) 
  
  mediation_result_12 <- summary(med_out_12)
  # mediation_result_13 <- summary(med_out_13)
  
  output <- tibble( mediation_proportion = mediation_result_12$n.avg, P = mediation_result_12$n.avg.p, direct_proportion = 1 - mediation_proportion)
  
  return(output)
  
}
Medication_proporation <- 
  map_df( c(unique(Final_code_list$organ_system), "any") %>% purrr::set_names(),
          exposure = "Tidy_total_score_BI",
          Causal_medication_analysis_func,
          df = M_Model_df_overall,
          .id = "target_outcome") 

# Individual lifestyle modelling ------------------------------------------
Model_binary_individual_factor_subgroup_naive_func <- function( df = Model_df_overall, 
                                                                outcomes = "any", 
                                                                exposure = "Tidy_PA"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0) %>% 
    mutate(  !!exposure := fct_rev(as.factor(.data[[exposure]])))
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( exposure, "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=1)
  
  return(output_model)
  
}

tidy_string <- c("Tidy_smoking", "Tidy_alcohol_intake", "Tidy_PA", "Tidy_BMI", "Tidy_TV_view", "Tidy_sleep_time", "Tidy_oil_fish", "Tidy_process_meat", "Tidy_fruit_vegetable", "Tidy_red_meat")

binary_individual_any_HR <- 
  map_df( tidy_string %>% purrr::set_names(),
          Model_binary_individual_factor_subgroup_naive_func,
          df = Model_df_overall,
          outcomes = "any",
          .id = "target_exposure") 

binary_individual_hospital_HR <- 
  map_df( tidy_string %>% purrr::set_names(),
          Model_binary_individual_factor_subgroup_naive_func,
          df = Model_df_overall,
          outcomes = "hospital",
          .id = "target_exposure") 

binary_individual_death_HR <- 
  map_df( tidy_string %>% purrr::set_names(),
          Model_binary_individual_factor_subgroup_naive_func,
          df = Model_df_overall,
          outcomes = "death",
          .id = "target_exposure") 

individual_lifestyle_semi_adjust_HR <- 
  bind_rows( list(any = binary_individual_any_HR,
                  hospital = binary_individual_hospital_HR,
                  death = binary_individual_death_HR),
             .id = "target_outcome")



Model_binary_individual_factor_subgroup_func <- function( df = Model_df_overall, 
                                                          outcomes = "any"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0) %>% 
    mutate( across( .cols = all_of(adjust_factor), ~ case_when(. == 1 ~ 0, . == 0 ~ 1)))
  
  adjust_factor <- c( "Tidy_smoking", "Tidy_alcohol_intake", "Tidy_PA", "Tidy_BMI", "Tidy_TV_view", "Tidy_sleep_time", "Tidy_oil_fish", "Tidy_process_meat", "Tidy_fruit_vegetable", "Tidy_red_meat")
  
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( adjust_factor, "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=10)
  
  
  
  return( output_model)
  
}

individual_lifestyle_full_adjust_HR <- 
  map_df( c("any", "hospital", "death") %>% purrr::set_names(),
          Model_binary_individual_factor_subgroup_func,
          df = Model_df_overall,
          .id = "target_outcome") 


Model_multiple_func <- function( df = Model_df_overall, 
                                 outcomes = "any"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0) %>% 
    mutate( across( .cols = all_of(adjust_factor), ~ case_when(. == 1 ~ 0, . == 0 ~ 1)))
  
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( "Tidy_total_score_multi", "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) %>% temp_func() %>% filter( row_number() <=6)
  
  output_statistics <- 
    df %>% 
    group_by( Tidy_total_score_multi) %>% 
    summarise( group_n = n(),
               group_event = sum(.data[[paste("incident_outcome_occur_date_", outcomes, sep = "")]])) %>% 
    ungroup() %>% 
    mutate( vars  = paste( "Tidy_total_score_multi", Tidy_total_score_multi, sep = "")) %>% 
    select( -Tidy_total_score_multi)
  
  output_all <- 
    output_statistics %>% 
    left_join( output_model, by = "vars") 
  
  return( output_all)
  
}

individual_number_of_lifestyle <- 
  map_df( c("any", "hospital", "death") %>% purrr::set_names(),
          Model_multiple_func,
          df = Model_df_overall,
          .id = "target_outcome") 

# Assumption_test ---------------------------------------------------------
Model_lifestyle_PH_assumption_func <- function( df = Model_df_overall, 
                                                outcomes = "Pulmonary", 
                                                exposure = "Tidy_total_score_CAT"){
  
  complete_df <- 
    df %>% 
    filter( .data[[paste("history_", outcomes, sep = "")]] == 0)
  
  # run model
  formula_input <- 
    as.formula( paste( 
      paste( "Surv( follow_up_days_occur_date_", outcomes, ",", "incident_outcome_occur_date_", outcomes, ")", "~", sep = ""), 
      paste( c( exposure, "age_year", "sex", "IMD", "Education"), collapse="+"),
      sep = ""))
  
  output_model <-  survival::coxph( formula = formula_input, data = complete_df) 
  test.ph <- survival::cox.zph(output_model)
  
  output_all <- test.ph$table %>% as.data.frame() %>% rownames_to_column() %>% filter(rowname == "Tidy_total_score_CAT")

  
  return(output_all)
  
}

lifestyle_PH_assumption <- 
  map_df( c(unique(Final_code_list$organ_system), "any", "hospital", "death") %>% purrr::set_names(),
          exposure = "Tidy_total_score_CAT",
          Model_lifestyle_PH_assumption_func,
          df = Model_df_postacute ,
          .id = "target_outcome") 

