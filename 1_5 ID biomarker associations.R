# Associations of ID (defined by Peeling et al 2008) and NCAA qualifier status
# on hematological parameters 


## Run all models, stratified by sex ------------------------------------
# Pivot_longer: Reshapes data from wide to long format, with a single
# column for each outcome and a single column for the 
# concentration (named "outcome_value")
# Second pivot longer:  Reshapes data from wide to long format, with a single
# column for each proteine and a single column for the 
# concentration (named "cytokine contcentration")
# group_by: groups data by sex, outcome, and cytokine to fit a 
# separate model for each outcome/cytokine/sex combination 
# mutate/map: runs the actual lme model 


mod_output <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(c("hgb", "ferritin","hepcidin", "stfr")), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  tidylog::pivot_longer(cols = all_of(c("id_peeling", "nationals")), 
                        names_to = "ind_var", 
                        values_to = "ind_var_value") %>% 
  filter(!is.na(outcome_value), 
         !is.na(ind_var_value)) %>% 
  droplevels() %>% 
  group_by(sex, outcome, ind_var) %>% 
  nest() %>% 
  mutate(
    # Run full model
    model_full = map(data, 
                     ~lme(outcome_value ~ as.factor(date)+ ind_var_value, 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          data = ., 
                          na.action = na.exclude )))

# Model Estimates
mod_ests <- mod_output %>% 
  mutate(mod_summary = map(model_full, 
                           ~summary(.)$tTable %>% 
                             as_tibble(rownames = "effect")), 
         mod_ests = map(model_full, 
                        ~intervals(., which = "fixed")$fixed %>% 
                          as_tibble(rownames = "effect")), 
         coef_ci = map2(mod_summary, mod_ests, 
                        ~full_join(.x, .y))) %>% 
  select(-c(data, model_full, mod_summary, mod_ests)) %>%
  unnest(coef_ci) %>% 
  ungroup() %>%
  clean_names() %>% 
  filter(str_detect(effect, "ind_var"))
