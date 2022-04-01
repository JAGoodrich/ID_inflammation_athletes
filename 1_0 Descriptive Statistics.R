# Tables and descriptive


# Get number of visits 
xc_data %>% 
  group_by(id) %>% 
  summarise(n_visits = length(id)) %>% 
  ungroup() %>% 
  summarise(n_visits_median = median(n_visits), 
            n_visits_iqr = IQR(n_visits),
            min = min(n_visits), 
            max = max(n_visits))
  



# Get maximum severuty of ID
id_once <- xc_data %>% 
  mutate(id_peeling = as.ordered(id_peeling)) %>%
  group_by(id, sex, nationals) %>% 
  summarise(id = max(id_peeling, na.rm = TRUE)) %>% 
  ungroup()


# Get percent ID at each stage
id_once %>% 
  group_by(sex, nationals) %>%
  summarise(n = length(id), 
            iron_replete = jag2::npct(id, "Iron Replete", n.digits = 2),  
            stage_1 = jag2::npct(id, "Stage 1", n.digits = 2), 
            stage_2 = jag2::npct(id, "Stage 2", n.digits = 2))  

# Test for differences between males and females
chisq.test(table(id_once$sex, id_once$id))



# Table 2
summary_data <- xc_data %>% 
  mutate(id_peeling = as.ordered(id_peeling)) %>%
  group_by(id, sex, nationals) %>% 
  summarise(n_visits = length(id), 
            across(where(is.numeric), ~mean(., na.rm = TRUE)), 
            id_peeling = max(id_peeling)) %>% 
  ungroup()


#Overall statistics
table_2_overall <- summary_data %>%
  group_by(sex) %>%
  summarise(`[Hb] (g/dL)` = jag2::avg_sd_fxn(hgb, n.digits = 1), 
            `Ferritin (ng/ml)` = jag2::median_iqr_fxn(ferritin, n.digits = 1), 
            `sTfR (nmol/L)` = jag2::median_iqr_fxn(stfr, n.digits = 1), 
            `Hepcidin (nmol/L)` = jag2::median_iqr_fxn(hepcidin, n.digits = 1), 
            `IL-1a (pg/ml)` = jag2::median_iqr_fxn(IL1a, n.digits = 1),  
            `IL-1b (pg/ml)` = jag2::median_iqr_fxn(IL1b, n.digits = 1),  
            `IL-2 (pg/ml)` =  jag2::median_iqr_fxn(IL2, n.digits = 1),  
            `IL-4 (pg/ml)` =  jag2::median_iqr_fxn(IL4, n.digits = 1),  
            `IL-6 (pg/ml)`  = jag2::median_iqr_fxn(IL6, n.digits = 1), 
            `IL-8 (pg/ml)` =  jag2::median_iqr_fxn(IL8, n.digits = 1),  
            `IL-10 (pg/ml)` = jag2::median_iqr_fxn(IL10, n.digits = 1),  
            `IL-12p70 (pg/ml)` = jag2::median_iqr_fxn(IL12p70, n.digits = 1),  
            `IFN-g (pg/ml)` = jag2::median_iqr_fxn(IFNg, n.digits = 1), 
            `TNF-a (pg/ml)` = jag2::median_iqr_fxn(TNFa, n.digits = 1), 
            `ERFE (ng/ml)`  = jag2::median_iqr_fxn(ery, n.digits = 1)) %>% 
  ungroup() %>% 
  column_to_rownames("sex") %>% 
  t(.) %>% 
  as_tibble(rownames = "Outcome")


#by nationals statistics
table_2_nats <- summary_data %>%
  group_by(sex, nationals) %>%
  summarise(`[Hb] (g/dL)` = jag2::avg_sd_fxn(hgb, n.digits = 1), 
            `Ferritin (ng/ml)` = jag2::median_iqr_fxn(ferritin, n.digits = 1), 
            `sTfR (nmol/L)` = jag2::median_iqr_fxn(stfr, n.digits = 1), 
            `Hepcidin (nmol/L)` = jag2::median_iqr_fxn(hepcidin, n.digits = 1), 
            `IL-1a (pg/ml)` = jag2::median_iqr_fxn(IL1a, n.digits = 1),  
            `IL-1b (pg/ml)` = jag2::median_iqr_fxn(IL1b, n.digits = 1),  
            `IL-2 (pg/ml)` =  jag2::median_iqr_fxn(IL2, n.digits = 1),  
            `IL-4 (pg/ml)` =  jag2::median_iqr_fxn(IL4, n.digits = 1),  
            `IL-6 (pg/ml)`  = jag2::median_iqr_fxn(IL6, n.digits = 1), 
            `IL-8 (pg/ml)` =  jag2::median_iqr_fxn(IL8, n.digits = 1),  
            `IL-10 (pg/ml)` = jag2::median_iqr_fxn(IL10, n.digits = 1),  
            `IL-12p70 (pg/ml)` = jag2::median_iqr_fxn(IL12p70, n.digits = 1),  
            `IFN-g (pg/ml)` = jag2::median_iqr_fxn(IFNg, n.digits = 1), 
            `TNF-a (pg/ml)` = jag2::median_iqr_fxn(TNFa, n.digits = 1), 
            `ERFE (ng/ml)`  = jag2::median_iqr_fxn(ery, n.digits = 1)) %>% 
  ungroup() %>%
  mutate(sex_nat = str_c(sex, " ", nationals)) %>% 
  select(-sex, -nationals) %>% 
  column_to_rownames("sex_nat") %>% 
  t(.) %>% 
  as_tibble(rownames = "Outcome")


# Join all 
table_2 <- tidylog::full_join(table_2_nats, table_2_overall) %>% 
  mutate(`_` = "") %>% 
  select(Outcome, 
         `Male NCAA Qualifiers`, 
         `Male Non NCAA Qualifiers`, 
         `Male`,
         `_`,
         `Female NCAA Qualifiers`,
         `Female Non NCAA Qualifiers`, 
         `Female`)

# Save 
write_csv(table_2, 
          file = fs::path(dir_xctab, 
                          "Table 2.csv"))
