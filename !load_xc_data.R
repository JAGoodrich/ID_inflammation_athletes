# Read in xc_data

# Read Full Dataset
xc_data <- read_rds(fs::path(dir_home,"FBV_data_for_analysis 1119.RDS")) %>% 
  filter(group == "Cross Country") %>% 
  droplevels() %>% 
  janitor::remove_empty(which = "cols")  %>% 
  select(-stfr_fer, 
         -stfr_fer_quant)


# Arrange by date
xc_data <- xc_data %>% 
  group_by(id) %>% 
  arrange(date) %>%
  ungroup()  %>% 
  mutate(year = if_else(ay == "AY1718", "Year 1", "Year 2"), 
         visit_ordered = str_c(time_point, year, sep = ", ") %>% 
           as.ordered() %>% 
           fct_reorder(., .x = as.numeric(date)))


# Calculate new stfr_logfer index  --------------------
# Ref:  https://doi.org/10.1002/ajh.22108
# Notes: In this dataset: stfr is in units of nmol/L
# "to manually convert nmol/L to mg/L, multiply by 0.0738" (from ref)

# Peeling 2019: stfr threshold of IDA: 2.5 mg/L
# Peeling 2019: stfr threshold of IDA: 29.5 nmol/L (actual ref: Koulaouzidis 2009)
xc_data <- xc_data %>% 
  mutate(stfr_logfer = (0.0738*stfr)/log10(ferritin), 
         id_peeling = case_when(ferritin < 20 | stfr > 29.5 ~ "Stage 2", 
                                ferritin < 35 ~ "Stage 1",
                                ferritin >= 35 & stfr <= 29.5 ~ "Iron Replete") %>% 
           as.factor(), 
         date_for_plotting = as.numeric(visit_ordered) %>% as.factor()) %>% 
  tidylog::select(-contains("iron_def"))

