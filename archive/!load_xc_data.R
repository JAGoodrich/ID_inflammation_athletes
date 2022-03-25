# Read in xc_data
proteins <- c("IFNg","IL1a","IL1b","IL2","IL4",
              "IL6","IL8","IL10","IL12p70","TNFa", "ery")

# Full xc_data
xc_data <- read_rds(here::here("FBV_data_for_analysis 1119.RDS")) %>% 
  filter(group == "Cross Country") %>% 
  droplevels() %>% 
  janitor::remove_empty(which = "cols")  %>% 
  select(-stfr_fer, 
         -stfr_fer_quant)


# Summary Data
fsd <- read_rds(here::here("FBV_summary_data.RDS")) %>% 
  filter(group == "Cross Country") %>% 
  droplevels() %>% 
  mutate(id_group = str_c(nationals, iron_def_avg, sep = ", ")%>% 
           as.factor %>% 
           fct_relevel("NCAA Qualifiers, Not ID", 
                       "Non NCAA Qualifiers, Not ID", 
                       "Non NCAA Qualifiers, ID")) %>%
  janitor::remove_empty(which = "cols") %>% 
  droplevels()  %>% 
  select(-stfr_fer, 
         -stfr_fer_quant)



# Calculate MCV, acute chronic id
xc_data <- xc_data %>% 
  mutate(mcv = (hct*10)/rbc, 
         acute_chronic_id = case_when(sex == "Male" & ferritin < 25 ~ "ID",
                                      TRUE ~ as.character(acute_chronic_id)) %>% 
           as.factor() %>% 
           fct_relevel("Not ID"))

# Create Dataframe for acute and chronic 
chronic_id <- fsd %>% filter(ferritin < 25)
acute_id <- xc_data %>% filter(ferritin < 25)


########  Group by acute ID
xc_data <- xc_data %>%
  mutate(acute_chronic_id_group = case_when(id %in% chronic_id$id ~ "Chronic ID", 
                                            id %in% acute_id$id ~ "Acute ID",
                                            TRUE ~ "No ID") %>% 
           as.factor() %>%  
           fct_relevel("No ID", "Acute ID","Chronic ID"))

# Clean environment
rm(chronic_id, acute_id)


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
         date_for_plotting = as.numeric(visit_ordered) %>% as.factor()) 

# Create quantcut vars 
xc_data <- xc_data %>% 
  mutate(across(all_of(proteins), 
                ~quantcut(., q = 4),
                .names = "{.col}_pro_quant"), 
         across(all_of(proteins), 
                ~quantcut(., q = 4) %>% as.numeric(),
                .names = "{.col}_pro_quantnum"), 
         across(all_of(proteins), 
                ~scale(.),
                # ~log(.+0.1),
                .names = "{.col}_pro_lg"), 
         across(all_of(proteins), 
                function(x){x+.1}))


