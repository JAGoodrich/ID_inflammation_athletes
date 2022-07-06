# Hepcidin versus ferritin and stfr 

# Statistics ----------------------
# Run all models, stratified by sex 
# Pivot_longer: Reshapes data from wide to long format, with a single
# column for each outcome and a single column for the 
# concentration (named "outcome_value")
# group_by: groups data by sex and outcome (ie, each protein) to fit a 
# separate model for each outcome/sex combination 
# mutate/map: runs the actual lme model 
mod_output_hep_id <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(c("ferritin", "stfr")), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  filter(!is.na(outcome_value)) %>% 
  droplevels() %>% 
  group_by(sex, outcome) %>% 
  nest() %>% 
  mutate(
    # Run full model
    model_full = map(data, 
                     ~lme(log(outcome_value) ~ as.factor(date)+
                            log2(hepcidin), 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          data = ., 
                          na.action = na.exclude)), 
    mod_summary = map(model_full, 
                      ~summary(.)$tTable %>% 
                        as_tibble(rownames = "effect")), 
    mod_ests = map(model_full, 
                   ~intervals(., which = "fixed")$fixed ), 
    mod_ests = map(mod_ests, 
                   ~((exp(.)-1)*100) %>% 
                     round(digits = 1) %>% 
      as_tibble(rownames = "effect")),
    coef_ci = map2(mod_summary, mod_ests, 
                   ~full_join(.x, .y)), 
    mod_summary = map(model_full, 
                      ~sjPlot::tab_model(.x)), 
    mod_summary = map(mod_summary, 
                      ~data.frame(readHTMLTable(htmlParse(.x))[1]) %>% 
                        row_to_names(row_number = 1)))


mod_output_hep_id_effect_est <- mod_output_hep_id %>% 
  select(-c(data, model_full, mod_summary, mod_ests)) %>%
  unnest(coef_ci) %>% 
  ungroup() %>%
  clean_names() %>% 
  filter(str_detect(effect, "log2"))


# Get all model results
hepcidin_outcomes_res <- mod_output_hep_id %>%
  select(-data, -model_full, -mod_ests, -coef_ci) %>%
  unnest(mod_summary) %>%
  mutate(description = "Associations of hepcidin with iron homeostasis biomarkers") %>%
  select(description, everything())

#View info about the one outlier participant
# View(xc_data %>% select(id, visit,sex, hepcidin, ferritin, stfr, id_peeling))

# Ferritin vs. Hepcidin -------------------------
(fer_hep <- xc_data %>% 
   mutate(normrange = if_else(ferritin < 30, 
                              "Out of range", 
                              "Within range")) %>% 
   filter(visit_ordered != "Early Fall, Year 1") %>%
   ggplot(aes(x = hepcidin, 
              y = ferritin, 
              group = id, 
              shape = normrange, 
              alpha = normrange)) + 
   geom_hline(yintercept = 30,  linetype = 2, color = "grey50") +
   geom_point(alpha = .6) + 
   geom_path(alpha = .1) +
   scale_x_log10() +
   scale_y_log10() +
   scale_shape_manual(values = c(1,19)) +
   scale_alpha_manual(values = c(.75, .5)) +
   xlab("Hepcidin (ng/mL)") +
   ylab("Ferritin (ng/mL)") +
   facet_wrap(~sex) +
   # facet_grid(sex ~ visit_ordered, drop = TRUE) + 
   theme(legend.position = "none"))
fer_hep


# Hepcidin vs. stfr -------------------------
(stfr_hep <- xc_data %>% 
   mutate(`Normal Range` = if_else(stfr > (29.5), 
                                   "Out of range", 
                                   "Within range"
   )) %>% 
   filter(!is.na(hepcidin), !is.na(stfr)) %>% 
   ggplot(aes(y = stfr, 
              x = hepcidin, 
              group = id, 
              shape = `Normal Range`, 
              alpha = `Normal Range`)) + 
   geom_hline(yintercept = 29.5,
              linetype = 2,
              color = "grey50") +
   geom_point() + 
   geom_path(alpha = .25) +
   scale_x_log10() +
   scale_y_log10() +
   scale_shape_manual(values = c(1,19)) +
   scale_alpha_manual(values = c(.75, .5)) +
   xlab("Hepcidin (ng/mL)") +
   ylab("sTfR (nmol/L)") +
   facet_wrap(~sex) + 
   theme(legend.position = "bottom"))
stfr_hep


# Hepcidin IL6 
# Hepcidin vs. stfr -------------------------
(il6_hep <- xc_data  %>% 
   filter(!is.na(hepcidin), !is.na(IL6), !is.na(id_peeling)) %>% 
   rename(`Iron Deficiency:` = id_peeling) %>% 
   ggplot(aes(x = IL6, 
              y = hepcidin, 
              group = id, 
              shape = `Iron Deficiency:`, 
              alpha = `Iron Deficiency:`)) + 
   geom_point() + 
   geom_path(alpha = .25) +
   scale_x_log10() +
   scale_y_log10() +
   scale_shape_manual(values = c(19, 1, 2)) +
   scale_alpha_manual(values = c(.5, .75, .75)) +
   ylab("Hepcidin (ng/mL)") +
   xlab("IL-6 (pg/mL)") +
   facet_wrap(~sex) + 
   theme(legend.position = "bottom"))
il6_hep



# Combine figures -------------------------------------------------
(figure_3 <- plot_grid(NULL, fer_hep, 
                       NULL, stfr_hep, 
                       ncol = 1, 
                       labels = c("A. Hepcidin and Ferritin", "", 
                                  "B. Hepcidin and sTfR", ""),
                       hjust = 0,
                       align = "v",
                       rel_heights = c(0.1, 1.0,
                                       0.1, 1.3)))

# Save Figure 
ggsave(figure_3, 
       filename = fs::path(dir_xcfig, 
                           "Figure 3 hepcidin ID scatterplots.jpg"), 
       width = 7, height = 6)



# Save IL-6 Figure 
ggsave(il6_hep, 
       filename = fs::path(dir_xcfig, 
                           "Supplemental Figure 1 hepcidin IL6 scatterplot.jpg"), 
       width = 7, height = 4)
