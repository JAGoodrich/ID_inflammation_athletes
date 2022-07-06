# XC Cytokine Outcome Analysis

outcomes <- c("ferritin", "stfr", "hepcidin")


proteins <- c("IFNg","IL1a","IL1b","IL2","IL4",
              "IL6","IL8","IL10","IL12p70","TNFa", "ery")

# Statistics -----------------------------------------

## Run all models, stratified by sex ------------------------------------
# Pivot_longer: Reshapes data from wide to long format, with a single
# column for each outcome and a single column for the 
# concentration (named "outcome_value")
# Second pivot longer:  Reshapes data from wide to long format, with a single
# column for each protein and a single column for the 
# concentration (named "cytokine concentration")
# group_by: groups data by sex, outcome, and cytokine to fit a 
# separate model for each outcome/cytokine/sex combination 
# mutate/map: runs the actual lme model 

mod_output <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(outcomes), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  tidylog::pivot_longer(cols = all_of(str_c(proteins)), 
                        names_to = "explanatory_variable", 
                        values_to = "cytokine_concentration") %>% 
  filter(!is.na(outcome_value), 
         !is.na(cytokine_concentration)) %>% 
  droplevels() %>% 
  group_by(sex, outcome, explanatory_variable) %>% 
  nest() %>% 
  mutate(
    # Run full model
    model_full = map(data, 
                     ~lme(log(outcome_value) ~ as.factor(date)+
                            log2(cytokine_concentration), 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          data = ., 
                          na.action = na.exclude )), 
    mod_summary = map(model_full, 
                      ~sjPlot::tab_model(.x)), 
    mod_summary = map(mod_summary, 
                      ~data.frame(readHTMLTable(htmlParse(.x))[1]) %>% 
                        row_to_names(row_number = 1)))


# Get all model results
cytokine_outcome_res <- mod_output %>% 
  select(-data, -model_full) %>% 
  unnest(mod_summary) %>% 
  mutate(description = "Associations of cytokine concentrations with iron homeostasis biomarkers") %>% 
  select(description, everything())


## Obtain model estimates ------------------------------
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
  filter(str_detect(effect, "log2"))




# Format outcomes and cytokines 
mod_ests_final <- mod_ests %>% 
  mutate(sig = if_else(p_value < 0.05, "*", ""), 
         est_pct = (exp(est) - 1) * 100, 
         lower_pct = (exp(lower) - 1) * 100, 
         upper_pct = (exp(upper) - 1) * 100, 
         cytokine_fmt = cytokine %>% 
           str_replace("ery", "Erythroferrone") %>% 
           str_replace("hepcidin", "Hepcidin") %>% 
           fct_relevel(#"Hepcidin", 
             "IL1a", 
             "IL1b", 
             "IL2", 
             "IL4", 
             "IL6", 
             "IL8", 
             "IL10", 
             "IL12p70", 
             "IFNg", 
             "TNFa", 
             "Erythroferrone"))

View(temp <- mod_ests_final %>% 
  filter( sex == "Female"))



# Plot Effects ---------------------------------------
# hepcidin
(hepcidin_coefplot <- ggplot(mod_ests_final %>% 
                               filter(outcome == "hepcidin"), 
                             aes(x = cytokine_fmt,
                                 y=est_pct, 
                                 ymin = lower_pct,
                                 ymax = upper_pct, 
                                 # color= sig
                                 )) + 
   geom_hline(yintercept = 0, linetype = 3, color= "grey50") + 
   geom_pointrange() + 
   scale_color_manual(values = c("black", "grey40")) +
   geom_text(aes(label = sig), y = 18, size = 7) +
   facet_wrap( ~  sex) + 
   labs(x = NULL, y= "Pct. Change\n(95% CI)") +
   theme(axis.text.x = element_blank(), 
         legend.position = "none"))

# ferritin
(fer_coefplot <- ggplot(mod_ests_final %>% 
                          filter(outcome == "ferritin"), 
                        aes(x = cytokine_fmt,
                            y=est_pct, 
                            ymin = lower_pct,
                            ymax = upper_pct)) + 
    geom_hline(yintercept = 0, linetype = 3, color= "grey50") + 
    geom_pointrange() + 
    geom_text(aes(label = sig), y = 19, size = 7) +
    scale_color_manual(values = c("black", "grey50")) +
    facet_wrap( ~  sex) + 
    labs(x = NULL, y= "Pct. Change\n(95% CI)") +
    theme(axis.text.x = element_blank(), 
          legend.position = "none"))

# sTfR
(stfr_coefplot <- ggplot(mod_ests_final %>% 
                           filter(outcome == "stfr"), 
                         aes(x = cytokine_fmt,
                             y=est_pct, 
                             ymin = lower_pct,
                             ymax = upper_pct)) + 
    geom_hline(yintercept = 0, linetype = 3, color= "grey50") + 
    geom_pointrange() + 
    scale_color_manual(values = c("grey50")) +
    facet_wrap( ~  sex) + 
    labs(x = NULL, y= "Pct. Change\n(95% CI)") +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1, 
                                     vjust = 1), 
          legend.position = "none"))



# Figures --------------------------------------------

# Combine figures -------------------------------------------------
(figure_4 <- plot_grid(NULL, hepcidin_coefplot, 
                       NULL, fer_coefplot, 
                       NULL, stfr_coefplot,
                       ncol = 1, 
                       labels = c("A. Hepcidin", "", 
                                  "B. Ferritin", "", 
                                  "C. sTfR", ""),
                       hjust = 0,
                       rel_heights = c(0.2, 1.0,
                                       0.2, 1.0, 
                                       0.2, 1.5)))

# Save Figure 
ggsave(figure_4, 
       filename = fs::path(dir_xcfig, "Figre 4 coef plot.jpg"), 
       width = 7, height = 8)
