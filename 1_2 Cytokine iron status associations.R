# XC Cytokine Outcome Analysis
outcomes <- c("ferritin", "stfr_logfer", "stfr")


proteins <- c("IFNg","IL1a","IL1b","IL2","IL4",
              "IL6","IL8","IL10","IL12p70","TNFa",
              "hepcidin", "ery")

# Statistics -----------------------------------------

## Run all models ------------------------------------
mod_output <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(outcomes), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  tidylog::pivot_longer(cols = all_of(str_c(proteins, "_pro_quantnum")),
                        names_to = "cytokine", 
                        values_to = "cytokine_concentration") %>% 
  filter(!is.na(outcome_value), 
         !is.na(cytokine_concentration)) %>% 
  droplevels() %>% 
  group_by(sex, outcome, cytokine) %>% 
  nest() %>% 
  mutate(
    # Run null model
    model_null = map(data, 
                     ~lme(scale(outcome_value) ~ as.factor(date), 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          data = . )), 
    # Run full model
    model_full = map(data, 
                     ~lme(scale(outcome_value) ~ as.factor(date)+cytokine_concentration, 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          data = ., 
                          na.action = na.exclude )), 
    # Compare models
    anova_results = map2(model_full, 
                         model_null, 
                         ~lmtest::lrtest(.x,.y )))


(xxx <- mod_output %>% 
    mutate(est = map(model_full,
                     ~tidy(.))) %>% 
  unnest(est) %>% 
    select(-c(data, model_null, model_full, anova_results)) %>%
  filter(str_detect(term, "cytokine"), 
         p.value < 0.5))

## Get p values for effect estimates -------------------
(anova_results <- mod_output %>% 
   mutate(term = map(anova_results, ~as_tibble(., rownames = "term"))) %>% 
   dplyr::select(-data, -model_null, -model_full, -anova_results) %>% 
   unnest(c(term)) %>% 
   clean_names() %>% 
   filter(!is.na(pr_chisq)) %>% 
   ungroup() %>%
   mutate(fdr = p.adjust(pr_chisq))) %>% 
  filter(pr_chisq < 0.1)


## Obtain model estimates ------------------------------
mod_ests <- mod_output %>% 
  mutate(mod_ests = map(model_full, 
                        ~emmeans(.,~cytokine_concentration, 
                                 type = "response") %>% 
                          tidy(conf.int = TRUE, 
                               conf.level = 0.95))) %>% 
  dplyr::select(-data, -model_null, -model_full, -anova_results) %>% 
  unnest(mod_ests) 

# Get p_values and date for plotting
mod_ests_final <- mod_ests %>%
  tidylog::left_join(anova_results, by = c("cytokine", "sex", "outcome"))
#   tidylog::left_join(xc_data %>% 
#                        select(date, visit_ordered, date_for_plotting) %>% 
#                        unique(), 
#                      by = "date")


# Plot data ---------------------------------------
# ggplot(xc_data %>% filter(!is.na(IL6)), 
#        aes(x = quantcut(IL6), y=log(ferritin))) + 
#   geom_boxplot() + 
#   facet_wrap(~sex)
#
#
ggplot(mod_ests_final %>% filter(sex == "Female"), 
       aes(x = cytokine_concentration,
           y=response, ymin = conf.low, ymax = conf.high,
           color = pr_chisq < 0.05)) + 
  geom_pointrange() + 
  geom_line(aes(group = 1)) +
  facet_grid(outcome ~ cytokine, scales = "free")
#

ggplot(mod_ests_final %>% filter(sex == "Male"), 
       aes(x = cytokine_concentration,
           y=response, ymin = conf.low, ymax = conf.high,
           color = pr_chisq < 0.05)) + 
  geom_pointrange() + 
  geom_line(aes(group = 1)) +
  facet_grid(outcome ~ cytokine, scales = "free")
#

















# Figures --------------------------------------------

## Ferritin ------------------------------------------
(fer_fig <- xc_data %>% 
   mutate(norm_value = if_else(ferritin < 25, 
                               "Low", "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = ferritin, 
            color = norm_value)) + 
   geom_hline(yintercept = 25, linetype = 2) + 
   # Individual points/lines
   geom_point(aes(group = id), alpha = .5, size = .75) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   scale_color_manual(values = c("red", "grey50")) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "ferritin")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "ferritin")) +
   geom_text(aes(x = date_for_plotting, label = sig),
             y = 90, size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "ferritin")) +
   labs(x = NULL, y="Ferritin (ng/mL)") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))

## sTfR ------------------------------------------
(stfr_fig <- xc_data %>% 
   mutate(norm_value = if_else(stfr > 21, 
                               "Elevated",
                               "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = stfr, 
            color = norm_value)) + 
   geom_hline(yintercept = 21, linetype = 2) + 
   # Individual points/lines
   geom_point(aes(group = id), alpha = .5, size = .75) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   scale_color_manual(values = c("red", "grey50")) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "stfr")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "stfr")) +
   geom_text(aes(x = date_for_plotting, label = sig),
             y = 35, size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "stfr")) +
   labs(x = NULL, y="sTfR (nmol/L)") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))


## sTfR-F index ------------------------------------------
(stfr_f_fig <- xc_data %>% 
   mutate(norm_value = if_else(stfr_logfer > 1.03, 
                               "Elevated", "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = stfr_logfer, 
            color = norm_value)) + 
   geom_hline(yintercept = 1.03, linetype = 2) + 
   # Individual points/lines
   geom_point(aes(group = id), alpha = .5, size = .75) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   scale_color_manual(values = c("red", "grey50")) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "stfr_logfer")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "stfr_logfer")) +
   # Asterisks for significance
   geom_text(aes(x = date_for_plotting, label = sig),
             y = 1.8, size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "stfr_logfer")) +
   labs(x = "Time point", 
        y="sTfR-F index") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none"))


# Combine figures -------------------------------------------------
(figure_1 <- plot_grid(NULL, fer_fig, 
                       NULL, stfr_fig, 
                       NULL, stfr_f_fig,
                       ncol = 1, 
                       labels = c("A. Ferritin", "", 
                                  "B. sTfR", "", 
                                  "C. sTfR Ferritin Index", ""),
                       hjust = 0,
                       rel_heights = c(0.2, 1.0,
                                       0.2, 1.0, 
                                       0.2, 1.0)))

# Save Figure 
ggsave(figure_1, 
       filename = fs::path(dir_xcfig, "Changes in iron status over time.jpg"), 
       width = 6, height = 9)
