# Iron Status Changes across time
colnames(xc_data)

outcomes <- c("ferritin", "stfr_logfer", "stfr", "hgb")

# Statistics -----------------------------------------

## Run all models ------------------------------------
mod_output <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(outcomes), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  filter(!is.na(outcome_value)) %>% 
  mutate(date = relevel(as.factor(date),
                 ref = "Oct 2017")) %>%
  group_by(sex, outcome) %>% 
  nest() %>% 
  mutate(
    # Run iron Models
    model = map(data, 
                ~lme((outcome_value) ~ date, 
                     random = ~1|id, 
                     correlation = corExp(form=~test_number|id),
                     data = .)))


## Get p values for effect estimates -------------------
effect_est_p <- mod_output %>% 
  mutate(model_tidy = map(model, ~tidy(.))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(model_tidy) %>% 
  filter(str_detect(term, "date")) %>% # select "date" effects
  mutate(sig = if_else(p.value < 0.05, "*", ""), 
         date_chr = str_remove(term, "date")) %>% 
  select(sex, outcome, date_chr, sig)


## Calculate type-III/marginal F test -------------------
anova_results <- mod_output %>% 
  mutate(anova_results = map(model, 
                             ~anova(.,type='marginal'))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(anova_results) %>% 
  filter(numDF != 1) # Remove numDF = 1 (the intercept)

## Obtain model estimates ------------------------------
mod_ests <- mod_output %>% 
  mutate(mod_ests = map(model, 
                        ~emmeans(.,~date, type = "response") %>% 
                          tidy(conf.int = TRUE, 
                               conf.level = 0.95))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(mod_ests) %>% 
  mutate(date_chr = as.character(date))

# Get p_values and date for plotting
mod_ests_final <- mod_ests %>% 
  tidylog::left_join(effect_est_p) %>% 
  tidylog::left_join(xc_data %>% 
                       select(date, visit_ordered, date_for_plotting) %>% 
                       unique() %>% 
                       mutate(date = as.character(date)), 
                     by = "date") %>% 
  mutate(sig = replace_na(sig, ""))

# Figures --------------------------------------------

## Ferritin ------------------------------------------
(fer_fig <- xc_data %>% 
   mutate(norm_value = if_else(ferritin < 30, 
                               "Low", "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = ferritin, 
            color = norm_value)) + 
   geom_hline(yintercept = 30, linetype = 2) + 
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
             y = 180, size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "ferritin")) +
   scale_y_continuous(breaks = c(30, 60, 90, 120, 150, 180), 
                      limits = c(15,190)) +
   labs(x = NULL, y="Ferritin (ng/mL)") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))

## sTfR ------------------------------------------
(stfr_fig <- xc_data %>% 
   mutate(norm_value = if_else(stfr > (1.8/0.0738), 
                               "Elevated", "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = stfr, 
            color = norm_value)) + 
   geom_hline(yintercept = (1.8/0.0738), linetype = 2) + 
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
   mutate(norm_value = if_else(stfr_logfer > 1.4, 
                               "Elevated", "Normal")) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = stfr_logfer, 
            color = norm_value)) + 
   geom_hline(yintercept = 1.4, linetype = 2) + 
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

## Hemoglobin ------------------------------------------
(hgb_fig <- xc_data %>% 
 ggplot(aes(x = date_for_plotting, 
            y = hgb)) + 
   # geom_hline(yintercept = 30, linetype = 2) + 
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
                     filter(outcome == "hgb")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "hgb")) +
   geom_text(aes(x = date_for_plotting, label = sig),
             y = 18, size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "hgb")) +
   # scale_y_continuous(breaks = c(30, 60, 90, 120, 150, 180), 
   #                    limits = c(15,190)) +
   labs(x = NULL, y="Hemoglobin (g/dL)") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))

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


