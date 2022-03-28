# Changes in ferritin and sTfR across seasons

# Set outcomes
outcomes <- c("ferritin", "stfr")

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
  filter(numDF != 1) %>% # Remove numDF = 1 (the intercept)
  janitor::clean_names() %>%
  mutate(sig = if_else(p_value < 0.05, "*", ""))
  
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
            shape = norm_value)) + 
   geom_hline(yintercept = 30, linetype = 2) + 
   # Individual points/lines
   geom_point(aes(group = id), alpha = .5, size = 1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   scale_shape_manual(values = c(1,19)) +
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
               filter(outcome == "ferritin", sex == "Male")) +
   scale_y_continuous(breaks = c(30, 60, 90, 120, 150, 180), 
                      limits = c(15,190)) +
   labs(x = NULL, y="Ferritin (ng/mL)") +
   facet_wrap(~sex, nrow = 1) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))

## sTfR ------------------------------------------
(stfr_fig <- xc_data %>% 
   mutate(`Normal Range` = if_else(stfr > (29.5), 
                               "Out of range", 
                               "Within range                                                      "
                               )) %>% 
   ### Plot -----------
 ggplot(aes(x = date_for_plotting, 
            y = stfr, 
            shape = `Normal Range`)) + 
   geom_hline(yintercept = 29.5, linetype = 2) + 
   # Individual points/lines
   geom_point(aes(group = id), alpha = .5, size = 1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   scale_shape_manual(values = c(1,19)) +
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
   labs(x = "Competitive Season and Timepoint",
        y="sTfR (nmol/L)") +
   facet_wrap(~sex, nrow = 1) +
   scale_x_discrete(breaks = as.character(1:8),
                    labels = c("1, EF", "1, LF", "1, ES", "1, LS",
                               "2, EF", "2, LF", "2, ES", "2, LS")) +
   theme(legend.position = "bottom",
         axis.text.x = element_text(angle = 90, vjust = .5)))



# Combine figures -------------------------------------------------
(figure_2 <- plot_grid(NULL, fer_fig, 
                       NULL, stfr_fig, 
                       ncol = 1, 
                       labels = c("A. Ferritin", "", 
                                  "B. Soluble Transferrin Receptor", ""),
                       hjust = 0,
                       align = "v",
                       rel_heights = c(0.1, 1.0,
                                       0.1, 1.5)))

# Save Figure 
ggsave(figure_2, 
       filename = fs::path(dir_xcfig, "Figure 2 Changes in iron status over time.jpg"), 
       width = 7, height = 7)


