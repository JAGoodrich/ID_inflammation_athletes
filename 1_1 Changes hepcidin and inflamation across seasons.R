# Changes in cytokines, hepcidin, and erythroferrone across study visits
proteins <- c("IFNg","IL1a","IL1b","IL2","IL4",
              "IL6","IL8","IL10","IL12p70","TNFa",
              "hepcidin", "ery")

# Statistics -----------------------------------------

## Run all models ------------------------------------
mod_output <- xc_data %>% 
  tidylog::pivot_longer(cols = all_of(proteins), 
                        names_to = "outcome", 
                        values_to = "outcome_value") %>% 
  filter(!is.na(outcome_value)) %>% 
  group_by(sex, outcome) %>% 
  nest() %>% 
  mutate(
    model = map(data, 
                ~lme(log(outcome_value) ~ as.factor(date), 
                     random = ~1|id, 
                     correlation = corExp(form=~test_number|id),
                     data = . )))

## Examine residual plots by hand (change "1" for different models) -------------
i = 1
plot(mod_output$model[[i]])
# Normal plot of standardized residuals by gender
qqnorm(mod_output$model[[i]], ~resid(., type = "p"))
# Normal plot of standardized residuals by gender
qqnorm(mod_output$model[[i]], ~ranef(.))

## Get p values for effect estimates -------------------
effect_est_p <- mod_output %>% 
  mutate(model_tidy = map(model, ~tidy(.))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(model_tidy) %>% 
  filter(str_detect(term, "as.factor\\(date\\)")) %>% # select "date" effects
  mutate(sig = if_else(p.value < 0.05, "*", ""), 
         date_chr = str_remove(term, "as.factor\\(date\\)")) %>% 
  select(sex, outcome, date_chr, sig)


## Calculate type-III/marginal F test -------------------
anova_results <- mod_output %>% 
  mutate(anova_results = map(model, 
                             ~car::Anova(.))) %>% 
  dplyr::select(-data, -model) %>% 
  unnest(anova_results) %>% 
  janitor::clean_names() %>% 
  rename(p_value = pr_chisq) %>%
  mutate(sig = if_else(p_value < 0.05, "Significant", ""))


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
                       unique(), 
                     by = "date") 


# Get location of asterisks 
mod_ests_final <- mod_ests_final %>% 
  mutate(sig_locaion = conf.high*1.1)


# Figures -------------------------------------------

## Hepcidin ------------------------------------------------
(hepcidin_fig <- ggplot(xc_data, 
                        aes(x = date_for_plotting, 
                            y = hepcidin, 
                            group = id)) +  
   geom_point(aes(group = id), 
              alpha = .5, 
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "hepcidin")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "hepcidin")) +
   scale_y_log10() +
   labs(x = NULL, y="Hepcidin (ng/mL)") +
   facet_wrap(~sex) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))


## IL-6 ------------------------------------------------
(il6_fig <- ggplot(xc_data, 
                   aes(x = date_for_plotting, 
                       y = IL6, 
                       group = id)) +  
   geom_point(aes(group = id), 
              alpha = .5, 
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "IL6")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL6")) + 
   geom_text(aes(x = date_for_plotting,
                 y = 70,
                 label = sig),
             size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL6")) +
   scale_y_log10() +
   labs(x = NULL, y="IL-6 (pg/ml)") +
   facet_wrap(~sex) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))


## IL-1a ------------------------------------------------
(il1a_fig <- ggplot(xc_data, 
                    aes(x = date_for_plotting, 
                        y = IL1a, 
                        group = id)) +  
   geom_point(aes(group = id), 
              alpha = .5, 
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "IL1a")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL1a")) + 
   geom_text(aes(x = date_for_plotting,
                 y = sig_locaion + 300,
                 label = sig),
             size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL1a")) + 
   scale_y_log10() +
   labs(x = NULL, y="IL-1\u03B1 (pg/ml)") +
   facet_wrap(~sex) +
   theme(legend.position = "none",  
         axis.text.x = element_blank()))


## IL-1b ------------------------------------------------
(il1b_fig <- ggplot(xc_data, 
                    aes(x = date_for_plotting, 
                        y = IL1b, 
                        group = id)) +  
   geom_point(aes(group = id), 
              alpha = .5, 
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting, 
                       ymin = conf.low,
                       ymax = conf.high, 
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>% 
                     filter(outcome == "IL1b")) +
   geom_line(aes(x = date_for_plotting, 
                 y = response, 
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL1b")) + 
   geom_text(aes(x = date_for_plotting,
                 y = sig_locaion + 80,
                 label = sig),
             size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>% 
               filter(outcome == "IL1b")  %>% 
               mutate(sig = if_else(sex == "Male", "", sig) )) + 
   scale_y_log10() +
   labs(x = "Competitive Season and Timepoint",
        y="IL-1\u03B2 (pg/ml)") +
   scale_x_discrete(breaks = as.character(1:8),
                    labels = c("1, EF", "1, LF", "1, ES", "1, LS",
                               "2, EF", "2, LF", "2, ES", "2, LS")) +
   facet_wrap(~sex) +
   theme(legend.position = "none",  
         axis.text.x = element_text(angle = 90, vjust = .5)))


## Erythroferrone ------------------------------------------------
(ery_fig <- ggplot(xc_data,
                   aes(x = date_for_plotting,
                       y = ery,
                       group = id)) +
   geom_point(aes(group = id),
              alpha = .5,
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting,
                       ymin = conf.low,
                       ymax = conf.high,
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_final %>%
                     filter(outcome == "ery")) +
   geom_line(aes(x = date_for_plotting,
                 y = response,
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_final %>%
               filter(outcome == "ery")) +
   geom_text(aes(x = date_for_plotting,
                 y = 800,
                 label = sig),
             size = 9,
             inherit.aes = FALSE,
             data = mod_ests_final %>%
               filter(outcome == "ery", !is.na(sig))) +
   scale_y_log10() +
   labs(x = "Competitive Season and Timepoint",
        y="ERFE (ng/ml)") +
   facet_wrap(~sex) +
   scale_x_discrete(breaks = as.character(1:8),
                    labels = c("1, EF", "1, LF", "1, ES", "1, LS",
                               "2, EF", "2, LF", "2, ES", "2, LS")) +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 90, vjust = .5)
   ))


# Combine figures -------------------------------------------------
figure_1 <- plot_grid(NULL, hepcidin_fig,
                      NULL, il6_fig, 
                      NULL, il1a_fig,
                      NULL, il1b_fig,
                      ncol = 1, 
                      labels = c("A. Hepcidin", "", 
                                 "B. IL-6", "", 
                                 "C. IL-1\u03B1", "", 
                                 "D. IL-1\u03B2", ""
                                 ),
                      hjust = 0,
                      align = "v",
                      rel_heights = c(0.2, 1.0,
                                      0.2, 1.0,
                                      0.2, 1.0,
                                      0.2, 1.5))

# Save Figure ----------------
ggsave(figure_1, 
       filename = fs::path(dir_xcfig, 
                           "Figure 1 Changes in cytokines across time.jpg"), 
       width = 6, height = 9)









# Supplemental figure 1: All cytokines ------------------------------------

mod_ests_reduced <- mod_ests_final  %>% 
  filter(!(outcome %in% c("ery", "hepcidin"))) %>% 
  mutate(sig.location = 2.5*conf.high, 
         sig = if_else(outcome == "IL1b" & sex == "Male", "", sig))
  
# pivot outcome data longer on outcomes 
xc_data_l <- xc_data %>% 
  pivot_longer(cols = all_of(proteins), 
               names_to = "outcome", 
               values_to = "conc") %>% 
  filter(!(outcome %in% c("ery", "hepcidin")))

## all cytokines ------------------------------------------------
(supplemental_figure_1 <- ggplot(xc_data_l,
                   aes(x = date_for_plotting,
                       y = conc,
                       group = id)) +
   geom_point(aes(group = id),
              alpha = .5,
              size = .75, shape =1) +
   geom_line(aes(group = id), color = "grey50", alpha = .5) +
   # Summary points/lines
   geom_pointrange(aes(x = date_for_plotting,
                       ymin = conf.low,
                       ymax = conf.high,
                       y = response),
                   inherit.aes = FALSE,
                   data = mod_ests_reduced) +
   geom_line(aes(x = date_for_plotting,
                 y = response,
                 group = 1),
             inherit.aes = FALSE,
             data = mod_ests_reduced) +
   geom_text(aes(x = date_for_plotting,
                 y = sig.location,
                 label = sig),
             size = 9,
             inherit.aes = FALSE,
             data = mod_ests_reduced %>%
               filter(!is.na(sig))) +
   scale_y_log10() +
   labs(x = "Competitive Season and Timepoint",
        y="Concentration (pg/ml)") +
   facet_grid(outcome~sex, scales = "free_y") +
   scale_x_discrete(breaks = as.character(1:8),
                    labels = c("1, EF", "1, LF", "1, ES", "1, LS",
                               "2, EF", "2, LF", "2, ES", "2, LS")) +
   theme(legend.position = "none",
         axis.text.x = element_text(angle = 90, vjust = .5)))


# Save Figure ----------------
ggsave(supplemental_figure_1, 
       filename = fs::path(dir_xcfig, 
                           "Figure S1 Changes in all cytokines across time.jpg"), 
       width = 6, height = 14)

