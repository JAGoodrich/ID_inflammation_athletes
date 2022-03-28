# Scatter Plots


# Statistics ----------------------
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
                   ~full_join(.x, .y))) %>% 
  select(-c(data, model_full, mod_summary, mod_ests)) %>%
  unnest(coef_ci) %>% 
  ungroup() %>%
  clean_names() %>% 
  filter(str_detect(effect, "log2"))


#View info about the one outlier participant
# View(xc_data %>% select(id, visit,sex, hepcidin, ferritin, stfr, id_peeling))


# Ferritin vs. Hepcidin -------------------------
(fer_hep <- xc_data %>% 
   mutate(normrange = if_else(ferritin < 30, 
                              "Out of range", 
                              "Within range")) %>% 
   ggplot(aes(x = hepcidin, 
              y = ferritin, 
              group = id, 
              shape = normrange, 
              alpha = normrange)) + 
   geom_hline(yintercept = 30,  linetype = 2, color = "grey50") +
   geom_point() + 
   geom_path(alpha = .25) +
   scale_x_log10() +
   scale_y_log10() +
   scale_shape_manual(values = c(1,19)) +
   scale_alpha_manual(values = c(.75, .5)) +
   xlab("Hepcidin (ng/mL)") +
   ylab("Ferritin (ng/mL)") +
   facet_wrap(~sex) + 
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