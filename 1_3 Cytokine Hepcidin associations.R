# Cytokine --> hepcidin models


# Outcome: Hgb, 
outcomes <- c("hepcidin")

proteins <- c("IFNg","IL1a","IL1b","IL2","IL4",
              "IL6","IL8","IL10","IL12p70","TNFa", "ery")


# Statistics -----------------------------------------
## Run all models ------------------------------------
cytokine_hepcidin_associations <- xc_data %>% 
  tidylog::pivot_longer(
    cols = all_of(str_c(proteins,"_pro_quantnum")), #_pro_quantnum
    names_to = "cytokine", 
    values_to = "cytokine_concentration") %>% 
  filter(!is.na(hepcidin), 
         !is.na(ferritin),
         !is.na(cytokine_concentration)) %>% 
  droplevels() %>% 
  group_by(cytokine) %>% 
  nest() %>% 
  mutate(
    # Run null model
    model_null = map(data, 
                     ~lme(log(hepcidin) ~ as.factor(date)+sex,#+log(ferritin), 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          weights=varIdent(form=~1|sex),
                          data = . )), 
    # Run full model
    model_full = map(data, 
                     ~lme(log(hepcidin) ~ as.factor(date)+sex+cytokine_concentration, #log(ferritin)+ 
                          random = ~1|id, 
                          correlation = corExp(form=~test_number|id),
                          weights=varIdent(form=~1|sex),
                          data = ., 
                          na.action = na.exclude )), 
    # Compare models
    anova_results = map2(model_full, 
                         model_null, 
                         ~lmtest::lrtest(.x,.y )))


(xxx <- cytokine_hepcidin_associations %>% 
    mutate(est = map(model_full,
                     ~tidy(.))) %>% 
    unnest(est) %>% 
    select(-c(data, model_null, model_full, anova_results)) %>%
    filter(str_detect(term, "cytokine"), 
           p.value < 0.5))


## Get p values for effect estimates -------------------
anova_results <- cytokine_hepcidin_associations %>%
   mutate(term = map(anova_results, ~as_tibble(., rownames = "term"))) %>%
   dplyr::select(-data, -model_null, -model_full, -anova_results) %>%
   unnest(c(term)) %>%
   clean_names() %>%
   filter(!is.na(pr_chisq)) %>%
   ungroup() %>%
   mutate(fdr = p.adjust(pr_chisq))


## Obtain model estimates ------------------------------
## For categorical: ------------------------------------
# mod_ests <- cytokine_hepcidin_associations %>% 
#   mutate(mod_ests = map(model_full, 
#                         ~emmeans(.,~cytokine_concentration, 
#                                  type = "response") %>% 
#                           tidy(conf.int = TRUE, 
#                                conf.level = 0.95))) %>% 
#   dplyr::select(-data, -model_null, -model_full, -anova_results) %>% 
#   unnest(mod_ests) 


## For continuous: ------------------------------------
mod_ests <- cytokine_hepcidin_associations %>%
  mutate(mod_ests = map(model_full,
                        ~intervals(.,
                              which = "fixed")$fixed %>% 
                          as_tibble(rownames = "term")), 
         p_value  = map(model_full,
                        ~coef(.,
                                   which = "fixed") %>% 
                          as_tibble(rownames = "term")), 
         output = map2(mod_ests, p_value, ~full_join(.x, .y))
         ) %>%
  dplyr::select(-data, -model_null, -model_full, -anova_results, 
                -mod_ests, -p_value) %>%
  unnest(output) %>% 
  janitor::clean_names() %>% 
  filter(term == "cytokine_concentration")

xxx<-summary(cytokine_hepcidin_associations$model_full[[1]])
coef(xxx, which = "fixed")

# Get p_values and date for plotting
# mod_ests_final <- mod_ests %>%
#   tidylog::left_join(anova_results %>% select(-term), by = c("cytokine"))

# format names of cytokines
mod_ests_final <- mod_ests %>% 
  mutate(cytokine_fmt = str_remove_all(cytokine, "_pro_quantnum") %>% 
           str_replace("ery", "ERFE"))


# Plot data ---------------------------------------
# ggplot(xc_data %>% filter(!is.na(IL6)), 
#        aes(x = quantcut(IL6), y=log(ferritin))) + 
#   geom_boxplot() + 
#   facet_wrap(~sex)
#
#
ggplot(mod_ests_final , 
       aes(x = cytokine_fmt,
           y = est,
           ymin = lower,
           ymax = upper)) + 
  geom_pointrange() + 
  geom_hline(yintercept = 0, linetype = 2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

  

