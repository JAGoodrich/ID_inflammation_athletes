library(MuMIn)
table(data$iron_def_acute)

# Changes in Ferritin vs. hepcidin
data <- data %>% 
  group_by(id, ay) %>% 
  arrange(date) %>%
  mutate(delta_ferritin = lead(ferritin) - ferritin, 
         delta_fer_pct = 100*(delta_ferritin)/ferritin) %>% 
  ungroup() %>% 
  mutate(change_in_id_next_visit = 
           case_when(
             iron_def_acute == "ID" & lead(iron_def_acute) == "ID" ~ "Stayed ID", 
             iron_def_acute == "ID" & lead(iron_def_acute) == "Not ID" ~ "Resolved ID",
             iron_def_acute == "Not ID" & lead(iron_def_acute) == "ID" ~ "Became ID",
             iron_def_acute == "Not ID" & lead(iron_def_acute) == "Not ID" ~ "Stayed non ID"),
         delta_ferritin_cat = 
           case_when(
             delta_ferritin < -20 ~ "Decreased >20 ng/dL", 
             delta_ferritin < -10 ~ "Decreased >10 ng/dL", 
             delta_ferritin <  0  ~ "Decreased > 1 ng/dL", 
             delta_ferritin >= 0  ~ "Increased") %>% 
           fct_relevel("Decreased >20 ng/dL",
                       "Decreased >10 ng/dL", 
                       "Decreased > 1 ng/dL", 
                       "Increased"), 
         delta_fer_pct_cat = 
           case_when(
             delta_fer_pct < -20 ~ "Decreased > 25%", 
             delta_fer_pct < -10 ~ "Decreased > 10%", 
             delta_fer_pct <  0  ~ "Decreased > 1%", 
             delta_fer_pct >= 0  ~ "Increased") %>% 
           fct_relevel("Decreased > 25%",
                       "Decreased > 10%", 
                       "Decreased > 1%", 
                       "Increased"))

quantile(data$delta_fer_pct, probs = c(.25, 1/3, .5, 2/3, .75), na.rm = TRUE)
table(data$delta_ferritin_cat, data$sex)

# Plot Change in Ferritin vs. Hepcidin --------------------------
(deltafer_vs_hep <- ggplot(data %>% filter(!is.na(delta_ferritin_cat)), 
       aes(x=delta_fer_pct_cat, y = hepcidin)) + 
  geom_boxplot() + 
  facet_wrap(~sex, scales = "free") +
  scale_y_log10()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Change in Ferritin at next visit") + 
  ylab("Hepcidin at Current Visit"))

ggsave(deltafer_vs_hep, filename = fs::path(dir_xcfig, "deltafer_vs_hep.jpeg"))



# IL6 and Hepcidin ----------------------------
summary(lmer(hepcidin ~ log2(IL6) +   visit*sex + (1|id), data = data))



# Hepcidin and delta Ferritin ------------------------------
hep_vs_deltafer <- lmer(delta_ferritin ~ log2(hepcidin)*scale(ferritin) +
                          sex + 
                          splines2::bSpline(scale(as.numeric(date))) + 
                          (1|id), 
                        data = data %>% filter(!is.na(hepcidin)))

fer_vs_deltafer <- lmer(delta_ferritin ~ scale(ferritin) +
                          sex + 
                          splines2::bSpline(scale(as.numeric(date))) + 
                          (1|id), 
                        data = data %>% filter(!is.na(hepcidin)))


anova(hep_vs_deltafer, fer_vs_deltafer)
summary(hep_vs_deltafer)
MuMIn::r.squaredGLMM(hep_vs_deltafer)



# IL6 and delta Ferritin ------------------------------
hep_vs_deltafer <- lmer(delta_ferritin ~ log2(IL6)*scale(ferritin) +
                          sex + 
                          splines2::bSpline(scale(as.numeric(date))) + 
                          (1|id), 
                        data = data %>% filter(!is.na(hepcidin)))

fer_vs_deltafer <- lmer(delta_ferritin ~ scale(ferritin) +
                          sex + 
                          splines2::bSpline(scale(as.numeric(date))) + 
                          (1|id), 
                        data = data %>% filter(!is.na(hepcidin)))


anova(hep_vs_deltafer, fer_vs_deltafer)
summary(hep_vs_deltafer)
MuMIn::r.squaredGLMM(hep_vs_deltafer)



ggplot(data, aes(x = hepcidin +.01, y = delta_ferritin, color = sex)) + 
  geom_point() + 
  stat_smooth(method = "lm") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  scale_x_log10() +
  facet_wrap(~sex) + 
  labs(y = "Change in Ferritin at Next Visit" 
  )




# il6 and delta Ferritin ------------------------------------
il6_vs_deltafer <- lmer(delta_ferritin ~ hepcidin + log(IL6) + sex + visit + (1|id), 
                        data = data )
summary(il6_vs_deltafer)

