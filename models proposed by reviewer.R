fsd2 <- data %>% 
  group_by(id) %>% 
  summarise(mean_ferritin = mean(ferritin,
                                 na.rm = TRUE), 
            sex = sex[1])


table(fsd$iron_def_once, fsd2$sex)



data <- data %>% 
  group_by(id) %>% 
  mutate(delta_ferritin = lead(ferritin)-ferritin)

summary(lmer(delta_ferritin ~ log2(hepcidin) + (1|id),
             data = data))



library(nlme)

data2 <- data %>% 
  filter(!is.na(time_point), 
         !is.na(ferritin), 
         !is.na(stfr_fer)) %>% 
  mutate(stfr_mg_l = stfr*0.0738, 
         stfr_fer = (stfr_mg_l)/log10(ferritin))


# Model
mod_out <- lme((stfr_fer_new) ~ time_point + ay, # + log(hepcidin),
      random = ~1|id, 
      weights = varIdent(form=~1|sex),
      correlation = corExp(form=~test_number|id),
    data = data2)
summary(mod_out)

BIC(mod_out)
sjPlot::tab_model(mod_out)

hist(data2$stfr_fer)

table(data2$stfr_fer>1.5, data2$ferritin < 25)

plot(data2$stfr_mg_l, log(data2$hepcidin))

ggplot(data, aes(x = stfr_fer)) +
  geom_histogram(fill = "grey90", color = "black") + facet_wrap(~sex)

plot(mod_out)

summary(data$stfr*0.0738)




hist(data2$hepcidin)

data2$hgb

summary(mod_out <- lme((hgb) ~ time_point + ay + sex + log(hepcidin),
               random = ~1|id, 
               weights = varIdent(form=~1|sex),
               correlation = corExp(form=~test_number|id),
               data = data2 %>% filter(!is.na(hgb))))
