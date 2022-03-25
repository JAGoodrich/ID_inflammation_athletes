# Statistical Analysis
library(effectsize)
#Males:
mymod <- lm(log(hepcidin) ~ iron_def_once * sex, data = id )
  
effectsize(mymod)
anova_table <- anova(mymod)
effectsize(anova_table)

Anova(mymod)
# summary(mymod)
confint(emmeans(mymod, pairwise ~ iron_def_once, type = "response"))

#Females:

table(id$acute_chronic_id, id$sex)

mymod <- lm(log(stfr) ~  as.ordered(acute_chronic_id), data = id %>% filter(sex == "Female") %>% droplevels())
Anova(mymod)
summary(mymod)
emmeans(mymod, pairwise ~ acute_chronic_id, type = "response")
confint(emmeans(mymod, pairwise ~ acute_chronic_id, type = "response"))



#Hepcidin
hist(fsd$stfr)
mymod <- lm(log(ery+1) ~  (iron_def_acute), data = id %>% filter(sex == "Female", !is.nan(gamma(ery_cor+1))) %>% droplevels())
Anova(mymod)
mymod <- update(ref_grid(mymod), tran = make.tran("power", -.5))

# summary(mymod)
emmeans(mymod, pairwise ~ iron_def_acute, type = "response")
# 

table(id$iron_def_once, id$ferritin < 25)


temp <- lm((percent_fm) ~  nationals, data = fsd %>% filter(sex == "Female"))
emmeans(temp, pairwise ~ nationals, type = "response")

table(id$hepcidin < 20,  id$iron_def_once, id$sex)

npoints <- data %>% 
  group_by(id)%>% 
  summarise(n = n_groups(.)) %>% 
  ungroup()


summary(lmer(hepcidin ~ IL6 + (1|id), data = data))
