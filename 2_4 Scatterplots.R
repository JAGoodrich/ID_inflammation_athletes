# Scatter Plots

# Ferritin vs. Hepcidin
ggplot(xc_data, 
       aes(x = ferritin, y = hepcidin, 
           group = id, color = id)) + 
  geom_vline(xintercept = 30,  linetype = 2, color = "grey50") +
  geom_point(aes(color = time_point)) + 
  geom_path(alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sex) + 
  theme(legend.position = "none")


# Ferritin vs. hgb ------------------------------------------------------------
lme(hgb ~ as.factor(date) + log(ferritin), 
    random = ~1|id, 
    correlation = corExp(form=~test_number|id),
    data = xc_data, 
    na.action = "na.exclude") %>% 
  tidy()

ggplot(xc_data, 
       aes(x = ferritin, y = hgb, 
           group = id, 
           color = sex, 
           shape = sex)) + 
  geom_vline(xintercept = 30,  linetype = 2, color = "grey50") +
  geom_point() + 
  geom_path(alpha = .15) +
  scale_color_manual(values = c("black", "grey50")) +
  scale_shape_manual(values = c(21, 18)) +
  scale_x_log10() +
  facet_wrap(~sex, scales = "free_x") +
  theme(legend.position = "none")


# Ferritin vs. IL6
ggplot(xc_data, 
       aes(x = ferritin, y = IL6, 
           group = id, color = id)) + 
  geom_vline(xintercept = 30,  linetype = 2, color = "grey50") +
  geom_point(aes(color = time_point)) + 
  geom_path(alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~sex) + 
  theme(legend.position = "none")



# Hepcidin vs. all inflammatory cytokines
ggplot(xc_data %>% 
         pivot_longer(cols = proteins) %>% 
         filter(hepcidin < 100), 
       aes(x = hepcidin, y = value+.01)) +
  geom_point(size = .5) + 
  geom_smooth() +
  geom_path(aes( group = id), alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(name~sex, scales =  "free") + 
  theme(legend.position = "none")


# ferritin vs. all inflammatory cytokines
ggplot(xc_data %>% 
         pivot_longer(cols = proteins), 
       aes(y = stfr,
           x = value+.01)) +
  geom_point(size = .5) + 
  geom_smooth() +
  geom_path(aes( group = id), alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(sex ~ name, scales =  "free") + 
  theme(legend.position = "none")
