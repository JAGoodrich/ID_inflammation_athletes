

library(ggh4x)


# Figure of ferritin over time
# ggplot(data=xc_data, 
#        aes(x = time_point,
#            group = id, 
#            y = ferritin, 
#            shape = sex, color = sex))+
#   geom_point(alpha = .5) + 
#   geom_line( alpha = .5) + 
#   labs(x = NULL, y="Ferritin (ng/mL)") +
#   geom_ribbon(data = create_norm_range(30, 100, 4, .5, 
#                                        offset = .75, sex = "Male"), 
#               aes(x = tp, ymin = min, ymax = max),
#               inherit.aes = FALSE)+
#   geom_ribbon(data = create_norm_range(15, 100, 4, .5, 
#                                        offset = 1.25, sex = "Female"), 
#               aes(x = tp, ymin = min, ymax = max),
#               inherit.aes = FALSE, color = 'grey', fill = 'grey')+
#   scale_colour_manual(values=c("black", "grey30")) + 
#   theme(plot.title = element_text(size = 25), 
#         # axis.text.x = element_blank()) +
#         axis.text.x = element_text(angle = 45, 
#                                    hjust = 1, 
#                                    size = 15)) +
#   # scale_y_log10() +
#   facet_nested(id ~ ay)


#  Ferritin -----------------------------------------------------
## Females ---------------------
(females_ferritin <- ggplot(data=xc_data %>% 
                              filter(sex == "Female",  
                                     !is.na(ferritin)), 
                            aes(x = date,
                                group = ay, 
                                y = ferritin, 
                                shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   geom_hline(yintercept = 25, linetype = 2, color = "grey80") +
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    hjust = 1,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   labs(x = NULL, y="Ferritin (ng/mL)", title = "A. Females, Ferritin") +
   facet_wrap(~id))

# ggsave(females_ferritin, 
#        filename = fs::path(dir_xcfig, "Female Individual Ferritin Changes.jpeg"), 
#        width = 7, height = 8.5)

## Males ---------------------
(males_ferritin <- ggplot(data=xc_data %>% 
                            filter(sex == "Male",
                                   !is.na(ferritin)), 
                          aes(x = date,
                              group = ay, 
                              y = ferritin, 
                              shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   geom_hline(yintercept = 25, linetype = 2, color = "grey80") +
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    hjust = 1,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   labs(x = NULL, y="Ferritin (ng/mL)", title = "A. Males, Ferritin") +
   facet_wrap(~id))

ggsave(males_ferritin, 
       filename = fs::path(dir_xcfig, "Males Individual Ferritin Changes.jpeg"), 
       width = 7, height = 8.5)



#  Hepcidin -----------------------------------------------------
## Females ---------------------
(females_hepcidin <- ggplot(data=xc_data %>%
                              filter(sex == "Female",
                                     !is.na(hepcidin)), 
                            aes(x = date,
                                group = ay, 
                                y = hepcidin, 
                                shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    hjust = 1,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   scale_y_log10() + 
   labs(x = NULL, y="Hepcidin (ng/mL)", title = "B. Females, Hepcidin") +
   facet_wrap(~id))

ggsave(females_hepcidin, 
       filename = fs::path(dir_xcfig, "Female Individual hepcidin Changes.jpeg"), 
       width = 7, height = 8.5)

## Males ---------------------
(males_hepcidin <- ggplot(data=xc_data %>% 
                            filter(sex == "Male", 
                                   !is.na(hepcidin)), 
                          aes(x = date,
                              group = ay, 
                              y = hepcidin, 
                              shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   labs(x = NULL, 
        y ="Hepcidin (ng/mL)", 
        title = "B. Males, Hepcidin") +
   # scale_y_log10() +
   facet_wrap(~id))

ggsave(males_hepcidin, 
       filename = fs::path(
         dir_xcfig, 
         "Males Individual hepcidin Changes.jpeg"), 
       width = 7, height = 8.5)



# IL6 ---------------------------------
## Females ---------------------

(females_il6 <- ggplot(data=xc_data %>%
                         filter(sex == "Female", 
                                !is.na(IL6)), 
                       aes(x = date,
                           group = ay, 
                           y = IL6, 
                           shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    hjust = 1,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   scale_y_log10() + 
   labs(x = NULL, y="IL-6 (pg/mL)", 
        title = "C. Females, IL-6") +
   facet_wrap(~id))

ggsave(females_il6, 
       filename = fs::path(dir_xcfig, "Female Individual IL-6 Changes.jpeg"), 
       width = 7, height = 8.5)

## Males ---------------------
(males_il6 <- ggplot(data=xc_data %>% 
                       filter(sex == "Male"), 
                     aes(x = date,
                         group = ay, 
                         y = IL6, 
                         shape = iron_def_acute))+
   geom_point(aes(color=iron_def_acute, group = ay)) + 
   geom_line() + 
   scale_colour_brewer(palette = "Set1", direction = -1) +
   theme(plot.title = element_text(size = 25), 
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 90,
                                    vjust = .5,
                                    size = 15), 
         panel.background = element_rect(fill = 'grey95'), 
         legend.title = element_blank()) +
   labs(x = NULL, y="IL-6 (pg/mL)", 
        title = "C. Males, IL-6") +
   scale_y_log10() +
   facet_wrap(~id))

ggsave(males_il6, 
       filename = fs::path(
         dir_xcfig, 
         "Males Individual il6 Changes.jpeg"),
       width = 7, height = 8.5)




