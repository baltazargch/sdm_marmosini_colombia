library(tidyverse)
library(wesanderson)
dir.create('output/supp_materials/', showWarnings = F)
plots <- list.files('output/supp_materials/', pattern = '.png')
# if(length(plots) == 8) stop('Ya están todos los gráficos realizados.') 
fit_table <- read.csv('output/models/marmosini_fitting_results.csv')

less80 <- c('L', 'LQ', 'LQP')
more80 <- c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')

fit_table %>% 
  na.omit() %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>% 
  mutate(features = factor(features, levels = more80)) %>% 
  mutate(rm = factor(as.character(rm), levels = unique(as.character(rm)))) %>% 
  group_by(case, area, features, cross.validation, rm) -> fit_table

fit_table%>% 
  mutate(mean.auc = mean(AICc)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean AICc") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.meanAICc.gg

fit_table %>% 
  mutate(mean.auc = mean(avg.test.AUC)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.test.AUC") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.test.AUC

fit_table %>% 
  mutate(mean.auc = mean(avg.diff.AUC)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.diff.AUC") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.diff.AUC

fit_table %>% 
  mutate(mean.auc = mean(avg.test.orMTP)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.test.orMTP") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.test.orMTP

gg.list <- list(Average_diff_AUC = gg.avg.diff.AUC, 
                Average_test_AUC = gg.avg.test.AUC, 
                Average_test_orMTP = gg.avg.test.orMTP, 
                Average_AICc = gg.meanAICc.gg)

for(i in seq_along(gg.list)){
  ggsave(filename = paste0('output/supp_materials/full_', names(gg.list[i]), '.png'),
         plot = gg.list[[i]], device = 'png', width = 28, height = 32, units = 'cm')
}

fit_table <- read.csv('output/models/marmosini_fitting_results.csv')

filterted_table <- fit_table %>% 
  na.omit() %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>% 
  mutate(features = factor(features, levels = more80)) %>% 
  mutate(
    type = case_when(
      species == 'Marmosa_robinsoni' & features  %in% more80 ~ "in",
      species != 'Marmosa_robinsoni' & features  %in% less80 ~ "in",
      species == 'Marmosa_robinsoni' & !features  %in% more80 ~ "out",
      species != 'Marmosa_robinsoni' & !features  %in% less80 ~ "out"
    )
  ) %>% 
  filter(type == 'in') %>% 
  mutate(rm = factor(as.character(rm), levels = unique(as.character(rm)))) %>% 
  group_by(case, area, features, cross.validation, rm) 

filterted_table%>% 
  mutate(mean.auc = mean(AICc)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean AICc") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.meanAICc.gg

filterted_table %>% 
  mutate(mean.auc = mean(avg.test.AUC)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.test.AUC") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.test.AUC

filterted_table %>% 
  mutate(mean.auc = mean(avg.diff.AUC)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.diff.AUC") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.diff.AUC

filterted_table %>% 
  mutate(mean.auc = mean(avg.test.orMTP)) %>% 
  arrange(mean.auc) %>% 
  ggplot(aes(shape=area, color=case, 
             group=method, y=mean.auc,
             x=rm)) + 
  geom_line(show.legend = T, lwd = 0.2) +
  geom_point(cex=1.8) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(features ~ cross.validation, scales = 'free_y', ncol=3, 
             labeller = label_value,  strip.position =  "right", drop = F)  +
  labs(title = '', x='', 
       y= "Mean avg.test.orMTP") +
  theme_light() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        strip.background = element_blank(), 
        strip.placement = "outside", 
        strip.text = element_text(colour = 'gray50')) -> gg.avg.test.orMTP

gg.list <- list(Average_diff_AUC = gg.avg.diff.AUC, 
                Average_test_AUC = gg.avg.test.AUC, 
                Average_test_orMTP = gg.avg.test.orMTP, 
                Average_AICc = gg.meanAICc.gg)

for(i in seq_along(gg.list)){
  ggsave(filename = paste0('output/supp_materials/filtered_', names(gg.list[i]), '.png'),
         plot = gg.list[[i]], device = 'png', width = 28, height = 32, units = 'cm')
}


