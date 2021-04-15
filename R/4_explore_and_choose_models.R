library(tidyverse)
source('R/3_fitmodels.R')

if(!file.exists('output/models/marmosini_fitting_results.csv')){
  
  fit_table <- data.frame()
  
  for(cs in names(cases)){
    case.table <- data.frame()
    for(sp in names(OCCS)){
      dir.sp <- paste0('output/models/', cs, '/fitting/',sp)
      
      sp.table <- data.frame()
      for(d in dir.sp){
        t <- list.files(d, full.names = T, pattern = '.csv')
        t_results <- lapply(t, read.csv)
        names(t_results) <- gsub("_results.csv", "", list.files(dir.sp, full.names = F, 
                                                                pattern = '.csv'))
        
        t_results <- lapply(seq_along(t_results), function(i){
          x <-t_results[[i]]
          area <- gsub('random|block|jackk', '', names(t_results[i]))
          cross.validation <- gsub('M1|M2', '', names(t_results[i]))
          x <- cbind(area = area, cross.validation = cross.validation, x)
        })
        tr <- cbind(species = sp, case = cs, do.call(rbind, t_results))
        sp.table <- rbind(sp.table, tr)
      }
      case.table <- rbind(case.table, sp.table)
    }
    fit_table <- rbind(fit_table, case.table)
  }
  
  write_csv(fit_table, 'output/models/marmosini_fitting_results.csv')
} else if(file.exists('output/models/marmosini_fitting_results.csv')){
  fit_table <- read.csv('output/models/marmosini_fitting_results.csv')
  # chunk added for new species: M. germana and jansae ----------------------
  if(any(!names(OCCS)  %in% unique(fit_table$species))){
    new_spp <-names(OCCS)[ !names(OCCS)  %in% unique(fit_table$species) ]
    
    fit_table_np <- data.frame()
    
    for(cs in names(cases)){
      case.table <- data.frame()
      for(sp in new_spp){
        dir.sp <- paste0('output/models/', cs, '/fitting/',sp)
        
        sp.table <- data.frame()
        for(d in dir.sp){
          t <- list.files(d, full.names = T, pattern = '.csv')
          t_results <- lapply(t, read.csv)
          names(t_results) <- gsub("_results.csv", "", list.files(dir.sp, full.names = F, 
                                                                  pattern = '.csv'))
          
          t_results <- lapply(seq_along(t_results), function(i){
            x <-t_results[[i]]
            area <- gsub('random|block|jackk', '', names(t_results[i]))
            cross.validation <- gsub('M1|M2', '', names(t_results[i]))
            x <- cbind(area = area, cross.validation = cross.validation, x)
          })
          tr <- cbind(species = sp, case = cs, do.call(rbind, t_results))
          sp.table <- rbind(sp.table, tr)
        }
        case.table <- rbind(case.table, sp.table)
      }
      fit_table_np <- rbind(fit_table_np, case.table)
    }
    
    fit_table <- rbind(fit_table, fit_table_np)
    fit_table <- split(fit_table, fit_table$species)
    fit_table <- fit_table[ names(OCCS) ]
    fit_table <- do.call(rbind, fit_table)
    
    stopifnot(identical(unique(fit_table$species), names(OCCS)))
    
    write_csv(fit_table, 'output/models/marmosini_fitting_results.csv')
  }
}

less80 <- c('L', 'LQ', 'LQP')
more80 <- c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')

fit_table %>% 
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
  filter(cross.validation != 'jackk')-> filtered_table #importante para manuscrito

if(!file.exists('output/models/marmosini_fitting_results_filtered.csv')){ 
  write_csv(filtered_table, 'output/models/marmosini_fitting_results_filtered.csv')
} else if(file.exists('output/models/marmosini_fitting_results_filtered.csv')){
  filtered_wr <- read.csv('output/models/marmosini_fitting_results_filtered.csv')
  
  if(any(!unique(filtered_table$species)  %in% unique(filtered_wr$species))){
    write_csv(filtered_table, 'output/models/marmosini_fitting_results_filtered.csv')
  }
}

source('R/4-1_make_supp_mat_plots.R')

filtered_table %>% 
  group_by(cross.validation) %>% 
  mutate(trainAUC = mean(train.AUC)) %>% 
  mutate(testAUC = mean(avg.test.AUC)) %>% 
  mutate(Akaike = mean(AICc)) %>% 
  mutate(Omission = mean(avg.test.orMTP)) %>% 
  ggplot(aes(x=cross.validation, y = AICc, fill=case)) + 
  geom_boxplot() 

df.grp <- filtered_table %>% 
  group_by(species, area, case) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>% 
  filter(avg.diff.AUC == min(avg.diff.AUC)) %>% 
  filter(avg.test.orMTP == min(avg.test.orMTP)) %>% 
  filter(AICc == min(AICc)) %>% 
  arrange(species, area, case, cross.validation) %>% 
  ungroup()

df.grp <- df.grp %>% 
  mutate(method = str_c(area, cross.validation)) %>% 
  mutate(Species = str_replace(species, 'Marmosa_|Monodelphis_', 'M. '))

df.grp %>% 
  ggplot(aes(x=Species, y = avg.diff.AUC, fill=case, 
             alpha=avg.test.AUC)) +
  scale_alpha(range = c(0.2,0.9)) +
  geom_bar(stat = 'identity', position = 'dodge', colour='gray60', lwd=0.2,  width = 0.7) + 
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  labs(title = 'Average difference AUC')+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic'), 
    strip.text = element_text(size = 6),
  )

df.grp %>% 
  ggplot(aes(x=Species, y = avg.test.orMTP, fill=case, 
             alpha=avg.test.AUC)) +
  scale_alpha(range = c(0.2,0.9)) +
  geom_bar(stat = 'identity', position = 'dodge', colour='gray60', lwd=0.2,  width = 0.7) + 
  facet_wrap(~ area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  labs(title = 'Average test ommision rate at MTP')+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic'), 
    strip.text = element_text(size = 6),
  )

df.grp %>% 
  ggplot(aes(x=Species, y = train.AUC, fill=case, 
             alpha=avg.test.AUC)) +
  scale_alpha(range = c(0.2,0.9)) +
  geom_bar(stat = 'identity', position = 'dodge', colour='gray60', lwd=0.2, width = 0.7) + 
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  labs(title = 'Training AUC')+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic'), 
    strip.text = element_text(size = 6),
  )

df.grp %>% 
  ggplot(aes(x=Species, y = AICc, fill=case, 
             alpha=avg.test.AUC)) +
  scale_alpha(range = c(0.4,0.9)) +
  geom_bar(stat = 'identity', position = 'dodge', colour='gray60', lwd=0.2, 
           width = 0.7) + 
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  labs(title = 'Akaike Information Criteria')+
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, face = 'italic'), 
    strip.text = element_text(size = 6),
  )

df.final <- filtered_table %>% 
  group_by(species, area) %>% 
  filter(avg.test.AUC == max(avg.test.AUC)) %>% 
  filter(avg.diff.AUC == min(avg.diff.AUC)) %>% 
  filter(avg.test.orMTP == min(avg.test.orMTP)) %>% 
  filter(AICc == min(AICc)) %>% 
  filter(rm == min(rm)) %>% 
  filter(train.AUC == max(train.AUC)) %>% 
  filter(parameters == min(parameters)) %>% 
  filter(avg.test.or10pct == min(avg.test.or10pct)) 

df.final <- split(df.final, df.final$species)

dff <- lapply(df.final, function(x) {
  if(length(x$case) > 2){
    
    xm1 <- subset(x, area == 'M1')
    xm2 <- subset(x, area == 'M2')
    
    if(NROW(xm1) > 1 & any(xm1$case %in% c('ud.all', 'ud.noplants'))){
      xm1 <- subset(xm1, case %in% c('ud.all', 'ud.noplants'))
      if(NROW(xm1) > 1) xm1 <- xm1[sample(1:NROW(xm1), 1),]
    } else if (NROW(xm1) > 1) {
      xm1 <- xm1[sample(1:NROW(xm1), 1),]
    }
    
    if(NROW(xm2) > 1 & any(xm2$case %in% c('ud.all', 'ud.noplants'))){
      xm2 <- subset(xm2, case %in% c('ud.all', 'ud.noplants'))
      if(NROW(xm2) > 1) xm2 <- xm2[sample(1:NROW(xm2), 1),]
    } else if (NROW(xm2) > 1) {
      xm2 <- xm2[sample(1:NROW(xm2), 1),]
    }
    
    xfin <- rbind(xm1, xm2)
    return(xfin)
  } else {
    return(x)
  }
})  

dff <- do.call(rbind, dff)

dir.create('output/models/final_models')

write_csv(dff, 'output/models/final_models/Choosen_models_marmosini_m1&m2.csv')

source('R/4a_explore_and _choose_suboptimal.R')
