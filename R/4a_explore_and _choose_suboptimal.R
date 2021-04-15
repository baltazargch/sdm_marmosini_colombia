library(dplyr)
library(stringr)
library(readr)

#### Load table ####
filtered_table <- read.csv('output/models/marmosini_fitting_results_filtered.csv')
optimal <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')

head(filtered_table)

species <- unique(filtered_table$species)

sp_filtered <- list()
for (sp in species){
  # sp=species[2]
  df.sp <- filtered_table %>% 
    group_by(area, case) %>% 
    filter(species == sp) %>% 
    filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>% 
    filter(avg.diff.AUC == min(avg.diff.AUC)) %>% 
    filter(avg.test.orMTP == min(avg.test.orMTP)) %>% 
    arrange(area, case, cross.validation) %>% 
    ungroup()
  
  opt <- subset(optimal, species == sp & area == 'M1')
  
  xm1 <- subset(df.sp, area == 'M1')
  xm1 <- xm1 %>% filter(settings != opt$settings | case != opt$case)
  xm1 <- subset(xm1, AICc == xm1$AICc[which(order(xm1$AICc) == 2)])
  
  opt <- subset(optimal, species == sp & area == 'M2')
  
  xm2 <- subset(df.sp, area == 'M2')
  xm2 <- xm2 %>% filter(settings != opt$settings | case != opt$case)
  xm2 <- subset(xm2, AICc == xm2$AICc[which(order(xm2$AICc) == 2)])
  
  
  dff <- rbind(xm1, xm2)
  
  sp_filtered[[sp]] <- dff
}

df.final <- do.call(rbind, sp_filtered)

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

dir.create('output/models/final_subopt_models')

write_csv(dff, 'output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')
