#### MODEL FILTERING, EVALUATION AND CHOOSNG ####
library(tidyverse)
source('R/3_fit_models.R')

#Access each species model directory and results and make one table. 
#If new species were added, this script add them automatically to the table.
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
  # chunk added for new species
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

#Clean the data: here we exlucded jackknife cross-validation
#since the prediction from this models were considered by the authors as 
#inadequate. 
filtered_table <- fit_table %>% 
  na.omit() %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>% 
  mutate(features = factor(features, 
                           levels = c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')
                           )
         ) %>% 
  filter(cross.validation != 'jackk') 

if(!file.exists('output/models/marmosini_fitting_results_filtered.csv')){ 
  write_csv(filtered_table, 'output/models/marmosini_fitting_results_filtered.csv')
} else if(file.exists('output/models/marmosini_fitting_results_filtered.csv')){
  filtered_wr <- read.csv('output/models/marmosini_fitting_results_filtered.csv')
  
  if(any(!unique(filtered_table$species)  %in% unique(filtered_wr$species))){
    write_csv(filtered_table, 'output/models/marmosini_fitting_results_filtered.csv')
  }
}

#### OPTIMAL MODELS PER SPECIES ####

#We used primarily the avg.test.AUC to filter the results. See details in the
#manuscript associated with this script (pending to be published).
df.final <- filtered_table %>% 
  group_by(species, area) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>% #main criterium
  filter(avg.diff.AUC == min(avg.diff.AUC)) %>% 
  filter(avg.test.orMTP == min(avg.test.orMTP)) %>% 
  filter(AICc == min(AICc)) %>% 
  filter(rm == min(rm)) %>% 
  filter(train.AUC == max(train.AUC)) %>% 
  filter(parameters == min(parameters)) %>% 
  filter(avg.test.or10pct == min(avg.test.or10pct)) 

df.final <- split(df.final, df.final$species)

#For the rare cases where more than one model for each area were chose, 
#assign one at random.
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

#Combine and write results
dff <- do.call(rbind, dff)
dir.create('output/models/final_models')
write_csv(dff, 'output/models/final_models/Choosen_models_marmosini_m1&m2.csv')

#### SUBOPTIMAL MODELS PER SPECIES ####
optimal <- dff

species <- unique(filtered_table$species)

#For each species chose the second best model
sp_filtered <- list()
for (sp in species){
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

#Combine and write results
dff <- do.call(rbind, dff)
dir.create('output/models/final_subopt_models')
write_csv(dff, 'output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')