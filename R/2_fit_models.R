#### MODEL FITTING ####
#assign ram for java usage (for MaxEnt)
options(java.parameters = "-Xmx4g" )
library(terra)
library(ENMeval)
library(tidyverse)

#Source previous step to get all data in environment
source('R/1_prepare_areas_and_cases.R')

#Define general flags for modeling
algorithm   <- 'maxent.jar' # 'maxnet' 'maxent.jar'
RMlist      <- seq(0.5, 5, 0.5)

#Re-load the predictors layer to reduced problems of memory 
#representation of the data (note that some objects are from previous script)
envPredics <- read.envs(dir.rast)
envPredics <- envPredics[[ vars.chosen ]]

#A simple progress bar for each species process. It doesn't work when parallel
#computing is used.
pb <- txtProgressBar(min = 0, max = length(OCCS), style = 2)
list.mask <- list(M1bgmask = M1bgmask, M2bgmask = M2bgmask)
list.bg <- list(M1bgxy = M1bgxy, M2bgxy = M2bgxy)

for (i in seq_along(OCCS)){
  # i = 1
  #Conditional features based on number of records
  if(NROW(OCCS[[i]]) > 80)
    FClist <- c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')
  else 
    FClist <- c('L', 'LQ', 'LQP')

  
  saveDirs <- paste0('output/fitting/', c('eval.tables', 'other'))
  walk(saveDirs, ~ if(!dir.exists(.x)){ dir.create(.x, recursive = T) })
  #This is a resource- and time-consuming computation (relative to the computer
  #power available). This tries to reduced by not calculating species already
  #assessed. This way, the process can be divided in different sessions. 
  sp.files <- paste0(saveDirs[1], '/', names(OCCS[i]), '.csv')
  if(file.exists(sp.files)){
    sp.files <- read.csv(sp.files)
    
    if(nrow(sp.files) == (length(FClist) * length(RMlist)) * 2 * 2 *  4) {
      cat('All models for', names(OCCS[i]), 'are saved. Skipping species \n')
      next
    }
  }
  
  cat(paste0("Runing model for ", names(OCCS[i]), " at:  ", Sys.time(), '\n'))
  
  # predictors mask and crop -------------------------------------------------
  for(m in c('M1', 'M2')){
    # m = 'M1'
    bg <- list.bg[[ switch(m, M1 = 'M1bgxy', M2 = 'M2bgxy') ]][[ names(OCCS[i]) ]]
    occ <- OCCS[[i]][,2:3]
    colnames(occ) <- colnames(bg) <- c('x', 'y')
    
    bg <- cbind(bg, terra::extract(envPredics, bg, cell=T))
    bg <- bg[ !duplicated(bg$cell), ] %>% dplyr::select(!cell) %>% na.omit()
    
    occ <- cbind(occ, terra::extract(envPredics, as.matrix(occ), cell=T))
    occ <- occ[ !duplicated(occ$cell), ] %>%  dplyr::select(!cell) %>% 
      na.omit()
    
    m.mask  <- st_sf(
      data.frame(id=1,
                 geometry= list.mask[[ switch(m, M1 = 'M1bgmask', 
                                              M2 = 'M2bgmask') ]] [[i]]))
    
    area <- as.numeric(sum(st_area(m.mask))) * 1e-06
    
    nCores <- ifelse(area < 5e05, 8,
                     ifelse(between(area, 5e05, 1e06), 7,  
                            ifelse(between(area, 1e06 + 1, 5e06), 6, 
                                   ifelse(area < 10e06, 4, 2))))
    
    setParallel <- ifelse(nCores > 1, T, F)
    
    # spM <- terra::mask(terra::crop(envPredics, m.mask), vect(m.mask))
    # spM <- terra::app(spM, fun = \(x) {if(sum(is.na(x)) > 0) x * NA else x})
    # tmpFiles(current = F, orphan = T, old = T, remove = T)
    
    for(j in seq_along(cases)){
      # j=1
      cat('Analysing model:   ', names(cases[j]), '\n')
      
      if(j != 4){
        vars.case <- cases[[j]]
      } else if(j==4){
        vars.case <- cases[[j]][[ tolower(m) ]][[ names(OCCS[i]) ]]
      }
      
      bg.case <- bg[ , c('x', 'y', vars.case) ]
      occ.case <- occ[, c('x', 'y',vars.case) ]
      
      cross.validations <- c('block', ifelse(NROW(OCCS[[i]]) < 25, 'jackknife', 'randomkfold'))
      list.dts <- list.vars <- list()
      
      out.names <- paste0(saveDirs, '/', names(OCCS[i]), '.csv')
      
      for(cv in cross.validations){
        # cv = cross.validations[1]
        
        model <- myenm_eval(occs = occ.case, 
                            bg = bg.case,
                            tune.args	= list(fc = FClist, rm = RMlist),
                            partitions = cv, 
                            algorithm = algorithm,
                            numCores = nCores, 
                            doClamp = FALSE, 
                            parallel = setParallel, 
                            quiet=F,
                            calcular.naMismatch = F)
        
        list.dts[[ cv ]] <- cbind(
          species= names(OCCS[i]),
          aream = m, 
          case = names(cases[j]),
          cv = cv, 
          model@results)
        
        list.vars[[ cv ]] <- cbind(
          settings = map(names(model@variable.importance), ~rep(.x, length(vars.case))) %>% 
            do.call(c, .),
          aream = m,
          case = names(cases[j]),
          cv = cv,
          do.call(rbind, model@variable.importance))
      }
      
      do.call(rbind, list.dts) %>% 
        write_csv(., out.names[1], append = ifelse(m == 'M1' && j == 1, FALSE, TRUE))
      do.call(rbind, list.vars) %>% 
        write_csv(., out.names[2], append = T,  col_names = ifelse(j == 1, T, F))
      
    }
  }
  # rm(list = c('spM', 'Predics', 'model', 'm.mask'))
  setTxtProgressBar(pb, i)
}
