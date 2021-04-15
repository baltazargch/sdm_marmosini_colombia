options(java.parameters = "-Xmx4g" )
source('R/2_get_modelCases.R')
library(dismo)
library(raster)
library(ENMeval)
library(dplyr)

algorithm   <- 'maxent.jar' # 'maxnet' 'maxent.jar'
rasterPreds <- TRUE
RMlists     <- seq(0.5, 5, 0.5)
FClist      <- c("LQP")
rm(envPredics)

envsFiles <- list.files('LayersBank/', pattern = '.grd$', full.names = T)
namesEnv <- gsub('.grd', '', list.files('LayersBank/', pattern = '.grd$', full.names = F))

names(envsFiles) <- namesEnv
preds <- c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi','EVI')

for (i in seq_along(OCCS)){
  
    cat(paste('For ', names(OCCS[i]), ' running extra models LQP\n'))
    
    all.files <- list.files('output/models/', recursive = T, pattern = '.csv$', full.names = T)
    all.files <- all.files[!str_detect(string = all.files, pattern = 'final_models')]
    sp.files <- all.files[str_detect(all.files, names(OCCS[i]))]
    sp.csv <- lapply(sp.files, function(x) read.csv(x))
    
    flag.csv <- vector()
    for(csv in seq_along(sp.csv)){
      flag.csv[csv] <- NROW(sp.csv[[csv]]) > 60
    }
    if(all(flag.csv)){
      cat('All models for', names(OCCS[i]), 'are saved. Skipping species \n')
      next
    }
      
    # if(length(sp.files) == 16) {
    #   cat('All models for', names(OCCS[i]), 'are saved. Skipping species \n')
    #   next
    # }
    cat(paste0("Runing model for ", names(OCCS[i]), " at:  ", Sys.time(), '\n'))
    
    # predictors mask and crop --------------------------------------------------------------------------------------------------
    m1mask <- data.frame(pol='pol')
    st_geometry(m1mask) <- M1bgmask[[i]]
    m2mask <- data.frame(pol='pol')
    st_geometry(m2mask) <- M2bgmask[[i]]
    
    area <- as.numeric(sum(st_area(m2mask))) * 1e-06
    print(area)
    nCores <- ifelse(area < 5e05, 5,
                     ifelse(between(area, 5e05, 1e06), 2,  
                            ifelse(between(area, 1e06 + 1, 5e06), 2, 
                                   ifelse(area<10e06, 2, 1))))
    
    setParallel <- ifelse(nCores > 1, T, F)
    
    spM1 <- terra::mask(crop(rast(envsFiles[preds]), vect(m1mask)), vect(m1mask))
    spM2 <- terra::mask(crop(rast(envsFiles[preds]), vect(m2mask)), vect(m2mask))
    
    for(j in seq_along(cases)){
      cat('Analysing model:   ', names(cases[j]), '\n')
      if(j != 4){
        m1Predics <- stack(spM1[[ cases[[j]] ]])
        m2Predics <- stack(spM2[[ cases[[j]] ]])
      } else if(j==4){
        m1Predics <- stack(spM1[[ cases[[j]]$m1[[names(OCCS[i])]] ]])
        m2Predics <- stack(spM2[[ cases[[j]]$m1[[names(OCCS[i])]] ]])
      }
      
      if (!dir.exists(paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i])))){
        dir.create(paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i])), 
                   recursive = T)
      }
      saveDir <- paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i]))
      
      # block for all -------------------------------------------------------------------------------------------------------------
      test <- read.csv(paste0(saveDir, "/M1block_results.csv"))
      if(nrow(test) == 60){
        dt_results <- lapply(RMlists, function(k) {
          ENMevaluate(OCCS[[i]][,2:3], m1Predics, bg.coords = M1bgxy[[i]], 
                      RMvalues = k,
                      fc = FClist,
                      method = "block", 
                      algorithm = algorithm,
                      numCores = nCores, 
                      parallel = setParallel,
                      bin.output = F,
                      rasterPreds = rasterPreds)@results
        })
        dt_results <- do.call(plyr::rbind.fill, dt_results)
        test <- rbind(test, dt_results)
        write.csv(test, paste0(saveDir, '/M1block_results.csv'), row.names = F)
        rm(list=c('dt_results', 'test'))
      } else {
        cat(paste0(saveDir, "/M1block_results.csv already exists. Skipping to next model\n"))
      }
      
      test <- read.csv(paste0(saveDir, "/M2block_results.csv"))
      if(nrow(test) == 60){
        dt_results <- lapply(RMlists, function(k) {
          ENMevaluate(OCCS[[i]][,2:3], m2Predics, bg.coords = M2bgxy[[i]], 
                      RMvalues = k,
                      fc = FClist,
                      method = "block", 
                      algorithm = algorithm,
                      numCores = nCores, 
                      parallel = setParallel,
                      bin.output = F,
                      rasterPreds = rasterPreds)@results
        })
        dt_results <- do.call(plyr::rbind.fill, dt_results)
        test <- rbind(test, dt_results)
        write.csv(test, paste0(saveDir, '/M2block_results.csv'), row.names = F)
        rm(list=c('dt_results', 'test'))
      } else {
        cat(paste0(saveDir, "/M2block_results.csv already exists. Skipping to next model\n"))
      }
      
      # conditional of n records --------------------------------------------------------------------------------------------------
      if (NROW(OCCS[[i]]) < 25){
        cat(names(OCCS[i]), 'has less than 25 records. Aplying models accordingly \n')
        
        # jackks --------------------------------------------------------------------------------------------------------------------
        test <- read.csv(paste0(saveDir, "/M1jackk_results.csv"))
        if(nrow(test) == 60){
          dt_results <- lapply(RMlists, function(k) {
            ENMevaluate(OCCS[[i]][,2:3], m1Predics, bg.coords = M1bgxy[[i]], 
                        RMvalues = k,
                        fc = FClist,
                        method = "jackknife", 
                        algorithm = algorithm,
                        numCores = nCores, 
                        parallel = setParallel,
                        bin.output = F,
                        rasterPreds = rasterPreds)@results
          })
          dt_results <- do.call(plyr::rbind.fill, dt_results)
          test <- rbind(test, dt_results)
          write.csv(test, paste0(saveDir, '/M1jackk_results.csv'), row.names = F)
          rm(list=c('dt_results', 'test'))
        } else {
          cat(paste0(saveDir, "/M1jackk_results.csv already exists. Skipping to next model\n"))
        }
        
        test <- read.csv(paste0(saveDir, "/M2jackk_results.csv"))
        if(nrow(test) == 60){
          dt_results <- lapply(RMlists, function(k) {
            ENMevaluate(OCCS[[i]][,2:3], m2Predics, bg.coords = M2bgxy[[i]], 
                        RMvalues = k,
                        fc = FClist,
                        method = "jackknife", 
                        algorithm = algorithm,
                        numCores = 4, 
                        parallel = setParallel,
                        bin.output = F,
                        rasterPreds = rasterPreds)@results
          })
          dt_results <- do.call(plyr::rbind.fill, dt_results)
          test <- rbind(test, dt_results)
          write.csv(test, paste0(saveDir, '/M2jackk_results.csv'), row.names = F)
          rm(list=c('dt_results', 'test'))
        } else {
          cat(paste0(saveDir, "/M2jackk_results.csv already exists. Skipping to next model\n"))
        }
        
      } else {
        # randoms -------------------------------------------------------------------------------------------------------------------
        cat(names(OCCS[i]), 'has more than 25 records. Aplying models accordingly \n')
        
        test <- read.csv(paste0(saveDir, "/M1random_results.csv"))
        if(nrow(test) == 60){
          dt_results <- lapply(RMlists, function(k) {
            ENMevaluate(OCCS[[i]][,2:3], m1Predics, bg.coords = M1bgxy[[i]], 
                        RMvalues = k,
                        fc = FClist,
                        method = "randomkfold", 
                        kfolds = 5, 
                        algorithm = algorithm,
                        numCores = nCores, 
                        parallel = setParallel,
                        bin.output = F,
                        rasterPreds = rasterPreds)@results
          })
          dt_results <- do.call(plyr::rbind.fill, dt_results)
          test <- rbind(test, dt_results)
          write.csv(test, paste0(saveDir, '/M1random_results.csv'), row.names = F)
          rm(list=c('dt_results', 'test'))
        } else {
          cat(paste0(saveDir, "/M1random_results.csv already exists. Skipping to next model\n"))
        }
        
        test <- read.csv(paste0(saveDir, "/M2random_results.csv"))
        if(nrow(test) == 60){
          dt_results <- lapply(RMlists, function(k) {
            ENMevaluate(OCCS[[i]][,2:3], m2Predics, bg.coords = M2bgxy[[i]], 
                        RMvalues = k,
                        fc = FClist,
                        method = "randomkfold", 
                        kfolds = 5, 
                        algorithm = algorithm,
                        numCores = nCores, 
                        parallel = setParallel,
                        bin.output = F,
                        rasterPreds = rasterPreds)@results
          })
          dt_results <- do.call(plyr::rbind.fill, dt_results)
          test <- rbind(test, dt_results)
          write.csv(test, paste0(saveDir, '/M2random_results.csv'), row.names = F)
          rm(list=c('dt_results', 'test'))
        } else {
          cat(paste0(saveDir, "/M2random_results.csv already exists. Skipping to next model\n"))
        }
      }
      rm(list = c('m1Predics', 'm2Predics'))
    }
    rm(list = c('spM1', 'spM2'))
  # }
}