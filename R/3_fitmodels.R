options(java.parameters = "-Xmx4g" )
source('R/2_get_modelCases.R')
library(dismo)
library(raster)
library(ENMeval)
library(dplyr)

algorithm   <- 'maxent.jar' # 'maxnet' 'maxent.jar'
rasterPreds <- TRUE
RMlists     <- seq(0.5, 5, 0.5)
rm(envPredics)

envsFiles <- list.files('LayersBank/', pattern = '.grd$', full.names = T)
namesEnv <- gsub('.grd', '', list.files('LayersBank/', pattern = '.grd$', full.names = F))

names(envsFiles) <- namesEnv
preds <- c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi')

pb <- txtProgressBar(min = 0, max = length(OCCS), style = 2)

for (i in seq_along(OCCS)){
  
  if(NROW(OCCS[[i]]) > 80){
    FClist <- c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')
  } else {
    FClist <- c('L', 'LQ', 'LQP')
  }
  
  all.files <- list.files('output/models/', recursive = T, pattern = '.csv$')
  all.files <- all.files[ !str_detect(all.files, 'final_models|final_subopt_models') ]
  sp.files <- all.files[str_detect(all.files, names(OCCS[i]))]
  if(length(sp.files) == 16) {
    cat('All models for', names(OCCS[i]), 'are saved. Skipping species \n')
    next
  }
  cat(paste0("Runing model for ", names(OCCS[i]), " at:  ", Sys.time(), '\n'))
  
  # predictors mask and crop --------------------------------------------------------------------------------------------------
  m1mask <- data.frame(pol='pol')
  st_geometry(m1mask) <- M1bgmask[[i]]
  m2mask <- data.frame(pol='pol')
  st_geometry(m2mask) <- M2bgmask[[i]]
  
  area <- as.numeric(sum(st_area(m2mask))) * 1e-06
  # print(area)
  nCores <- ifelse(area < 5e05, 6,
                   ifelse(between(area, 5e05, 1e06), 5,  
                          ifelse(between(area, 1e06 + 1, 5e06), 4, 
                                 ifelse(area < 10e06, 3, 1))))
  
  setParallel <- ifelse(nCores > 1, T, F)
  
  spM1 <- raster::mask(raster::crop(stack(rast(envsFiles[preds])), m1mask), m1mask)
  spM2 <- raster::mask(raster::crop(stack(rast(envsFiles[preds])), m2mask), m2mask)
  # plot(spM2$bio_2)
  for(j in seq_along(cases)){
    cat('Analysing model:   ', names(cases[j]), '\n')
    if(j != 4){
      m1Predics <- stack(spM1[[ cases[[j]] ]])
      m2Predics <- stack(spM2[[ cases[[j]] ]])
    } else if(j==4){
      m1Predics <- stack(spM1[[ cases[[j]]$m1[[ names(OCCS[i]) ]] ]])
      m2Predics <- stack(spM2[[ cases[[j]]$m1[[ names(OCCS[i]) ]] ]])
    }
    
    if (!dir.exists(paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i])))){
      dir.create(paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i])), 
                 recursive = T)
    }
    saveDir <- paste0('output/models/', names(cases[j]), '/fitting/', names(OCCS[i]))
    
    # block for all -------------------------------------------------------------------------------------------------------------
    if(!file.exists(paste0(saveDir, "/M1block_results.csv"))){
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
      write.csv(dt_results, paste0(saveDir, '/M1block_results.csv'), row.names = F)
      rm(dt_results)
    } else {
      cat(paste0(saveDir, "/M1block_results.csv already exists. Skipping to next model\n"))
    }
    
    if(!file.exists(paste0(saveDir, "/M2block_results.csv"))){
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
      write.csv(dt_results, paste0(saveDir, '/M2block_results.csv'), row.names = F)
      rm(dt_results)
    } else {
      cat(paste0(saveDir, "/M2block_results.csv already exists. Skipping to next model\n"))
    }
    
    # conditional of n records --------------------------------------------------------------------------------------------------
    if (NROW(OCCS[[i]]) < 25){
      cat(names(OCCS[i]), 'has less than 25 records. Aplying models accordingly \n')
      
      # jackks --------------------------------------------------------------------------------------------------------------------
      if(!file.exists(paste0(saveDir, "/M1jackk_results.csv"))){
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
        write.csv(dt_results, paste0(saveDir, '/M1jackk_results.csv'), row.names = F)
        rm(dt_results)
      } else {
        cat(paste0(saveDir, "/M1jackk_results.csv already exists. Skipping to next model\n"))
      }
      
      if(!file.exists(paste0(saveDir, "/M2jackk_results.csv"))){
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
        write.csv(dt_results, paste0(saveDir, '/M2jackk_results.csv'), row.names = F)
        rm(dt_results)
      } else {
        cat(paste0(saveDir, "/M2jackk_results.csv already exists. Skipping to next model\n"))
      }
      
    } else {
      # randoms -------------------------------------------------------------------------------------------------------------------
      cat(names(OCCS[i]), 'has more than 25 records. Aplying models accordingly \n')
      
      if(!file.exists(paste0(saveDir, "/M1random_results.csv"))){
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
        write.csv(dt_results, paste0(saveDir, '/M1random_results.csv'), row.names = F)
        rm(dt_results)
      } else {
        cat(paste0(saveDir, "/M1random_results.csv already exists. Skipping to next model\n"))
      }
      
      if(!file.exists(paste0(saveDir, "/M2random_results.csv"))){
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
        write.csv(dt_results, paste0(saveDir, '/M2random_results.csv'), row.names = F)
        rm(dt_results)
      } else {
        cat(paste0(saveDir, "/M2random_results.csv already exists. Skipping to next model\n"))
      }
    }
    rm(list = c('m1Predics', 'm2Predics'))
  }
  rm(list = c('spM1', 'spM2'))
  setTxtProgressBar(pb, i)
}
