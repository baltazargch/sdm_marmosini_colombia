#### MAKE MODELS FOR VISUAL INSPECTION AND EVALUATION ####
#Assign ram limit to java process (MaxEnt)
options(java.parameters = "-Xmx8g" )

#Source previous script to include data into the environment
source('R/3_fit_models.R')

#Remove what isn't necessary for this step
rm(list=c('namesEnv', 'NombresOCCS', 'preds', 'rasterPreds', 'RMlists', 'sp.files', 
          'algorithm', 'all.files', 'csv', 'DATA_WD', 'envsFiles', 'FClist', 'filesR',
          'flag.csv', 'i', 'sp.csv'))

#Load, name, and subset the predictors layers
envsFiles <- list.files('LayersBank', pattern = '.grd$', full.names = T)
names(envsFiles) <- gsub('.grd', '', list.files('LayersBank/', pattern = '.grd$', full.names = F))

preds <- c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi')

predictors <- rast(envsFiles[preds])
names(predictors) <- preds

  #### OPTIMAL MODELS ####

#Load table of chosen models
table_models <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv') 

# helper functions ----------------------------------------------------------
getPartitions <- function(x, bg, method){
  if (method == 'random'){
    y <- ENMeval::get.randomkfold(x[,c('longitude', 'latitude')], bg[,c('x','y')], 5)
  } else if (method == 'jackk') {
    y <- ENMeval::get.jackknife(x[,c('longitude', 'latitude')], bg[,c('x','y')])
  } else if(method == 'block'){
    y <- ENMeval::get.block(x[,c('longitude', 'latitude')], bg[,c('x','y')])
  }
  return(y)
}

getFeaturesArgs <- function(x){
  switch (x,
          L = c('-h','-q','-p', 'threshold=false'), 
          LQ = c('-h','-p', 'threshold=false'),
          LQP = c( '-h', 'threshold=false'),
          H = c('-l', '-q','-p', 'threshold=false'),
          LQH = c('-p', 'threshold=false'),
          LQHP = c('threshold=false'), 
          LQHPT = c('threshold=true')
  )
}

# loop for species -----------------------------------------------------------
library(doParallel)
options(cores=6) #adjust according to computer hardware
registerDoParallel()
getDoParWorkers()

#Following is the process for final optimal models (and it applies for suboptimal too):
#According to the model chosen for each species, we ran MaxEnt models again to 
#extract the prediction rasters and boyce metrics, not previusly obtained by the 
#ENMevaluate function in '3_fit_models.R' script. Note that we used ENMeval 1.0
#so Boyce metrics was calculated manually outside the ENMevaluate function. 

#These following steps create a tables, html, figures and maps fo each species, 
#that were inspected visually one by one to evaluate the models performance and fitting.
#since partitions are applied at random, some metrics varied slightly. Other data where
#the exact same as for the initial models. 

#Importantly, in this step we used the maxSSS threshold to explore the results 
#in a absent/present map for further application. That is, to transform the 
#continuous prediction to discrete maps: range distribution maps. 

#Shut on or off the parrallel processing by commenting or uncommenting the first
#line of the code
for (sp in unique(table_models$species)){
# foreach(sp=unique(table_models$species)) %dopar% {
  flag <- list.files(paste0('output/models/final_models/', sp), recursive = T)
  if(length(flag) == 24){
    message(sp, ' está completo. Pasando a siguiente especie')
    next
  }
  tbl_sp <- subset(table_models, species == sp)
  
  m1mask <- data.frame(pol='pol')
  st_geometry(m1mask) <- M1bgmask[[sp]]
  m2mask <- data.frame(pol='pol')
  st_geometry(m2mask) <- M2bgmask[[sp]]
  
  masks <- list(M1 = m1mask, M2=m2mask)
  bg.points <- list(M1 = M1bgxy, M2= M2bgxy)
  
  # extract settings ---------------------------------------------------------
  settings <- list(
    sp = sp,
    case = tbl_sp$case,
    fea = tbl_sp$features,
    rm = tbl_sp$rm,
    area = tbl_sp$area,
    cv = tbl_sp$cross.validation, 
    auc = tbl_sp$train.AUC
  )
  
  dir.sp <- paste0('output/models/final_models/', sp)
  dir.create(dir.sp)
  
  # loop for area M ---------------------------------------------------------
  for(i in 1:NROW(tbl_sp)){
    occ.sp <- OCCS[[ sp ]]
    
    fold <- getPartitions(occ.sp, bg.points[[ settings$area[i] ]][[sp]], settings$cv[i])
    
    nk <- length(unique(fold$occ.grp))
    
    occtest <- occ.sp[ fold$occ.grp == 1 , ]
    occtrain <- occ.sp[ fold$occ.grp != 1 , ]
    
    p <- c(rep(1, NROW(occ.sp)), rep(0, NROW(bg.points[[ settings$area[i] ]][[sp]])))
    
    # get data --------------------------------------------------------------
    if(settings$case[i] != 'uncorr'){
      data <- as.data.frame(
        rbind(extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                      occ.sp[,c('longitude', 'latitude')])[,-1],
              extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                      bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1])
      )
      
      pres.data <- extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                           occ.sp[,c('longitude', 'latitude')])[,-1]
      bg.data <- extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                         bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1]
      
      case.vars <- cases[[ settings$case[i] ]] 
    } else {
      data <- as.data.frame(
        rbind(
          extract(predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]  ]],
                  occ.sp[,c('longitude', 'latitude')])[,-1],
          extract(predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]] ]],
                  bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1])
        )
      
      pres.data <- extract(
        predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]  ]],
                           occ.sp[,c('longitude', 'latitude')])[,-1]
      bg.data <- extract(
        predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]] ]],
                         bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1]
      
      case.vars <- cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]
    }
    
    AUC.TEST <- double()
    AUC.DIFF <- double()
    maxSSS <- double()
    no_omi <- double()
    for(j in 1:nk){
      train.val <- pres.data[fold$occ.grp != j, , drop = FALSE]
      test.val <- pres.data[fold$occ.grp == j, , drop = FALSE]
      bg.val <- bg.data[fold$bg.grp != j, , drop = FALSE]
      x <- as.data.frame(rbind(train.val, bg.val))
      p.fit <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
      mod <- maxent(x=x, p=p.fit, args = c('autofeature=false',
                                           paste0('betamultiplier', '=', settings$rm[i]), 
                                           getFeaturesArgs(settings$fea[i])))
      
      AUC.TEST[j] <- dismo::evaluate(test.val, bg.data, mod)@auc
      AUC.DIFF[j] <- max(0, dismo::evaluate(train.val, bg.data, 
                                            mod)@auc - AUC.TEST[j])
       
      testp <- dismo::predict(mod, data.frame(test.val))
      testa <- dismo::predict(mod, data.frame(bg.data))
       
      e <- evaluate(p=testp, a=testa)
       
      tab <- base::cbind(e@t, e@TPR + e@TNR)
      maxSSS[j] <- (base::subset(tab, tab[, 2] == max(tab[, 2])))[1, 1]
      no_omi[j] <- as.numeric(threshold(e)[3])
    }
    
    pdat <- cbind(p, data)
    
    data <- as.data.frame(na.omit(pdat)[,-1])
    p <- as.data.frame(na.omit(pdat)[,1])
    
    # make model -------------------------------------------------------------
    mod_sp <- maxent(x=data, p=p, removeDuplicates=T, 
                     args=c('autofeature=false',
                            paste0('betamultiplier', '=', settings$rm[i]), 
                            getFeaturesArgs(settings$fea[i])
                     )
    )
    
    r <- dismo::predict(mod_sp, stack(mask(crop(predictors[[ case.vars  ]],
                                                vect(masks[[ settings$area[i] ]])), 
                                           vect(as(masks[[ settings$area[i] ]], 'Spatial')))), 
                        args=c("outputformat=cloglog", 'threads=4'))
    
    # evaluate ---------------------------------------------------------------
    e3 <- evaluate(p=data[p == 1,], 
                   a=data[p == 0,], 
                   mod_sp)

    tab <- base::cbind(e3@t, e3@TPR + e3@TNR)
    SSS <- (base::subset(tab, tab[, 2] == max(tab[, 2])))[1, 1]
    NO <- as.numeric(threshold(e3)[3])
    # boyce ------------------------------------------------------------------
    boyce.fit <- r
    boyce.obs <- extract(r, occ.sp[c('longitude', 'latitude')])
    boyce.analysis <- ecospat::ecospat.boyce(boyce.fit, boyce.obs, PEplot = F)
    
    # write info -------------------------------------------------------------
    subdirs <- paste0(dir.sp, '/', c('plots', 'tables', 'rasters', 'models'))
    walk(subdirs, dir.create)
    
    basename <- paste(settings$area[i], settings$case[i], 
                      settings$fea[i], settings$rm[i], settings$cv[i],sep='_',
                      collapse = '')
    
    writeRaster(r, paste0(subdirs[3], '/', basename, '.tif'))
    
    png(paste0(subdirs[1], '/', 'thresholded_maps_', basename, '.png'), 
        width = 30, height = 18,
        units = 'cm', res = 150)
    par(mfrow=c(1,2))
    plot(r >= SSS, main=paste('Prediction threholded \n at maxSSS = ', 
                                       round(SSS,3)), legend=F)
    points(occ.sp[,2:3], pch='+')
    plot(r > as.numeric(threshold(e3)[3]), main=paste('Prediction threholded \n at no_omission = ', 
                    round(as.numeric(threshold(e3)[3]), 3)), legend=F)
    points(occ.sp[,2:3], pch='+')
    dev.off()
    
    png(paste0(subdirs[1], '/', 'variable_response_', basename, '.png'), width = 20, 
        height = 18,units = 'cm', res = 150)
    response(mod_sp)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'variable_contribution_', basename, '.png'), width = 20, 
        height = 18,units = 'cm', res = 150)
    plot(mod_sp)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'ROC_TPR_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    par(mfrow=c(1,2))
    plot(e3, 'ROC')
    plot(e3, 'TPR')
    dev.off()
    
    png(paste0(subdirs[1], '/', 'miscellaneus_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    par(mfrow=c(1,2))
    boxplot(e3, notch=F)
    density(e3)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'boyce_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    r <- c(1:length(boyce.analysis$F.ratio))[boyce.analysis$F.ratio != 
                                               c(boyce.analysis$F.ratio[-1], FALSE)]
    plot(boyce.analysis$HS, boyce.analysis$F.ratio, xlab = "Habitat suitability", 
         ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75, 
         main= paste('Boyce analysis', 'Spearman cor:', boyce.analysis$Spearman.cor))
    points(boyce.analysis$HS[r], boyce.analysis$F.ratio[r], pch = 19, cex = 0.75)
    dev.off()
    
    # tables -----------------------------------------------------------------
    results_out <- data.frame(variables = dimnames(mod_sp@results)[[1]], 
                              result = unname(mod_sp@results))
    write_csv(results_out, paste0(subdirs[2], '/contri_permu_', basename, '.csv'))
    lambdas_out <- do.call(rbind, as.list(mod_sp@lambdas))
    write.table(lambdas_out, 
                paste0(subdirs[2], '/lambdas_', basename, '.txt'), 
                sep='\t', row.names = F, col.names = F)
    
    values_out <- data.frame(
      pres_abs = c(rep(1, NROW(mod_sp@presence)), rep(0, NROW(mod_sp@absence))),
      rbind(mod_sp@presence, mod_sp@absence)
    )
    
    fit_values <- data.frame(
      fullmod.AUC = e3@auc,
      avg.test.AUC = mean(AUC.TEST), 
      avg.diff.AUC = mean(AUC.DIFF), 
      avg.test.maxSSS = mean(maxSSS), 
      maxSSS = SSS,
      avg.test.no_omission = mean(no_omi), 
      no_omi = NO,
      cv = settings$cv[i], 
      groups = nk
    )
    
    write_csv(fit_values, paste0(subdirs[2], '/fit_metrics_', basename, '.csv'))
    write_csv(values_out, paste0(subdirs[2], '/values_presence_absence_', basename, '.csv'))
    file.copy(mod_sp@html, paste0(subdirs[4], '/model_page_', basename, '.html'))
  }
}

  #### SUBOPTIMAL MODELS ####
#Load suboptimal model table
table_models <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')  

#Shut on or off the parrallel processing by commenting or uncommenting the first
#line of the code

for (sp in unique(table_models$species)){
  # foreach(sp=unique(table_models$species)) %dopar% {
  # sp = 'Marmosa_alstoni'
  flag <- list.files(paste0('output/models/final_subopt_models/', sp), recursive = T)
  if(length(flag) == 24){
    message(sp, ' está completo. Pasando a siguiente especie')
    next
  }
  tbl_sp <- subset(table_models, species == sp)
  
  m1mask <- data.frame(pol='pol')
  st_geometry(m1mask) <- M1bgmask[[sp]]
  m2mask <- data.frame(pol='pol')
  st_geometry(m2mask) <- M2bgmask[[sp]]

  masks <- list(M1 = m1mask, M2=m2mask)
  bg.points <- list(M1 = M1bgxy, M2= M2bgxy)
  
  # extract settings ---------------------------------------------------------
  settings <- list(
    sp = sp,
    case = tbl_sp$case,
    fea = tbl_sp$features,
    rm = tbl_sp$rm,
    area = tbl_sp$area,
    cv = tbl_sp$cross.validation, 
    auc = tbl_sp$train.AUC
  )
  
  dir.sp <- paste0('output/models/final_subopt_models/', sp)
  dir.create(dir.sp)
  
  # loop for area M ----------------------------------------------------------
  library(doParallel)
  options(cores=4)
  registerDoParallel()
  getDoParWorkers()
  
  foreach(i=1:NROW(tbl_sp)) %dopar% {
    # for(i in 1:NROW(tbl_sp)){
    occ.sp <- OCCS[[ sp ]]
    
    fold <- getPartitions(occ.sp, bg.points[[ settings$area[i] ]][[sp]], settings$cv[i])
    
    nk <- length(unique(fold$occ.grp))
    
    occtest <- occ.sp[ fold$occ.grp == 1 , ]
    occtrain <- occ.sp[ fold$occ.grp != 1 , ]
    
    p <- c(rep(1, NROW(occ.sp)), rep(0, NROW(bg.points[[ settings$area[i] ]][[sp]])))
    
    # get data -----------------------------------------------------------------------------------------
    if(settings$case[i] != 'uncorr'){
      
      data <- as.data.frame(
        rbind(extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                      occ.sp[,c('longitude', 'latitude')])[,-1],
              extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                      bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1])
      )
      
      pres.data <- extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                           occ.sp[,c('longitude', 'latitude')])[,-1]
      bg.data <- extract(predictors[[ cases[[ settings$case[i] ]]  ]], 
                         bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1]
      
      case.vars <- cases[[ settings$case[i] ]] 
    } else {
      
      data <- as.data.frame(
        rbind(
          extract(
            predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]  ]],
            occ.sp[,c('longitude', 'latitude')])[,-1],
          extract(
            predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]] ]],
            bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1])
      )
      
      pres.data <- extract(
        predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]  ]],
        occ.sp[,c('longitude', 'latitude')])[,-1]
      bg.data <- extract(
        predictors[[ cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]] ]],
        bg.points[[ settings$area[i] ]][[sp]][,c('x', 'y')])[,-1]
      
      case.vars <- cases[[ settings$case[i] ]][[ tolower(settings$area[i]) ]][[sp]]
    }
    
    AUC.TEST <- double()
    AUC.DIFF <- double()
    maxSSS <- double()
    no_omi <- double()
    for(j in 1:nk){
      train.val <- pres.data[fold$occ.grp != j, , drop = FALSE]
      test.val <- pres.data[fold$occ.grp == j, , drop = FALSE]
      bg.val <- bg.data[fold$bg.grp != j, , drop = FALSE]
      x <- as.data.frame(rbind(train.val, bg.val))
      p.fit <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))

      mod <- maxent(x=x, p=p.fit, args = c('autofeature=false',
                                           paste0('betamultiplier', '=', settings$rm[i]), 
                                           getFeaturesArgs(settings$fea[i])))
      
      AUC.TEST[j] <- dismo::evaluate(test.val, bg.data, mod)@auc
      AUC.DIFF[j] <- max(0, dismo::evaluate(train.val, bg.data, 
                                            mod)@auc - AUC.TEST[j])
      
      testp <- dismo::predict(mod, data.frame(test.val))
      testa <- dismo::predict(mod, data.frame(bg.data))
      
      e <- evaluate(p=testp, a=testa)
      
      tab <- base::cbind(e@t, e@TPR + e@TNR)
      maxSSS[j] <- (base::subset(tab, tab[, 2] == max(tab[, 2])))[1, 1]
      no_omi[j] <- as.numeric(threshold(e)[3])
    }
    
    pdat <- cbind(p, data)
    
    data <- as.data.frame(na.omit(pdat)[,-1])
    p <- as.data.frame(na.omit(pdat)[,1])
    
    # make model -------------------------------------------------------------
    mod_sp <- maxent(x=data, p=p, removeDuplicates=T, 
                     args=c('autofeature=false',
                            paste0('betamultiplier', '=', settings$rm[i]), 
                            getFeaturesArgs(settings$fea[i])
                     )
    )
    
    r <- dismo::predict(mod_sp, stack(mask(crop(predictors[[ case.vars  ]],
                                                vect(as(masks[[ settings$area[i] ]], 'Spatial'))), 
                                           vect(as(masks[[ settings$area[i] ]], 'Spatial')))), 
                        args=c("outputformat=cloglog", 'threads=4'))
    # evaluate ---------------------------------------------------------------
    e3 <- evaluate(p=data[p == 1,], 
                   a=data[p == 0,], 
                   mod_sp)
    
    tab <- base::cbind(e3@t, e3@TPR + e3@TNR)
    SSS <- (base::subset(tab, tab[, 2] == max(tab[, 2])))[1, 1]
    NO <- as.numeric(threshold(e3)[3])
    
    # boyce ------------------------------------------------------------------
    boyce.fit <- r
    boyce.obs <- extract(r, occ.sp[c('longitude', 'latitude')])
    boyce.analysis <- ecospat::ecospat.boyce(boyce.fit, boyce.obs, PEplot = F)
    
    # write info ----------------------------------------------------------------------------------------------------------------
    subdirs <- paste0(dir.sp, '/', c('plots', 'tables', 'rasters', 'models'))
    walk(subdirs, dir.create)
    
    basename <- paste(settings$area[i], settings$case[i], 
                      settings$fea[i], settings$rm[i], settings$cv[i],sep='_',
                      collapse = '')
    
    writeRaster(r, paste0(subdirs[3], '/', basename, '.tif'))
    
    png(paste0(subdirs[1], '/', 'thresholded_maps_', basename, '.png'), 
        width = 30, height = 18,
        units = 'cm', res = 150)
    par(mfrow=c(1,2))
    plot(r >= SSS, main=paste('Prediction threholded \n at maxSSS = ', 
                              round(SSS,3)), legend=F)
    points(occ.sp[,2:3], pch='+')
    plot(r > as.numeric(threshold(e3)[3]), main=paste('Prediction threholded \n at no_omission = ', 
                                                      round(as.numeric(threshold(e3)[3]), 3)), legend=F)
    points(occ.sp[,2:3], pch='+')
    dev.off()
    
    png(paste0(subdirs[1], '/', 'variable_response_', basename, '.png'), width = 20, 
        height = 18,units = 'cm', res = 150)
    response(mod_sp)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'variable_contribution_', basename, '.png'), width = 20, 
        height = 18,units = 'cm', res = 150)
    plot(mod_sp)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'ROC_TPR_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    par(mfrow=c(1,2))
    plot(e3, 'ROC')
    plot(e3, 'TPR')
    dev.off()
    
    png(paste0(subdirs[1], '/', 'miscellaneus_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    par(mfrow=c(1,2))
    boxplot(e3, notch=F)
    density(e3)
    dev.off()
    
    png(paste0(subdirs[1], '/', 'boyce_', basename, '.png'), width = 30, 
        height = 18,units = 'cm', res = 150)
    r <- c(1:length(boyce.analysis$F.ratio))[boyce.analysis$F.ratio != 
                                               c(boyce.analysis$F.ratio[-1], FALSE)]
    plot(boyce.analysis$HS, boyce.analysis$F.ratio, xlab = "Habitat suitability", 
         ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75, 
         main= paste('Boyce analysis', 'Spearman cor:', boyce.analysis$Spearman.cor))
    points(boyce.analysis$HS[r], boyce.analysis$F.ratio[r], pch = 19, cex = 0.75)
    dev.off()
    
    # tables -----------------------------------------------------------------
    results_out <- data.frame(variables = dimnames(mod_sp@results)[[1]], 
                              result = unname(mod_sp@results))
    write_csv(results_out, paste0(subdirs[2], '/contri_permu_', basename, '.csv'))
    lambdas_out <- do.call(rbind, as.list(mod_sp@lambdas))
    write.table(lambdas_out, 
                paste0(subdirs[2], '/lambdas_', basename, '.txt'), 
                sep='\t', row.names = F, col.names = F)
    
    values_out <- data.frame(
      pres_abs = c(rep(1, NROW(mod_sp@presence)), rep(0, NROW(mod_sp@absence))),
      rbind(mod_sp@presence, mod_sp@absence)
    )
    
    fit_values <- data.frame(
      fullmod.AUC = e3@auc,
      avg.test.AUC = mean(AUC.TEST), 
      avg.diff.AUC = mean(AUC.DIFF), 
      avg.test.maxSSS = mean(maxSSS), 
      maxSSS = SSS,
      avg.test.no_omission = mean(no_omi), 
      no_omi = NO,
      cv = settings$cv[i], 
      groups = nk
    )
    
    write_csv(fit_values, paste0(subdirs[2], '/fit_metrics_', basename, '.csv'))
    write_csv(values_out, paste0(subdirs[2], '/values_presence_absence_', basename, '.csv'))
    file.copy(mod_sp@html, paste0(subdirs[4], '/model_page_', basename, '.html'))
  }
}

