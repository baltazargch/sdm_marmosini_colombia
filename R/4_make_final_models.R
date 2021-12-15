#### MAKE MODELS FOR VISUAL INSPECTION AND EVALUATION ####
#Assign ram limit to java process (MaxEnt)
options(java.parameters = "-Xmx4g" )

library(sf)
library(tidyverse)
library(terra)
library(ENMeval)
library(doParallel)
walk(list.files('R/UDF/', '.R$', full.names = T), source)

dir.rast <- '/mnt/2TB/GIS/Rasters/Clima/South-Central Ame Climate/'

#Load the data
OCCS <- read.occs()
cases <- readRDS('Cases/cases.rds')

envsPredics <- read.envs(dir.rast)
envsPredics <- envsPredics[[ cases$ud.all ]]

flag <- list.files('output/models/', recursive = T) %>% 
  .[ !str_detect(., 'maps/')]
if(length(flag) >= 48 * length(OCCS)){
  message('All species complete')
  stop()
}

envsPredics <- terra::app(envsPredics, fun = \(x) {if(sum(is.na(x)) > 0) x * NA else x})

masks <- lapply(list.files('AreaM/', 'RData', full.names = T), \(x) get(load(x)))

bg.points <- map(list.files('records/XYMs/', pattern = 'RDS', full.names = T), readRDS)
names(bg.points) <- names(masks) <- c('M1', 'M2')

#### OPTIMAL MODELS ####

#Load table of chosen models
table_models <- rbind(
  read.csv('output/fitting/Choosen_models_marmosini_m1&m2.csv'), 
  read.csv('output/fitting/Subopt_models_marmosini_m1&m2.csv')[,-26] 
)


# loop for species -----------------------------------------------------------
options(cores=2) #adjust according to computer hardware
registerDoParallel()
getDoParWorkers()

#Following is the process for final optimal models (and it applies for suboptimal too):
#According to the model chosen for each species, we ran MaxEnt models again to 
#extract the prediction rasters and boyce metrics, not previously obtained by the 
#ENMevaluate function in '2_fit_models.R' script. 

#These following steps create a tables, html, figures and maps fo each species, 
#that were inspected visually one by one to evaluate the models performance.
#Since partitions are applied at random, some metrics varied slightly. Other data where
#the exact same as for the initial models. 

#Importantly, in this step we used the maxSSS threshold to explore the results 
#in a absent/present map for further application. That is, to transform the 
#continuous prediction to discrete maps: range distribution maps. 

#Shut on or off the parallel processing by commenting or uncommenting the first
#line of the code
for(sp in unique(table_models$species)){
  # foreach(sp=unique(table_models$species), 
  #         .packages = c('tidyverse', 'terra', 'raster', 'ENMeval', 'sf')) %dopar% {
  # sp = unique(table_models$species)[1]
  
  flag <- list.files(paste0('output/models/', sp), recursive = T)
  if(length(flag) == 48){
    message(sp, ' está completo. Pasando a siguiente especie')
    next
  }
  
  # extract settings ---------------------------------------------------------
  tbl_sp <- subset(table_models, species == sp)
  
  dir.sp <- paste0('output/models/', sp)
  dir.create(dir.sp, recursive = T)
  
  subdirs <- paste0(dir.sp, '/', c('plots', 'tables', 'rasters', 'models'))
  walk(subdirs, dir.create)
  
  occ.sp <- OCCS[[ sp ]][, 2:3]
  
  # loop for area M ---------------------------------------------------------
  # for(i in 1:NROW(tbl_sp)){
  foreach(i=1:NROW(tbl_sp),
          .packages = c('sf', 'tidyverse', 'terra', 'ENMeval', 'raster')) %dopar% {
            # i=1
            
            bg.sp <- bg.points[[ tbl_sp$aream[i] ]][[ sp ]]
            colnames(bg.sp) <- colnames(occ.sp) <- c('x', 'y')
            
            aream <- masks[[ tbl_sp$aream[i] ]][[ sp ]]
            
            basename <- paste(tbl_sp$aream[i], tbl_sp$case[i], 
                              tbl_sp$fc[i], tbl_sp$rm[i], tbl_sp$cv[i], sep='_',
                              collapse = '')
            if(length(list.files(paste0('output/models/',sp), pattern = basename, recursive = T)) ==12){
              return(NULL)
            } else {
              
              if(tbl_sp$case[i] == 'uncorr'){ 
                vars.case <- cases$uncorr[[ tolower(tbl_sp$aream[i] ) ]][[ sp ]]
              } else { 
                vars.case <- cases[[ tbl_sp$case[i] ]]
              }
              
              preds <- mask(crop(envsPredics[[ vars.case ]], aream), vect(aream))
              
              model <- myenm_eval(
                occs = occ.sp, 
                envs = raster:::stack(preds), 
                bg = bg.sp, 
                tune.args = list(fc=tbl_sp$fc[i], rm = tbl_sp$rm[i]), 
                partitions = tbl_sp$cv[i], algorithm = 'maxent.jar', 
                doClamp = TRUE,
                taxon.name = sp, calcular.naMismatch = FALSE)
              
              
              
              # variable responses and contribution ------------------------------------------------
              png(paste0(subdirs[1], '/', basename, 'variable_response.png'), width = 20, 
                  height = 18,units = 'cm', res = 150)
              dismo::response(eval.models(model)[[1]])
              title(basename)
              dev.off()
              
              png(paste0(subdirs[1], '/', basename, 'variable_contribution.png'), width = 20, 
                  height = 18,units = 'cm', res = 150)
              plot(model@models[[1]])
              dev.off()
              
              # boyce ------------------------------------------------------------------
              boyce <- ENMeval:::boyce.cm(fit = model@predictions[[1]], 
                                          obs = terra::extract(model@predictions[[1]], occ.sp) %>%
                                            na.omit(), PEplot = F)
              r <- c(1:length(boyce$F.ratio))[boyce$F.ratio != c(boyce$F.ratio[-1], FALSE)]
              
              png(paste0(subdirs[1], '/', basename, 'boyce_index.png'), width = 30, 
                  height = 18,units = 'cm', res = 150)
              plot(boyce$HS, boyce$F.ratio, xlab = "Habitat suitability", 
                   ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75, 
                   main= paste('Boyce analysis', 'Spearman cor:', boyce$Spearman.cor))
              points(boyce$HS[r], boyce$F.ratio[r], pch = 19, cex = 0.75)
              dev.off()
              
              
              # Predictions maps -------------------------------------------------------------------
              writeRaster(model@predictions[[1]], paste0(subdirs[3], '/', basename, '.tif'))
              
              ths <- c(
                'Maximum.training.sensitivity.plus.specificity.Cloglog.threshold', 
                'X10.percentile.training.presence.Cloglog.threshold'
              )
              v <- model@models[[1]]@results
              names(v) <- dimnames(v)[[1]]
              
              png(paste0(subdirs[1], '/',basename, 'thresholded_maps.png'), 
                  width = 30, height = 18,
                  units = 'cm', res = 150)
              par(mfrow=c(1,2))
              plot(model@predictions[[1]] >= v[ ths[1] ], main=paste('Prediction threholded \n at maxSSS = ', 
                                                                     round(v[ ths[1]], 3)), legend=F)
              points(occ.sp, pch='+')
              plot(model@predictions[[1]] >= v[ ths[2] ], 
                   main=paste('Prediction threholded \n at X10 percentil training = ', 
                              round(v[ ths[2] ], 3)), legend=F)
              points(occ.sp, pch='+')
              dev.off()
              
              # ROC and TPS plot -------------------------------------------------------------------
              png(paste0(subdirs[1], '/', basename, 'ROC_TPR.png'), width = 30, 
                  height = 18,units = 'cm', res = 150)
              par(mfrow=c(1,2))
              dismo::evaluate(model@models[[1]]@presence, model@models[[1]]@absence, 
                              model@models[[1]]) %>% plot(., 'ROC')
              dismo::evaluate(model@models[[1]]@presence, model@models[[1]]@absence, 
                              model@models[[1]]) %>% plot(., 'TPR')
              dev.off()
              
              # Other plots  -----------------------------------------------------------------------
              png(paste0(subdirs[1], '/', basename, 'miscellaneus.png'), width = 30, 
                  height = 18,units = 'cm', res = 150)
              par(mfrow=c(1,2))
              dismo::evaluate(model@models[[1]]@presence, model@models[[1]]@absence, 
                              model@models[[1]]) %>% boxplot(., notch=F)
              dismo::evaluate(model@models[[1]]@presence, model@models[[1]]@absence, 
                              model@models[[1]]) %>% density(.)
              dev.off()
              
              
              # tables -----------------------------------------------------------------
              results_out <- enframe(model@models[[1]]@results) %>% as.data.frame()
              colnames(results_out) <- c('variable', 'value')
              
              write.csv(results_out, paste0(subdirs[2], '/', basename, '_contri_permu.csv'), 
                        row.names = F)
              
              lambdas_out <- rmaxent::parse_lambdas(model@models[[1]])
              
              saveRDS(lambdas_out,paste0(subdirs[2], '/', basename, '_lambdas.RDS'))
              
              values_out <- data.frame(
                pres_abs = c(rep(1, NROW(model@models[[1]]@presence)), 
                             rep(0, NROW(model@models[[1]]@absence))),
                rbind(model@models[[1]]@presence, model@models[[1]]@absence)
              )
              
              write_csv(values_out, paste0(subdirs[2], '/', basename, '_values_presence_absence.csv'))
              write_csv(model@results, paste0(subdirs[2], '/', basename, '_fit_metrics.csv'))
              saveRDS(model@models[[1]], paste0(subdirs[4], '/', basename, '_model.RDS'))
              rm(model)
            }
          }
  tmpFiles(FALSE, TRUE, TRUE, TRUE)
  gc()
}
