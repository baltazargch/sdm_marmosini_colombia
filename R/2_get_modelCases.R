library(terra)
library(parallel)
source('R/UDF/terra_Remove_Cor.R')

source('R/1_get_write_areas_m.R')

if (file.exists(paste0(DATA_WD, "/XYMs/M1bgxy.RData"))){
  # chunk added for new species: M. germana and jansae ----------------------
  load(file=paste0(DATA_WD, "/XYMs/", "M1bgxy.RData"))
  
  if(any(!names(OCCS)  %in% names(M1bgxy))){
    sp_new <- names(OCCS)[ which(!names(OCCS)  %in% names(M1bgxy))]
    sp_new <- OCCS[ sp_new ]
    
    masks_sp_new <- M1bgmask[ names(sp_new) ]
    
    new_bgxy  <- lapply(masks_sp_new, function(x){
      p <- st_sample(x, 10000)
      p <- as.data.frame(st_coordinates(p)); colnames(p) <- c('x', 'y')
      p
    })
    
    M1bgxy[ names(sp_new) ] <- new_bgxy
    M1bgxy <- M1bgxy[ names(OCCS) ]
    
    #error control
    stopifnot(identical(names(M1bgxy), names(OCCS)))
    save(M1bgxy, file=paste0(DATA_WD, "/XYMs/M1bgxy.RData"))
  }
  
} else {
  
  M1bgxy  <- lapply(M1bgmask, function(x){
    p <- st_sample(x, 10000)
    p <- as.data.frame(st_coordinates(p)); colnames(p) <- c('x', 'y')
    p
  })
  save(M1bgxy, file=paste0(DATA_WD, "/XYMs/M1bgxy.RData"))
}

if (file.exists(paste0(DATA_WD, "/XYMs/M2bgxy.RData"))){
  # chunk added for new species: M. germana and jansae ----------------------
  load(file=paste0(DATA_WD, "/XYMs/", "M2bgxy.RData"))
  
  if(any(!names(OCCS)  %in% names(M2bgxy))){
    sp_new <- names(OCCS)[ which(!names(OCCS)  %in% names(M2bgxy))]
    sp_new <- OCCS[ sp_new ]
    
    masks_sp_new <- M2bgmask[ names(sp_new) ]
    
    new_bgxy  <- lapply(masks_sp_new, function(x){
      p <- st_sample(x, 10000)
      p <- as.data.frame(st_coordinates(p)); colnames(p) <- c('x', 'y')
      p
    })
    
    M2bgxy[ names(sp_new) ] <- new_bgxy
    M2bgxy <- M2bgxy[ names(OCCS) ]
    
    #error control
    stopifnot(identical(names(M2bgxy), names(OCCS)))
    save(M2bgxy, file=paste0(DATA_WD, "/XYMs/M2bgxy.RData"))
  }
  
} else {
  M2bgxy  <- lapply(M2bgmask, function(x){
    p <- st_sample(x, 10000)
    p <- as.data.frame(st_coordinates(p)); colnames(p) <- c('x', 'y')
    p
  })
  save(M2bgxy, file=paste0(DATA_WD, "/XYMs/M2bgxy.RData"))
}

envsFiles <- list.files('LayersBank/', pattern = '.grd$', full.names = T)
namesEnv <- gsub('.grd', '', list.files('LayersBank/', pattern = '.grd$', full.names = F))
envPredics <- rast(envsFiles)
names(envPredics) <- namesEnv

envPredics <- envPredics[c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi')]

if(file.exists('Cases/cases.rds')){
  
  cases <- readRDS('Cases/cases.rds')

# chunk added for new species: M. germana and jansae ----------------------
  if(any(!names(OCCS)  %in% names(cases$uncorr$m1))){
    new_spp <- names(OCCS)[!names(OCCS)  %in% names(cases$uncorr$m1)]
    
    m1_cases <- lapply(M1bgmask[ new_spp ], function(x){
      terraRemoveCor(envPredics, x, n.points = 500000)
    })
    
    m2_cases <- lapply(M2bgmask[ new_spp ], function(x){
      terraRemoveCor(envPredics, x, n.points = 500000)
    })
    
    cases$uncorr$m1[ new_spp ] <- m1_cases
    cases$uncorr$m2[ new_spp ] <- m2_cases
    
    cases$uncorr$m1 <- cases$uncorr$m1[ names(OCCS) ]
    cases$uncorr$m2 <- cases$uncorr$m2[ names(OCCS) ]
    
    stopifnot(
      identical(names(OCCS), names(cases$uncorr$m1)), 
      identical(names(OCCS), names(cases$uncorr$m2))
    )
    
    saveRDS(cases, 'Cases/cases.rds')
  }
  
  
} else {
  dir.create('Cases')
  cases <- list(onlywc = c(paste0('bio_', c(2,4,6,10,11,15,16,17))), 
                ud.all = c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi'),
                ud.noplants = c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri'),
                uncorr = list(m1 = NA, m2=NA))
  
  cases$uncorr$m1 <- lapply(M1bgmask, function(x){
    terraRemoveCor(envPredics, x, n.points = 500000)
  })
  cases$uncorr$m2 <- lapply(M2bgmask, function(x){
    terraRemoveCor(envPredics, x, n.points = 500000)
  })
  saveRDS(cases, 'Cases/cases.rds')
}
