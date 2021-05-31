#### MAKE MODEL AREAS ####
library(sf)
library(here)
library(parallel)
library(tidyverse)

#Source area M function
source('R/UDF/make_area_m.R')
source('R/0_prepare_data.R')

#Species' records directory
DATA_WD   <- here("records") #path where the records are

#Load and name species records
OCCS    <- lapply(list.files(DATA_WD, pattern = ".csv$", full.names = T), read.csv) 
names(OCCS) <- gsub("\\.csv$", "", list.files(DATA_WD, pattern = ".csv$", full.names = F)) 

#Exclude duplicated records
OCCS <- lapply(OCCS, function(x) { subset(x, rep != TRUE) })

#This is a conditional approach. If data exists, it does not calculated it again
#but instead load it to the environment. If new species are added, erase the file
#"AreaM/M*bgmask.RData" and run the code. Note that M2 uses parallel computing. 
#Adjust cores according to your computer capacity.

#Read data if areas are to be calculated. Otherwise do not read the data.
if(!file.exists("AreaM/M1bgmask.RData") & !file.exists("AreaM/M2bgmask.RData")) {
  world <- read_sf('wmap/ne_10m_admin_0_countries.shp')
  ecoregions <- read_sf('BiomeSHP/wwf_terr_ecos.shp')
  ecoregions <- ecoregions %>% filter(REALM == 'NT')
  ecoregions <- st_buffer(ecoregions, 0) #make sure data is valid
}

dir.create('AreaM', showWarnings = F)
#Calculate M1 area: buffer-derived
if(file.exists("AreaM/M1bgmask.RData")){  
  load(file="AreaM/M1bgmask.RData")
  names(M1bgmask) <- names(OCCS)
} else {
  M1bgmask <- lapply(OCCS, function(x){
    make_area_m(x, 'simple', ecoregion = ecoregions, cutline = world,
                save.area.m.as.shapes = F)})
  save(M1bgmask, file = 'AreaM/M1bgmask.RData')
}

#Calculate M2 area: ecoregion-derived
if(file.exists("AreaM/M2bgmask.RData")){  
  load(file="AreaM/M2bgmask.RData")
  names(M2bgmask) <- names(OCCS)
} else {
  M2bgmask <- mclapply(OCCS, function(x){
    make_area_m(x, 'ecoregion', ecoregion = ecoregions,
                save.area.m.as.shapes = F) 
  }, mc.cores = 8) #adjust number of cores
  save(M2bgmask, file = 'AreaM/M2bgmask.RData')
}

#### BACKGROUND POINTS ####
library(terra)
library(parallel)

#This is and automated process: calculate background points for each species
#and each area M. If new species are added it checks and append the new ones
#in alphabetical order. If data exist and no new species, only load the data
#without calculating it.

dir.create('records/XYMs', showWarnings = FALSE)
#Background points for M1
if (file.exists(paste0(DATA_WD, "/XYMs/M1bgxy.RData"))){
  # chunk added for new species
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

#Background points for M2
if (file.exists(paste0(DATA_WD, "/XYMs/M2bgxy.RData"))){
  # chunk added for new species
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

#### MAKE MODEL CASES ####

#source user-define function to remove correlation among predictors with
#package terra, for reduced multicollinearity case
source('R/UDF/terra_Remove_Cor.R')


#load and name predictors layers
envsFiles <- list.files('LayersBank/', pattern = '.grd$', full.names = T)
namesEnv <- gsub('.grd', '', list.files('LayersBank/', pattern = '.grd$', full.names = F))
envPredics <- rast(envsFiles)
names(envPredics) <- namesEnv

#choose those predictors known to be important for marsupial distribution
envPredics <- envPredics[c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 
                           'topoWet', 'tri', 'msavi')]

#As above, if new species are added, this script automatically add the predictors
#cases to the data. If already all species have data, only read it and it does not 
#calculate again. 
dir.create('Cases', showWarnings = F)

if(file.exists('Cases/cases.rds')){
  cases <- readRDS('Cases/cases.rds')
  
  # chunk added for new species
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
  dir.create('Cases', showWarnings = F)
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
