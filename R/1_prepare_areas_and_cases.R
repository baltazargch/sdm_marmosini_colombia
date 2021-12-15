#### MAKE MODEL AREAS ####
library(sf)
library(parallel)
library(tidyverse)

# source user-defined functions
walk(list.files('R/UDF/', '.R$', full.names = T), source)
source('R/0_prepare_data.R')

#Load and name species records
OCCS <- read.occs()

#This is a conditional approach. If data exists, it does not calculate it again
#but instead load it to the environment. If new species are added, erase the file
#"AreaM/M*bgmask.RData" and run the code. Note that M2 uses parallel computing. 
#Adjust cores according to your computer capacity.

# Turn off spherical geometries
sf_use_s2(F)

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

env_rast <- rast('/mnt/2TB/GIS/Rasters/Clima/South-Central Ame Climate/bio_1.grd')
dir.create('records/XYMs', showWarnings = FALSE)

area.types <- c('M1', 'M2')
mmask <- list(M1bgmask = M1bgmask, M2bgmask = M2bgmask)

for(type in area.types){
  # type='M1'
  fl <- paste0("records/XYMs/", type, "bgxy.RDS")
  case.mask <- switch(type, M1 = 'M1bgmask', M2 = 'M2bgmask')
  case.mask <- mmask[[ case.mask ]]
  
  if(file.exists(fl)){
    
    # chunk added for new species
    mxy <- readRDS(file=fl)
    
    if(any(!names(OCCS)  %in% names(mxy))){
      sp_new <- names(OCCS)[ which(!names(OCCS)  %in% names(mxy))]
      sp_new <- OCCS[ sp_new ]
      
      masks_sp_new <- case.mask[ names(sp_new) ]
      
      new_bgxy  <- lapply(masks_sp_new, function(x){
        p <- mask(crop(env_rast, x), vect(x))
        p <- terra::spatSample(p, 10000, na.rm=T, xy=T, values=F, as.df=T)
        p
      })
      
      mxy[ names(sp_new) ] <- new_bgxy
      mxy <- mxy[ names(OCCS) ]
      
      #error control
      stopifnot(identical(names(mxy), names(OCCS)))
      saveRDS(mxy, file=fl)
    } else {
      
      case.mask <- switch(type, M1 = 'M1bgxy', M2 = 'M2bgxy')
      
      assign(case.mask, readRDS(fl))
      rm('mxy')
    }
  } else {
    mxy  <- lapply(case.mask, function(x){
      p <- mask(crop(env_rast, x), vect(x))
      p <- terra::spatSample(p, 10000, na.rm=T, xy=T, values=F, as.df=T)
      p
    })
    saveRDS(mxy, file=fl)
    
    case.mask <- switch(type, M1 = 'M1bgxy', M2 = 'M2bgxy')
    
    assign(case.mask, readRDS(fl))
    rm('mxy')
  }
}

#### MAKE MODEL CASES ####

#load and name predictors layers
dir.rast ='/mnt/2TB/GIS/Rasters/Clima/South-Central Ame Climate/'

#choose those predictors known to be important for marsupial distribution
envPredics <- read.envs(dir.rast)
vars.chosen <- c(paste0('bio_', c(2,4,6,10,11,15,16,17)), 'topoWet', 'tri', 'msavi')

envPredics <- envPredics[[ vars.chosen ]]

#error control
stopifnot(all(names(envPredics) == vars.chosen))

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

rm(list=c('dt_rast', 'env_rast', 'mmask', 'X', 'area.types', 'case.mask', 'fl', 
             'i', 'type'))
