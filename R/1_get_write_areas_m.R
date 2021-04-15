library(sf)
library(here)
library(parallel)
library(tidyverse)
source('R/UDF/make_area_m.R')

DATA_WD   <- here("records") #path where the records are

filesR  <- list.files(DATA_WD, pattern = ".csv$", full.names = T) #leer directorio
OCCS    <- lapply(filesR, read.csv) #leer .csv en objeto tipo lista

NombresOCCS <- list.files(DATA_WD, pattern = ".csv$", full.names = F) #para nombrar elementos de la lista
names(OCCS) <- gsub("\\.csv$", "", NombresOCCS) #eliminar nombre de la extensión en el nombre y asignar a lista OCCS
# names(OCCS) #corroborar nombres son cor

OCCS <- lapply(OCCS, function(x) { subset(x, rep != TRUE) })

if(!file.exists("AreaM/M1bgmask.RData") & !file.exists("AreaM/M2bgmask.RData")) {
  world <- read_sf('/mnt/2TB/GIS_datos/Vectores/Mundo/ne_10m_admin_0_countries.shp')
  
  ecoregions <- read_sf('/mnt/2TB/GIS_datos/Vectores/Biologicos/Bioregiones/Ecoregions2017/America_ecoregions2017.shp')
  ecoregions <- st_buffer(ecoregions, 0)
}

if(file.exists("AreaM/M1bgmask.RData")){  
  load(file="AreaM/M1bgmask.RData")
  names(M1bgmask) <- names(OCCS)
  
} else {
  
  M1bgmask <- lapply(OCCS, function(x){
    make_area_m(x, 'simple', ecoregion = ecoregions, cutline = world,
                save.area.m.as.shapes = F) 
  })
  save(M1bgmask, file = 'AreaM/M1bgmask.RData')
}

if(file.exists("AreaM/M2bgmask.RData")){  
  load(file="AreaM/M2bgmask.RData")
  names(M2bgmask) <- names(OCCS)
} else {
  
  M2bgmask <- mclapply(OCCS, function(x){
    make_area_m(x, 'ecoregion', ecoregion = ecoregions,
                save.area.m.as.shapes = F) 
  }, mc.cores = 8)
  save(M2bgmask, file = 'AreaM/M2bgmask.RData')
}

# plot(M2bgmask$Marmosa_germana)
