#### CLEAN RECORDS ####
library(terra)
library(dplyr)
library(purrr)

# source user-defined functions
walk(list.files('R/UDF/', '.R$', full.names = T), source)

#read and name record's list
OCCS <- read.occs(rm.dups = F)

#load one predictor to match records and predictors resolution (~1 km)
env_rast <- rast('/mnt/2TB/GIS/Rasters/Clima/South-Central Ame Climate/bio_1.grd')

#loop for each species' records: extract cell number and flag duplicated
for(i in names(OCCS)){
  # i = 'Monodelphis_brevicaudata'
  X <- OCCS[[i]]
  
  dt_rast <- terra::extract(env_rast, X[c(2, 3)], cells=T, df=T)
  dt_rast$rep <- duplicated(dt_rast$cell)
  
  OCCS[[i]]$rep <- dt_rast$rep 
}

#write the flagged data.frames into the original records
walk2(OCCS, list.files('records', pattern = '.csv$', full.names = T), 
             readr::write_csv)
