#### CLEAN RECORDS ####
library(raster)
library(dplyr)

#read and name record's list
OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)
names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

#load one predictor to match records and predictors resolution (~1 km)
env_rast <- raster('LayersBank/bio_1.grd')

#loop for each species' records: extract cell number and flag duplicated
for(i in names(OCCS)){
  # i = 'Monodelphis_brevicaudata'
  X <- OCCS[[i]]
  
  dt_rast <- raster::extract(env_rast, X[c(2, 3)], cellnumbers=T, df=T)
  dt_rast$rep <- duplicated(dt_rast$cells)
  
  OCCS[[i]]$rep <- dt_rast$rep 
}

#write the flagged data.frames into the original records
purrr::walk2(OCCS, list.files('records', pattern = '.csv$', full.names = T), 
             readr::write_csv)

#### MAKE MSAVI PREDICTOR LAYER ####
library(terra)
library(stringr)
library(tidyverse)
library(rlang)

#load and name MODIS Terra data for the year 2000
envfiles <- list.files('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/GIS_datos/Rasters/MODIS/Terra/Vegetation_Index/', 
                       all.files = T, full.names = T)
names <- list.files('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/GIS_datos/Rasters/MODIS/Terra/Vegetation_Index/', 
                    all.files = T, full.names = F)

#extract and name red band for every raster
reds <- lapply(envfiles[str_detect(envfiles, 'red')][-1], function(x) {
  r <- rast(x)/10000
  r
})
names(reds) <- gsub("MOD13A3.006__1_km_monthly_|_aid0001.tif","", 
                    names[str_detect(names, 'red')][-1])

#extract and name nir band for every raster
nirs <- lapply(envfiles[str_detect(envfiles, 'NIR')][-1], function(x) {
  r <- rast(x)/10000
  r
})
names(nirs) <- gsub("MOD13A3.006__1_km_monthly_|_aid0001.tif","", 
                    names[str_detect(names, 'NIR')][-1])

#Qi et al. (1994) function to calculate MSAVI
msavi <- function(NIR, R){
  f <- expr({(2 * expr(!!NIR) + 1 - sqrt((2 * expr(!!NIR) + 1)^2 - 8 * (expr(!!NIR) - expr(!!R)))) / 2})
}

#create formula's promises for each month of 2000
formulasMSAVI <- map2(nirs, reds, msavi) 

#evaluate each formula
rastMSAVI <- map(formulasMSAVI, eval)

#combine in one stack
brickMSAVI <- rast(rastMSAVI)

#create and write mean value of MSAVI for the year 2000
meanMSAVI <- mean(brickMSAVI)
writeRaster(meanMSAVI, 'LayersBank/msavi.grd', overwrite=T)
