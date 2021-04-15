library(terra)
library(stringr)
library(tidyverse)
library(rlang)

envfiles <- list.files('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/GIS_datos/Rasters/MODIS/Terra/Vegetation_Index/', 
                       all.files = T, full.names = T)
names <- list.files('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/GIS_datos/Rasters/MODIS/Terra/Vegetation_Index/', 
                    all.files = T, full.names = F)
reds <- lapply(envfiles[str_detect(envfiles, 'red')][-1], function(x) {
  r <- rast(x)/10000
  r
})
names(reds) <- gsub("MOD13A3.006__1_km_monthly_|_aid0001.tif","", 
                    names[str_detect(names, 'red')][-1])
nirs <- lapply(envfiles[str_detect(envfiles, 'NIR')][-1], function(x) {
  r <- rast(x)/10000
  r
})
names(nirs) <- gsub("MOD13A3.006__1_km_monthly_|_aid0001.tif","", 
                    names[str_detect(names, 'NIR')][-1])
msavi <- function(NIR, R){
  f <- expr({(2 * expr(!!NIR) + 1 - sqrt((2 * expr(!!NIR) + 1)^2 - 8 * (expr(!!NIR) - expr(!!R)))) / 2})
}

formulasMSAVI <- map2(nirs, reds, msavi) 

rastMSAVI <- map(formulasMSAVI, eval)

brickMSAVI <- rast(rastMSAVI)

meanMSAVI <- mean(brickMSAVI)

writeRaster(meanMSAVI, 'LayersBank/msavi.grd', overwrite=T)
