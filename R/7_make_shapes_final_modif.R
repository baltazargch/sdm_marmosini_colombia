library(sf)
library(dplyr)
library(raster)
library(smoothr)
source('R/UDF/getmaxSSS.R')
source('/mnt/2TB/UDF/polygonizer.R')

dtmodels <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')
dtmodels$model <- 'o'
dtmodels1 <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')
dtmodels1$model <- 's'

CHOS <- rbind(dtmodels, dtmodels1)

OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)

names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

goodmodels <- data.frame(
  species = names(OCCS), 
  area = c('M2', 'M2', 'M1', 'M2', 'M2',
           'M2', 'M2', 'M2', 'M2', 'M2', 
           'M2', 'M2', 'M2', 'M1', 'M2', 
           'M2'), 
  model = c('o', 's', 's', 'o', 'o', 
            'o', 'o', 'o', 'o', 's', 
            'o', 's', 's', 'o', 'o', 
            'o')
)

wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
col <- subset(wmap, NAME == 'Colombia')

col <- col %>%  st_cast('POLYGON') %>% 
  mutate(area = st_area(.)) %>% 
  filter(area == max(area))

for(i in 1:NROW(goodmodels)){
  cc <- CHOS[ CHOS$species == goodmodels$species[i], ]
  
  cc <- cc %>% filter(area == goodmodels$area[i] & model == goodmodels$model[i])
  
  goodmodels[ i , colnames(cc)[ -c(1, 3, 22, 23) ]] <- cc[ ,-c(1, 3, 22, 23) ]
}

dir.create('output/final_shapes', showWarnings = F)

shapes_species <- vector('list', length(OCCS))

library(doParallel)
options(cores=4)
registerDoParallel()
getDoParWorkers()

# foreach(i=1:NROW(goodmodels)) %dopar% {
for(i in 1:NROW(goodmodels)){
  # i = 1
  
  if(file.exists(paste0('output/final_shapes/', goodmodels$species[i], '.shp'))){
    print( paste0('output/final_shapes/', goodmodels$species[i], '.shp', ' already exists'))
    next
  } else {
    basedir <- ifelse(goodmodels$model[i] == 'o', 'output/models/final_models/', 'output/models/final_subopt_models/')
    
    rfile <- list.files(paste0(basedir, goodmodels$species[i], '/rasters'), 
                        pattern = '.tif$', full.names = T) 
    r <- raster(rfile[grep(goodmodels$area[i], rfile)])
    
    maxSSS <- getmaxSSS(goodmodels$species[i], goodmodels$area[i], goodmodels$model[i])
    
    r[r < maxSSS] <- NA
    r[!is.na(r)] <- 1
    
    r1 <- raster::aggregate(r, fact=2, fun=mean)
    # area(r1)
    r_modif <- rasterToPolygons(r1, dissolve = T)
    r_modif <- st_as_sf(r_modif)
    
    r_modif <- st_transform(r_modif, 
                            crs='+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs') 
    
    th <- units::set_units(200, km^2)
    
    r_modif <- fill_holes(r_modif, threshold = th)
    # plot(r_modif)
    r_modif <- st_cast(r_modif, 'POLYGON')
    
    th <- units::set_units(50, km^2)
    r_modif <- drop_crumbs(r_modif, threshold = th)
    
    r_modif <- smooth(r_modif, method = "ksmooth", smoothness = 2)
    r_modif <- st_transform(r_modif, 
                            crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ') 
    
    range_col <- st_intersection(r_modif, col)
    
    # plot(col['geometry'], col='gray90')
    # plot(range_col['geometry'], col='forest green', add=T)
    # points(OCCS[[i]][ , c(2, 3)], pch='+', cex=1.5)
    
    write_sf(range_col, 
             paste0('output/final_shapes/', goodmodels$species[i], '.shp'))
  }
}



