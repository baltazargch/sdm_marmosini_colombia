# After visual inspection of optimal and sub optimal models, this script chose the 
# final selected models and create a shapefile according to the data.frame 'goodmodels'.
# Final shapes where then modify according to geographic barriers that can be consulted 
# at the manuscript supporting information (pending acceptance).
library(sf)
library(tidyverse)
library(terra)
library(smoothr)
sf_use_s2(F) #turn off spherical geometries

walk(list.files('R/UDF/', '.R$', full.names = T), source)

OCCS <- read.occs()

#This were chosen manually. So the table should be modified according to species number
# and selection of area M and optimal or suboptimal models. 
final_mods <- read.csv('output/final_chosen_models.csv')

all.files <- list.files('output/models/', recursive = T, full.names = T)

wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
col <- subset(wmap, NAME == 'Colombia')

col <- col %>%  st_cast('POLYGON') %>% 
  mutate(area = st_area(.)) %>% 
  filter(area == max(area))

dir.create('output/ranges/original', showWarnings = F, recursive = T)

shapes_species <- vector('list', length(OCCS))

library(doParallel)
options(cores=10) #adjust according to computer hardware
registerDoParallel()
getDoParWorkers()

#Shut on or off parallel processing by commenting or uncommenting the next lines
foreach(i=seq_along(OCCS)) %dopar% {
# for(i in seq_along(OCCS)){
  # i = 1
  sp <- names(OCCS[i])
  save.dir <- paste0('output/ranges/original/', sp, '.gpkg')
  if(file.exists(save.dir)){
    print(paste('Final range for', sp, 'already exists'))
    # next
  } else {
    sp.data <- filter(final_mods, Species == sp)
    rast.sp <- all.files[str_detect(all.files, paste0(sp, '/rasters'))] %>% 
      .[str_detect(., sp.data$final.model)]
    
    sp.prediction <- rast(rast.sp)
    
    maxSSS <- all.files[str_detect(all.files, paste0(sp, '/tables'))] %>% 
      .[str_detect(., 'contri_permu')] %>% 
      .[str_detect(., sp.data$final.model)] %>% 
      read.csv() %>% 
      filter(variable == 'Maximum.training.sensitivity.plus.specificity.Cloglog.threshold') %>% 
      select(value) %>% as.numeric()
    
    r1 <- ifel(sp.prediction < maxSSS, NA, sp.prediction)
    r1 <- ifel(!is.na(r1), 1, r1)
    
    #reduce resolution to avoid excesive presenc of single pixels 
    r1 <- aggregate(r1, fact=2, fun="max", na.rm=T)
    
    #polygonize the prediction map
    r1 <- as.polygons(r1, dissolve = T)
    r1 <- st_as_sf(r1)
    
    #project to equal area
    r1 <- st_transform(r1, crs='+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs') 
    
    #erase crumbs and fill holes
    th <- units::set_units(200, km^2)
    r1 <- fill_holes(r1, threshold = th)
    r1 <- st_cast(r1, 'POLYGON')
    
    th <- units::set_units(50, km^2)
    r1 <- drop_crumbs(r1, threshold = th)
    
    #smooth edges of polygon and reprojecto to WGS84
    r1 <- smooth(r1, method = "ksmooth", smoothness = 2)
    r1 <- st_transform(r1, crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') 
    r1 <- st_buffer(r1, 0) %>% st_make_valid()
    r1$species <- sp
    r1$area.km2 <- st_area(r1) * 1e-06
    r1$area.km2 <- round(as.numeric(r1$area.km2), 3)
    r1$model <- sp.data$final.model
    r1$datum <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
    r1$method <- 'Maxent v. 3.4.3'
    r1$reference <- 'Gonzalez et al.'
    r1 <- r1 %>% select(species:last_col())
    write_sf(r1, save.dir)
  }
}

#### CLIP TO COLOMBIA FOR MANUAL FILTERING ####
library(sf)
library(tidyverse)
walk(list.files('R/UDF/', '.R$', full.names = T), source)

OCCS <- read.occs()


ranges <- lapply(list.files('output/ranges/original/', pattern = 'gpkg', full.names = T), read_sf)
names(ranges) <- names(OCCS)
COL <- read_sf('/mnt/2TB/GIS/Vectores/Paises/Colombia/COL_adm/COL_adm0.shp')

for(sp in names(ranges)){
  # sp='Marmosa_zeledoni'
  sp.range <- ranges[[ sp ]]
  
  sp.def <- st_intersection(sp.range, COL)
  
  sp.def <- sp.def %>% select(species:reference)
  write_sf(sp.def, paste0('output/ranges/definitive/', sp, '.gpkg'))
}

