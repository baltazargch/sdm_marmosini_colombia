#### CALCULATE CONSERVATION METRICS ####
library(sf)
library(terra)
library(tidyverse)
library(wdpar)
library(parallel)
sf_use_s2(FALSE)
walk(list.files('R/UDF/', '.R$', full.names = T), source)

dir.create('results')
#Base world and Colombia map
wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
col <- subset(wmap, NAME == 'Colombia')

#Clean WDPA data set or read it if already cleaned
if(!file.exists('WDPA/wdpa_colombia_cleaned.shp')) {
  wdpa0 <- read_sf('WDPA/WDPA_WDOECM_Feb2021_Public_COL_shp/WDPA_WDOECM_Feb2021_Public_COL_shp_0/WDPA_WDOECM_Feb2021_Public_COL_shp-polygons.shp')
  
  wdpa1 <- read_sf('WDPA/WDPA_WDOECM_Feb2021_Public_COL_shp/WDPA_WDOECM_Feb2021_Public_COL_shp_1/WDPA_WDOECM_Feb2021_Public_COL_shp-polygons.shp')
  
  wdpa2 <- read_sf('WDPA/WDPA_WDOECM_Feb2021_Public_COL_shp/WDPA_WDOECM_Feb2021_Public_COL_shp_2/WDPA_WDOECM_Feb2021_Public_COL_shp-polygons.shp')
  
  NPA <- rbind(wdpa0, wdpa1, wdpa2)
  unique(NPA$MARINE)
  plot(NPA[ NPA$MARINE == '1', 'geometry'])
  
  NPA <- NPA %>% 
    filter(MARINE != 3)
  
  NPA_cleaned <- wdpa_clean(NPA)
  
  NPA_cleaned <- st_buffer(NPA_cleaned, 0) %>% st_make_valid() %>% st_transform(4326)
  
  write_sf(NPA_cleaned, 'WDPA/wdpa_colombia_cleaned.shp')
} else {
  NPA_cleaned <- read_sf('WDPA/wdpa_colombia_cleaned.shp')
}

#Apply further filters, project, buffer and clean. Extract IUCN and GOV type
NPA_iucn <- NPA_cleaned %>% 
  filter(!IUCN_CA  %in% c('Not Reported', 'Not Applicable')) %>% 
  mutate(type = ifelse(IUCN_CA  %in% c('Ia', 'II'), 
                       'strict-conservation', 'managed resources')) %>%
  st_transform(crs=3116) %>% st_buffer(0) %>%  st_transform(crs=4326)

NPA_gov <- NPA_cleaned %>% 
  mutate(type = ifelse(GOV_TYP  %in% c('Federal or national ministry or agency'), 
                       'national', ifelse(GOV_TYP  %in% c("Sub-national ministry or agency"), 
                                          'sub-national', 
                                          ifelse(GOV_TYP == 'Individual landowners', 'private', 
                                                 'not reported')))) %>% 
  st_transform(crs=3116) %>% st_buffer(0) %>%  st_transform(crs=4326)

#Read ranges manually modified by expert. See first commento of "7:make_shapes_final_modif.R"
ranges <- lapply(list.files('output/ranges/definitive/', pattern = '.gpkg$', full.names = T), read_sf)
names(ranges) <- gsub('.gpkg$', '', list.files('output/ranges/definitive/', '.gpkg$'))

#Prepare Human Index data. Make discrete scale
if(!file.exists('HFP/pressure_Colombia.tif')){
  dir.create('HFP')
  hfp <- rast('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/GIS_datos/Rasters/Footprint/LHFI_Colombia_Correa2020/IHEH_2015.tif')
  
  hfp_p <- hfp
  
  ref <- 40
  hfp_p[hfp_p <= ref] <- 0
  hfp_p[hfp_p > ref] <- 1
  
  hfp_p <- terra::project(hfp_p, method='near', 'epsg:4326')
  
  writeRaster(hfp_p, 'HFP/pressure_Colombia.tif', overwrite=T)
  
} else {
  hfp <- rast('HFP/pressure_Colombia.tif')
}

#Prepare data.frame to fill with data
dt_results <- data.frame( species = names(ranges) )

if(file.exists('results/species_overlap_pressure_conservation.csv')){
  tbl_ref <-  read.csv('results/species_overlap_pressure_conservation.csv')
}

#Loop for every species calculating coverage areas of conservation and pressure. 
#If new species are added, the script automatically add them and do not calculate
#all over again. 
for(i in 1:NROW(dt_results)){
  # i=2
  sp <- names(ranges)[i]
  cat('Processing data of ', sp, '\n')
  
  if (exists("tbl_ref") && sp %in% tbl_ref$species) {
    cat('Species ', sp, 'already calculated \n')
    next
  }
  rg_sp <- ranges[[sp]]
  rg_sp <- st_union(rg_sp)
  rg_sp <- st_sf(
    data.frame(
      species = sp, 
      area.km2=st_area(rg_sp) * 1e-06, 
      datum='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0',
      method= 'Maxent v. 3.4.3',
      reference= 'Gonzalez et al.',
      geometry=rg_sp), sf_column_name = 'geometry')
  rg_sp$area.km2 <- round(as.numeric(rg_sp$area.km2), 3)
  # rg_sp <- rg_sp %>% st_buffer(0) %>% st_make_valid()
  
  hfp_crop <- terra::mask(terra::crop(hfp, rg_sp), vect(as(rg_sp, 'Spatial')), touches=T)
  
  area_r <- cellSize(hfp_crop)
  area_p <- zonal(area_r, hfp_crop, 'sum', na.rm=T)
  
  cols <- c('total_area', 'low_pressure', 'high_pressure', 'no_data_by_mask')
  dt_results[i, cols ] <- c(
    as.numeric(sum(st_area(rg_sp) * 1e-06)), 
    area_p[ , 2 ] * 1e-06, 
    as.numeric(sum(st_area(rg_sp) * 1e-06)) - sum(area_p[ , 2 ] * 1e-06) 
  )
  
  #### IUCN type overlap ####
  pa_types <- list(
    iucn = st_intersection(rg_sp, NPA_iucn), 
    gov = st_intersection(rg_sp, NPA_gov)
  )
  
  for(k in seq_along(pa_types)){
    # k=2
    for(cat in unique(pa_types[[ k ]]$type)){
      # cat = 'private'
      pa_target <- pa_types[[ k ]] %>% filter(type == cat)
      
      pa_hfp <- crop(hfp_crop, pa_target)
      
      suppressMessages({
        flag <- map(1:NROW(pa_target),
                    ~ st_intersects(st_as_sfc(st_bbox(pa_target[.x,])),
                                    st_as_sfc(st_bbox(pa_hfp)), sparse = F)[,1]) %>% 
          unlist()})
      
      pa_target <- pa_target[flag, ]
      
      dts_out <- mclapply(1:NROW(pa_target), function(j){
        # j=1
        
        pp <- rasterize(vect(pa_target[j,]), crop(pa_hfp, pa_target[j,]), cover=T)
        
        pols_rast <- ifel(pp == 0 ,NA, pp)
        
        area_pa <- cellSize(pols_rast)
        
        v <- values(area_pa)
        v[ is.na(values(pols_rast)) ] <- NA
        v <- (v * 1e-06) * values(pols_rast)
        values(area_pa) <- v
        
        area_z <- zonal(area_pa, terra::crop(pa_hfp, pols_rast), 'sum', na.rm=T)
        if(NROW(area_z) == 0){
          area_z <- data.frame(
            zone = c('0', '1'), 
            pol_proportions = c(0,0)
          )
        } else if(NROW(area_z[,1]) == 1 && area_z[,1] == 0){
          area_z[2,] <- c('1', 0)
          # area_z[2,2] <- as.numeric(area_z[1,2])
          
        } else if(NROW(area_z[,1]) == 1 && area_z[,1] == 1) {
          surr <- area_z
          area_z <- data.frame(
            zone= c('0','1'), 
            pol_proportions = c(0, surr[1,2])
          )
        }
        
        area_z[,2] <- as.numeric(area_z[,2])
        cols <- paste0(cat, c('_total_area', '_low_pres_area', '_high_pres_area', '_no_data'))
        
        dt_out <- data.frame(
          sp, 
          sum(st_area(pa_target[j,]) * 1e-06) %>% as.numeric(), 
          area_z[1,2], area_z[2,2],  
          sum(st_area(pa_target[j,]) * 1e-06) %>% as.numeric() - sum(area_z[ , 2]) 
        )
        
        colnames(dt_out) <- c('species', cols)
        
        return(dt_out)
        
      }, mc.cores=8) #modify cores according to computer hardware
      
      tmpFiles(FALSE, TRUE,FALSE, TRUE)
      
      tbl_out <- do.call(rbind, dts_out)
      tbl_out <- apply(tbl_out[,-1], 2, 'sum')
      
      cols <- paste0(cat, c('_total_area', '_low_pres_area', '_high_pres_area', '_no_data'))
      
      dt_results[ i , cols ] <- c(
        tbl_out
      )
    }
  }
}

# chunk added for new species
if(exists('tbl_ref')){
  dt_results <- dt_results %>% filter(!is.na(total_area))
  colnames(dt_results) <- gsub('-', '.', colnames(dt_results))
  colnames(dt_results) <- gsub('[ ]', '.', colnames(dt_results))
  
  dt_results <- dt_results[ , colnames(tbl_ref)]
  
  stopifnot(identical(colnames(dt_results), colnames(tbl_ref)))
  
  final_dt <- rbind(tbl_ref, dt_results) %>% arrange(species)
  
  stopifnot(identical(final_dt$species, names(ranges)))
  
  write.csv(final_dt, 'results/species_overlap_pressure_conservation.csv', row.names = F, na='')
  
} else {
  
  write.csv(dt_results, 'results/species_overlap_pressure_conservation.csv', row.names = F, na='')
}

#### SPECIES RICHNESS ####

library(BIEN)
library(maptools)    
library(ggplot2)
library(classInt)
library(plotly)
library(raster)
library(letsR)
library(rgdal)
library(sf)

#Load and name range maps
ranges <- lapply(list.files('output/ranges/definitive/', pattern = '.gpkg$', full.names = T), read_sf)
names(ranges) <- gsub('.gpkg$', '', list.files('output/ranges/definitive/', '.gpkg$'))

ranges <- lapply(names(ranges), function(x) { 
  sp <- ranges[[x]]
  
  geom <- st_combine(sp)
  sp <- st_sf(binomial = x, geometry = geom)
})

#Combine in one sf object
SppMap <- do.call(rbind, ranges)
SppMap <- st_transform(SppMap, 
                       '+proj=tmerc +lat_0=4.596200416666666 +lon_0=-74.07750791666666 +k=1 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
e <- extent(SppMap)
PAM <- lets.presab(as(SppMap, 'Spatial'), 
                   xmn = e[1],
                   xmx = e[2],
                   ymn = e[3], 
                   ymx = e[4], 
                   crs = crs(SppMap),
                   resol = 5000, #resolution of output raster
                   remove.sp = T)

#Optional: visualize results
plot(PAM$Richness_Raster, main = "Species richness of the Marmosini from \n Colombia",
     col=c(NA, heat.colors(15, rev=T)))

#Write richness map
richness <- PAM$Richness_Raster
writeRaster(richness, 'results/richness_map_MAGNA-SIRGAS_EPSG3116.tif', 
            overwrite=TRUE)

