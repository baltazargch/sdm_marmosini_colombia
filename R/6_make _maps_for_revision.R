#### OPTIMAL MAPS ####

# read and load libraries and info
library(sf)
library(raster)
library(rangeBuilder)
library(smoothr)
library(ggplot2)
library(scales)
library(tidyverse)

dir.create(save.dir <- 'output/models/final_models/maps/images', recursive = T)

#Read and name the species records
OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)
OCCS <- lapply(OCCS, function(x) { subset(x, rep != TRUE) })
OCCS <- lapply(seq_along(OCCS), function(x) {
  X <- OCCS[[x]]
  X <- X[,1:3]
  colnames(X)  <- c('name', 'longitude', 'latitude')
  X
})

names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

#A base map for plotting
wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')

#Lists the final models files and table
all.files <- list.files('output/models/final_models/', recursive = T, all.files = T)
tuned.table <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')

# make plots
library(doParallel)
options(cores=7) #adjust according to computer hardware
registerDoParallel()
getDoParWorkers()

#Shut on or off parallel processing by commenting or uncommenting the next lines
# foreach(i=seq_along(OCCS)) %dopar% {
  for( i in seq_along(OCCS)){
  sp <- names(OCCS[i])
  sp.tuned <- subset(tuned.table, species == sp)
  
  for(row in 1:nrow(sp.tuned)){
    table <- read.csv( paste0(
      'output/models/final_models/',
      sp.tuned$species[row], '/tables/fit_metrics_',
      sp.tuned$area[row], '_', sp.tuned$case[row], '_', 
      sp.tuned$settings[row], '_', sp.tuned$cross.validation[row],
      '.csv') )
    
    sp.prediction <- raster( paste0(
      'output/models/final_models/',
      sp.tuned$species[row], '/rasters/',
      sp.tuned$area[row], '_', sp.tuned$case[row], '_', 
      sp.tuned$settings[row], '_', sp.tuned$cross.validation[row],
      '.tif') )
    
    maxSSS <- table$maxSSS
    
    if(file.exists(paste0(save.dir, '/', sp, '_', 
                          paste(sp.tuned$area[row], sp.tuned$case[row], 
                                sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                                sep='_'), '_maxSSS_', maxSSS,'.pdf'))) {
      next
      
    } else {
      
      #use maxSSS threshold
      rm <- r <- sp.prediction
      r[r < maxSSS] <- NA
      
      r1 <- r
      r1[!is.na(r1)] <- 1
      
      #reduce resolution to avoid excesive presenc of single pixels 
      r1 <- raster::aggregate(r1, fact=2, fun=mean)
      
      #polygonize the prediction map
      r_modif <- rasterToPolygons(r1, dissolve = T)
      r_modif <- st_as_sf(r_modif)
      
      #project to equal area
      r_modif <- st_transform(r_modif, crs='+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs') 
      
      #erase crumbs and fill holes
      th <- units::set_units(200, km^2)
      
      r_modif <- fill_holes(r_modif, threshold = th)
      r_modif <- st_cast(r_modif, 'POLYGON')
      
      th <- units::set_units(50, km^2)
      r_modif <- drop_crumbs(r_modif, threshold = th)
      
      #smooth edges of polygon and reprojecto to WGS84
      r_modif <- smooth(r_modif, method = "ksmooth", smoothness = 2)
      r_modif <- st_transform(r_modif, 
                              crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ') 
      
      #prepare data to plot with ggplot
      rm[rm < maxSSS] <- NA
      rm[is.na(rm)] <- min(getValues(rm), na.rm = T)
      rm <- raster::mask(crop(rm, as(r_modif, 'Spatial')), r_modif)
      
      wmap_crop <- st_crop(wmap, st_bbox(extend(extent(rm), 1)))
      rma <- aggregate(rm, 4)
      
      toPlot <- raster::as.data.frame(rma, xy=T)
      
      #re-scale prediction to threshold area
      aa <- rescale(toPlot[,3], to=c(0,1))
      rma <- setValues(rma, aa)
      toPlot[,3] <- aa
      
      #plot and save
      g <- ggplot(wmap_crop) + 
        geom_sf(aes(fill=REGION_UN), fill='#BCC3CF', colour='#D8E0ED') +
        geom_raster(data= toPlot, aes(x, y, fill= toPlot[,3]), alpha = 0.8, 
                    na.rm = T, interpolate = T) + 
        scale_fill_gradientn(colours = c('#0274ba','#1ddd00', '#ffff00', 
                                         '#fd7a00', '#d70003'), 
                             na.value = NA, guide = 'colourbar', name='Suitability', 
                             breaks=c(0,0.5,1),
                             labels=c('lowest','medium', 'highest')) + 
        geom_sf(data=r_modif, aes(fill=layer),fill=NA, colour='black', size=0.2) +
        scale_size_identity()+
        coord_sf() +
        geom_point(data = OCCS[[sp]], 
                   aes(x = longitude, y = latitude,
                       shape=factor(name, labels = 'known localities')), 
                   size = 1.8, show.legend = T) + 
        scale_shape_manual(values=4) +
        labs(shape = "")+
        labs(title = gsub("_", " ", sp), x = 'Longitude', y = 'Latitude', 
             subtitle = paste(sp.tuned$area[row], sp.tuned$case[row], 
                              sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                              sep=' ')) +
        theme_light() +
        theme(title = element_text(face = "bold.italic"), 
              plot.subtitle = element_text(face = 'plain'), 
              axis.title = element_text(face = 'plain'), 
              legend.title = element_text(face = 'bold'))
      
      ggsave(paste0(save.dir, '/', sp, '_', 
                    paste(sp.tuned$area[row], sp.tuned$case[row], 
                          sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                          sep='_'), '_maxSSS_', maxSSS,'.pdf'), 
             plot = g, dpi=320, units= 'cm', width = 30, height = 20)
    }
  }
}

# create final PDF ----------------------------------------------------------
library(pdftools)

files <- list.files(save.dir, full.names = T, pattern = '.pdf$')

pdftools::pdf_combine(files, 
                      output = 'output/models/final_models/maps/all_maps_Marmosini.pdf')

#### SUBOPTIMAL MAPS ####

dir.create(save.dir <- 'output/models/final_subopt_models/maps/images', recursive = T)

all.files <- list.files('output/models/final_subopt_models/', recursive = T, all.files = T)
tuned.table <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')

library(doParallel)
options(cores=7)
registerDoParallel()
getDoParWorkers()

foreach(i=seq_along(OCCS)) %dopar% {
  # for( i in seq_along(OCCS)){
  sp <- names(OCCS[i])
  sp.tuned <- subset(tuned.table, species == sp)

  for(row in 1:nrow(sp.tuned)){
    table <- read.csv( paste0(
      'output/models/final_subopt_models/',
      sp.tuned$species[row], '/tables/fit_metrics_',
      sp.tuned$area[row], '_', sp.tuned$case[row], '_', 
      sp.tuned$settings[row], '_', sp.tuned$cross.validation[row],
      '.csv') )
    
    sp.prediction <- raster( paste0(
      'output/models/final_subopt_models/',
      sp.tuned$species[row], '/rasters/',
      sp.tuned$area[row], '_', sp.tuned$case[row], '_', 
      sp.tuned$settings[row], '_', sp.tuned$cross.validation[row],
      '.tif') )
    
    maxSSS <- table$maxSSS
    
    if(file.exists(paste0(save.dir, '/', sp, '_', 
                          paste(sp.tuned$area[row], sp.tuned$case[row], 
                                sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                                sep='_'), '_maxSSS_', maxSSS,'.pdf'))) {
      print(sp)
      next
      
    } else {
      
      rm <- r <- sp.prediction
      r[r < maxSSS] <- NA
      
      r1 <- r
      r1[!is.na(r1)] <- 1
      
      r1 <- raster::aggregate(r1, fact=2, fun=mean)
      r_modif <- rasterToPolygons(r1, dissolve = T)
      r_modif <- st_as_sf(r_modif)
      
      r_modif <- st_transform(r_modif, crs='+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs') 
      
      th <- units::set_units(200, km^2)
      
      r_modif <- fill_holes(r_modif, threshold = th)
      r_modif <- st_cast(r_modif, 'POLYGON')
      
      th <- units::set_units(50, km^2)
      r_modif <- drop_crumbs(r_modif, threshold = th)
      
      r_modif <- smooth(r_modif, method = "ksmooth", smoothness = 2)
      r_modif <- st_transform(r_modif, 
                              crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ') 
      
      rm[rm < maxSSS] <- NA
      rm[is.na(rm)] <- min(getValues(rm), na.rm = T)
      rm <- raster::mask(crop(rm, as(r_modif, 'Spatial')), r_modif)
      
      wmap_crop <- st_crop(wmap, st_bbox(extend(extent(rm), 1)))
      rma <- aggregate(rm, 4)
      
      toPlot <- raster::as.data.frame(rma, xy=T)
      
      aa <- rescale(toPlot[,3], to=c(0,1))
      rma <- setValues(rma, aa)
      toPlot[,3] <- aa
      
      g <- ggplot(wmap_crop) + 
        geom_sf(aes(fill=REGION_UN), fill='#BCC3CF', colour='#D8E0ED') +
        geom_raster(data= toPlot, aes(x, y, fill= toPlot[,3]), alpha = 0.8, 
                    na.rm = T, interpolate = T) + 
        scale_fill_gradientn(colours = c('#0274ba','#1ddd00', '#ffff00', 
                                         '#fd7a00', '#d70003'), 
                             na.value = NA, guide = 'colourbar', name='Suitability', 
                             breaks=c(0,0.5,1),
                             labels=c('lowest','medium', 'highest')) + 
        geom_sf(data=r_modif, aes(fill=layer),fill=NA, colour='black', size=0.2) +
        scale_size_identity()+
        coord_sf() +
        geom_point(data = OCCS[[sp]], 
                   aes(x = longitude, y = latitude,
                       shape=factor(name, labels = 'known localities')), 
                   size = 1.8, show.legend = T) + 
        scale_shape_manual(values=4) +
        labs(shape = "")+
        labs(title = gsub("_", " ", sp), x = 'Longitude', y = 'Latitude', 
             subtitle = paste(sp.tuned$area[row], sp.tuned$case[row], 
                              sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                              sep=' ')) +
        theme_light() +
        theme(title = element_text(face = "bold.italic"), 
              plot.subtitle = element_text(face = 'plain'), 
              axis.title = element_text(face = 'plain'), 
              legend.title = element_text(face = 'bold'))

      ggsave(paste0(save.dir, '/', sp, '_', 
                    paste(sp.tuned$area[row], sp.tuned$case[row], 
                          sp.tuned$settings[row], sp.tuned$cross.validation[row], 
                          sep='_'), '_maxSSS_', maxSSS,'.pdf'), 
             plot = g, dpi=320, units= 'cm', width = 30, height = 20)
    }
  }
}

# create final PDF ----------------------------------------------------------
library(pdftools)

files <- list.files(save.dir, full.names = T, pattern = '.pdf$')

pdftools::pdf_combine(files, 
                      output = 'output/models/final_subopt_models/maps/subopt_all_maps_Marmosini.pdf')

