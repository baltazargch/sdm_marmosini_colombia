#### OPTIMAL MAPS ####
source('R/4_make_final_models.R')
# read and load libraries and info
library(sf)
library(terra)
library(rangeBuilder)
library(smoothr)
library(ggplot2)
library(scales)
library(tidyverse)
walk(list.files('R/UDF/', '.R$', full.names = T), source)

OCCS <- read.occs()
dir.create(save.dir <- 'output/models/maps/images', recursive = T)

#A base map for plotting
wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')

#Lists the final models files and table
all.files <- list.files('output/models/', recursive = T, full.names = T)
# tuned.table <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')

# make plots
library(doParallel)
options(cores=5) #adjust according to computer hardware
registerDoParallel()
getDoParWorkers()

#Shut on or off parallel processing by commenting or uncommenting the next lines
# foreach(i=seq_along(OCCS)) %dopar% {
for(i in seq_along(OCCS)){
    # i=2
  sp <- names(OCCS[i])
  # sp <- "Marmosa_phaea"
  
  colnames(OCCS[[sp]])[1] <- 'name'
  occs.sp <- OCCS[[sp]]
  sp::coordinates(occs.sp) <- ~longitude+latitude
  
  rast.sp <- all.files[str_detect(all.files, paste0(sp, '/rasters'))]
  tbls.sp <- all.files[str_detect(all.files, paste0(sp, '/tables'))] %>% 
  .[str_detect(., 'contri_permu')]

  for(row in seq_along(rast.sp)){
    # row=2
    table <- read.csv(tbls.sp[row])
    
    sp.prediction <- rast(rast.sp[row])
    
    maxSSS <- table %>% 
      filter(variable == 'Maximum.training.sensitivity.plus.specificity.Cloglog.threshold') %>% 
      dplyr::select(value) %>% as.numeric()
    
    basename <- paste0('output/models/maps/images/', sp, '_',
                       gsub(paste0('output/models//', sp, '/rasters/|.tif$'),'', rast.sp[row]), 
                       '.pdf')
    
    if(file.exists(basename)){
      next
    } else {
      
      #use maxSSS threshold
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
      
      #prepare data to plot with ggplot
      rm <- ifel(sp.prediction < maxSSS, NA, sp.prediction)
      rm <- ifel(is.na(rm), maxSSS, rm)
      rm <- mask(crop(rm, r1), vect(r1))
      
      wmap_crop <- st_crop(wmap, st_bbox(extend(ext(occs.sp) + ext(rm), 0.5)))
      rm <- aggregate(rm, fact=2, fun='max', na.rm=T)
      
      rm <- as.data.frame(rm, xy=T, na.rm=F)
      
      #re-scale prediction to threshold area
      rm$vals <- rescale(rm[,3], to=c(0,1))
      
      #plot and save
      g <- ggplot(wmap_crop) + 
        geom_sf(aes(fill=REGION_UN), fill='#BCC3CF', colour='#D8E0ED') +
        geom_raster(data= rm, aes(x, y, fill= vals), alpha = 0.8, 
                    na.rm = T, interpolate = T) + 
        scale_fill_gradientn(colours = c('#0274ba','#1ddd00', '#ffff00', 
                                         '#fd7a00', '#d70003'), 
                             na.value = NA, guide = 'colourbar', name='Suitability', 
                             breaks=c(0,0.5,1),
                             labels=c('lowest','medium', 'highest')) + 
        geom_sf(data=r1, aes(fill=layer),fill=NA, colour='black', size=0.2) +
        scale_size_identity()+
        coord_sf() +
        geom_point(data = OCCS[[sp]], 
                   aes(x = longitude, y = latitude,
                       shape=factor(name, labels = 'known localities')), 
                   size = 1.8, show.legend = T) + 
        scale_shape_manual(values=4) +
        labs(shape = "")+
        labs(title = gsub("_", " ", sp), x = 'Longitude', y = 'Latitude', 
             subtitle = basename %>% gsub(paste0('output/models/maps/images/', sp, '_|.pdf'), '', .)) +
        theme_light() +
        theme(title = element_text(face = "bold.italic"), 
              plot.subtitle = element_text(face = 'plain'), 
              axis.title = element_text(face = 'plain'), 
              legend.title = element_text(face = 'bold'))
      
      ggsave(basename, plot = g, dpi=320, units= 'cm', width = 30, height = 20)
      # gc()
      # tmpFiles(FALSE, TRUE, TRUE, TRUE)
    }
  }
}

 # create final PDF ----------------------------------------------------------
library(pdftools)

files <- list.files(save.dir, full.names = T, pattern = '.pdf$')

pdftools::pdf_combine(files, 
                      output = 'output/models/maps/all_maps_Marmosini.pdf')
