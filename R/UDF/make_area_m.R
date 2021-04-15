make_area_m <- function(sp.occ=NULL, method=c('simple', 'ecoregion'), ecoregion = NULL, 
                        save.area.m.as.shapes=T, file=NULL,
                        buffer.eco=0.5, buffer.simple=3, cut.terrestrial=T, cutline=NULL,
                        crs.proj = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs",
                        crs.out ="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") {
  require(sp)
  require(sf)
  require(adehabitatHR)
  # sp.occ <- read.csv('records/Marmosa_alstoni.csv')
  coordinates(sp.occ) <-  cbind(sp.occ[ ,2], sp.occ[ ,3])
  proj4string(sp.occ) <-  CRS(crs.out)
  
  sp.occ <- sp::spTransform(sp.occ, CRS(crs.proj))
  
  bg.mcp <- as(mcp(as(sp.occ, 'SpatialPoints'), percent = 100), 'sf')
  
  if(method == 'simple'){
    
    area.m <- st_buffer(bg.mcp, buffer.simple * 1e05) %>% st_transform(crs.out)
    
    if(cut.terrestrial & !is.null(cutline)) {
      bbox <- st_bbox(area.m)
      
      cutted <- st_intersection(cutline, st_as_sfc(bbox))
      
      area.m <- area.m %>% st_intersection(cutted) %>% st_union()
    }
    
    if(save.area.m.as.shapes) sf::write_sf(area.m, file)
    return(area.m)
    
  } else if(method == 'ecoregion' & !is.null(ecoregion)) {
    
    area.m <- st_buffer(bg.mcp, buffer.eco * 1e05) %>% 
      st_transform(crs.out)
    
    # bbox <- st_bbox(area.m)
    ecos.in <- st_intersection(ecoregion, area.m)
    area.m <- subset(ecoregion, ECO_NAME  %in% ecos.in$ECO_NAME) %>% st_union() 
    
    if(save.area.m.as.shapes) sf::write_sf(area.m, file)
    return(area.m)
  }
}


