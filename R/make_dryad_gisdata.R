library(sf)
library(dplyr)

dir.shapes <- '/mnt/2TB/Research/Submitted/Marmosini Paper/Marmosini SDM project/output/definitive_rages'

files <- lapply(list.files(dir.shapes, pattern = '.shp$', full.names = T), read_sf)
names <- list.files(dir.shapes, pattern = '.shp$') %>% 
  gsub('.shp', '', .) %>% 
  gsub('_', ' ', .)
names(files) <- names

shapes <- lapply(seq_along(files), function(x){
  SHP <- st_combine(files[[x]])
  SHP <- data.frame(area.km2 = sum(as.numeric(st_area(files[[x]]) * 1e-06) %>% round(., 3)), 
                    geometry = SHP) %>% st_as_sf()
  SHP$species <- names(files[x])
  SHP$author <- 'González et al. 2021'
  SHP$datum <- 'WGS84'
  SHP$data.link <- 'BIOC-D-21-00385'
  SHP$data.source <- 'maxent models'
  SHP <-  SHP %>% 
    summarise(species, area.km2 = sum(area.km2), 
              author, datum, data.source, geometry)
  
})

shapes <- do.call(rbind, shapes)

write_sf(shapes, 
         '/mnt/2TB/Research/Submitted/Marmosini Paper/Marmosini SDM project/marmosini_ranges_Colombia_WGS84_20210526.geojson')
