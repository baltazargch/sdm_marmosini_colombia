library(raster)
library(dplyr)

OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)

names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))


env_rast <- raster('LayersBank/bio_1.grd')

for(i in names(OCCS)){
  # i = 'Monodelphis_brevicaudata'
  X <- OCCS[[i]]
  
  dt_rast <- raster::extract(env_rast, X[c(2, 3)], cellnumbers=T, df=T)
  
  dt_rast$rep <- duplicated(dt_rast$cells)
  
  OCCS[[i]]$rep <- dt_rast$rep 
  
}

for(i in names(OCCS)){
  X <- OCCS[[i]]
  
  if(any(X$rep == T)){
    cat('Species', i, 'has duplicated coordinates \n\n')
  } else {
    cat('Species', i, 'does not have duplicated coordinates\n\n')
  }
}

purrr::walk2(OCCS, list.files('records', pattern = '.csv$', full.names = T), 
             readr::write_csv)

OCCS$Marmosa_germana %>% View()
