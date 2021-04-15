speciesAreas <- function(SpeciesDist, th, Country, ConservationArea){
  
  SpeciesDist[SpeciesDist >= as.numeric(th)] <- 1
  SpeciesDist[SpeciesDist < as.numeric(th)]  <- NA
  
  cell_size        <- raster::area(SpeciesDist, na.rm=TRUE, weights=FALSE)
  cell_size        <- cell_size[!is.na(cell_size)]
  TotalDist        <- length(cell_size)*median(cell_size)
  
  M1masked         <- raster::crop(SpeciesDist, Country)
  M1masked         <- raster::mask(SpeciesDist, Country)
  
  flagCol <- M1masked
  flagCol[is.na(flagCol)] <- 0
  flagCol <- unique(flagCol@data@values)
  
  if (length(flagCol) > 1){
    cell_size        <- raster::area(M1masked, na.rm=TRUE, weights=FALSE)
    cell_size        <- cell_size[!is.na(cell_size)]
    CountryDist      <- length(cell_size)*median(cell_size)
    
    M1xPNN           <- raster::crop(M1masked, ConservationArea) 
    M1xPNN           <- raster::mask(M1masked, ConservationArea)
    
    cell_size        <- raster::area(M1xPNN, na.rm=TRUE, weights=FALSE)
    cell_size        <- cell_size[!is.na(cell_size)]
    ConservationDist <- length(cell_size)*median(cell_size)
    
    XY <- as.data.frame(rasterToPoints(M1xPNN))
    sp::coordinates(XY) <- cbind(XY$x, XY$y)
    crs(XY) <- crs(ConservationArea)
    number <- XY %over% ConservationArea
    
    numberAreas <- length(unique(number$nombre))
    
  } else{
    CountryDist <- "No overlap w/ Country"
    ConservationDist <- "No overlap w/ Conservation Area"
    numberAreas <- ConservationDist
  }
  
  DistributionInfo <- cbind(Total.Area = TotalDist, 
                            Dist.Area.Col = CountryDist, 
                            Dist.Area.PNN = ConservationDist,
                            Number.PNN = numberAreas)
  
  return(DistributionInfo)
}