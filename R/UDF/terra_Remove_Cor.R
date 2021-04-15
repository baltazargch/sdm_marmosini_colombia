terraRemoveCor <- function(terraStack, cutline=NULL, n.points=100000, plot=F, 
                           cutoff=0.75, method='pearson', select.vars = T) {
  stopifnot(!is.null(cutline), !missing(terraStack), class(terraStack)[1] == "SpatRaster", 
            class(cutline)[1]  %in% c('sf', 'sfc_POLYGON', "sfc_MULTIPOLYGON"))
  
  require(terra)
  require(sf)
  
  if(any(class(cutline)[1]  %in% c('sfc_POLYGON', "sfc_MULTIPOLYGON"))){
    vsf <- data.frame(pol='pol')
    st_geometry(vsf) <- cutline
    # plot(vsf['geometry'])
  } else {
    vsf <- cutline
  }
  
  stackMask   <- terra::mask(crop(terraStack, vsf), vect(vsf))
  
  if(n.points=='all') {
    
    dt <- as.data.frame(na.omit(values(stackMask)))
    
  } else if(is.numeric(n.points)) {
    
    # set.seed(30)
    flag <- which(!is.na(values(stackMask[[1]])))
    if(length(flag) <= n.points) {
      
      message(paste('Número de pixeles con datos es menor a n.points.', 
                    'Se estimará con el máximo posible', 
                    'ncell: ', length(flag), collapse = ''))
      
      dt <- as.data.frame(na.omit(values(stackMask)))
      
    } else {
      
      dt <- sample(flag, n.points)
      dt <- stackMask[dt]
      dt <- as.data.frame(na.omit(dt))
      
      if(NROW(dt) < n.points) {
        # warning(paste('Los valores obtenidos son menos que el n.points (',n.points, ') : \n', 
        #             '\t Numero de valores: ', NROW(dt), '.\n', 
        #             'se intentará completar automaticamente \n', collapse = ''), call. = F)
        plus <- n.points - NROW(dt)
        dt <- sample(flag, n.points + plus)
        dt <- stackMask[dt]
        dt <- as.data.frame(na.omit(dt))
        
      } else {
        dt <- dt[1:n.points,]
      }
    }
  } else {
    stop('n.points debe ser numerico')
  }
  
  cor_dt <- 1 - abs(stats::cor(dt, method = method))
  dist.matrix <- stats::as.dist(cor_dt)
  ahc <- stats::hclust(dist.matrix, method = "complete")
  groups <- stats::cutree(ahc, h = 1 - cutoff)
  
  if (length(groups) == max(groups)) {
    message(paste("  - No multicollinearity detected in your data at threshold ", 
                  cutoff, "\n", sep = ""))
    return(names(stackMask))
  } else {
    
    if(plot){
      op <- par(no.readonly = TRUE)
      graphics::par(mar = c(5.1, 5.1, 4.1, 3.1))
      plot(ahc, hang = -1, xlab = "", ylab = "Distance (1 - Pearson's r)", 
           main = "", las = 1, sub = "", axes = F)
      graphics::axis(2, at = seq(0, 1, length = 6), las = 1)
      graphics::title(paste("Groups of intercorrelated variables at cutoff", 
                            cutoff))
      par(xpd = T)
      rect.hclust(ahc, h = 1 - cutoff)
      par(op)
    }
    
    if(select.vars) {
      sel.vars <- NULL
      for (i in 1:max(groups)) {
        sel.vars <- c(sel.vars, sample(names(groups[groups == i]), 1))
      }
      return(sel.vars)
    }
  }
}


