AIC.predict <- function(TableAIC, OCCS, modelsDir = ''){
  if (class(tableAIC) == "data.frame" && class(OCCS) == "list") {
    cat(paste("Procesing results for", NROW(TableAIC), "different models and for", length(OCCS), "species", "\n"))
    
    AIC.per.species <- split(tableAIC, tableAIC[,1])
    
    for (i in seq(AIC.per.species)){
      
        if (NROW(AIC.per.species[[i]]) > 2){
        cat("More than 1 option for either M1 or M2. Choosing by diff.AUC metric \n")
        
        save.path <- paste0(modelsDir, "/Finalmodels/", names(OCCS[i]))
        
        load.path <- paste0(modelsDir, "/Second Tuning/", names(OCCS[i]))
        
        K <- AIC.per.species[[i]] 
        
        toChoose <- K[,2]
        
        filterChoose <- gsub("[[:lower:]]*._.*", "", toChoose)
        
        toChoose <- cbind(Area = filterChoose, K[,-1])
        
        W <- split(toChoose, toChoose$Area)
        
        if (NROW(W[[1]]) == 1){
          M1is <- row.names(W[[1]])
        } else {
          M1is <- W[[1]]
          M1is <- row.names(M1is[which(M1is$avg.diff.AUC == min(M1is$avg.diff.AUC)),])
        }
          
        if (NROW(W[[2]]) == 1) {
          Z <- row.names(W[[2]])
        } else {
          Z <- W[[2]]
          Z <- row.names(Z[which(Z$avg.diff.AUC == min(Z$avg.diff.AUC)),])
        }
        
        toChoose <- rbind(K[M1is,], K[Z,])
        
        e <- extent(M1bgxy[[i]])
        print(paste("Cutting predictors to M1 extent of", names(OCCS[i])))
        envPredics1 <- raster::crop(envPredics, e)
        envPredics1 <- raster::mask(envPredics1, M1bgmask[[i]])
        
        e1 <- extent(M2bgxy[[i]])
        print(paste("Cutting predictors to M2 extent of", names(OCCS[i])))
        envPredics2 <- raster::crop(envPredics, e1)
        envPredics2 <- raster::mask(envPredics2, M2bgmask[[i]])
        
        M1 <- get(load(paste0(load.path, "/", toChoose[1,2], ".RData")))
        
        dismo::predict(M1@models[[1]], envPredics1, args=c("outputformat=logistic"), progress="text", 
                       filename = paste0(save.path, "/", toChoose[1,2], "_log.grd"))
        ModM1 <- M1@models
        
        a <- as.numeric((7+length(envPredics@layers))) #seven is where variables importances are in ENMeval
        
        var.contri <- ModM1[[1]]@results[7:as.numeric(a-1),]
        
        var.permu  <- ModM1[[1]]@results[as.numeric(a):as.numeric((a+length(envPredics@layers)-1)),]
        
        vars.names <- gsub(".contribution", "", names(var.contri))
        
        df <- as.data.frame(cbind(Variables = vars.names, 
                                  Contribution = var.contri, 
                                  Permutation = as.numeric(var.permu)))
        
        write.table(df, file=paste0(save.path, "/", toChoose[1,2], "var.importance.csv"), 
                    sep=";", row.names = F)
        
        M2 <- get(load(paste0(load.path, "/", toChoose[2,2], ".RData")))
        
        dismo::predict(M2@models[[1]], envPredics2, args=c("outputformat=logistic"), progress="text", 
                       filename = paste0(save.path, "/", toChoose[2,2], "_log.grd"))
        ModM2 <- M2@models
        
        a <- as.numeric((7+length(envPredics@layers))) #seven is where variables importances are in ENMeval
        
        var.contri <- ModM2[[1]]@results[7:as.numeric(a-1),]
        
        var.permu  <- ModM2[[1]]@results[as.numeric(a):as.numeric((a+length(envPredics@layers)-1)),]
        
        vars.names <- gsub(".contribution", "", names(var.contri))
        
        df <- as.data.frame(cbind(Variables = vars.names, 
                                  Contribution = var.contri, 
                                  Permutation = as.numeric(var.permu)))
        
        write.table(df, file=paste0(save.path, "/", toChoose[2,2], "var.importance.csv"), 
                    sep=";", row.names = F)
        
      } else {
        
        save.path <- paste0(modelsDir, "/Finalmodels/", names(OCCS[i]))
        
        load.path <- paste0(modelsDir, "/Second Tuning/", names(OCCS[i]))
        
        e <- extent(M1bgxy[[i]])
        print(paste("Cutting predictors to M1 extent of", names(OCCS[i])))
        envPredics1 <- raster::crop(envPredics, e)
        envPredics1 <- raster::mask(envPredics1, M1bgmask[[i]])
        
        e1 <- extent(M2bgxy[[i]])
        print(paste("Cutting predictors to M2 extent of", names(OCCS[i])))
        envPredics2 <- raster::crop(envPredics, e1)
        envPredics2 <- raster::mask(envPredics2, M2bgmask[[i]])
        
        X <- AIC.per.species[[i]]
        
        toChoose <- as.character(X[,2])
        
        M1 <- get(load(paste0(load.path, "/", toChoose[1], ".RData")))
        
        dismo::predict(M1@models[[1]], envPredics1, args=c("outputformat=logistic"), progress="text", 
                       filename = paste0(save.path, "/", toChoose[1], "_log.grd"))
        ModM1 <- M1@models
        
        a <- as.numeric((7+length(envPredics@layers))) #seven is where variables importances are in ENMeval
        
        var.contri <- ModM1[[1]]@results[7:as.numeric(a-1),]
        
        var.permu  <- ModM1[[1]]@results[as.numeric(a):as.numeric((a+length(envPredics@layers)-1)),]
        
        vars.names <- gsub(".contribution", "", names(var.contri))
        
        df <- as.data.frame(cbind(Variables = vars.names, 
                                  Contribution = var.contri, 
                                  Permutation = as.numeric(var.permu)))
        
        write.table(df, file=paste0(save.path, "/", toChoose[1], "var.importance.csv"), 
                    sep=";", row.names = F)
        
        M2 <- get(load(paste0(load.path, "/", toChoose[2], ".RData")))
        
        dismo::predict(M2@models[[1]], envPredics2, args=c("outputformat=logistic"), progress="text", 
                       filename = paste0(save.path, "/", toChoose[2], "_log.grd"))
        ModM2 <- M2@models
        
        a <- as.numeric((7+length(envPredics@layers))) #seven is where variables importances are in ENMeval
        
        var.contri <- ModM2[[1]]@results[7:as.numeric(a-1),]
        
        var.permu  <- ModM2[[1]]@results[as.numeric(a):as.numeric((a+length(envPredics@layers)-1)),]
        
        vars.names <- gsub(".contribution", "", names(var.contri))
        
        df <- as.data.frame(cbind(Variables = vars.names, 
                                  Contribution = var.contri, 
                                  Permutation = as.numeric(var.permu)))
        
        write.table(df, file=paste0(save.path, "/", toChoose[2], "var.importance.csv"), 
                    sep=";", row.names = F)
        
      }
    }
  }
}

