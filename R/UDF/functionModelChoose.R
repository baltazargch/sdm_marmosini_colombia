choose.Models <- function(x){
  if (typeof(x) == "list"){
    
    tableMetrics <- vector("list", length(x))
    names(tableMetrics) <- names(x)
    
    for (i in seq(x)){
      
      X <- x[[i]]
      
      minAUCdiff <- as.numeric(row.names(X[order(X$avg.diff.AUC),][1,]))
      
      bestMTP <- X[order(X$avg.test.orMTP),][1:5,]
      
      if (any(bestMTP$avg.test.orMTP[2:5] == bestMTP$avg.test.orMTP[1])){
        tableMTP <- X[which(X$avg.test.orMTP == min(bestMTP$avg.test.orMTP)),]
        tableMTP <- tableMTP[which(tableMTP$avg.diff.AUC == min(tableMTP$avg.diff.AUC)),][1,]
        minMTP <- as.numeric(row.names(tableMTP))
      } else {
        minMTP <- as.numeric(row.names(bestMTP[1,]))
        }
      
      if (minAUCdiff != minMTP){
        
        cat("For", names(x[i]), "\n","minimum AUCdiff and orMTP are different", "\n")
        
        DiffnMTP <- as.data.frame(rbind(X[minAUCdiff,], 
                                        X[minMTP,]))
        
        DiffnMTP <- cbind(Criteria=c("lwstAUCdiff", "lwstorMTP"), 
                          DiffnMTP)
        
        tableMetrics[[i]] <-  DiffnMTP
        
      } else {
        
        cat("For", names(x[i]), 
            "Minimum AUCdiff and orMTP are the same", "\n", "picking second best orMTP", "\n")
        
        Dfnn <- X[-minAUCdiff,]
        Dfnn <- Dfnn[order(Dfnn$avg.test.orMTP),][1:5,]
        Dfnn <- Dfnn[1,]
        
        minMTPn <- as.numeric(row.names(Dfnn))
        
        M1blockProcessed <- as.data.frame(rbind(X[minAUCdiff,], 
                                                 X[minMTPn,]))
        
        M1blockProcessed <- cbind(Criteria=c("lwstAUCdiff", "lwstorMTP"), 
                                  M1blockProcessed)

        tableMetrics[[i]] <-  M1blockProcessed
      }
    }
    
    tableFinal <- tableMetrics
    
    return(tableFinal)
    
  } else {
    cat("Wrong type of data. Object x must be a list", "\n")
    }
  }
