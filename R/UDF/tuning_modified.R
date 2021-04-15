tuning.modified <- function (occ, env, bg.coords, occ.grp, bg.grp, method, algorithm, 
                             args, args.lab, categoricals, aggregation.factor, kfolds, 
                             bin.output, clamp, alg, rasterPreds, parallel, numCores, 
                             progbar, updateProgress, userArgs, occs.p, bg.coords.p) 
{
  # pres <- as.data.frame(extract(env, occ))
  # bg <- as.data.frame(extract(env, bg.coords))
  pres <- occs.p
  bg <- bg.coords.p
  numNA.pres <- sum(rowSums(is.na(pres)) > 0)
  numNA.bg <- sum(rowSums(is.na(bg)) > 0)
  if (numNA.pres > 0) {
    message(paste("There are", numNA.pres, "occurrence records with NA for at least\n                  one predictor variable. Removing these records from analysis,\n                  resulting in", 
                  nrow(pres) - numNA.pres, "records..."))
    keepInds <- !apply(pres, 1, anyNA)
    occ <- occ[keepInds, ]
    occ.grp <- occ.grp[keepInds]
    pres <- pres[keepInds, ]
  }
  if (numNA.bg > 0) {
    message(paste("There are", numNA.bg, "background records with NA for at least one predictor variable.\n                  Removing these records from analysis, resulting in", 
                  nrow(bg) - numNA.bg, "records..."))
    keepInds <- !apply(bg, 1, anyNA)
    bg.coords <- bg.coords[keepInds, ]
    bg.grp <- bg.grp[keepInds]
    bg <- bg[keepInds, ]
  }
  if (!is.null(categoricals)) {
    for (i in 1:length(categoricals)) {
      pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
      bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])
    }
  }
  if ("checkerboard1" %in% method) {
    method <- c(method = "checkerboard1", aggregation.factor = aggregation.factor)
    group.data <- get.checkerboard1(occ, env, bg.coords, 
                                    aggregation.factor)
  }
  if ("checkerboard2" %in% method) {
    method <- c(method = "checkerboard2", aggregation.factor = aggregation.factor)
    group.data <- get.checkerboard2(occ, env, bg.coords, 
                                    aggregation.factor)
  }
  if ("block" %in% method) 
    group.data <- get.block(occ, bg.coords)
  if ("jackknife" %in% method) 
    group.data <- get.jackknife(occ, bg.coords)
  if ("randomkfold" %in% method) {
    method <- c(method = "randomkfold", number.folds = kfolds)
    group.data <- get.randomkfold(occ, bg.coords, kfolds)
  }
  if ("user" %in% method) {
    method <- c(method = "user", number.folds = length(unique(occ.grp)))
    group.data <- get.user(occ.grp, bg.grp)
  }
  nk <- length(unique(group.data$occ.grp))
  if (parallel == TRUE) {
    allCores <- detectCores()
    if (is.null(numCores)) {
      numCores <- allCores
    }
    c1 <- makeCluster(numCores)
    registerDoParallel(c1)
    numCoresUsed <- getDoParWorkers()
    message(paste("Of", allCores, "total cores using", numCoresUsed))
    message("Running in parallel...")
    if (algorithm == "maxnet") {
      out <- foreach(i = seq_len(length(args)), .packages = c("dismo", 
                                                              "raster", "ENMeval", "maxnet")) %dopar% {
                                                                modelTune.maxnet(pres, bg, env, nk, group.data, 
                                                                                 args[[i]], rasterPreds, clamp)
                                                              }
    }
    else if (algorithm == "maxent.jar") {
      out <- foreach(i = seq_len(length(args)), .packages = c("dismo", 
                                                              "raster", "ENMeval", "rJava")) %dopar% {
                                                                modelTune.maxentJar(pres, bg, env, nk, group.data, 
                                                                                    args[[i]], userArgs, rasterPreds, clamp)
                                                              }
    }
    stopCluster(c1)
  }
  else {
    out <- list()
    if (progbar == TRUE & !is.function(updateProgress)) {
      pb <- txtProgressBar(0, length(args), style = 3)
    }
    for (i in 1:length(args)) {
      if (length(args) > 1) {
        if (is.function(updateProgress)) {
          text <- paste0("Running ", args.lab[[1]][i], 
                         args.lab[[2]][i], "...")
          updateProgress(detail = text)
        }
        else if (progbar == TRUE) {
          setTxtProgressBar(pb, i)
        }
      }
      if (algorithm == "maxnet") {
        out[[i]] <- modelTune.maxnet(pres, bg, env, 
                                     nk, group.data, args[[i]], rasterPreds, clamp)
      }
      else if (algorithm == "maxent.jar") {
        out[[i]] <- modelTune.maxentJar(pres, bg, env, 
                                        nk, group.data, args[[i]], userArgs, rasterPreds, 
                                        clamp)
      }
    }
    if (progbar == TRUE) 
      close(pb)
  }
  full.mods <- lapply(out, function(x) x[[1]])
  statsTbl <- as.data.frame(t(sapply(out, function(x) x[[2]])))
  if (rasterPreds) {
    predictive.maps <- stack(sapply(out, function(x) x[[3]]))
  }
  else {
    predictive.maps <- stack()
  }
  AUC.DIFF <- statsTbl[, 1:nk]
  AUC.TEST <- statsTbl[, (nk + 1):(2 * nk)]
  OR10 <- statsTbl[, ((2 * nk) + 1):(3 * nk)]
  ORmin <- statsTbl[, ((3 * nk) + 1):(4 * nk)]
  names(AUC.DIFF) <- paste("diff.AUC_bin", 1:nk, sep = ".")
  Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
  Var.AUC.DIFF <- corrected.var(AUC.DIFF, nk)
  names(AUC.TEST) <- paste("AUC_bin", 1:nk, sep = ".")
  Mean.AUC <- rowMeans(AUC.TEST)
  Var.AUC <- corrected.var(AUC.TEST, nk)
  names(OR10) <- paste("test.or10pct_bin", 1:nk, sep = ".")
  Mean.OR10 <- rowMeans(OR10)
  Var.OR10 <- apply(OR10, 1, var)
  names(ORmin) <- paste("test.orMTP_bin", 1:nk, sep = ".")
  Mean.ORmin <- rowMeans(ORmin)
  Var.ORmin <- apply(ORmin, 1, var)
  full.AUC <- double()
  for (i in 1:length(full.mods)) {
    if (algorithm == "maxnet") {
      full.AUC[i] <- dismo::evaluate(pres, bg, full.mods[[i]])@auc
    }
    else if (algorithm == "maxent.jar") {
      full.AUC[i] <- full.mods[[i]]@results[5]
    }
  }
  nparam <- numeric()
  for (i in 1:length(full.mods)) {
    if (algorithm == "maxnet") {
      nparam[i] <- length(full.mods[[i]]$betas)
    }
    else if (algorithm == "maxent.jar") {
      nparam[i] <- get.params(full.mods[[i]])
    }
  }
  aicc <- calc.aicc(nparam, occ, predictive.maps)
  features <- args.lab[[1]]
  rm <- args.lab[[2]]
  settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")
  res <- data.frame(settings, features, rm, train.AUC = full.AUC, 
                    avg.test.AUC = Mean.AUC, var.test.AUC = Var.AUC, avg.diff.AUC = Mean.AUC.DIFF, 
                    var.diff.AUC = Var.AUC.DIFF, avg.test.orMTP = Mean.ORmin, 
                    var.test.orMTP = Var.ORmin, avg.test.or10pct = Mean.OR10, 
                    var.test.or10pct = Var.OR10, aicc)
  if (bin.output == TRUE) {
    res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, 
                               OR10, ORmin))
  }
  if (rasterPreds == TRUE) {
    names(predictive.maps) <- settings
  }
  results <- ENMevaluation(algorithm = alg, results = res, 
                           predictions = predictive.maps, models = full.mods, partition.method = method, 
                           occ.pts = occ, occ.grp = group.data[[1]], bg.pts = bg.coords, 
                           bg.grp = group.data[[2]])
  return(results)
}
