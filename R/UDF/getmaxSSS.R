getmaxSSS <- function(sp, aream, model){
  require(dplyr)
  # model='s'
  # aream = 'M2'
  # sp='Marmosa_alstoni'
  
  basedir <- ifelse(model == 'o', 'output/models/final_models/', 'output/models/final_subopt_models/')
  
  f <- list.files(paste0(basedir, sp,'/tables/'), 'csv', full.names = T)
  
  dt <- read.csv(f[ grep(paste0('fit_metrics_', aream), f)])
  
  maxSSS <- dt$maxSSS
  
  return(maxSSS)
}
