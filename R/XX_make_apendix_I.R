files <- list.files('records/', pattern = '.csv', full.names = T)

df <- data.frame()


for(i in files){
  dt <- read.csv(i)
  
  colnames(dt) <- c("species", "longitude", "latitude", "reference", "GBIF ID", "full reference", "comment")
  
  df <- rbind(df, dt)

}

df$comment[df$comment == ""] <- NA

write.csv(df, 'output/Marmosini_Colombia_appendix_I.csv', na='', row.names = F)


