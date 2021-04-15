library(dplyr)
library(ggplot2)
library(here)
library(sf)
library(stringr)

OCCS <- read.csv('output/Marmosini_Colombia_appendix_I.csv')

OCCS %>% 
  count(species) %>% arrange(n)

CHOS <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')
CHOS1 <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')

CHOS <- rbind(CHOS, CHOS1)

CHOS %>% 
  group_by(species) %>% 
  summary(avg.test.AUC)

wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
COL <- wmap[ wmap$NAME == 'Colombia',]

pp <- OCCS

sp::coordinates(pp) <- ~longitude+latitude
pp <- as(pp, 'sf') 
st_crs(pp) <- 4326
pp_col <- st_intersection(pp, wmap)

pp_col %>% group_by(species) %>% 
  select(species, NAME) %>% 
  mutate(In_COL = NAME == 'Colombia') %>% 
  filter(all(In_COL == F)) %>% 
  count()

pp_col %>% as.data.frame() %>% 
  select(NAME) %>% 
  count(NAME) #%>% arrange(-n)

pp_col %>% as.data.frame() %>% 
  select(full.reference) %>% distinct() %>% View()

#### plots AUC ####
library(ggsci)

all_results <- read.csv('output/models/marmosini_fitting_results.csv')

all_results %>% group_by(species) %>% count(case) #%>% View()

results <- all_results %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>% 
  mutate(rm = factor(as.character(rm), levels = unique(as.character(rm)))) %>% 
  group_by(case) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4])

dir.create('/mnt/2TB/Articulos/En_proceso/Marmosini Paper/Manuscript/MS_Marmosini_2021/Fig2' -> dir.save)
#Original plot proposed, then was switched to the Local Polynomial Fit
# results %>%
#   ggplot(aes(x=rm, y=train.AUC, fill=case)) +
#   ggsci::scale_fill_aaas(alpha = 0.8) +
#   geom_boxplot(outlier.size = 0.8) +
#   theme_light() +
#   xlab('regularization multiplier') +
#   ylab('training AUC') + 
#   theme(
#     legend.position = 'bottom', 
#     axis.text.x = element_text(size=14, family = 'sans'),
#     axis.text.y = element_text(size=14, family = 'sans'),
#     axis.title.x = element_text(size=14, family = 'sans'), 
#     axis.title.y = element_text(size=14, family = 'sans'), 
#     panel.grid = element_blank(),
#   ) +facet_grid(area ~ cross.validation)

results %>% 
  filter(cross.validation != 'jackk', 
         # avg.diff.AUC < 0.25, 
         ) %>%
  ggplot(aes(x=AICc, color=case)) + 
  geom_freqpoly() +
  # geom_vline(xintercept = 0.87, linetype='dashed')+
  ggsci::scale_color_d3(alpha = 0.7)+
  facet_grid(cross.validation~area, scales = 'free') +
  theme_light()+
  theme(
    panel.grid = element_blank()
  )

results %>% 
  ggplot(aes(x=rm, y=train.AUC, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(method = NULL, se=T) +
  theme_light() +
  xlab('regularization multiplier') +
  ylab('training AUC') + 
  # facet_grid( ~ cross.validation, scales = 'free')+
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'), 
    panel.grid = element_blank(),
  ) 

ggsave(paste0(dir.save, '/trainAUC.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=avg.test.AUC, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(method = NULL, se=T) +
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average test AUC') + 
  # facet_grid( ~ cross.validation, scales = 'free')+
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/testAUC.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=avg.test.orMTP, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(method = NULL, se=T) +
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average test orMTP') + 
  # facet_grid( ~ cross.validation, scales = 'free')+
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/testorMTP.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=AICc, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(method = NULL, se=T) +
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average AICc') + 
  # facet_grid( ~ cross.validation, scales = 'free')+
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/AICc.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

ggsave(paste0(dir.save, '/legend.png' ), dpi = 150)


# results for M1 and M2 ---------------------------------------------------
m.area <- all_results %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>%
  # mutate(area_cross = str_c(area, '_', cross.validation)) %>% 
  group_by(area, case, rm, cross.validation) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>% 
  mutate(avg_TAUC = mean(train.AUC), 
         avg_TestAUC = mean(avg.test.AUC), 
         avg_orMTP = mean(avg.test.orMTP),
         avg_AIC = mean(AICc)) %>% 
  distinct()

dir.create('/mnt/2TB/Articulos/En_proceso/Marmosini Paper/Manuscript/MS_Marmosini_2021/Fig3' -> dir.save)

m.area %>%
  ggplot(aes(x=rm, shape=cross.validation, y = avg_TAUC, group=method, color=case)) +
  ggsci::scale_color_d3(alpha = 0.8)+
  geom_point(size=2) +
  geom_line()+
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'),     
    strip.text = element_text(size = 14),
    legend.position = 'bottom'
  ) +
  ylab('average train AUC') + xlab('regularization multiplier') 

ggsave(paste0(dir.save, '/trainAUC.png' ),dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3.5,
       units = 'mm')


m.area %>% 
  ggplot(aes(x=rm, shape=cross.validation, y = avg_TestAUC, group=method, color=case)) +
  ggsci::scale_color_d3(alpha = 0.8)+
  geom_point() + 
  geom_line()+
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'),     
    strip.text = element_text(size = 14),
    legend.position = 'bottom'
  ) +
  ylab('average test AUC') + xlab('regularization multiplier') 

ggsave(paste0(dir.save, '/testAUC.png' ),dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3.5,
       units = 'mm')


m.area %>% 
  ggplot(aes(x=rm, shape=cross.validation, y = avg_orMTP, group=method, color=case)) +
  ggsci::scale_color_d3(alpha = 0.8)+
  geom_point() + 
  geom_line()+
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'),     
    strip.text = element_text(size = 14),
    legend.position = 'bottom'
  ) +
  ylab('average orMTP') + xlab('regularization multiplier') 

ggsave(paste0(dir.save, '/orMTP.png' ),dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3.5,
       units = 'mm')

m.area %>% 
  ggplot(aes(x=rm, shape=cross.validation, y = avg_AIC, group=method, color=case)) +
  ggsci::scale_color_d3(alpha = 0.8)+
  geom_point() + 
  geom_line()+
  facet_wrap(~area, ncol=2, strip.position = 'top', scales = 'fixed') +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size=14, family = 'sans'),
    axis.text.y = element_text(size=14, family = 'sans'),
    axis.title.x = element_text(size=14, family = 'sans'), 
    axis.title.y = element_text(size=14, family = 'sans'),     
    strip.text = element_text(size = 14),
    legend.position = 'bottom'
  ) +
  ylab('average AICc') + xlab('regularization multiplier') 

ggsave(paste0(dir.save, '/AICc.png' ),dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3.5,
       units = 'mm')
ggsave(paste0(dir.save, '/legend.png' ), dpi = 300)


# -------------------------------------------------------------------------


# table_areas <- 
#   all_results %>% 
#   group_by(rm, area) %>% 
#   filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>%
#   mutate(mean.TrainAUC = mean(train.AUC), 
#          SD.TrainAUC = sd(train.AUC), 
#          mean.TestAUC = mean(avg.test.AUC), 
#          SD.TestAUC = sd(avg.test.AUC), 
#          mean.orMTP = mean(avg.test.orMTP), 
#          SD.orMTP = sd(avg.test.orMTP), 
#          mean.AICc = mean(AICc), 
#          SD.AICc = sd(AICc)) %>% 
#   add_count() %>% 
#   select(rm, area, 
#          mean.TrainAUC, SD.TrainAUC, 
#          mean.TestAUC, SD.TestAUC, 
#          mean.orMTP, SD.orMTP, 
#          mean.AICc, SD.AICc, n) %>% 
#   distinct() %>% arrange(rm, area) %>% 
#   print()
# 
# write.csv(table_areas, 'output/table_areas_m_metrics.csv', row.names = F)


####model results and range maps####
all_results %>% 
  count(case, species)

optimal <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')
optimal <- optimal %>% mutate(model = 'o')

subopti <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')
subopti <- subopti %>% mutate(model= 's')

NROW(rbind(optimal, subopti))

tbl_final <- rbind(optimal, subopti)

OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)

names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

goodmodels <- data.frame(
  species = names(OCCS), 
  area = c('M2', 'M2', 'M1', 'M2', 'M2',
           'M2', 'M2', 'M2', 'M2', 'M2', 
           'M2', 'M2', 'M2', 'M1', 'M2', 
           'M2'), 
  model = c('o', 's', 's', 'o', 'o', 
            'o', 'o', 'o', 'o', 's', 
            'o', 's', 's', 'o', 'o', 
            'o')
)

for(i in 1:NROW(goodmodels)){
  cc <- tbl_final[ tbl_final$species == goodmodels$species[i], ]
  
  cc <- cc %>% filter(area == goodmodels$area[i] & model == goodmodels$model[i])
  
  goodmodels[ i , colnames(cc)[ -c(1, 3, 22, 23) ]] <- cc[ ,-c(1, 3, 22, 23) ]
}

goodmodels %>% count(case) %>% arrange(n)
goodmodels %>% count(cross.validation)

range(goodmodels$train.AUC)

write.csv(goodmodels, 'output/final_models_chosen.csv', row.names = F)

####predictors importance####
library(dplyr)
library(ggplot2)

tbl_final <- read.csv('output/final_models_chosen.csv')

dt_conper <- data.frame(species=NULL)
dt_toplot <- dt_conper
for(i in 1:NROW(tbl_final)){
  # i=1
  
  sp <- tbl_final[ i , ]
  
  dir <- ifelse(sp$model == 'o', 
                paste0('output/models/final_models/', sp$species, 
                       '/tables/'), 
                paste0('output/models/final_subopt_models/', sp$species, 
                       '/tables/'))
  
  file.chose <- list.files(dir, pattern = paste0('contri_permu_', sp$area), full.names = T)
  
  contri <- read.csv(file.chose)
  
  rows.keep <- grep('contribution|permutation' ,contri$variables)
  
  filtered <- contri[rows.keep,]
  
  filtered$species <- sp$species
  
  filtered <- filtered[ , c(3,1,2)]
  dt_toplot <- plyr::rbind.fill(dt_toplot, filtered)
  tfilt <- filtered[, 3] %>% t() %>%  cbind(sp$species, .) 
  
  colnames(tfilt) <- c('species', filtered[ , 2])
  tfilt <- as.data.frame(tfilt)
  
  dt_conper <- plyr::rbind.fill(dt_conper, tfilt)
}

contri <- c(1, grep('contribution', colnames(dt_conper)))

dt_conper[,-1] <- as.numeric(as.matrix(dt_conper[,-1]))

dt_conper[ , contri] -> cont

apply(cont, 1, which.min) -> col.max 
colnames(cont)[sort(col.max)]

dt_toplot %>% 
  mutate(type = ifelse(stringr::str_detect(variables, 'contribution'), 
                       'contribution', 'permutation')) %>% 
  # filter(type == 'contribution') %>% 
  mutate(variables = factor(variables, levels = sort(unique(variables)))) %>% 
  mutate(variables = stringr::str_replace(variables, '.contribution', '')) %>% 
  mutate(variables = stringr::str_replace(variables, '.permutation.importance', '')) %>% 
  ggplot(aes(y=result, fill=variables, x=species))+
  scale_fill_manual(values = RColorBrewer::brewer.pal(10, 'BrBG'))+
  geom_bar(stat = 'identity',position = 'stack', na.rm = F, show.legend = T) +
  facet_wrap(~type, nrow = 2)+
  coord_flip()+
  theme_light() +
  theme(axis.text.x = element_text(angle = 90)
  )

dt_toplot %>% 
  mutate(type = ifelse(stringr::str_detect(variables, 'contribution'), 
                       'contribution', 'permutation')) %>% 
  # filter(type == 'contribution') %>% 
  mutate(variables = stringr::str_replace(variables, '.permutation.importance|.contribution', '')) %>% 
  arrange(result) %>% 
  mutate(variables = factor(variables, levels = c(paste0('bio_', c(4,6,10,11,15,16,17)), 
                                                  'topoWet','tri', 'msavi'))) %>% 
  ggplot(aes(y=result, fill=variables, x=species), shape='black')+
  scale_fill_manual(name='Variables',
                    labels=c(paste0('bio_', c(4,6,10,11,15,16,17)), 
                             'topoWet', 'tri', 'msavi'),
                    values = c(
                      sample(RColorBrewer::brewer.pal(n = 9, name = "OrRd"), size = 4), 
                      sample(RColorBrewer::brewer.pal(n = 9, name = "Blues"), size = 4), 
                      'hotpink', 'forest green'
                    ))+
  geom_bar(stat = 'identity', na.rm = F) +
  facet_wrap(~type)+
  theme_gray()+
  coord_flip()

write.csv(dt_conper, 'output/contribution_permutation_marmosini.csv', row.names = F)

#### area shapes ####
library(sf)

wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
COL <- wmap[ wmap$NAME == 'Colombia',]

file.shp <- list.files('output/definitive_rages/', pattern = '.shp', full.names = T)
nm.shp <- list.files('output/definitive_rages/', pattern = '.shp', full.names = F)
ranges <- lapply(file.shp, read_sf)
names(ranges) <- gsub('.shp', '', nm.shp)

ranges <- lapply(seq_along(ranges), function(i) {
  x <- ranges[[i]]
  x$species <- names(ranges[i])
  x[, c('species', 'geometry')]
})

joint <- st_union(do.call(rbind, ranges))

plot(joint)

sapply(ranges, function(x){
  paste0(x[1,1], ' ', 
         st_area(st_union(x)) * 1e-06, 
         '\n')
})
(st_area(joint) * 1e-06) / (st_area(COL) * 1e-06 ) *100

#### make appendix 1 S2 ####

library(tidyverse)

cases <- readRDS('Cases/cases.rds')

cases_simple <- lapply(1:3, function(x){
  X <- cases[[x]]
  
  dt <- data.frame(
    Predictors_scenario = rep(names(cases[x]), NROW(X)),
    Variables = X, 
    Species = 'all species', 
    Area = 'M1/M2'
  )
}) %>% do.call(rbind, .)

dt_complicated <- lapply(cases[4], function(x){
  dt_all <- data.frame(
    Predictors_scenario= NULL, 
    Variables =NULL, 
    Species = NULL, 
    Area = NULL
  )
  for(area in names(x)){
    df <- data.frame(
      Predictors_scenario= NULL, 
      Variables =NULL, 
      Species = NULL, 
      Area = NULL
    )
    for(sp in names(x[[area]])){
      dt <- data.frame(
        Predictors_scenario = rep('uncorrelated', NROW(x[[area]][[sp]])),
        Variables = x[[area]][[sp]], 
        Species = names(x[[area]][sp]), 
        Area =str_to_upper(area)
      )
      
      df <- rbind(df, dt)
    }
      dt_all <- rbind(dt_all, df)
  }
  
  dt_all
}) 

dt_final <- rbind(cases_simple, dt_complicated$uncorr)

homolog <- data.frame(
  short.name = unique(dt_final$Variables), 
  long.name = c('Mean Diurnal Range', 'Temperature Seasonality', 
            'Min Temperature of Coldest Month', 
            'Mean Temperature of Warmest Quarter', 
            'Mean Temperature of Coldest Quarter', 
            'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 
            'Precipitation of Driest Quarter', 
            'SAGA-GIS topographic wetness index', 
            'terrain roughness index', 
            'Modified Soil-Adjusted Vegetation Index'
                )
)

dt_final$long.name <- homolog$long.name[ match(dt_final$Variables, homolog$short.name) ]
dt_final[order(dt_final$Species),]

#### make supp mat S1 ####
library(dplyr)
library(raster)
library(sf)
library(ggplot2)
library(ggpubr)

OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)

names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

OCCS %>% names()

wmap <- read_sf('wmap/ne_10m_admin_0_countries.shp')
COL <- wmap[ wmap$NAME == 'Colombia',]

file.shp <- list.files('output/definitive_rages/', pattern = '.shp', full.names = T)
nm.shp <- list.files('output/definitive_rages/', pattern = '.shp', full.names = F)
ranges.def <- lapply(file.shp, read_sf)
names(ranges.def) <- gsub('.shp', '', nm.shp)

file.shp <- list.files('output/final_shapes/', pattern = '.shp', full.names = T)
nm.shp <- list.files('output/final_shapes/', pattern = '.shp', full.names = F)
ranges.fin <- lapply(file.shp, read_sf)
names(ranges.fin) <- gsub('.shp', '', nm.shp)

all.equal(names(OCCS), names(ranges.def), names(ranges.fin))

for(i in seq_along(OCCS)){

ggplot(COL) + 
  geom_sf(aes(fill=CONTINENT), fill='#BCC3CF', colour='#D8E0ED') +
  geom_sf(data=ranges.fin[[ names(OCCS[i]) ]], aes(fill=SOVEREI),fill='orange', colour=NA, size=0.2) +
  coord_sf() +
  labs(title = "Model's map") +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    ) -> ori

ggplot(COL) + 
  geom_sf(aes(fill=CONTINENT), fill='#BCC3CF', colour='#D8E0ED') +
  geom_sf(data=ranges.def[[ names(OCCS[i]) ]], aes(fill=SOVEREI),fill='forest green', colour=NA, size=0.2) +
  coord_sf() +
  labs(title = 'Definitive map') +
  theme_light() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    ) -> def

ggarrange(ori, def) -> gg
ggsave(filename = paste0('/mnt/b9dcef9e-ce57-4276-b96b-9c6402387b2f/Articulos/En_proceso/Marmosini Paper/Manuscript/MS_Marmosini_2021/SuppMat/', names(OCCS[i]), '.png'),
       plot = gg,  dpi = 300
      )
}
