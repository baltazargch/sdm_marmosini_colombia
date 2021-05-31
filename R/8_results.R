library(here)
library(sf)
library(tidyverse)
library(ggsci)
library(raster)
library(ggpubr)

#### GENERAL RESULTS ####

# Make Appendix I - Sheet 1 
files <- list.files('records/', pattern = '.csv', full.names = T)

df <- data.frame()
for(i in files){
  dt <- read.csv(i)
  colnames(dt) <- c("species", "longitude", "latitude", "reference", 
                    "GBIF ID", "full reference", "comment")
  df <- rbind(df, dt)
}

df$comment[df$comment == ""] <- NA

write.csv(df, 'output/Marmosini_Colombia_appendix_I.csv', na='', row.names = F)

# Make Appendix I - Sheet 2

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
dt_final[order(dt_final$Species),] #Copied and pasted manually into excel sheet

#------------------------------///------------------------------#

##### Other results ####
OCCS <- read.csv('output/Marmosini_Colombia_appendix_I.csv')

# Records per species
OCCS %>% 
  count(species) %>% arrange(n)

# Models chosen generalities
CHOS <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')
CHOS1 <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')

CHOS <- rbind(CHOS, CHOS1)

CHOS %>% 
  group_by(species) %>% 
  summary(avg.test.AUC)

# Records per country
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
  count(NAME)

# Number of different references in Appendix I - Sheet 1
pp_col %>% as.data.frame() %>% 
  select(full.reference) %>% distinct() %>% View()

#### EVALUATION METRICS ####
all_results <- read.csv('output/models/marmosini_fitting_results.csv')

all_results %>% group_by(species) %>% count(case) #%>% View()

results <- all_results %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>% 
  mutate(rm = factor(as.character(rm), levels = unique(as.character(rm)))) %>% 
  group_by(case) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4])

dir.create(dir.save <- 'results/Fig2', recursive = T)

##### Fig 2 plots #####
results %>% 
  ggplot(aes(x=rm, y=train.AUC, colour=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(aes(linetype=case),method = NULL, se=T, show.legend = T) +
  scale_linetype_manual(values = c('dotted', 'twodash', 'solid', 'dotdash'))+
  theme_light() +
  xlab('regularization multiplier') +
  ylab('training AUC') + 
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=16, family = 'sans'),
    axis.text.y = element_text(size=16, family = 'sans'),
    axis.title.x = element_text(size=16, family = 'sans'), 
    axis.title.y = element_text(size=16, family = 'sans'), 
    panel.grid = element_blank(),
  ) 

ggsave(paste0(dir.save, '/trainAUC.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=avg.test.AUC, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(aes(linetype=case),method = NULL, se=T, show.legend = T) +
  scale_linetype_manual(values = c('dotted', 'twodash', 'solid', 'dotdash'))+
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average test AUC') + 
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=16, family = 'sans'),
    axis.text.y = element_text(size=16, family = 'sans'),
    axis.title.x = element_text(size=16, family = 'sans'), 
    axis.title.y = element_text(size=16, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/testAUC.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=avg.test.orMTP, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(aes(linetype=case),method = NULL, se=T, show.legend = T) +
  scale_linetype_manual(values = c('dotted', 'twodash', 'solid', 'dotdash'))+
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average test orMTP') + 
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=16, family = 'sans'),
    axis.text.y = element_text(size=16, family = 'sans'),
    axis.title.x = element_text(size=16, family = 'sans'), 
    axis.title.y = element_text(size=16, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/testorMTP.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

results %>% 
  ggplot(aes(x=rm, y=AICc, color=case, group=case))+
  ggsci::scale_color_aaas(alpha = 0.7)+
  geom_smooth(aes(linetype=case),method = NULL, se=T, show.legend = T) +
  scale_linetype_manual(values = c('dotted', 'twodash', 'solid', 'dotdash'))+
  theme_light() +
  xlab('regularization multiplier') +
  ylab('average AICc') + 
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(size=16, family = 'sans'),
    axis.text.y = element_text(size=16, family = 'sans'),
    axis.title.x = element_text(size=16, family = 'sans'), 
    axis.title.y = element_text(size=16, family = 'sans'), 
    panel.grid = element_blank(),
  ) 
ggsave(paste0(dir.save, '/AICc.png' ), dpi = 300, 
       width = 70.27 *3, height = 48.77 * 3,
       units = 'mm')

##### Other results #####
# Count predictors scenarios frequencies overall
all_results %>% 
  count(case, species)

optimal <- read.csv('output/models/final_models/Choosen_models_marmosini_m1&m2.csv')
optimal <- optimal %>% mutate(model = 'o')

subopti <- read.csv('output/models/final_subopt_models/Choosen_models_marmosini_m1&m2.csv')
subopti <- subopti %>% mutate(model= 's')

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

# Count predictors scenarios frequencies models chosen
goodmodels %>% count(case) %>% arrange(n)

# Count cross-validation frecuencies models chosen
goodmodels %>% count(cross.validation)

# Range of training AUC
range(goodmodels$train.AUC)

#Write results
write.csv(goodmodels, 'output/final_models_chosen.csv', row.names = F)


#### PREDICTORS IMPORTANCE IN MODELS ####
tbl_final <- read.csv('output/final_models_chosen.csv')

#Prepare data.frames to fill with data
dt_conper <- data.frame(species=NULL)
dt_toplot <- dt_conper

#Extract the data for each species from optimal or supoptimal directory
for(i in 1:NROW(tbl_final)){
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

#Which rows are contribution type
contri <- c(1, grep('contribution', colnames(dt_conper)))

dt_conper[,-1] <- as.numeric(as.matrix(dt_conper[,-1]))

#Extract contribution rows
cont <- dt_conper[ , contri] 

col.max <- apply(cont, 1, which.min) 
colnames(cont)[sort(col.max)]

##### Fig 3 plot #####
dir.create('results/Fig3/')
gg.impt <- dt_toplot %>% 
  mutate(type = ifelse(stringr::str_detect(variables, 'contribution'), 
                       'contribution', 'permutation')) %>% 
  mutate(species = gsub('_', ' ', species)) %>% 
  mutate(species = as.factor(species)) %>% 
  mutate(variables = stringr::str_replace(variables, '.permutation.importance|.contribution', '')) %>% 
  arrange(result) %>% 
  mutate(variables = factor(variables, levels = c(paste0('bio_', c(4,6,10,11,15,16,17)), 
                                                  'topoWet','tri', 'msavi'))) %>% 
  ggplot(aes(y=result, fill=variables, x=species), shape='black')+
  scale_fill_manual(name='Variables',
                    labels=c(paste0('bio_', c(4,6,10,11,15,16,17)), 
                             'topoWet', 'tri', 'msavi'),
                    values = c(
                      RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(7,6,5,4)], 
                      RColorBrewer::brewer.pal(n = 9, name = "Blues")[4:7], 
                      'hotpink', 'forest green'
                    ))+
  scale_x_discrete(limits = rev)+
  geom_bar(stat = 'identity', na.rm = F) +
  facet_wrap(~type)+
  ylab('%')+ xlab('')+
  theme_light()+
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = 'italic')
  )+
  coord_flip()

#Write plot and table
ggsave(gg.impt, 'results/Fig3/fig3.tiff', dpi=300, units = 'cm')
write.csv(dt_conper, 'output/contribution_permutation_marmosini.csv', row.names = F)

#### CONSERVATION AND PRESSURE ####
dt_con <- read.csv('output/species_overlap_pressure_conservation.csv')

##### Conservation generalities #### 
dt_con %>% 
  group_by(species) %>% 
  mutate(total_con_sp  = sum(strict.conservation_total_area, managed.resources_total_area,
                             na.rm=T)) %>% 
  mutate(percent_total_con = (total_con_sp / total_area) * 100) %>% 
  select(species, percent_total_con) %>% ungroup() %>% 
  mutate(median_val = median(percent_total_con)) %>% 
  mutate(overall_wo_con = 100 - median_val) %>% 
  arrange(percent_total_con) 

##### Pressure generalities #### 
dt_con %>% 
  group_by(species) %>% 
  mutate(percent_pressure_high = (high_pressure / total_area) * 100) %>% 
  mutate(percent_pressure_low = (low_pressure / total_area) * 100) %>% 
  mutate(percent_pressure_nd = (no_data_by_mask / total_area) * 100) %>% 
  select(species, percent_pressure_high, percent_pressure_low, 
         percent_pressure_nd) %>% 
  arrange(percent_pressure_nd)

##### Type of conservation #### 
  ###### IUCN #####
dt_con %>% 
  mutate(percent_strict  = (strict.conservation_total_area / total_area) * 100) %>% 
  mutate(percent_managed = (managed.resources_total_area / total_area) * 100) %>% 
  select(species, percent_strict, percent_managed) %>% 
  summarise(strict_med = median(percent_strict, na.rm=T), 
            managed_med = median(percent_managed, na.rm=T))

  ###### Governance #####
dt_con %>% 
  mutate(percent_national = (national_total_area / total_area) * 100) %>% 
  mutate(percent_subnatio = (sub.national_total_area / total_area) * 100) %>% 
  mutate(percent_private  = (private_total_area / total_area) * 100) %>% 
  select(species, percent_national, percent_subnatio, percent_private) %>% 
  summarise(national_med = median(percent_national, na.rm=T), 
            subnational_med = median(percent_subnatio, na.rm=T), 
            private_med = median(percent_private, na.rm=T))

##### Conservation pressure ####
dt_con %>% 
  group_by(species) %>% 
  mutate(
    percent_high_strict = 
      (strict.conservation_high_pres_area / strict.conservation_total_area) * 100, 
    percent_high_managed = 
  (managed.resources_high_pres_area / managed.resources_total_area) * 100,
    percent_high_national = 
      (national_high_pres_area / national_total_area) * 100,
    percent_high_subnational = 
      (sub.national_high_pres_area / sub.national_total_area) * 100,
    percent_high_private = 
      (private_high_pres_area/ private_total_area) * 100
  ) %>% #ungroup() %>%
  select(contains('percent')) %>%  ungroup() %>%
  apply(., 2, function(i) {median(i, na.rm = T)})

##### Table 2 ####
dt_con %>% 
  mutate('Range area' = round(total_area, 2),
         Strict  = round(strict.conservation_total_area, 2), 
         Managed = round(managed.resources_total_area, 2), 
         High    = round(high_pressure, 2), 
         Low     = round(low_pressure, 2)) %>% 
  select(species, 'Range area', Strict, Managed, High, Low)

##### Table 3 ####
dt_con %>% 
  mutate(Species      = gsub('_', ' ', species),
         Strict       = round(strict.conservation_total_area / total_area, 4) * 100, 
         Managed      = round(managed.resources_total_area / total_area, 4) * 100, 
         National     = round(national_total_area / total_area, 4) * 100, 
         'Sub-national' = round(sub.national_total_area / total_area, 4) * 100,
         Private      = round(private_total_area / total_area, 4) * 100) %>% 
  select(Species, Strict, Managed, National, 'Sub-national', Private)

##### Table 4 ####
dt_con %>% 
  mutate(Species      = gsub('_', ' ', species),
         Strict.h       = round(strict.conservation_high_pres_area / strict.conservation_total_area, 4) * 100, 
         Managed.h      = round(managed.resources_high_pres_area / managed.resources_total_area, 4) * 100, 
         National.h     = round(national_high_pres_area / national_total_area, 4) * 100, 
         'Sub-national.h' = round(sub.national_high_pres_area / sub.national_total_area, 4) * 100,
         Private.h      = round(private_high_pres_area / private_total_area, 4) * 100, 
         
         Strict.l       = round(strict.conservation_low_pres_area / strict.conservation_total_area, 4) * 100, 
         Managed.l      = round(managed.resources_low_pres_area / managed.resources_total_area, 4) * 100, 
         National.l     = round(national_low_pres_area / national_total_area, 4) * 100, 
         'Sub-national.l' = round(sub.national_low_pres_area / sub.national_total_area, 4) * 100,
         Private.l      = round(private_low_pres_area / private_total_area, 4) * 100, 
         ) %>% 
  mutate(Strict = paste0(Strict.h, ' [', Strict.l, ']'), 
         Managed = paste0(Managed.h, ' [', Managed.l, ']'),
         National = paste0(National.h, ' [', National.l, ']'), 
         'Sub-national' = paste0(`Sub-national.h`, ' [', `Sub-national.l`, ']'), 
         Private = paste0(Private.h, ' [', Private.l, ']')) %>% 
  select(Species, Strict, Managed, National, 'Sub-national', Private)

#### RANGES AREAS ####
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

#Areas of each species (km2)
sapply(ranges, function(x){
  paste0(x[1,1], ' ', 
         st_area(st_union(x)) * 1e-06, 
         '\n')
})

#Proportion of in-country distribution 
(st_area(joint) * 1e-06) / (st_area(COL) * 1e-06 ) * 100

#### SUPPLEMENTARY INFORMATION ####
###### Fig S2 of supplementary information S2 #####
m.area <- all_results %>% 
  mutate(method = str_c(area, '_', cross.validation, '_', case)) %>%
  group_by(area, case, rm, cross.validation) %>% 
  filter(avg.test.AUC >= quantile(avg.test.AUC)[4]) %>% 
  mutate(avg_TAUC = mean(train.AUC), 
         avg_TestAUC = mean(avg.test.AUC), 
         avg_orMTP = mean(avg.test.orMTP),
         avg_AIC = mean(AICc)) %>% 
  distinct()

dir.create('results/FigS2' -> dir.save)

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



###### Figures of supplementary information S3 #####
dir.create('results/FigS3' -> dir.save)

OCCS <- lapply(list.files('records', pattern = '.csv$', full.names = T), read.csv)
names(OCCS) <- gsub(".csv", "", list.files('records', pattern = '.csv$', full.names = F))

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

stopifnot(all.equal(names(OCCS), names(ranges.def), names(ranges.fin)))

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
  ggsave(filename = paste0('/results/FigS3/', names(OCCS[i]), '.png'),
         plot = gg,  dpi = 300
  )
}

#### MAKE DATASET FOR ZENODO ####

dir.shapes <- '/output/definitive_rages' #these are the manually modified ranges

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

write_sf(shapes, 'results/marmosini_ranges_Colombia_WGS84_20210526.geojson')
