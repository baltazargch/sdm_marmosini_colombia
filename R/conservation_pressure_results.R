library(tidyverse)

dt_con <- read.csv('output/species_overlap_pressure_conservation.csv')

head(dt_con)
colnames(dt_con)

#### Conservation generalities #### 
dt_con %>% 
  group_by(species) %>% 
  mutate(total_con_sp  = sum(strict.conservation_total_area, managed.resources_total_area,
                             na.rm=T)) %>% 
  mutate(percent_total_con = (total_con_sp / total_area) * 100) %>% 
  select(species, percent_total_con) %>% ungroup() %>% 
  mutate(median_val = median(percent_total_con)) %>% 
  mutate(overall_wo_con = 100 - median_val) %>% 
  arrange(percent_total_con) 

#### Pressure generalities #### 
dt_con %>% 
  group_by(species) %>% 
  mutate(percent_pressure_high = (high_pressure / total_area) * 100) %>% 
  mutate(percent_pressure_low = (low_pressure / total_area) * 100) %>% 
  mutate(percent_pressure_nd = (no_data_by_mask / total_area) * 100) %>% 
  select(species, percent_pressure_high, percent_pressure_low, 
         percent_pressure_nd) %>% 
  arrange(percent_pressure_nd)

#### Type of conservation #### 
  ##### IUCN #####
dt_con %>% 
  mutate(percent_strict  = (strict.conservation_total_area / total_area) * 100) %>% 
  mutate(percent_managed = (managed.resources_total_area / total_area) * 100) %>% 
  select(species, percent_strict, percent_managed) %>% 
  summarise(strict_med = median(percent_strict, na.rm=T), 
            managed_med = median(percent_managed, na.rm=T))

  ##### Governance #####
dt_con %>% 
  mutate(percent_national = (national_total_area / total_area) * 100) %>% 
  mutate(percent_subnatio = (sub.national_total_area / total_area) * 100) %>% 
  mutate(percent_private  = (private_total_area / total_area) * 100) %>% 
  select(species, percent_national, percent_subnatio, percent_private) %>% 
  summarise(national_med = median(percent_national, na.rm=T), 
            subnational_med = median(percent_subnatio, na.rm=T), 
            private_med = median(percent_private, na.rm=T))

#### Conservation pressure ####
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

#### Table 2 ####
dt_con %>% 
  mutate('Range area' = round(total_area, 2),
         Strict  = round(strict.conservation_total_area, 2), 
         Managed = round(managed.resources_total_area, 2), 
         High    = round(high_pressure, 2), 
         Low     = round(low_pressure, 2)) %>% 
  select(species, 'Range area', Strict, Managed, High, Low)

#### Table 3 ####
dt_con %>% 
  mutate(Species      = gsub('_', ' ', species),
         Strict       = round(strict.conservation_total_area / total_area, 4) * 100, 
         Managed      = round(managed.resources_total_area / total_area, 4) * 100, 
         National     = round(national_total_area / total_area, 4) * 100, 
         'Sub-national' = round(sub.national_total_area / total_area, 4) * 100,
         Private      = round(private_total_area / total_area, 4) * 100) %>% 
  select(Species, Strict, Managed, National, 'Sub-national', Private)

#### Table 4 ####
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

