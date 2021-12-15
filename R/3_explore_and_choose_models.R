#### MODEL FILTERING, EVALUATION AND CHOOSNG ####
library(tidyverse)
source('R/2_fit_models.R')

tbl_fl <- 'output/fitting/marmosini_fitting_results.csv'

#Access each species model directory and results and make one table. 
#If new species were added, this script add them automatically to the table.
fit_table <- lapply(list.files('output/fitting/eval.tables/', '.csv', full.names = T), 
                      read_csv) %>% do.call(rbind, .)
write_csv(fit_table, tbl_fl)
#Clean the data: here we exlucded jackknife cross-validation
#since the prediction from this models were considered by the authors as 
#inadequate. 
unique(fit_table$cv)
filtered_table <- fit_table %>% 
  na.omit() %>% 
  mutate(method = str_c(aream, '_', cv, '_', case)) %>% 
  mutate(features = factor(fc, 
                           levels = c('L', 'LQ', 'LQP', 'H', 'LQH', 'LQHP', 'LQHPT')
                           )
         ) %>% 
  filter(cv != 'jackknife') 

#### OPTIMAL MODELS PER SPECIES ####

#We used primarily the auc.val.avg to filter the results. See details in the
#manuscript associated with this script (pending to be published).
df.final <- filtered_table %>% 
  group_by(species, aream) %>% 
  filter(auc.val.avg >= quantile(auc.val.avg)[4]) %>% #main criterium
  filter(auc.diff.avg == min(auc.diff.avg)) %>% 
  filter(or.mtp.avg == min(or.mtp.avg)) %>% 
  filter(AICc == min(AICc)) %>% 
  filter(rm == max(rm)) %>% 
  filter(auc.train == max(auc.train)) %>% 
  filter(ncoef == min(ncoef)) %>% 
  filter(or.10p.avg == min(or.10p.avg)) 

#For the rare cases where more than one model for each area were chose, 
#assign one at random.

df.final <- df.final %>% 
  group_by(species, aream) %>% 
  slice_sample(n = 1)

#Write results
write_csv(df.final, 'output/fitting/Choosen_models_marmosini_m1&m2.csv')

#### SUBOPTIMAL MODELS PER SPECIES ####
optimal <- df.final %>% 
  mutate(id = str_c(species, method, tune.args))

df_subopt <- filtered_table %>% 
  mutate(id = str_c(species, method, tune.args)) %>% 
  filter(!id  %in% optimal$id)

df_subopt <- df_subopt %>% 
  group_by(species, aream) %>% 
  filter(auc.val.avg >= quantile(auc.val.avg)[4]) %>% #main criterium
  filter(auc.diff.avg == min(auc.diff.avg)) %>% 
  filter(or.mtp.avg == min(or.mtp.avg)) %>% 
  filter(AICc == min(AICc)) %>% 
  filter(rm == max(rm)) %>% 
  filter(auc.train == max(auc.train)) %>% 
  filter(ncoef == min(ncoef)) %>% 
  filter(or.10p.avg == min(or.10p.avg)) 

#For the rare cases where more than one model for each area were chose, 
#assign one at random.

df_subopt <- df_subopt %>% 
  group_by(species, aream) %>% 
  slice_sample(n = 1)

#Write results
write_csv(df_subopt, 'output/fitting/Subopt_models_marmosini_m1&m2.csv')

