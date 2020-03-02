library(tidyverse)
library(glue)
library(progress)
library(RMySQL)


pairwise_df = data.table::fread('pairwise.csv')
empty_df = data.table::fread('data/df_combined_gene_all_1009_bestids.tsv')
pairwise_df_best = data.table::data.table()
pairwise_df_best_ = pairwise_df_best
pairwise_df_best <- pairwise_df_best %>% filter(species2!='Spar')
pr = progress_bar$new(total = nrow(empty_df))
for(i in 1:nrow(empty_df)){
  o_name = empty_df$orf_name[i]
  
  pr$tick()
  #pairwise_df[orf_name==orf_name]
  if(empty_df[i,]$Spar_orf_exists){
    spar_id = empty_df[i,]$Spar_ids
    pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Spar' & id == spar_id))]))
  }
  
  # if(empty_df[i,]$Smik_orf_exists){
  #   smik_id = empty_df[i,]$Smik_ids
  #   pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Smik' & id == smik_id))]))
  #   
  # }
  # 
  # if(empty_df[i,]$Skud_orf_exists){
  #   skud_id = empty_df[i,]$Skud_ids
  #   pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Skud' & id == skud_id))]))
  #   
  # }
  # 
  # if(empty_df[i,]$Seub_orf_exists){
  #   seub_id = empty_df[i,]$Seub_ids
  #   pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Seub' & id == seub_id))]))
  #   
  # }
  # 
  # if(empty_df[i,]$Suva_orf_exists){
  #   suva_id = empty_df[i,]$Suva_ids
  #   
  #   pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Suva' & id == suva_id))]))
  #   
  # }
  
  if(empty_df[i,]$Sarb_orf_exists){
    sarb_id = empty_df[i,]$Sarb_ids
    pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Sarb' & id == sarb_id))]))
    
  }
  
  # if(empty_df[i,]$Sjur_orf_exists){
  #   sjur_id = empty_df[i,]$Sjur_ids
  #   pairwise_df_best = data.table::rbindlist(list(pairwise_df_best,pairwise_df[orf_name==o_name & ((species2 == 'Sjur' & id == sjur_id))]))
  #   
  # }
  
}
#
pairwise_df_best %>% write_delim('pairwise_best_corrected.tsv','\t')

#empty_df %>% as_tibble() %>% mutate_all(as.character) %>% pivot_longer(-orf_name)

empty_df%>% 
  as_tibble() %>% 
  rename(Scer_length=orf_len) %>% 
  mutate_at(
    .vars = vars(contains("common")),
    .funs = list(~ ifelse(is.na(.), NA, round((. / Scer_length),digits = 3)))) %>% 
  mutate_at(
    .vars = vars(contains("dna_identity")),
    .funs = list(~ ifelse(is.na(.), NA, round(.,digits = 3)))) %>%
  mutate_all(as.character) %>% 
  pivot_longer(-orf_name) %>% 
  separate(name,sep='_',into = c('species','data'),extra='merge') %>% 
  mutate(data=replace(data, data=='common_aa', 'aa_identity')) %>%
  pivot_wider(id_cols = c('orf_name','species','data'),names_from='species',values_from = 'value') %>% 
  write_delim('data/df_combined_gene_all_1009_bestids_long_form.tsv',delim='\t')
