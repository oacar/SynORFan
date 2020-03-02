library(tidyverse)
library(glue)
library(progress)
library(RMySQL)

df <- data.frame()

fls = list.files('data/python_analysis/gene_all_1009/')
count=0
pb <- progress_bar$new(total = length(fls))


for(i in 1:length(fls)){
#  pb$tick()
  fname = fls[i]
  orf_name = str_split(fname,'_')[[1]][1]
  path = glue('data/python_analysis/gene_all_1009/{fname}/{orf_name}_data.csv')
  if(!file.exists(path)){
    count=count+1
    print(fname)
    #data = data.table::fread(path)
  }
  #df = data.table::rbindlist(list(df, data),use.names = TRUE,fill = TRUE)
}


for(i in 1:length(fls)){
  pb$tick()
  fname = fls[i]
  orf_name = str_split(fname,'_')[[1]][1]
  path = glue('data/python_analysis/gene_all_1009/{fname}/{orf_name}_data.csv')
  if(file.exists(path)){
    #count=count+1
    data = data.table::fread(path)
  }
  df = data.table::rbindlist(list(df, data),use.names = TRUE,fill = TRUE)
}
pb$terminate()
df_2 <- data.frame()
orf_names_list = unique(df$orf_name)
pb <- progress_bar$new(total = length(orf_names_list))
df_distinct = df %>% distinct()
for(i in 1:length(orf_names_list)){
  pb$tick()
  o_name = orf_names_list[i]
  df_sub_ = df_distinct %>% filter(orf_name==o_name)
  df_2_2 = df_sub_[1,] %>% 
    map(~.x) %>%
    discard(~all(is.na(.x))) %>%
    map_df(~.x) 
  if(nrow(df_sub_)>1){
    for( i_ in 2:nrow(df_sub_)){ 
      df_2_2 = df_sub_[i_,] %>%  map(~.x) %>%
        discard(~all(is.na(.x))) %>%
        map_df(~.x) %>% inner_join(df_2_2,by=c('V1', 'orf_name','orf_length'))
    }  
  }
  
  df_2 = data.table::rbindlist(list(df_2,df_2_2),use.names = TRUE,fill=TRUE)
}



df_ <- df_2 %>% 
  #data.table::fread('~/pgsNetwork/analysis/data/derived_data/synal_7spec_python_0923_comb.csv') %>%
  mutate(orf_len=orf_length/3) %>%
  select(-orf_length) %>%
  distinct(orf_name,.keep_all = T) %>%
  filter(str_starts(orf_name,'Y'))

extract_best <- function(df, species) {
  # common_str = paste0(species,'_common_aa')
  common= rlang::sym(paste0(quo_name(enquo(species)),'_common_aa'))
  length = rlang::sym(paste0(quo_name(enquo(species)),'_length'))
  df_sub = df %>%
    select(orf_name,orf_len,contains(quo_name(enquo(species)))) %>%
    pivot_longer(cols=c(contains('length'),contains('common')),names_to=c('colnames','ids'),names_pattern = '(.*)_(.*)') %>%
    pivot_wider(names_from=colnames) %>%
    mutate(!! common := replace_na(!! common,0), !! (length) := replace_na(!! (length),0)) %>%
    group_by(orf_name) %>%
    # #drop_na() %>%
    filter(!! common == max(!! common ,na.rm=T))%>%
    filter(!!length==max(!!length,na.rm=T)) %>%
    ungroup() %>%
    distinct(orf_name,.keep_all = TRUE)
  
  df_sub
}


species = c('Spar','Smik','Skud','Seub','Suva','Sarb','Sjur')
empty_df = tibble(orf_name = df_$orf_name,orf_len=df_$orf_len)
for(i in 1:7){
  spec = species[i]
  df_mut = df_ %>% extract_best(!! spec)  %>% mutate(!! (paste0(spec,'_orf_exists')) := ifelse(is.na(.[[5]])==F&.[[7]]!=0,TRUE,FALSE))
  #%>% select(-contains('ids'))
  colnames(df_mut)[colnames(df_mut)=='ids'] <- str_glue(spec,'_ids')
  empty_df = full_join(empty_df,df_mut,by=c('orf_name'='orf_name','orf_len'='orf_len'))
}

empty_df
mydb = dbConnect(MySQL(), user='oma21', password='dktgp2750', dbname='omer',host='127.0.0.1')
#db_list_tables(mydb)
#var='orf_name'
#table = 'orf'
#q='YBR196C-A'
#sql = glue("SELECT * from orf where {var}='{q}'")
#query <- DBI::dbSendQuery(mydb, sql)
#DBI::dbBind(query)

#fetch(query)
pb = progress_bar$new(total = nrow(empty_df))

for(i in 1:nrow(empty_df)){
  pb$tick()
  orf_name = empty_df$orf_name[i]
  if(empty_df[i,]$Spar_orf_exists){
    spar_id = empty_df[i,]$Spar_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Spar' and id!={spar_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Smik_orf_exists){
    smik_id = empty_df[i,]$Smik_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Smik' and id!={smik_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Skud_orf_exists){
    skud_id = empty_df[i,]$Skud_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Skud' and id!={skud_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Seub_orf_exists){
    seub_id = empty_df[i,]$Seub_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Seub' and id!={seub_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Suva_orf_exists){
    suva_id = empty_df[i,]$Suva_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Suva' and id!={suva_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Sarb_orf_exists){
    sarb_id = empty_df[i,]$Sarb_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Sarb' and id!={sarb_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
  
  if(empty_df[i,]$Sjur_orf_exists){
    sjur_id = empty_df[i,]$Sjur_ids
    #DELETE FROM `table_name` [WHERE condition];
    #tryCatch()
    sql = glue("DELETE from pairwise where orf_name='{orf_name}' and species2='Sjur' and id!={sjur_id}")
    query <- DBI::dbSendQuery(mydb, sql)
    dbClearResult(query)
  }
}

