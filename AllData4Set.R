library(readr)
library(tidyverse)
library(data.table)
library(sf)
AllData3 <- readRDS("~/Documents/Site_Selection_DK/AllData3.rds")

alledata_abiotiske <-  fread(file="Novana/alledata-abiotiske.csv", na.strings = "mv") %>% as.data.frame() %>% 
  dplyr::select(site, plot, sekhabtype, terhabtype, year) %>% 
  mutate(terhabtype = str_remove_all(terhabtype, "\\}"),
         terhabtype = str_remove_all(terhabtype, "\\{"),
         terhabtype = as.integer(terhabtype),
         equal = ifelse(terhabtype == sekhabtype, "Yes", "No")) %>% 
  distinct() %>% 
  dplyr::filter(!is.na(terhabtype) & !is.na(sekhabtype)) %>% 
  unite(col = "ID", site, plot) %>% 
  group_split(ID) %>% 
  purrr::map(~dplyr::filter(.x, year == max(year))) %>% 
  purrr::reduce(bind_rows) %>% 
  dplyr::select(ID, terhabtype) %>% 
  rename(habtype = terhabtype) %>% 
  dplyr::filter(ID %in% AllData3$ID)

AllData4 <- left_join(AllData3,alledata_abiotiske)

saveRDS(AllData4, "AllData4.rds")

