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

AllData4 <- merge(AllData3,alledata_abiotiske)

saveRDS(AllData4, "AllData4.rds")

Habs <- unique(AllData4$habtype)
IDs <- character()
Final <- list()
for(i in 1:500){
  if(i == 1){
    Temp <- AllData4
    Selected <- Temp %>% sample_n(size = 1)
    IDs <- pull(Selected, ID)
    Temp <- Temp %>% dplyr::filter(!(ID %in% IDs)) 
    Habs <- Habs[!(Habs %in% Selected$habtype)]
  }
  if(i > 1 & length(Habs) > 0){
    
    IDs <- c(IDs, pull(Selected, ID)) %>% unique()
    Temp <- Temp %>% dplyr::filter(habtype %in% Habs,!(ID %in% IDs))
    Selected <- rbind(Selected,sample_n(Temp, size = 1))
    Habs <- Habs[!(Habs %in% Selected$habtype)]
    if(length(Habs) == 0){
      next
    }
  }
  if(i > 1 & length(Habs) == 0){
    IDs <- c(IDs, pull(Selected, ID)) %>% unique()
    Temp <- AllData4 %>% dplyr::filter(!(ID %in% IDs))
    Habs <- unique(Temp$habtype)
    Final[[length(Final) + 1]] <- Selected
    Selected <- Temp %>% sample_n(size = 1)
    Habs <- Habs[!(Habs %in% Selected$habtype)]
    print(paste("Set!!!", Sys.time()))
    IDs <- c(IDs, pull(Selected, ID)) %>% unique()
  }
  message(paste(i, length(IDs)))
}

Final[[length(Final) + 1]] <- Selected
Final <- Final %>%  purrr::reduce(rbind)

Final$ID %>% length()

Final$habtype %>% table()
