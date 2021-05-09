library(tidyverse)
library(readxl)


c('PGM', 'ACYP', 'PGI', 'PGK','PYK', 'HEX1', 'DPGase', 'TPI', 'PFK', 'ENO', 'GAPD', 'DPGM', 'FBA', 'G3PD2m' ,
  'ACYP_Neuron', 'DPGM_Neuron', 'DPGase_Neuron', 'ENO_Neuron', 'FBA_Neuron', 'G3PD2m_Neuron', 'GAPD_Neuron', 'HEX1_Neuron', 'PFK_Neuron', 'PGI_Neuron', 'PGK_Neuron', 'PGM_Neuron', 'PYK_Neuron', 'TPI_Neuron' ,
  'ATPS4m_Neuron', 'CYOOm2_Neuron', 'CYOR-u10m_Neuron', 'NADH2-u10m_Neuron', 'PPA_Neuron', 'PPAm_Neuron' ,
  'PPAm', 'ATPS4m', 'CYOOm2', 'CYOR-u10m', 'NADH2-u10m', 'PPA', 'GLCt1r', 'GLCt1r_Neuron', 'GLCt2r') -> subsystems
filter_out_subsystems <- function(df){df[!(df$ID %in% subsystems),]}
get_genes_from_module <- function(hubs){     hubs %>% 
    dplyr::select(matches('gene')) %>% purrr::reduce(c) %>% na.exclude() %>% str_c(collapse = ", ") %>% 
    str_split(pattern=",") %>%str_extract_all('\\d+') %>% unlist() %>% str_trim %>% unique()}
#Functions
get_genes_from_each_hub <- function(file, cut_off_optimality = 0.01 , cut_off_centrality =  0.01){
  
  #central_hubs <- readxl::read_excel(file)%>% dplyr::select(-"Optimality")%>%filter(Centrality > cut_off_centrality)%>% filter_out_subsystems
  #optimal_hubs <- readxl::read_excel(file)%>% dplyr::select(-"Centrality")%>%filter(Optimality > cut_off_optimality)%>% filter_out_subsystems
  
  central_hubs <- readr::read_csv(file)%>% dplyr::select(-"Optimality")%>%filter(Centrality > cut_off_centrality)%>% filter_out_subsystems
  optimal_hubs <- readr::read_csv(file)%>% dplyr::select(-"Centrality")%>%filter(Optimality > cut_off_optimality)%>% filter_out_subsystems
  
  central_optimal_hubs <-  readr::read_csv(file)%>% filter_out_subsystems
  
  intersect(central_hubs$ID, optimal_hubs$ID) -> hyper_IDs
  setdiff(central_hubs$ID,hyper_IDs)          -> central_IDs
  setdiff(optimal_hubs$ID,hyper_IDs)         -> optimal_IDs
  
  central_hubs         %>% filter(ID %in% central_IDs) %>%  get_genes_from_module -> central_hub_genes
  optimal_hubs         %>% filter(ID %in% optimal_IDs) %>%  get_genes_from_module -> optimal_hub_genes
  central_optimal_hubs %>% filter(ID %in% hyper_IDs) %>% get_genes_from_module -> hyper_hub_genes
  
  list_of_genes_by_hubs <- list(optimal_hub_genes,central_hub_genes,hyper_hub_genes) %>% set_names(c('optimal_hub_genes','central_hub_genes','hyper_hub_genes'))
  
  return(list_of_genes_by_hubs)
  
}


get_pure_gene_groups <- function(list_of_genes_by_hubs){
  
  list_of_genes_by_hubs$optimal_hub_genes -> optimal_hub_genes
  list_of_genes_by_hubs$central_hub_genes ->central_hub_genes
  list_of_genes_by_hubs$hyper_hub_genes ->hyper_hub_genes
  
  
  Type_5                <- Reduce(intersect, list(central_hub_genes,optimal_hub_genes,  hyper_hub_genes))
  Type_4_hyper_optimal  <- setdiff( Reduce(intersect, list(optimal_hub_genes,  hyper_hub_genes)), Type_5)
  Typer_4_hyper_central <- setdiff( Reduce(intersect, list(central_hub_genes,  hyper_hub_genes)), Type_5)
  
  Type_3                <- setdiff(hyper_hub_genes, c(Type_5, Type_4_hyper_optimal, Typer_4_hyper_central) )
  
  
  Type_2                <-  setdiff(Reduce(intersect, list(central_hub_genes,optimal_hub_genes)), c(Type_5, Type_4_hyper_optimal, Typer_4_hyper_central, Type_3) )
  Type_1_optimal        <-  setdiff(optimal_hub_genes, c(Type_2, Type_3, Typer_4_hyper_central, Type_4_hyper_optimal, Type_5) )
  Type_1_central        <-  setdiff(central_hub_genes, c(Type_2, Type_3, Typer_4_hyper_central, Type_4_hyper_optimal, Type_5) )
  
  list(Type_1_optimal, Type_1_central, Type_2, Type_3, Type_4_hyper_optimal, Typer_4_hyper_central, Type_5) %>% 
    purrr::set_names('Type_1_optimal', 'Type_1_central', 'Type_2', 'Type_3', 'Type_4_hyper_optimal', 'Typer_4_hyper_central', 'Type_5') -> my_list
  # sanity check
  
  c(optimal_hub_genes, central_hub_genes, hyper_hub_genes) %>% unique() -> all_genes
  all((Reduce(c, my_list) %in% all_genes) & (all_genes %in% Reduce(c, my_list)) ) -> sanity_check
  
  if (sanity_check){message('all classified')}
  
  return(my_list)
}