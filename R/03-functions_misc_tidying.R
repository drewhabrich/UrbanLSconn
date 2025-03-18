## HEADER---------------------------
## Script name: misc functions for targets pipeline
##
## Purpose of script: A collection of functional tasks for 'targets' package
##
## Author: Andrew Habrich
##
## Date Created: 2025-01-22
## Date last Modified: 2025-01-22
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
## Only need to define functions here, no need to run them or require packages to be loaded

## Summarize counts of checklists and species per city
citychkl_sumtab <- function(city_zfdf) {
  ## read in the data
  df <- readRDS(city_zfdf$rds_tidyzf[1]) 
  ## remove rows where species observed is FALSE
  df <- df %>% filter(species_observed == TRUE)
  ## create a new dataframe to output data
  citysummary <- tibble(cityname = city_zfdf$cityname[1], 
                        n_chkl = n_distinct(df$checklist_id),
                        n_spp = n_distinct(df$scientific_name))
  return(citysummary)
}

## Create a summary dataframe to references file locations to be used in the pipeline
get_df_divLC <- function(divdata, chkldata, LCrasterfiles, BVrasterfiles) {
  # Join with the div dataframe with the LC files by cityname
  bycity_divlc <- divdata %>% left_join(LCrasterfiles, by = "cityname") %>% rename(LCpath = filepath)
  # Join with the BV dataframe by 
  bycity_divlc <- bycity_divlc %>% left_join(BVrasterfiles, by = "cityname") %>% rename(BVpath = filepath)
  # Join with the checklist data by cityname
  bycity_filessummary <- bycity_divlc %>% left_join(chkldata, by = "cityname") %>%
    select(!file) %>% 
    select(!geometry)
  return(bycity_filessummary)
}

## save results of lsm to csv; takes results from multi-buffer lists and saves them to csv
# save_results_csv <- function(results_list, city_name, outputdir) {
#   # Get list names
#   list_names <- names(results_list) 
#   # Save each tibble
#   for (name in list_names) {
#     file_name <- paste0(outputdir, city_name, "_", name, ".csv")
#     write_csv(x = results_list[[name]], file = file_name)
#   }
#   #message("Saved ", length(results_list), " files to: ", normalizePath(outputdir))
#   #plandfiles <- tibble(cityname = city_name, file = list.files(outputdir))
#   return(list.files(outputdir))
# }

save_results_csv <- function(results_list, city_name, outputdir, whatmetric) {
  #save the file to csv in the outputdir
  write_csv(x = results_list, file = paste0(outputdir, city_name, "_", whatmetric,"_results.csv"))
  
  #check what files are in the output directory and return the list
  plandfiles <- tibble(cityname = city_name, file = list.files(outputdir))
  return(list.files(outputdir))
}