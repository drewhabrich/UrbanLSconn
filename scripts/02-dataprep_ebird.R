## HEADER---------------------------
## Script name: dataprep_ebird.R
##
## Purpose of script: Extract the high quality, relevant ebird data for the selected urban centres
##
## Author: Andrew Habrich
##
## Date Created: 2024-08-15
## Date last Modified: 2024-08-19
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

## 1. Load relevant packages--------
#install.packages("pacman")
pacman::p_load(tidyverse, ggpubr, RColorBrewer, 
               ggspatial, sf, terra, rnaturalearth, mapview,
               auk, targets, tarchetypes)

## 2. Auk settings preparation ----
#auk::auk_set_ebd_path(path = "E:/DATA/ebird/") ## Run once and restart, should remember between sessions
auk::auk_get_ebd_path()

## Full dataset until jul 2023
ebd <- auk_ebd(file = "ebd_relJul-2023.txt", #this is the .txt with the data
               file_sampling = "ebd_sampling_relJul-2023.txt") #this is the .txt with info on the sampling

## Canada dataset from jul 2023 to jul 2024
ebd_can2024 <- auk_ebd(file = "ebd_CA_relJul-2024.txt",
                       file_sampling = "ebd_CA_relJul-2024_sampling.txt")

### 2.1 check filter options for the ebird dataset ----
# You can find more description for filter options here: auk_filter_options()

### 2.2 Filter the ebd dataset ----
### Define filters for the ebd (MORE filters, means smaller, more work-able dataset)
auk_filt_canada <- ebd %>% 
  auk_country(country = "CA") %>% 
  auk_protocol(protocol = c("Stationary", "Traveling", "Area")) %>% #stationary, transect, or area-based sampling
  auk_duration(duration = c(5, 240)) %>% 
  auk_distance(distance = c(2, 15), distance_units = "km") %>% 
  auk_date(date = c("*-03-01", "*-08-31")) %>% #breeding bird season in CANADA
  auk_year(year = c(2010:2022)) %>% #last 12 years of data
  auk_complete() #ONLY completed checklists

auk_filt_canada2024 <- ebd_can2024 %>% 
  auk_protocol(protocol = c("Stationary", "Traveling", "Area")) %>% #stationary, transect, or area-based sampling
  auk_duration(duration = c(5, 240)) %>% 
  auk_distance(distance = c(2, 15), distance_units = "km") %>% 
  auk_date(date = c("*-03-01", "*-08-31")) %>% #breeding bird season in CANADA
  auk_year(year = c(2023:2024)) %>%
  auk_complete() #ONLY completed checklists

### 2.3 Filter the FULL ebd dataset for Canada ----
start_time <- Sys.time()
if (!file.exists("./raw_data/ebd_can_obsdata.txt")) {
  auk_filter(auk_filt_canada, 
             file = "./raw_data/ebd_can_obsdata.txt", 
             file_sampling = "./raw_data/ebd_can_sampdata.txt", 
             filter_sampling = T, overwrite = F)
}
end_time <- Sys.time()
end_time - start_time ##Time difference of 1.734443 hours

### filter the newest Canadian data
start_time <- Sys.time()
if (!file.exists("./raw_data/ebd_can_obsdata_2024.txt")) {
  auk_filter(auk_filt_canada2024, 
             file = "./raw_data/ebd_can_obsdata_2024.txt", 
             file_sampling = "./raw_data/ebd_can_sampdata_2024.txt", 
             filter_sampling = T, overwrite = F)
}
end_time <- Sys.time()
end_time - start_time ##Time difference 

## import dataset and take a look; THIS DATA IS FED INTO THE TARGETS PIPELINE
# co <- read_ebd("./raw_data/ebd_can_obsdata.txt")
# str(co)
# cs <- read_sampling("./raw_data/ebd_can_sampdata.txt")
# str(cs)

## 3.0 Filter ebird by urban area BOUNDING BOX in sf object ----
fuacan <- read_sf("./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua")
## Create a list of strings to name the files after; ebd_CITYNAME_fua.txt
ebirdcity_out <- data.frame(checklist = str_c("./raw_data/ebd2022_", fuacan$eFUA_name, "_fua.txt"),
                            sampling = str_c("./raw_data/ebdsamp2022_", fuacan$eFUA_name, "_fua.txt"))

## Define the reduced ebd dataset as the new dataset to run the filters on
ebd_can <- auk_ebd(file = "./raw_data/ebd_can_obsdata.txt", #this is the .txt with the data
                   file_sampling = "./raw_data/ebd_can_sampdata.txt") 

### 3.1 filter to each city polygon ----
## generate a list of filters to run a loop on
city_filters <- list()
for (i in 1:nrow(fuacan)) {
  # define the filter using the polygon of the city
  filt <- ebd_can %>% 
    auk_bbox(fuacan %>% slice(i)) %>% 
    auk_complete()
  # append to the list
  city_filters[[i]] <- filt
}

## Export for EACH city
data_list <- list()
for (i in 1:length(city_filters)) {
  auk_filter(city_filters[[i]], 
             file = ebirdcity_out[i,1], 
             file_sampling = ebirdcity_out[i,2], 
             filter_sampling=T, overwrite = T)
}

####################################################################