## HEADER---------------------------
## Script name: 04-ebird_data
## 
## Purpose of script: Filter ebird data and sampling information to relevant regions of interest. Using 'auk_filter' to reduce dataset to research grade checklists for modelling urban bird communities.
##
## Author: Andrew Habrich
##
## Date Created: 2023-08-21
## Date last modified: 2023-11-23
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

## 1. Load relevant packages--------
library(auk)
library(sf)
library(tidyverse)
library(lubridate)
library(rnaturalearth)

### 1.1. Define path to ebird data and sampling info ----
#auk::auk_set_ebd_path(path = "F:/DATA/ebird/") ## Run once and restart, should remember between sessions
auk::auk_get_ebd_path()
### create the ebd object to keep track of defined filters and input data
ebd <- auk_ebd(file = "ebd_relJul-2023.txt", #this is the .txt with the data
             file_sampling = "ebd_sampling_relJul-2023.txt") #this is the .txt with info on the sampling

### 1.2 check filter options for the ebird dataset ----
### Define filters for the ebd (MORE filters, means smaller, more work-able dataset)
## FIRST; what countries and administrative zones can we filter by?
ebird_states %>% filter(country_code %in% c('CA')) %>% head(10)
valid_protocols #what protcols can we filter by?

## SECOND; define the filters to get relevant research quality data
### See options here; https://cornelllabofornithology.github.io/auk/articles/auk.html

### 1.3 Filter to just Canada data (should be smaller to extract data from) ----
auk_filt_canada <- ebd %>% 
  auk_country(country = "CA") %>% 
  auk_protocol(protocol = c("Stationary", "Traveling", "Area")) %>% #stationary, transect, or area-based sampling
  auk_duration(duration = c(5, 240)) %>% 
  auk_distance(distance = c(2, 15), distance_units = "km") %>% 
  auk_date(date = c("*-03-01", "*-08-31")) %>% #breeding bird season in CANADA
  auk_year(year = c(2010:2022)) %>% #last 12 years of data
  auk_complete() #ONLY completed checklists

start_time <- Sys.time()
if (!file.exists("./raw_data/ebd_can_obsdata.txt")) {
auk_filter(auk_filt_canada, 
           file = "./raw_data/ebd_can_obsdata.txt", 
           file_sampling = "./raw_data/ebd_can_sampdata.txt", 
           filter_sampling=T, overwrite = F)
}
end_time <- Sys.time()
end_time - start_time ##Time difference of 1.734443 hours

## import dataset and take a look; THIS DATA IS FED INTO THE TARGETS PIPELINE
# co <- read_ebd("./raw_data/ebd_can_obsdata.txt")
# cs <- read_sampling("./raw_data/ebd_can_sampdata.txt")

## 2.0 Filter ebird by urban area polygons in sf object ----
cancities <- read_sf("./data/01-can_cities_5kmbuffer.shp")
## Create a list of strings to name the files after; ebd_CITYNAME_5kmbuff.txt
ebirdcity_out <- data.frame(checklist = str_c("./raw_data/ebd_", cancities$PCNAME_x, "_5kmbuff.txt"),
                            sampling = str_c("./raw_data/ebdsamp_", cancities$PCNAME_x, "_5kmbuff.txt"))

## Define the reduced ebd dataset as the new dataset to run the filters on
ebd_can <- auk_ebd(file = "./raw_data/ebd_can_obsdata.txt", #this is the .txt with the data
                   file_sampling = "./raw_data/ebd_can_sampdata.txt") 

## generate a list of filters to run a loop on
city_filters <- list()
for (i in 1:nrow(cancities)) {
  # define the filter using the polygon of the city
  filt <- ebd_can %>% 
    auk_bbox(cancities %>% slice(i)) %>% 
    auk_complete()
  # append to the list
  city_filters[[i]] <- filt
}

## Export for EACH city
for (i in 1:length(city_filters)) {
  auk_filter(city_filters[[i]], file = ebirdcity_out[i,1], file_sampling = ebirdcity_out[i,2], filter_sampling=T, overwrite = T)
}

## 3.0. Zero-filling presence/absence checklists ----
## Clear the R environment
rm(list=ls())
## read in the filtered ebird subset
path <- "./raw_data/"
ebd_on <- read_ebd("./raw_data/ebd_ontario_obsdata.txt", unique=T) 
ebdsamp_on <- read_sampling("./raw_data/ebd_ontario_sampdata.txt")

## take a peek at the dataset
glimpse(ebd_on)

### Zero-filling and tidying
zf <- auk_zerofill(ebd_on, ebdsamp_on, collapse=T)
# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}
# clean up variables
zf <- zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )
# additional filtering
zf_filtered <- zf %>% 
  filter(effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)

### 3.1 Test-training split to evaluate models ----
zf_split <- zf_filtered %>% 
  mutate(type = if_else(runif(nrow(.)) <= 0.8, "train", "test"))
# confirm the proportion in each set is correct
table(zf_split$type) / nrow(zf_split)

# remove some redundant columns and save to .csv for analysis
checklists <- zf_split %>% 
  select(checklist_id, observer_id, type,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         hours_of_day, 
         effort_hours, effort_distance_km, effort_speed_kmph,
         number_observers)
write_csv(checklists, "./output/04-checklists-zf_can-on.csv", na = "")

### 3.2 Explore the zero-filled data
## setup map for Visualize
map_proj <- st_crs("EPSG:4326")

ne_land <- read_sf("./data/01-NA_gis_data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("./data/01-NA_gis_data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("./data/01-NA_gis_data.gpkg", "ne_adm1_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
on_boundary <- read_sf("./data/01-NA_gis_data.gpkg", "ne_adm1") %>% 
  filter(adm1_code == "CA-ON") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

## read in zero-filled checklists
mtl <- read_ebd("./raw_data/ebd_Montréal_5kmbuff.txt", unique=T)
mtl_s<- read_sampling("./raw_data/ebdsamp_Montréal_5kmbuff.txt")

mtl_sf <- mtl %>% 
  select(longitude, latitude) %>% 
  st_as_sf( coords = c("longitude", "latitude"), crs = 4326)
mtl_s_sf <- mtl_s %>% 
  select(longitude, latitude) %>% 
  st_as_sf( coords = c("longitude", "latitude"), crs = 4326)

## find points inside the polygon buffer
mtla <- st_within(mtl_sf, cancities %>% slice(16), sparse = F)
mtlb <- st_within(mtl_s_sf, cancities %>% slice(16), sparse = F)

# subset data frame
mtl_checklists <- mtl[mtla[, 1], ]
mtl_sampling <- mtl_s[mtlb[, 1], ]

# Make a quick map
par(mar = c(0, 0, 0, 0))
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
plot(cancities %>% slice(16) %>% st_geometry(), col = "grey40", border = NA, add = TRUE)
plot(mtl_sf, col = "black", pch = 19, cex = 0.5, add = TRUE)
plot(mtl_sf[mtla[, 1], ], 
     col = "forestgreen", pch = 19, cex = 0.5, 
     add = TRUE)

mtl_zf <- auk_zerofill(mtl_checklists, mtl_sampling, collapse=T)

