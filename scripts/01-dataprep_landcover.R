## HEADER---------------------------
## Script name: dataprep_landcover.R
##
## Purpose of script: match the ESA landcover data to the selected urban centres for extraction
##
## Author: Andrew Habrich
##
## Date Created: 2024-08-15
## Date last Modified: 2024-08-15
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
               targets, tarchetypes)
terraOptions(progress = 10, verbose = TRUE)
## 2. landcover data from Copernicus, European space agency ----
# Source: ESA WorldCover 10 m 2021 v200, downloadable here: https://zenodo.org/records/7254221

### 2.1. Identify ESA cells that overlap with urban centres ----
fuacan <- read_sf("./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua")
esagrid <- read_sf("./raw_data/esa_worldcover_grid_composites.fgb")

# identify the cells that overlap with the urban centres
esalist_fua <- st_intersects(esagrid, fuacan, sparse = F) %>% 
  #extract all elements that are not 0
  apply(., 1, any) %>% 
  #coerce to dataframe and rename the . column to overlap
  as_tibble() %>%
  rename(urbanoverlap = value) %>%
  #add the tile name to the dataframe
  mutate(tile = esagrid$tile) %>%
  #filter to only the cells that overlap
  filter(urbanoverlap == T)

## A couple of steps required to import the raster files that match the macrotiles
### 2.2. identify macrotiles from esalist_fua ----
## separate the tile column into N and W extracting only the numbers
esalist_fua <- esalist_fua %>% 
  separate(tile, into = c("N", "W"), sep = c(3, 7), remove = F) %>% 
  mutate(N = as.numeric(str_extract(N, "\\d{2}")),
         W = as.numeric(str_extract(W, "\\d{3}"))) %>% 
  # Round the N value up and the W value down to the closest macrotile, which are 3x3
  mutate(Nmacrotile = floor(N/3)*3, 
         Wmacrotile = ceiling(W/3)*3) %>% 
  #create a new column for the macrotile name
  mutate(macrotile = paste0("N", str_pad(Nmacrotile, 2, pad = "0"), "W", 
                            str_pad(Wmacrotile, 3, pad = "0")))

## 3.0 Crop layers to urban areas ----
### 3.1 Identify a list of the relevant macrotile rasters and import ----
datadir <- "E:/DATA/ESA2021LANDCOVER/"
esa_tilelist <- list.files(datadir) %>% 
  str_subset("tif$") %>% #only the raster images
  as_tibble() %>% rename(filename = value) %>% 
  #extract the tile name from the filename, and add it to the dataframe
  mutate(macrotile = str_extract(filename, "(N\\d{2}W\\d{3})")) %>% 
  #filter to only the tiles that overlap with urban centres
  filter(macrotile %in% esalist_fua$macrotile)

# import a list of raster files that match the macrotiles 
esa_rasters <- esa_tilelist %>% 
  mutate(landcover = map(filename, ~rast(paste0(datadir, .x)))) %>% 
  pull(landcover) %>% 
  sprc()

summary(esa_rasters)
# master <- terra::merge(esa_rasters, 
#                        filename = paste0(datadir, "esa_master_lc.tif"),
#                        overwrite = FALSE)

## for each of the cities in the fuacan, crop the landcover rastercollection so that each city has its own raster
TO <- vect(fuacan %>% filter(eFUA_name == "Toronto"))
MTL <- vect(fuacan %>% filter(eFUA_name == "Montreal"))

plot(esa_rasters[4])
plot(TO, add = T)
TOLC <- mask(esa_rasters[4], TO)
plot(TOLC)

###
## load packages
pacman::p_load(tidyverse, ggpubr, ggspatial, geodata, terra, sf, mapview, leaflet, stars)

## check what datasets are available on geodata
countries <- country_codes() #high resolution

## get the country code for Canada
can_adm <- gadm(country = "Canada", version = "latest", resolution = 2, path = "raw_data/") #low resolution
plot(can_adm)

## geodata, ESA landcover
esalc <- geodata::landcover(var = "built", path = "raw_data/")
plot(esalc)

## Crop the landcover data to Canada polygon
esalc_canada <- crop(esalc, can_adm)
plot(esalc_canada)

