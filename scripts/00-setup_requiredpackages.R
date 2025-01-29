## HEADER---------------------------
## Script name: 00-setup_requiredpackages.R
##
## Purpose of script: download and install required packages for the workflow.
## This script is meant to be run once at the beginning of the workflow.
## Author: Andrew Habrich
##
## Date Created: 2024-08-14
## Date last Modified: 2024-08-14
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

## 2. Administrative shapes data from natural earth ----
# Source: https://www.naturalearthdata.com/ 
## political boundaries - land border with lakes removed
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()
# state boundaries for COUNTRY
ne_adm1_canada <- ne_download(scale = 50, category = "cultural",
                       type = "admin_1_states_provinces",
                       returnclass = "sf") %>% 
  filter(adm0_a3 == "CAN") %>%  ## Change the filter here for different countries
  select(adm1 = name, adm1_code = iso_3166_2)
# country lines - global layer filtered to North America with st_intersect()
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
                    st_geometry()
ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%  
  as.logical() %>%  {ne_country_lines[.]}
# borderlines, North America
ne_adm1_lines <- ne_download(scale = 50, category = "cultural",
                             type = "admin_1_states_provinces_lines",
                             returnclass = "sf") %>%
  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
  select(country = ADM0_NAME, country_code = iso_a2)

## write all the layers into a geopackage in the data folder
st_write(ne_land, "./data/00-ne_admin_layers.gpkg", layer = "ne_land", driver = "GPKG")
st_write(ne_adm1_canada, "./data/00-ne_admin_layers.gpkg", layer = "ne_canada", append = TRUE)
st_write(ne_adm1_lines, "./data/00-ne_admin_layers.gpkg", layer = "ne_adm1borders", append = TRUE)
st_write(ne_country_lines, "./data/00-ne_admin_layers.gpkg", layer = "ne_countryborders", append = TRUE)

## 3. Urban centres from GHSL functional urban areas ----
# Source: https://human-settlement.emergency.copernicus.eu/index.php
## read in the ghs fua from the geopackage
ghs_fua <- read_sf(dsn = "./raw_data/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg") %>% 
# filter to Canada
  filter(Cntry_name == "Canada") %>% 
# filter to cities with a population greater than 100,000
  filter(FUA_p_2015 > 100000) %>% 
# convert crs to WG84
  st_transform(4326) %>% 
# create a 1km buffer around each city
  st_make_valid(snap_radius = -1, validate = T) %>% 
  st_buffer(dist = 1000)

## save to file in data folder
st_write(ghs_fua, "./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua", driver = "GPKG")

## read in the urban centres database file
ghs_fua_db <- read_sf("./raw_data/GHS_STAT_UCDB2015MT_GLOBE_R2019A/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg") %>% 
# filter to Canada
  filter(CTR_MN_NM == "Canada") %>% 
  # convert crs to WG84
  st_transform(4326) %>% 
  # create a 1km buffer around each city
  st_buffer(dist = 1000)
## save and append to the gpkg file in data folder
st_write(ghs_fua_db, "./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua_db", append = T)

## 4. Plot the data ----
# plot the cities in canada using mapview with a population greater than 100,000, 
# coloured by population, using Rcolourbrewer Yl-Rd, and
# label with the city name
mapview(ghs_fua, 
        layer.name = "Population size",
        zcol = "FUA_p_2015", 
        col.regions = viridisLite::magma(7), 
        label = "eFUA_name") +
mapview(ghs_fua_db, layer.name = "UrbanCentre Area(km2)", 
        zcol = "AREA", 
        col.regions = viridisLite::viridis(7), 
        label = "UC_NM_MN")