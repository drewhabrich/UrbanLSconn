## HEADER---------------------------
## Script name: 01-urban_landscape_extraction
##
## Purpose of script: To extract extents and land cover from urban areas in Canada using census data
##
## Author: Andrew Habrich
##
## Date Created: 2023-09-05
## Date last Modified: 2023-10-30
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

## 1.0. Load relevant packages--------
library(rnaturalearth)
library(terra)
library(tidyverse) #extract might conflict with terra, be sure to specify 'package'::extract
library(sf)
library(stars)

### 1.1 extract natural earth administrative maps for visualization ----
# create a file to save spatial data
gpkg_file <- "raw_data/01-NA_gis_data.gpkg"
dir.create(dirname(gpkg_file), showWarnings = FALSE, recursive = TRUE)

# political boundaries
# land border with lakes removed
ne_land <- ne_download(scale = 50, category = "cultural",
                       type = "admin_0_countries_lakes",
                       returnclass = "sf") %>%
  filter(CONTINENT == "North America") %>%
  st_set_precision(1e6) %>%
  st_union()
# state boundaries for COUNTRY
ne_adm1 <- ne_download(scale = 50, category = "cultural",
                       type = "admin_1_states_provinces",
                       returnclass = "sf") %>% 
  filter(adm0_a3 == "CAN") %>%  ## Change the filter here for different countries
  select(adm1 = name, adm1_code = iso_3166_2)
# country lines
# downloaded globally then filtered to north america with st_intersect()
ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                type = "admin_0_boundary_lines_land",
                                returnclass = "sf") %>% 
  st_geometry()
ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
  as.logical() %>%
  {ne_country_lines[.]}
# borderlines, north america
ne_adm1_lines <- ne_download(scale = 50, category = "cultural",
                             type = "admin_1_states_provinces_lines",
                             returnclass = "sf") %>%
  filter(ADM0_A3 %in% c("USA", "CAN")) %>%
  mutate(iso_a2 = recode(ADM0_A3, USA = "US", CAN = "CAN")) %>% 
  select(country = ADM0_NAME, country_code = iso_a2)

# save all layers to a geopackage
unlink(gpkg_file)
write_sf(ne_land, gpkg_file, "ne_land")
write_sf(ne_adm1, gpkg_file, "ne_adm1")
write_sf(ne_country_lines, gpkg_file, "ne_country_lines")
write_sf(ne_adm1_lines, gpkg_file, "ne_adm1_lines")

## 2.0. Read in relevant shape and attribute files ----
# population centres from 2021 Canadian census; https://www12.statcan.gc.ca/census-recensement/index-eng.cfm
can_popc <- st_read("./raw_data/lpc_000b21a_e.shp")
# relational database csv with attributes for each polygon
popctr <- read_csv("./raw_data/POPCTR.csv", col_types = "dc") %>% as_tibble() 
popctr <- popctr %>% select(DGUID = POPCTRRAdguid,
                            PCNAME = POPCTRRAname, 
                            pop2021 = POPCTRRApop_2021, 
                            urbarea = POPCTRRAarea, 
                            XPRuid = XPRuid) 
## Join polygons with attribute data, and filter to 'large' cities (PCCLASS 4 is >100000)
map_proj <- st_crs("EPSG:4326") #set crs 
urbcan <- can_popc %>% left_join(popctr, by = "DGUID") %>% filter(PCCLASS == 4) %>% st_transform(crs = map_proj)

## Visualize
ne_land <- read_sf("./raw_data/01-NA_gis_data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("./raw_data/01-NA_gis_data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("./raw_data/01-NA_gis_data.gpkg", "ne_adm1_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
on_boundary <- read_sf("./raw_data/01-NA_gis_data.gpkg", "ne_adm1") %>% 
  filter(adm1_code == "CA-ON") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# map
basemap <- ggplot() + 
  geom_sf(data = ne_land) + 
  geom_sf(data = on_boundary, fill="orange") + 
  geom_sf(data = ne_state_lines) + 
  geom_sf(data = ne_country_lines) + 
  geom_sf(data = urbcan, fill= "red") + 
  coord_sf(xlim = c(-76.22, -73.15), ylim = c(45, 46), expand = FALSE) +
    theme_bw()
basemap

# map
x <- list(ne_land = ne_land, ne_country_lines=ne_country_lines, ne_state_lines=ne_state_lines, on_boundary=on_boundary)
basemap <- ggplot() + 
  #geom_sf(data = x$ne_land) + 
  #geom_sf(data = x$on_boundary, fill="orange") + 
  geom_sf(data = NA_basemap$ne_state_lines) + 
  #geom_sf(data = x$ne_country_lines) + 
  geom_sf(data = urbcan, fill= "red") + 
  #coord_sf(xlim = c(-76.22, -73.15), ylim = c(45, 46), expand = FALSE) +
  theme_bw()
basemap

### 2.1. Buffer cities to include boundary land cover types ----
buff_urbcan <- st_buffer(urbcan, dist = 5000)
ggplot() + 
  geom_sf(data = ne_land) + 
  geom_sf(data = on_boundary, fill="orange") + 
  geom_sf(data = ne_state_lines) + 
  geom_sf(data = ne_country_lines) +
  geom_sf(data = buff_urbcan, fill= "red")+
  geom_sf(data = urbcan, fill= "blue") + 
  coord_sf(xlim = c(-76.22, -73.15), ylim = c(45, 46), expand = FALSE) +
  theme_bw()

### 2.2. Write to shapefile for use in other scripts ----
st_write(buff_urbcan, "./data/01-can_cities_5kmbuffer.shp")
