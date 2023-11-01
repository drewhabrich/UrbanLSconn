## HEADER---------------------------
## Script name: functions for urban centres identification
##
## Purpose of script: A collection of functional tasks for 'targets' package
##
## Author: Andrew Habrich
##
## Date Created: 2023-10-18
## Date last Modified: 2023-10-30
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
## Only need to define functions here, no need to run them or require packages to be loaded

####
get_naturalearth <- function (listname, continent, countrycode, countryvector) {
  # political boundaries
  # land border with lakes removed
  ne_land <- ne_download(scale = 50, category = "cultural",
                         type = "admin_0_countries_lakes",
                         returnclass = "sf") %>%
    filter(CONTINENT == continent) %>%
    st_set_precision(1e6) %>%
    st_union()
  # state boundaries for COUNTRY
  ne_adm1 <- ne_download(scale = 50, category = "cultural",
                         type = "admin_1_states_provinces",
                         returnclass = "sf") %>% 
    filter(adm0_a3 == countrycode) %>%  ## Change the filter here for different countries
    select(adm1 = name, adm1_code = iso_3166_2, country = adm0_a3)
  # country lines
  # downloaded globally then filtered to north america with st_intersect()
  ne_country_lines <- ne_download(scale = 50, category = "cultural",
                                  type = "admin_0_boundary_lines_land",
                                  returnclass = "sf") %>% st_geometry()
  ne_country_lines <- st_intersects(ne_country_lines, ne_land, sparse = FALSE) %>%
    as.logical() %>% {ne_country_lines[.]}
  # borderlines, north america
  ne_adm1_lines <- ne_download(scale = 50, category = "cultural",
                               type = "admin_1_states_provinces_lines",
                               returnclass = "sf") %>%
    filter(ADM0_A3 %in% as_vector(countryvector)) %>%
    select(country = ADM0_NAME, country_code = ADM0_A3)
  
  # save all layers to a list
  listname <- list(
    ne_land = ne_land,
    ne_country_lines = ne_country_lines,
    ne_state_lines = ne_adm1,
    ne_boundary = ne_adm1_lines)
}

####
get_populationcentres <- function (shapefile, censusdat, crs, minpopsize) {
  popc <- st_read(paste0("./raw_data/", shapefile))
  # relational database csv with attributes for each polygon
  popctr <- read_csv(paste0("./raw_data/", censusdat), col_types = "dc") %>% as_tibble() 
  popctr <- popctr %>% select(DGUID = POPCTRRAdguid,
                              PCNAME = POPCTRRAname, 
                              pop2021 = POPCTRRApop_2021, 
                              urbarea = POPCTRRAarea, 
                              XPRuid = XPRuid) 
  ## Join polygons with attribute data, and filter to 'large' cities 
  map_proj <- st_crs(crs) #set crs 
  urb <- popc %>% left_join(popctr, by = "DGUID") %>% filter(pop2021 >= minpopsize) %>% st_transform(crs = map_proj)
}

####
plot_urbanarea <- function (zoomtocity, zoomlevel, basemaplist, citylist) {
  city <- citylist %>% filter(PCNAME.x == zoomtocity) %>% st_geometry() %>% st_centroid()
  lon_span <- 360 / 2^zoomlevel
  lat_span <- 180 / 2^zoomlevel
  lon_bounds <- c(st_coordinates(city)[1] - lon_span / 2, st_coordinates(city)[1] + lon_span / 2)
  lat_bounds <- c(st_coordinates(city)[2] - lat_span / 2, st_coordinates(city)[2] + lat_span / 2)
  ggplot() + 
    geom_sf(data = basemaplist$ne_land) + 
    geom_sf(data = basemaplist$ne_state_lines) + 
    geom_sf(data = basemaplist$ne_country_lines) + 
    geom_sf(data = citylist, fill= "red") +
    coord_sf(xlim = lon_bounds, ylim = lat_bounds, expand = FALSE) +
    ggtitle(zoomtocity) +
    theme_bw()
}

