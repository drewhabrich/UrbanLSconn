## HEADER---------------------------
## Script name: functions for working with ebird data
##
## Purpose of script: A collection of functional tasks for 'targets' package working with ebird data
##
## Author: Andrew Habrich
##
## Date Created: 2023-10-30
## Date last Modified: 2023-10-30
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
## Only need to define functions here, no need to run them or require packages to be loaded

#### extract list of cities
ebird_by_city <- function (urb_bufferslist) {
  ebirdcity_out <- urb_bufferslist %>% select(PCNAME.x) %>% 
    mutate(cityname = PCNAME.x,
           geometry = geometry,
           checklist = str_c("./raw_data/ebd_", PCNAME.x, "_5kmbuff.txt"),
           sampling = str_c("./raw_data/ebdsamp_", PCNAME.x, "_5kmbuff.txt")) %>% select(-PCNAME.x) %>% 
    mutate(checklist_exists = file.exists(checklist),
           sampling_exists = file.exists(sampling))
}

#### generate ebird filters per city based on urban area bounding box      
generate_city_filters <- function(city_data, ebd_data) {
  city_filters <- list()
  for (i in 1:nrow(city_data)) {
    filt <- ebd_data %>%
      auk_bbox(city_data %>% slice(i)) %>%
      auk_complete()
    city_filters[[i]] <- filt
  }
  return(city_filters)
}

#### create ebird data subsets and put them into a list
ebird_data_bycity <- function(city_filters, ebird_citylist) {
 data_list <- list()
  for (i in 1:length(city_filters)) {
  citydat <- auk_filter(city_filters[[i]],
             file = ebird_citylist[i,2]$checklist, file_sampling = ebird_citylist[i,3]$sampling,
             filter_sampling=T, overwrite = T)
  data_list[[i]] <- citydat
  }
 return(data_list)
}

#### zerofill and tidy ebird datasets 
zf_tidy <- function(zf_dataframe){
  # function to convert time observation to hours since midnight
  time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
  }
  # tidy
  zf <- zf_dataframe %>% 
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
  zf_filtered <- zf %>% 
    filter(effort_hours <= 6,
           effort_distance_km <= 10,
           effort_speed_kmph <= 100,
           number_observers <= 10)
  checklists <- zf_filtered %>%
    select(checklist_id, observer_id,
           observation_count, species_observed, scientific_name, 
           state_code, locality_id, latitude, longitude,
           protocol_type, all_species_reported,
           observation_date, year, day_of_year,
           hours_of_day,
           effort_hours, effort_distance_km, effort_speed_kmph,
           number_observers)
}

#### summarize counts
city_summary_table <- function(zf_checklists, urb_popcentres) {
## create empty summary
  citynames <- urb_popcentres$PCNAME.x
  names(zf_checklists) <- citynames
  summary_list <- map(zf_checklists, ~ summarize(.x, checklist_count = n_distinct(checklist_id)))
  summary_df <- bind_rows(summary_list, .id = "Dataframe")
}

#### city summary plot
city_summary_plot <- function(summary_table)
summary_table %>% distinct() %>% 
  mutate(Dataframe = fct_reorder(Dataframe, checklist_count)) %>% 
  # Sort the data by IDcount in descending order
  ggplot(aes(y = Dataframe, x = checklist_count)) +
    geom_bar(stat = "identity") +
    labs(y = "City", x = "Checklist count") +
    theme_minimal()

#### pivot wide and calculate species diversity estimates for each checklist
calc_diversity <- function(checklist_dataframe) {
  zf <- checklist_dataframe
  zf_wide <- zf %>% select(c('checklist_id','scientific_name','observation_count')) %>% 
    pivot_wider(id_cols = checklist_id,
                names_from = scientific_name, 
                values_from = observation_count,
                values_fill = 1, ### fill observed, but abundance unknown with 1 so the calcs work
                names_sort = F)
    ## calculate diversity estimates
  zf_div <- zf_wide %>% mutate(checklist_id = checklist_id,
                                 shannondiv = diversity(zf_wide %>% select(!(checklist_id)), index = 'shannon'),
                                 specrich = specnumber(zf_wide), #species richness
                                 evenness = shannondiv / log(specnumber(zf_wide)))
  ## join with zf spatial features dataframe by checklist_id
  div <- zf_div %>% left_join(zf %>% distinct(checklist_id, .keep_all = TRUE) %>%
                                       select(checklist_id, latitude, longitude), by = "checklist_id") %>% 
                     relocate(c("shannondiv","specrich","evenness","latitude","longitude"), .after = checklist_id)
}

#### plot diversity
biodiv_plot <- function (div_datlist, zoomtocity, zoomlevel, basemaplist, citylist) {
  names(div_datlist) <- citylist$PCNAME.x
  city <- citylist %>% filter(PCNAME.x == zoomtocity) %>% st_geometry() %>% st_centroid()
  lon_span <- 360 / 2^zoomlevel
  lat_span <- 180 / 2^zoomlevel
  lon_bounds <- c(st_coordinates(city)[1] - lon_span / 2, st_coordinates(city)[1] + lon_span / 2)
  lat_bounds <- c(st_coordinates(city)[2] - lat_span / 2, st_coordinates(city)[2] + lat_span / 2)
  coi <- div_datlist[[zoomtocity]]
coi %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>% 
  ggplot() +
    geom_sf(data = basemaplist$ne_land) +
    geom_sf(data = basemaplist$ne_state_lines) +
    geom_sf(data = basemaplist$ne_country_lines) +
    geom_sf(data = citylist %>% filter(PCNAME.x == zoomtocity), fill= "pink") +
    geom_sf(mapping = aes(colour = specrich)) + viridis::scale_color_viridis(alpha=0.75, begin = 0.2, end = 1) +
    ggtitle(zoomtocity) +
    coord_sf(xlim = lon_bounds, ylim = lat_bounds, expand = FALSE, default_crs = sf::st_crs(4326)) +
    theme_bw()
}
# 
# biodiv_plot <- function (div_datlist, div_metric, zoomtocity, zoomlevel, basemaplist, citylist) {
#   names(div_datlist) <- citylist$PCNAME.x
#   cityofinterest <- div_datlist$"zoomtocity" %>% as_tibble(.name_repair = "minimal")
#   city <- citylist %>% filter(PCNAME.x == zoomtocity) %>% st_geometry() %>% st_centroid()
#   lon_span <- 360 / 2^zoomlevel
#   lat_span <- 180 / 2^zoomlevel
#   lon_bounds <- c(st_coordinates(city)[1] - lon_span / 2, st_coordinates(city)[1] + lon_span / 2)
#   lat_bounds <- c(st_coordinates(city)[2] - lat_span / 2, st_coordinates(city)[2] + lat_span / 2)
#   coi_sf <- cityofinterest %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326)
#   coi_sf %>% 
#     ggplot() + 
#     geom_sf(data = NA_basemaps$ne_land) + 
#     geom_sf(data = NA_basemaps$ne_state_lines) + 
#     geom_sf(data = NA_basemaps$ne_country_lines) + 
#     geom_sf(data = citylist %>% filter(PCNAME.x == zoomtocity), fill= "pink") +
#     geom_sf(aes(col = div_metric)) + viridis::scale_color_viridis(alpha=0.75, begin = 0.2, end = 1) +
#     ggtitle(zoomtocity) +
#     coord_sf(xlim = lon_bounds, ylim = lat_bounds, expand = FALSE) +
#     theme_bw()
# }