## HEADER---------------------------
## Script name: functions for working with ebird data
##
## Purpose of script: A collection of functional tasks for 'targets' package working with ebird data
##
## Author: Andrew Habrich
##
## Date Created: 2023-10-30
## Date last Modified: 2024-08-30
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
## Only need to define functions here, no need to run them or require packages to be loaded

#### check for data of extracted list of cities
checkdata_by_city <- function (urbanarea_vector) {
  urbarea_sf <- urbanarea_vector %>% st_as_sf()
  ebirdcity_out <- urbarea_sf %>% dplyr::select(eFUA_name, geometry) %>% 
    mutate(cityname = eFUA_name,
           checklist = str_c("./raw_data/ebd2022_", eFUA_name, "_fua.txt"),
           sampling = str_c("./raw_data/ebdsamp2022_", eFUA_name, "_fua.txt")) %>% 
    mutate(chklRDS = str_c("./output/ebird_data/rds_citydata/chkl_", cityname, ".rds"),
           smplRDS = str_c("./output/ebird_data/rds_citydata/smpl_", cityname, ".rds")) %>% 
    mutate(checklist_exists = file.exists(checklist),
           sampling_exists = file.exists(sampling)) %>% 
    #check for the existence of the RDS files
    mutate(chklRDS_exists = file.exists(chklRDS),
           smplRDS_exists = file.exists(smplRDS)) %>% 
    dplyr::select(-eFUA_name)
}

#### generate ebird filters per city based on urban area bounding box      
ebirdfilter_by_city <- function(urban_polygons, ebd_data, ebirdcity_out) {
  urbarea_sf <- st_as_sf(urban_polygons)
  city_filters <- list()
  for (i in 1:nrow(urbarea_sf)) {
    filt <- ebd_data %>%
      auk_bbox(urbarea_sf %>% slice(i)) %>%
      auk_complete()
    city_filters[[i]] <- filt
  }

  data_list <- list()
  for (i in 1:length(city_filters)) {
    citydat <- auk_filter(city_filters[[i]],
               file = ebirdcity_out[i,2]$checklist,
               file_sampling = ebirdcity_out[i,3]$sampling,
               filter_sampling = T, overwrite = T)
    data_list[[i]] <- citydat
  }
  names(data_list) <- ebirdcity_out$cityname
  return(data_list)
}

#### Import the checklist and sampling data for each city and save as RDS files
write_cityRDS <- function(citylist, ebird_cityfilters) {
  # Initialize a list to store results
  results <- list()
  for (i in 1:nrow(citylist)) {
    chkl_file <- str_c("./output/ebird_data/rds_citydata/chkl_", citylist$cityname[i], ".rds")
    smpl_file <- str_c("./output/ebird_data/rds_citydata/smpl_", citylist$cityname[i], ".rds")
    # Check if both files exist
    if (!file.exists(chkl_file) || !file.exists(smpl_file)) {
      # Read in the data
      chkl <- read_ebd(ebird_cityfilters[[i]]$output) 
      smpl <- read_sampling(ebird_cityfilters[[i]]$output_sampling)
      # Save the RDS files if they don't already exist
      saveRDS(chkl, file = chkl_file)
      saveRDS(smpl, file = smpl_file)
    }
    # Store the filenames in the results list
    results[[i]] <- list(
      cityname = citylist$cityname[i],
      chkl_filename = chkl_file,
      smpl_filename = smpl_file
    )
  }
  # Convert the results list to a dataframe
  results_df <- bind_rows(results)
  return(results_df)
}

#### filter ebird data by city polygon
spatial_filter <- function(citylist, ebirddat_RDS) {
  filtered_ebird_bycity <- list()
  # read in the filtered data
  for (i in 1:nrow(citylist)) {
    chkldat <- readRDS(ebirddat_RDS$chkl_filename[i])
    smpldat <- readRDS(ebirddat_RDS$smpl_filename[i])
    # convert to points geometries
    chkl_sf <- chkldat %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    smpl_sf <- smpldat %>% st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
    # spatial filtering
    chkl_region <- st_filter(chkl_sf, citylist$geometry)
    smpl_region <- st_filter(smpl_sf, citylist$geometry)
    # join by checklist_id
    chkl <- semi_join(chkldat, chkl_region, by = "checklist_id")
    smpl <- semi_join(smpldat, smpl_region, by = "checklist_id")
    # save to a list with the city name as the name of the list item
    datlist <- list(obs = chkl, smpl = smpl)
    saveRDS(datlist, file = str_c("./output/ebird_data/rds_ebird_filteredbycity/", citylist$cityname[i], ".rds"))
    filtered_ebird_bycity <- append(filtered_ebird_bycity, list(datlist))
    names(filtered_ebird_bycity)[i] <- citylist$cityname[i]
  }
  # create a tibble and create columns for the cityname and RDS file location
  files_df <- tibble(cityname = citylist$cityname, 
                     rds_file = str_c("./output/ebird_data/rds_ebird_filteredbycity/", citylist$cityname, ".rds"))
  return(files_df)
}

#### zerofill ebird datasets 
zerofill_citylist <- function(ebird_filtdata) {
  citydat <- readRDS(ebird_filtdata$rds_file[1])
  ## zero fill the data (0s for absence, 1s for presence) and collapse into 1 dataframe
  zf_citydat <- auk_zerofill(x = citydat$obs, sampling_events = citydat$smpl, collapse = T)
  output_dir <- "./output/ebird_data/rds_ebird_zf_bycity/"
  ## write the zero-filled data to a RDS file if it doesnt exist, otherwise do nothing
  if(!file.exists(paste0(output_dir, ebird_filtdata$cityname[1], "_zf.rds"))){
    saveRDS(zf_citydat, paste0(output_dir, ebird_filtdata$cityname[1], "_zf.rds"))
  }
  ## return the cityname and the RDS file location
  files_df <- tibble(cityname = ebird_filtdata$cityname, 
                     rds_zfdf = str_c(output_dir, ebird_filtdata$cityname, "_zf.rds"))
  return(files_df)
}

#### tidy zerofilled ebird datasets 
zf_tidy <- function(zf_bycity_dataframes){
# function to convert time observation to hours since midnight
  time_to_decimal <- function(x) {
    x <- hms(x, quiet = TRUE)
    hour(x) + minute(x) / 60 + second(x) / 3600
  }
  # tidy
  zf_dataframe <- readRDS(zf_bycity_dataframes$rds_zfdf[1])
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
  #save to RDS
  output_dir <- "./output/ebird_data/rds_ebird_tidyzf/"
  if(!file.exists(paste0(output_dir, zf_bycity_dataframes$cityname[1], "_tidyzf.rds"))){
    saveRDS(checklists, paste0(output_dir, zf_bycity_dataframes$cityname[1], "_tidyzf.rds"))
  }
  files_zfdf <- tibble(cityname = zf_bycity_dataframes$cityname, 
                       rds_tidyzf = str_c(output_dir, zf_bycity_dataframes$cityname, "_tidyzf.rds"))
  return(files_zfdf)
}

#### calculate species diversity estimates for each checklist
calc_biodiversity <- function(checklist_dataframe) {
  zf <- readRDS(checklist_dataframe$rds_tidyzf[1]) 
  zf <- zf %>% mutate(observation_count = ifelse(is.na(observation_count), 1, observation_count))
  zf_wide <- zf %>% select(c('checklist_id','scientific_name','observation_count')) %>%
    pivot_wider(id_cols = checklist_id,
                names_from = scientific_name,
                values_from = observation_count,
                values_fill = 1, ### fill observed, but abundance unknown with 1 so the calcs work
                names_sort = F)
  ## calculate diversity estimates
  zf_div <- zf_wide %>% mutate(checklist_id = checklist_id,
                               shannondiv = diversity(zf_wide %>% select(-(checklist_id)), index = 'shannon'),
                               specrich = specnumber(zf_wide%>% select(-(checklist_id))), #species richness
                               evenness = shannondiv / log(specrich),
                               abundance = rowSums(zf_wide %>% select(-checklist_id)))
  ## join with zf spatial features dataframe by checklist_id
  biodiv <- zf_div %>% left_join(zf %>% distinct(checklist_id, .keep_all = TRUE) %>%
                                   select(checklist_id, latitude, longitude), by = "checklist_id") %>%
    relocate(c("shannondiv","specrich","evenness","abundance", "latitude","longitude"), .after = checklist_id) 
  ## remove rows with 0 species richness (e.g. no species observed)
  biodiv <- biodiv %>% filter(specrich > 0)
  #save to RDS
  output_dir <- "./output/ebird_data/rds_div_bycity/"
  #if folder does not exist, create it
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  #if file does not exist, save it
  if(!file.exists(paste0(output_dir, checklist_dataframe$cityname[1], "_div.rds"))){
    saveRDS(biodiv, paste0(output_dir, checklist_dataframe$cityname[1], "_div.rds"))
  }
  files_divbycity <- tibble(cityname = checklist_dataframe$cityname, 
                       rds_div = str_c(output_dir, checklist_dataframe$cityname, "_div.rds"))
  return(files_divbycity)
}