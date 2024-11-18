## HEADER---------------------------
## Script name: functions for urban centres identification
##
## Purpose of script: A collection of functional tasks for 'targets' package
##
## Author: Andrew Habrich
##
## Date Created: 2023-10-18
## Date last Modified: 2024-10-24
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
## Only need to define functions here, no need to run them or require packages to be loaded

#### Function to identify urban centres from the GHSL functional urban areas
get_urbanareas <- function (country, minpopsize, buffersize) {
  # read in the GHSL functional urban areas
  ghsl_fua <- vect("./raw_data/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0/GHS_FUA_UCDB2015_GLOBE_R2019A_54009_1K_V1_0.gpkg") %>% 
    project(y = "EPSG:4326") %>% #convert to lat/long
    makeValid() #make sure the polygons are valid (no self-intersections or unconnected polylines))
  fua <- ghsl_fua %>% filter(Cntry_name == country, FUA_p_2015 > minpopsize)
  fua <- buffer(fua, width = buffersize) #create a buffer around the urban areas
}
    
#### Function to identify which ESA grid cells overlap with urban centres
get_urbanESA_grids <- function (urban_polygons, esa_datadir) {
  esagrid <- vect("./raw_data/esa_worldcover_grid_composites.fgb")
  # identify the cells that overlap with the urban centres
  esacells <- esagrid[urban_polygons] #subset the esagrid to only the cells that overlap with urban areas
  
  # identify the macrotiles with the data
  macrotiles <- esacells %>% as_tibble() %>% 
    #remove all the columns except the first
    select(1) %>%
    separate(tile, into = c("N", "W"), sep = c(3, 7), remove = F) %>% 
    mutate(N = as.numeric(str_extract(N, "\\d{2}")),
           W = as.numeric(str_extract(W, "\\d{3}"))) %>% 
    # Round the N value up and the W value down to the closest macrotile, which are 3x3
    mutate(Nmacrotile = floor(N/3)*3, 
           Wmacrotile = ceiling(W/3)*3) %>% 
    #create a new column for the macrotile name
    mutate(macrotile = paste0("N", str_pad(Nmacrotile, 2, pad = "0"), "W", 
                              str_pad(Wmacrotile, 3, pad = "0")))
  
  # list of macrotiles that overlap with urban centres
  esa_tilelist <- list.files(esa_datadir) %>% 
    str_subset("tif$") %>% #only the raster images
    as_tibble() %>% rename(filename = value) %>% 
    #extract the tile name from the filename, and add it to the dataframe
    mutate(macrotile = str_extract(filename, "(N\\d{2}W\\d{3})")) %>% 
    #filter to only the tiles that overlap with urban centres
    filter(macrotile %in% macrotiles$macrotile)
  
  # save to virtual raster dataset (collection of tiles)
  lc_vrt <- vrt(paste0(esa_datadir, esa_tilelist$filename), 
                filename = "./output/ESALC_tiles.vrt",
                overwrite = T)
  return(lc_vrt)
}

#### Function to merge and crop the landcover rasters to the urban area
get_urbanESA_LC <- function (rasterlist, outputdir, urban_polygons) {
  esa_city_output <- paste0(getwd(), outputdir)
  if(!dir.exists(esa_city_output)){
    dir.create(esa_city_output)
  }
  for (i in 1:nrow(urban_polygons)) {
    #get the city name
    cityname <- urban_polygons$eFUA_name[i] 
    #buffer by 5km to make sure we sample data from outside the urban area
    citypoly <- urban_polygons[i,] %>% buffer(width = 5000) 
    #crop and save the raster to the output directory, save the filename to a list
    #check if the file already exists, if it does skip this step
    if(!file.exists(paste0(esa_city_output, cityname, "_ESA.tif"))){
      crop(rasterlist, citypoly, filename = paste0(esa_city_output, cityname, "_ESA.tif"),
           overwrite = T)
    }
  }
  filenames <- list.files(esa_city_output, full.names = F)
  raster_df <- filenames %>% as_tibble() %>% rename(file = value) %>% mutate(datadir = outputdir)
  return(raster_df)
}

#### Function to get building heights and footprints raster layers
get_builtLCvolume <- function(datadir, outputdir, urban_polygons) {
  # import the raster layer
  ghsl_bh <- rast(paste0(getwd(), datadir, "GHS_BUILT_V_E2020_GLOBE_R2023A_54009_100_V1_0.tif"))
  # create the output list
  ghsl_city_output <- paste0(getwd(), outputdir)
  # create the folder for the output if it doesnt exist
  if(!dir.exists(ghsl_city_output)){
    dir.create(ghsl_city_output)
  }
  for (i in 1:nrow(urban_polygons)) {
    #get the city name
    cityname <- urban_polygons$eFUA_name[i]
    #buffer by 5km to make sure we sample data from outside the urban area
    citypoly <- urban_polygons[i,] %>% buffer(width = 5000) %>% project(y = ghsl_bh)
    #crop and save the raster to the output directory, save the filename to a list
    #check if the file already exists, if it does skip this step
    if(!file.exists(paste0(ghsl_city_output, cityname, "_builtV.tif"))){
      crop(ghsl_bh, citypoly, filename = paste0(ghsl_city_output, cityname, "_builtV.tif"),
           mask = T, overwrite = T)
      }
  }
  
  filenames <- list.files(ghsl_city_output, full.names = F)
  raster_df <- filenames %>% as_tibble() %>% 
    rename(file = value) %>% 
    mutate(datadir = outputdir) %>% 
    mutate(cityname = str_remove(file, "_builtV.tif")) %>% 
    mutate(filepath = paste0(datadir, file)) %>% 
    select(cityname, filepath)
  ## attach the polygon of the city to the dataframe
  urban_polysf <- st_as_sf(urban_polygons) %>% select(eFUA_name, geometry)
  raster_df <- left_join(urban_polysf, raster_df, by = c("eFUA_name" = "cityname"))
  raster_df <- raster_df %>% rename("cityname" = "eFUA_name")
  ## keep only the columns that are needed
  return(raster_df)
}

#### Function to calculate the proportion of each landcover class in each buffer
calc_landscape_metric <- function (citydiv, rasterlist, bufferradius, crs) {
  ## add an if statement to check if the output file already exists
  if(!file.exists(paste0("./output/ebird_data/lsm_div_bycity/", citydiv$cityname[1], bufferradius, "_lsm_div.rds"))){
  ## import data
  city_smpl <- readRDS(paste0("./output/ebird_data/rds_citydata/smpl_", citydiv$cityname[1], ".rds"))
  city <- readRDS(citydiv$rds_div[1])
  ## Remove the zero filled data, keep only the first 7 columns
  city <- city[,1:7]
  ## Join the sampling data to the diversity estimates using the checklist_id
  city_data <- left_join(city, city_smpl, by = c("checklist_id", "latitude", "longitude"), keep = FALSE) %>% 
    #remove useless columns
    select(-c("last_edited_date", "country", "country_code", "state", "state_code", "county", "county_code", 
              "iba_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type")) %>% 
    #create a year, month, and ordinal day column
    mutate(year = as.numeric(format(observation_date, "%Y")), 
           month = as.numeric(format(observation_date, "%m")), 
           ordinal_day = as.numeric(format(observation_date, "%j")))
  ## convert data to spatvector
  city_vect <- vect(city_data, geom=c("longitude", "latitude"), crs = "EPSG:4326", keepgeom=TRUE)
  ## generate circular buffers around each checklist location
  city_buffers <- buffer(city_vect, width = bufferradius)
  
  ## import the landcover raster
  rasterlist <- rasterlist %>% mutate(filepath = paste0(getwd(), rasterlist$datadir, rasterlist$file))
  raster <- rast(rasterlist %>% 
                 filter(str_detect(filepath, citydiv$cityname[1])) %>% 
                 pull(filepath))
  
  ## crop the raster to the city area
  city_lc <- crop(raster, city_vect, mask = T)
  
  ## project the city vector and raster to the same crs (a meter based reference system for lsm)
  city_buffers <- project(city_buffers, y = crs)
  city_lc <- project(city_lc, y = crs, method = "near") #need to include nearest neighbour method to keep categories
  
  ## calculate selected landscapemetrics in each buffer using landscapemetrics package
  lsm <- NULL
  for (i in seq_len(nrow(city_buffers))) {
    buffer_i <- city_buffers[i, ]
    
    # crop and mask landcover raster so all values outside buffer are missing
    lsm[[i]] <- crop(city_lc, buffer_i, mask = T) %>% 
      # calculate landscape metrics
      calculate_lsm(level = "class", metric = c("pland")) %>% 
      # add variables to uniquely identify each point
      mutate(checklist_id = buffer_i$checklist_id) %>% 
      select(checklist_id, class, metric, value)
  }
  lsm <- bind_rows(lsm)
  ### replace class values with the landcover class names
  lcclass_lookup <- tibble(class = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100), 
                           class_name = c("treecover", "shrubland", "grassland", "cropland", 
                                          "builtup", "bare", "snow", "water", "wetland", 
                                          "mangrove", "mosslichen"))
  
  ### convert to wide format and clean the data
  lsm_wide <- lsm %>% 
    #fill any missing classes with zeros
    complete(class, nesting(checklist_id), fill = list(value = 0)) %>%
    #class names
    inner_join(select(lcclass_lookup, class, class_name), by = "class") %>%  
    #transform to wide format
    pivot_wider(values_from = value, 
                names_from = c(class, class_name, metric), 
                names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{class_name}",
                names_sort = T) %>%
    arrange(checklist_id)
  
  #save to RDS
  output_dir <- "./output/ebird_data/lsm_div_bycity/"
  #if folder does not exist, create it
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
    }
  saveRDS(lsm_wide, paste0(output_dir, citydiv$cityname[1], bufferradius, "_lsm_div.rds"))
  }
  
  ## return a dataframe of the file names and location of the RDSfile
  files_lsmdiv_bycity <- tibble(cityname = citydiv$cityname, 
                                rds_lsm = str_c(output_dir, citydiv$cityname, bufferradius, "_lsm_div.rds"))
  return(files_lsmdiv_bycity)
}

#### Function to extract building volumes for the buffered checklists
extract_raster_inpolygon <- function (chkldata, rasterlist, bufferradius) {
  CN <- chkldata$cityname[1]
  city_rast <- rast(paste0(getwd(), rasterlist %>% filter(cityname == CN) %>% pull(filepath[1])))
  
  city_div <- readRDS(chkldata %>% filter(cityname == CN) %>% pull(rds_div[1])) 
  city_div <- city_div[,1:7]
  city_div <- vect(city_div, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  city_div <- project(city_div, y = city_rast)
  city_div_buff <- buffer(city_div, width = bufferradius)
  city_sfbuff <- st_as_sf(city_div_buff)
  
  ### extract the volume of built up area in each buffer
  city_extractedvalues <- exact_extract(city_rast, city_sfbuff, fun = c("mean","stdev","sum"),
                                        append_cols = "checklist_id")
  #save to csv
  write.csv(city_extractedvalues, paste0(getwd(), "/output/ebird_data/builtv_chkl_bycity/", CN, "_builtv_chkl.csv"))
  #save to RDS
  #create output directory if it doesnt exist
  output_dir <- "./output/ebird_data/builtv_chkl_bycity/"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(city_extractedvalues, paste0(output_dir, CN, "_builtv_chkl.rds"))

  ## return a dataframe of the file names and location of the RDSfile
  builtv_bycity_df <- tibble(cityname = CN,
                             rds_builtv = str_c(output_dir, CN, "_builtv_chkl.rds"))
  return(builtv_bycity_df)
}
# extract_raster_inpolygon <- function (chkldata, rasterlist, bufferradius, crs) {
#   city_smpl <- readRDS(chkldata$smplRDS[1])
#   city_chkl <- vect(city_smpl, geom = c("longitude", "latitude"), crs = "EPSG:4326", keepgeom = TRUE)
#   city_poly <- vect(chkldata$geometry[1])
#   city_chkl_sub <- city_chkl[city_poly] #subset to only the checklists in the polygon
#   
#   ##create a string of the filepath for the city by pasting the working directly, the datadir, and file
#   city_rast_filepath <- paste0(getwd(), rasterlist$filepath[1])
#   city_rast <- rast(city_rast_filepath) %>% project(y = crs) #import and project the raster
#   city_chklbuff <- buffer(city_chkl_sub, width = bufferradius) %>% project(y = crs) #buffer and reproject the checklists
#   
#   ### extract the volume of built up area in each buffer
#   city_extractedvalues <- exact_extract(city_rast, st_as_sf(city_chklbuff), fun = c("mean","stdev","sum"), 
#                                         append_cols = "checklist_id")
#   
#   #save to RDS
#   #create output directory if it doesnt exist
#   output_dir <- "./output/ebird_data/builtv_chkl_bycity/"
#   if (!dir.exists(output_dir)){
#     dir.create(output_dir)
#   }
#   saveRDS(city_extractedvalues, paste0(output_dir, chkldata$cityname[1], "_builtv_chkl.rds"))
# 
#   ## return a dataframe of the file names and location of the RDSfile
#   builtv_bycity_df <- tibble(cityname = chkldata$cityname[1], 
#                              rds_builtv = str_c(output_dir, chkldata$cityname[1], "_builtv_chkl.rds"))
#   return(builtv_bycity_df)
# }