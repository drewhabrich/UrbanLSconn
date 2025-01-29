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
get_urbanESA_LC <- function (rasterlist, outputdir, urban_polygons, targetcrs) {
  esa_city_output <- paste0(getwd(), outputdir)
  if(!dir.exists(esa_city_output)){
    dir.create(esa_city_output)
  }
  crs <- targetcrs
  for (i in 1:nrow(urban_polygons)) {
    #get the city name
    cityname <- urban_polygons$eFUA_name[i] 
    #buffer by 5km to make sure we sample data from outside the urban area
    citypoly <- urban_polygons[i,] %>% buffer(width = 5000) 
    #crop and save the raster to the output directory, save the filename to a list
    #check if the file already exists, if it does skip this step
    if(!file.exists(paste0(esa_city_output, cityname, ".tif"))){
      cityraster <- terra::crop(rasterlist, citypoly, mask = T) %>% 
                    terra::as.factor() %>% 
                    terra::project(y = crs, method = "near", threads = TRUE)
      writeRaster(cityraster, filename = paste0(esa_city_output, cityname, ".tif"), overwrite = TRUE)
    }
  }
  filenames <- list.files(esa_city_output, pattern = ".tif$", full.names = F, ignore.case = T)
  raster_df <- filenames %>% as_tibble() %>% 
    rename(file = value) %>% 
    mutate(filepath = paste0(outputdir, file)) %>% 
    mutate(cityname = str_remove(file, ".tif"))
  return(raster_df)
}

#### Function to calculate the distance from the checklist to the functional urban boundary
calc_urbandistance <- function(citychkls, urban_polygons, outputdir) {
  # read in urban area and centroids
  ua_centr <- vect("./data/00-ghs_fua_canada.gpkg", layer = "urbcentre_loc") %>% project(y = "EPSG:4326")
  ua <- urban_polygons #already a spatvector
  # Ensure the output directory exists
  if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }

  # Create a for loop to calculate the distance from each checklist to the city boundary
  for(i in 1:nrow(citychkls)){
    cityname <- pull(citychkls[i,1])
    chkldata <- readRDS(pull(citychkls[i,2])) %>% dplyr::select(1:7)
    #note: 1km buffer was added to the city polygons
    # coerce to spatvector
    pt_vect <- vect(chkldata, geom = c("longitude", "latitude"), crs = "EPSG:4326") 
    poly_vect <- ua %>% 
      filter(eFUA_name == cityname) %>% 
      as.lines() #need to coerce polygon to lines
    
    # Bind the data to the original data as a new column
    boundary_dist <- terra::distance(pt_vect, poly_vect, pairwise = F, unit = "m")
    nearest_ua  <-terra::nearest(pt_vect, ua_centr, pairs = F, method = "geo", lines = T)
    uacentr_dist <- terra::distance(pt_vect, ua_centr[nearest_ua,], pairwise = F, unit = "m")
    uacentr_min <- apply(uacentr_dist, MARGIN = 1, min) #margin = 1 for ROWWISE
    
    # Save to the orignal dataframe
    chkldistances <- chkldata %>% 
      mutate(boundarydist_m = as.vector(boundary_dist)) %>% 
      mutate(uacentrdist_m = as.vector(uacentr_min))
    
    # Save the data to csv
    write_csv(chkldistances, paste0(outputdir, cityname, "_chkl_dist.csv"))
  }
  
  # Create and return a summary dataframe
  filenames <- list.files(outputdir, full.names = F)
  raster_df <- filenames %>% as_tibble() %>% 
    rename(file = value) %>% 
    mutate(filepath = paste0(outputdir, file)) %>% 
    mutate(cityname = str_remove(file, "_chkl_dist.csv"))
    return(raster_df)
}

#### Function to get building heights and footprints raster layers
get_builtLCvolume <- function(datadir, outputdir, urban_polygons, targetcrs) {
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
      crop(ghsl_bh, citypoly, mask = T) %>% 
      project(y = targetcrs, method = "bilinear", threads = TRUE,
                filename = paste0(ghsl_city_output, cityname, "_builtV.tif"), 
                overwrite = TRUE)
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
  urban_polysf <- st_as_sf(urban_polygons) %>% 
    st_transform(crs = 3978) %>% 
    select(eFUA_name, geometry)
  raster_df <- left_join(urban_polysf, raster_df, by = c("eFUA_name" = "cityname"))
  raster_df <- raster_df %>% rename("cityname" = "eFUA_name")
  ## keep only the columns that are needed
  return(raster_df)
}

#### Function to calculate the proportion of each landcover class in each buffer
calc_landscape_metric <- function (citydiv, rasterlist, bufferradius, crs) {
  output_dir <- "./output/ebird_data/lsm_div_bycity/"
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

#### Function to calculate the proportion of each landcover class in each buffer####
## sample lsm from landscapemetrics
calc_pland_metric <- function(cityname, citydata, landcover, 
                              bufferradius, outputdir, targetcrs) {
  #if folder does not exist, create it
  if(!dir.exists(outputdir)){
    dir.create(outputdir)
  }
  ## get the specific city input
  if(!file.exists(paste0(outputdir, cityname, bufferradius, "_pland_div.csv"))){ ###IS IF STATEMENT RIGHT for FILENAME?
    ## import checklist data & convert data to spatvector
    citychkldata <- citydata 
    city_vect <- vect(citychkldata, geom = c("longitude", "latitude"), crs = "EPSG:4326") %>% 
      project(y = targetcrs)
    ## generate circular buffers around each checklist location
    city_buffers <- buffer(city_vect, width = bufferradius)
    
    ## calculate selected landscapemetrics in each buffer using landscapemetrics package
    chkl_pland <- sample_lsm(landscape = landcover, 
                             y = city_buffers, 
                             plot_id = city_buffers$checklist_id, 
                             shape = "circle", size = bufferradius, all_classes = T, 
                             what = "lsm_c_pland", 
                             directions = 8, neighbourhood = 8)
    
    ### replace class values with the landcover class names
    lcclass_lookup <- tibble(class = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100), 
                             class_name = c("treecover", "shrubland", "grassland", "cropland", 
                                            "builtup", "bare", "snow", "water", "wetland", 
                                            "mangrove", "mosslichen"))
    
    ### convert to wide format and clean the data
    chkl_pland %>% 
      rename(checklist_id = plot_id) %>% 
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
    }
}

## scale lsm from landscapemetrics
# calc_scalesample <- function(cityname, citydata, landcover, 
#                               bufferradius, outputdir, targetcrs) {
#   #if folder does not exist, create it
#   if(!dir.exists(outputdir)){
#     dir.create(outputdir)
#   }
#   ## get the specific city input
#   if(!file.exists(paste0(outputdir, cityname, "multiscale_pland.csv"))){
#     ## import checklist data & convert data to spatvector
#     citychkldata <- citydata 
#     city_vect <- vect(citychkldata, geom = c("longitude", "latitude"), crs = "EPSG:4326") %>% 
#       project(y = targetcrs)
#    
#     ## calculate selected landscapemetrics in each buffer using landscapemetrics package
#     chkl_pland <- scale_sample(landscape = landcover, 
#                              y = city_vect, 
#                              #plot_id = city_vect$checklist_id, 
#                              shape = "circle", size = bufferradius, 
#                              what = "lsm_c_pland", 
#                              directions = 8, neighbourhood = 8)
#     
    ### replace class values with the landcover class names
    # lcclass_lookup <- tibble(class = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100), 
    #                          class_name = c("treecover", "shrubland", "grassland", "cropland", 
    #                                         "builtup", "bare", "snow", "water", "wetland", 
    #                                         "mangrove", "mosslichen"))
    
    # ### convert to wide format and clean the data
    # lsm_wide <- chkl_pland %>% 
    #   rename(checklist_id = plot_id) %>% 
    #   #fill any missing classes with zeros
    #   complete(class, nesting(checklist_id), fill = list(value = 0)) %>%
    #   #class names
    #   inner_join(select(lcclass_lookup, class, class_name), by = "class") %>%  
    #   #transform to wide format
    #   pivot_wider(values_from = value, 
    #               names_from = c(class, class_name, metric), 
    #               names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{class_name}",
    #               names_sort = T) %>%
    #   arrange(checklist_id)
    
    #save to file
    #write_csv(chkl_pland, paste0(outputdir, cityname, "multiscale_pland.csv"), append = T)
    #return(chkl_pland)
  #}
  #files_planddiv_bycity <- tibble(rds_pland = list.files(output_dir, pattern = "_pland_div.rds"))
#}

# write results to a CSV file
# write_plandresults <- function(results, city_name, outputdir) {
#   file_name <- paste0(outputdir, city_name, "_pland_results.csv")
#   write_csv(results, file_name)
#   ## check what files are in the output directory and save to a tibble and return the tibble
#   plandfiles <- tibble(cityname = city_name, file = list.files(outputdir))
#   return(plandfiles)
# }

#### Function to sample landscape metrics around sample points (checklists)####
calc_aggrm_metric <- function(cityname, citydata, landcover, 
                             bufferradius, outputdir, targetcrs, what_metrics) {
  ## create the metric label
  agg_metrics <- what_metrics
  m_label <- agg_metrics %>% str_remove("lsm_l_") %>% paste(., collapse = "_")
  
  ## define output directory
  if(!dir.exists(outputdir)){
    dir.create(outputdir)
  }
  
  ## get the specific city input
  if(!file.exists(paste0(outputdir, cityname, bufferradius,"_", m_label, ".csv"))){
    citychkldata <- citydata 
    city_vect <- vect(citychkldata, geom = c("longitude", "latitude"), crs = "EPSG:4326") %>% 
      project(y = targetcrs)
    ## generate circular buffers around each checklist location
    city_buffers <- buffer(city_vect, width = bufferradius)
    ## sample landscape metrics around chkl locations
    city_aggrmetrics <- sample_lsm(landscape = landcover, y = city_buffers, 
                                 plot_id = city_buffers$checklist_id, 
                                 shape = "circle", size = bufferradius, all_classes = T, 
                                 what = agg_metrics, 
                                 directions = 8, neighbourhood = 8)
  
  ## convert to wide format
  city_aggrmetrics %>% 
    rename(checklist_id = plot_id) %>% 
    pivot_wider(values_from = value, 
                names_from = c(level, metric), 
                names_glue = "{level}_{metric}",
                names_sort = T) %>%
    arrange(checklist_id)
    }
}

#### Function to extract building volumes for the buffered checklists####
extract_bv_bychkl <- function (data, bufferradius, outputdir) {
  #create output directory if it doesnt exist
  output_dir <- outputdir
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  ## get the specific city input
  city <- data[1,]
  CN <- city$cityname
  
  if(!file.exists(paste0(output_dir, CN, bufferradius, "_builtv_chkl.rds"))){
  city_rast <- rast(city$LCpath)
  city_div <- readRDS(city$rds_div) %>% select(1:7)
  city_div <- vect(city_div, geom = c("longitude", "latitude"), crs = "EPSG:4326") %>% 
    project(y = city_rast)
  city_div_buff <- buffer(city_div, width = bufferradius)
  ## get the area of the buffer
  polyarea <- expanse(city_div_buff, unit = "m")
  
  ## vector needs to be sf for exact_extract
  city_sfbuff <- st_as_sf(city_div_buff)
  city_sfbuff$polyarea <- polyarea
  ### extract the volume of built up area in each buffer
  city_extractedvalues <- exact_extract(city_rast, city_sfbuff, fun = c("mean","stdev","sum"),
                                        append_cols = c("checklist_id", "polyarea"))
  #save to csv
  write.csv(city_extractedvalues, paste0(output_dir, CN, bufferradius, "_builtv_chkl.csv"))
  #save to RDS
  saveRDS(city_extractedvalues, paste0(output_dir, CN, bufferradius, "_builtv_chkl.rds"))
  }
  ## return a dataframe of the file names and location of the RDSfile
  builtv_bycity_df <- tibble(cityname = CN,
                             rds_builtv = str_c(output_dir, CN, bufferradius, "_builtv_chkl.rds"))
  return(builtv_bycity_df)
}

