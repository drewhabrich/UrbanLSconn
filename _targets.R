## Load required packages ##
library(targets)
library(tarchetypes)
library(tidyterra)
library(geotargets)
library(autometric)
library(crew)
library(qs2)
library(landscapemetrics)
library(dplyr)

# This is where you write source(\"R/functions.R\")
tar_source("R/")

# define the controller and the log using crew and autometrics
controller <- crew_controller_local(
  name = "ctrl",
  reset_globals = T,
  garbage_collection = T,
  workers = 2,
  seconds_idle = 30
  # options_metrics = crew_options_metrics(
  #   path = "worker_log_directory/",
  #   seconds_interval = 60)
)

# Start the controller when the target pipeline is active
# if (tar_active()) {
#   controller$start()
#   log_start(
#     path = "main_log_directory/log.txt",
#     seconds = 60,
#     pids = controller$pids()
#   )
# }

# Set target-specific options such as packages that are required:
tar_option_set(packages = c("targets", "tarchetypes", "geotargets","qs2","crew",
                            "ggspatial", "sf", "terra", "tidyterra", "rnaturalearth",
                            "auk", "vegan", "landscapemetrics", "exactextractr", #"gdistance", 
                            "tidyverse", "dplyr", "ggpubr", "RColorBrewer"),
               format = "qs", memory = "auto", error = "stop",
               deployment = "worker", storage = "main", retrieval = "main", 
               garbage_collection = 1,
               controller = controller, 
               resources = tar_resources(qs = tar_resources_qs())
)

list(
  #### Read in spatial data ####################################################
  tar_terra_vect(urbanareas, get_urbanareas(country = "Canada", 
                                            minpopsize = 100000, 
                                            buffersize = 1000)), 
  tar_terra_rast(urban_esarasters, get_urbanESA_grids(urban_polygons = urbanareas, 
                                                      esa_datadir = "F:/DATA/ESA2021LANDCOVER/"),
                 cue = tar_cue(mode = "never")),
  
  tar_target(urban_ESA_LC, get_urbanESA_LC(rasterlist = urban_esarasters,
                                           outputdir = "./output/ESArasters_factor/",
                                           urban_polygons = urbanareas,
                                           targetcrs = "EPSG:3978"),
             cue = tar_cue(mode = "never")),
  tar_target(ghsl_builtLC, get_builtLCvolume(datadir = "./raw_data/ghsl/",
                                             outputdir = "./output/GHSL_builtV/",
                                             urban_polygons = urbanareas,
                                             targetcrs = "EPSG:3978"),
             cue = tar_cue(mode = "never")),
  
  #### Read in the ebird data ##################################################
  tar_target(ebd_can, auk_ebd(file = "./raw_data/ebd_can_obsdata.txt", #this is the .txt with the data
                              file_sampling = "./raw_data/ebd_can_sampdata.txt"),
             cue = tar_cue("never")),
  # create and check list of cities with ebird data
  tar_target(ebird_citydatalist, checkdata_by_city(urbanarea_vector = urbanareas)), 
  # generate ebird filters per city
  tar_target(ebird_citydata, ebirdfilter_by_city(urban_polygons = urbanareas, 
                                                 ebd_data = ebd_can,
                                                 ebirdcity_out = ebird_citydatalist),
             cue = tar_cue("never")),
  #create and save RDS file of ebird obs/sampl data
  tar_target(ebird_RDS, write_cityRDS(citylist = ebird_citydatalist, 
                                      ebird_cityfilters = ebird_citydata), 
             pattern = map(ebird_citydatalist), iteration = "vector",
             cue = tar_cue("never")),
  #filter ebird data by city polygon
  tar_target(ebird_filteredbycity, spatial_filter(citylist = ebird_citydatalist, 
                                                  ebirddat_RDS = ebird_RDS),
             cue = tar_cue("never")),
  
  #zerofill the data to get presence/absence
  tar_target(zf_citydata, zerofill_citylist(ebird_filtdata = ebird_filteredbycity), 
             pattern = map(ebird_filteredbycity), iteration = "vector"),
  #create a tidy data frame
  tar_target(city_tidyzfdf, zf_tidy(zf_citydata), 
             pattern = map(zf_citydata), iteration = "vector"),
  
  #### calculate biodiversity indices for every checklist, per city
  tar_target(div_bycity, calc_biodiversity(city_tidyzfdf), 
             pattern = map(city_tidyzfdf), iteration = "vector"),
  
  #### create a summary table for all the cities combined 
  tar_target(city_ebirdsums, citychkl_sumtab(city_zfdf = city_tidyzfdf), 
             pattern = map(city_tidyzfdf), iteration = "vector"),
  
  #### Create a summary dataframe to references files to #######################
  tar_target(cityfiles_summary, get_df_divLC(divdata = div_bycity, 
                                             chkldata = city_ebirdsums,
                                             LCrasterfiles = urban_ESA_LC,
                                             BVrasterfiles = ghsl_builtLC)),
  
  #extract the raster values in a buffer around each checklist in a city
  ## 
  tar_target(builtv_bycity500, extract_bv_bychkl(data = cityfiles_summary,
                                                 bufferradius = 500,
                                                 outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity1000, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 1000,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity1500, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 1500,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity2000, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 2000,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  
  #calculate the distance of each checklist to the city boundary
  tar_target(chkldist_citybounds, calc_urbandistance(citychkls = cityfiles_summary, 
                                                     urban_polygons = urbanareas,
                                                     outputdir = "./output/ebird_data/dist_urbbound_bycity/")),

  #### Calculate landscape metrics for each checklist, per city #################
  tar_target(cityfiles_csv, read_csv("./cityfiles_summary.csv")),
  ## create a subset of just 2 cities to test
  tar_target(city_testsubset, filter(cityfiles_csv, 
                                     cityname == "Winnipeg" | 
                                     cityname == "Ottawa")),
  # Calculate landscape metrics 
  tar_map(
    values = tar_read(city_testsubset),
    names = cityname,
    delimiter = "_",
    
    ## Get the data for the city
    tar_target(city_data, {
      list(
        cityplacename = cityname,
        citydf = readRDS(rds_div) %>% select(1:7),
        cityLCraster = LCpath
      )
    }),
    ## Get the raster data for the city
    tar_terra_rast(cityraster, rast(city_data$cityLCraster), 
                   cue = tar_cue(mode = "never")),
    ## Split into groups by checklists; this avoids RAM limitations issues
    tar_group_size(grouped_chkl, city_data$citydf, size = 100), 
    # Define buffer radii
    tar_target(buffer_radii, c(500, 1000, 1500, 2000)),
    
    ## Calculate pland metrics for multiple buffer radii
    tar_target(pland_results, {
      calc_pland_metric(cityname = city_data$cityplacename,
                        citydata = grouped_chkl,
                        landcover = cityraster,
                        bufferradius = buffer_radii,
                        targetcrs = "EPSG:3978")
    }, pattern = cross(buffer_radii, grouped_chkl)),
    
    ###### NEED TO CONVERT TO WIDE AND SAVE TO CSV ###
    # Calculate euclidean nearest neighbour distances for multiple buffer radii
    tar_target(enn_results, {
      calc_aggrm_metric(cityname = city_data$cityplacename,
                        citydata = grouped_chkl,
                        landcover = cityraster,
                        bufferradius = buffer_radii,
                        what_metrics = c("lsm_c_enn_mn", "lsm_c_enn_sd"),
                        targetcrs = "EPSG:3978")
    }, pattern = cross(buffer_radii, grouped_chkl)),
    
    ## Calculate aggregation metrics for multiple buffer radii
    tar_target(aggrm_l_results, {
      calc_aggrm_metric(cityname = city_data$cityplacename,
                        citydata = grouped_chkl,
                        landcover = cityraster,
                        bufferradius = buffer_radii,
                        what_metrics = c("lsm_l_ai", "lsm_l_contag", "lsm_l_iji", "lsm_l_mesh"),
                        targetcrs = "EPSG:3978")
    }, pattern = cross(buffer_radii, grouped_chkl)),
    tar_target(aggrm_c_results, {
      calc_aggrm_metric(cityname = city_data$cityplacename, 
                        citydata = grouped_chkl, 
                        landcover = cityraster,
                        bufferradius = buffer_radii, 
                        what_metrics = c("lsm_c_ai", "lsm_c_contag", "lsm_c_iji", "lsm_c_mesh"),
                        targetcrs = "EPSG:3978")
    }, pattern = cross(buffer_radii, grouped_chkl))
    ## Checkup landcover classes and convert results to wide format
    # tar_target(pland_results_wide, {
    #   lsm_pivotwide(pland_results)
    # }),
    
    # Write results to CSV
    # tar_target(citycsv_pland, {
    #   save_results_csv(pland_results_wide, city_data$cityplacename, 
    #                    whatmetric = "pland",
    #                    outputdir = "./output/ebird_data/pland_chkl_bycity/")
    # }, cue = tar_cue(file = TRUE))
  )
)