## Load required packages ##
library(targets)
library(tarchetypes)
library(tidyterra)
library(geotargets)
library(crew)
library(qs2)
library(landscapemetrics)
landscapemetrics::options_landscapemetrics(to_disk = T)

# This is where you write source(\"R/functions.R\")
tar_source("R/")

# define the controller and the log using crew and autometrics
controller <- crew_controller_local(
  name = "ctrl",
  reset_globals = T,
  garbage_collection = T,
  workers = 2,
  seconds_idle = 5
)

# Set target-specific options such as packages that are required:
tar_option_set(packages = c("targets", "tarchetypes", "geotargets","qs2","crew",
                            "ggspatial", "sf", "terra", "tidyterra", "rnaturalearth",
                            "auk", "vegan", "landscapemetrics", "exactextractr", #"gdistance", 
                            "tidyverse", "dplyr", "ggpubr", "RColorBrewer"),
               format = "qs", memory = "transient", error = "stop",
               deployment = "worker", storage = "main", retrieval = "main", 
               garbage_collection = 1,
               controller = controller, 
               resources = tar_resources(qs = tar_resources_qs())
)

# Global variables
buffer_radii <- c(500, 1000, 1500, 2000)

# List of target objects and the functions to apply to each of them; 
# NOTE this does not RUN the pipeline, only defines it.
list(
  #### Read in spatial data ####
  tar_terra_vect(urbanareas, get_urbanareas(country = "Canada", 
                                            minpopsize = 100000, 
                                            buffersize = 1000)), 
  tar_terra_rast(urban_esarasters, get_urbanESA_grids(urban_polygons = urbanareas, 
                                                      esa_datadir = "F:/DATA/ESA2021LANDCOVER/"),
                 cue = tar_cue(depend = TRUE, mode = "never")),
  
  tar_target(urban_ESA_LC, get_urbanESA_LC(rasterlist = urban_esarasters,
                                           outputdir = "./output/ESArasters_factor/",
                                           urban_polygons = urbanareas,
                                           targetcrs = "EPSG:3978")),
  tar_target(ghsl_builtLC, get_builtLCvolume(datadir = "./raw_data/ghsl/",
                                             outputdir = "./output/GHSL_builtV/",
                                             urban_polygons = urbanareas,
                                             targetcrs = "EPSG:3978")),
  
  #### Read in the ebird data ####
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
  
  #### Create a summary dataframe to references files to ####
  tar_target(cityfiles_summary, get_df_divLC(divdata = div_bycity, 
                                             chkldata = city_ebirdsums,
                                             LCrasterfiles = urban_ESA_LC,
                                             BVrasterfiles = ghsl_builtLC)),
  
  #extract the raster values in a buffer around each checklist in a city
  tar_target(builtv_bycity500, extract_bv_bychkl(data = cityfiles_summary,
                                                 bufferradius = 250,
                                                 outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity1000, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 500,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity2000, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 1000,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  tar_target(builtv_bycity3000, extract_bv_bychkl(data = cityfiles_summary,
                                                  bufferradius = 1500,
                                                  outputdir = "./output/ebird_data/builtv_chkl_bycity/"),
             pattern = map(cityfiles_summary), iteration = "vector"),
  
  #calculate the distance of each checklist to the city boundary
  tar_target(chkldist_citybounds, calc_urbandistance(citychkls = cityfiles_summary, 
                                                     urban_polygons = urbanareas,
                                                     outputdir = "./output/ebird_data/dist_urbbound_bycity/")),
  
  # #calculate landscape metrics for each checklist, per city
  # tar_target(lsm_bycity750, calc_landscape_metric(citydiv = div_bycity,
  #                                              rasterlist = urban_ESA_LC,
  #                                              bufferradius = 750), #0.75km diameter buffer
  #            pattern = map(div_bycity), iteration = "vector", memory = "transient"),
  # tar_target(lsm_bycity1500, calc_landscape_metric(citydiv = div_bycity,
  #                                              rasterlist = urban_ESA_LC,
  #                                              bufferradius = 1500), #3km diameter buffer
  #            pattern = map(div_bycity), iteration = "vector", memory = "transient"),
  # tar_target(lsm_bycity2500, calc_landscape_metric(citydiv = div_bycity,
  #                                              rasterlist = urban_ESA_LC,
  #                                              bufferradius = 2500), #3km diameter buffer
  #            pattern = map(div_bycity), iteration = "vector", memory = "transient"),
  # 
  # #extract the aggregation metrics in a buffer around checklists
  # tar_target(agg_bycity750, sample_aggr_lsm(chkldata = div_bycity,
  #                                           rasterlist = urban_ESA_LC,
  #                                           bufferradius = 750,
  #                                           crs = "EPSG:3978",
  #                                           what_metrics = c("lsm_l_ai", "lsm_l_contag", 
  #                                                            "lsm_l_iji", "lsm_l_mesh", "lsm_l_pladj")),
  #            pattern = map(div_bycity), iteration = "vector", memory = "transient"),
  # ### 3km buffer
  # tar_target(ai_bycity1500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                           bufferradius = 1500, crs = "EPSG:3978",
  #                                           what_metrics = c("lsm_l_ai")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(iji_bycity1500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                           bufferradius = 1500, crs = "EPSG:3978",
  #                                           what_metrics = c("lsm_l_iji")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(contag_bycity1500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                            bufferradius = 1500, crs = "EPSG:3978",
  #                                            what_metrics = c("lsm_l_contag")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(mesh_bycity1500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                               bufferradius = 1500, crs = "EPSG:3978",
  #                                               what_metrics = c("lsm_l_mesh")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # ### 5km buffer
  # tar_target(ai_bycity2500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                           bufferradius = 2500, crs = "EPSG:3978",
  #                                           what_metrics = c("lsm_l_ai")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(iji_bycity2500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                            bufferradius = 2500, crs = "EPSG:3978",
  #                                            what_metrics = c("lsm_l_iji")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(contag_bycity2500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                               bufferradius = 2500, crs = "EPSG:3978",
  #                                               what_metrics = c("lsm_l_contag")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient"),
  # tar_target(mesh_bycity2500, sample_aggr_lsm(chkldata = div_bycity, rasterlist = urban_ESA_LC,
  #                                             bufferradius = 2500, crs = "EPSG:3978",
  #                                             what_metrics = c("lsm_l_mesh")),
  #            pattern = map(div_bycity), iteration = "vector", error = "null", memory = "transient")
  tar_target(city_testsubset, filter(cityfiles_summary, cityname == "Winnipeg" | 
                                                        cityname ==  "London")),
  tar_map(
    values = tar_read(city_testsubset),
    names = cityname,
    delimiter = "_",
    tar_target(city_data, {
      list(
        cityplacename = cityname,
        citydf = readRDS(rds_div) %>% select(1:7),
        cityLCraster = LCpath
      )
    }),
    tar_terra_rast(cityraster, rast(city_data$cityLCraster)),
    tar_group_size(grouped_chkl, city_data$citydf, size = 1000), #split into groups by # of checklists
    
    ## calculate pland metrics
    tar_target(pland500_results, {
      calc_pland_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster, 
                        bufferradius = 500, outputdir = "./output/ebird_data/pland_chkl_bycity/", targetcrs = "EPSG:3978")
    }, pattern = map(grouped_chkl)),
    tar_target(pland1000_results, {
      calc_pland_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster,
                        bufferradius = 1000, outputdir = "./output/ebird_data/pland_chkl_bycity/", targetcrs = "EPSG:3978")
    }, pattern = map(grouped_chkl)),
    tar_target(pland1500_results, {
      calc_pland_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster,
                        bufferradius = 1500, outputdir = "./output/ebird_data/pland_chkl_bycity/", targetcrs = "EPSG:3978")
    }, pattern = map(grouped_chkl)),
    tar_target(pland2000_results, {
      calc_pland_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster,
                        bufferradius = 2000, outputdir = "./output/ebird_data/pland_chkl_bycity/", targetcrs = "EPSG:3978")
    }, pattern = map(grouped_chkl)),
    #gather results
    tar_target(pland_results, {
      list(
        pland500 = pland500_results,
        pland1000 = pland100_results,
        pland1500 = pland1500_results,
        pland2000 = pland2000_results
      )
    }),
    #write results to csv
    # write results to csv
    tar_target(citycsv_pland, {
      save_results_csv(pland_results, city_data$cityplacename, outputdir = "./output/VM_pland/")},
      pattern = map(pland_results)
    
    ## calculate aggregation metrics
    # tar_target(aggrm500_results, {
    #   calc_aggrm_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster, 
    #                     bufferradius = 500, outputdir = "./output/ebird_data/aggrm_chkl_bycity/", targetcrs = "EPSG:3978")
    # }, pattern = map(grouped_chkl)),
    # tar_target(aggrm1000_results, {
    #   calc_aggrm_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster, 
    #                     bufferradius = 1000, outputdir = "./output/ebird_data/aggrm_chkl_bycity/", targetcrs = "EPSG:3978")
    # }, pattern = map(grouped_chkl)),
    # tar_target(aggrm1500_results, {
    #   calc_aggrm_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster, 
    #                     bufferradius = 1500, outputdir = "./output/ebird_data/aggrm_chkl_bycity/", targetcrs = "EPSG:3978")
    # }, pattern = map(grouped_chkl)),
    # tar_target(aggrm2000_results, {
    #   calc_aggrm_metric(cityname = city_data$cityplacename, citydata = grouped_chkl, landcover = cityraster, 
    #                     bufferradius = 2000, outputdir = "./output/ebird_data/aggrm_chkl_bycity/", targetcrs = "EPSG:3978")
    # }, pattern = map(grouped_chkl))
   )
)