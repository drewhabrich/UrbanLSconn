library(targets)
library(tarchetypes)
library(geotargets)

#run tar_make() to run the pipeline
#run tar_manifest() to see the steps in the pipeline
#run tar_visnetwork() to see the visual input-output steps for the pipeline
#tar_read() to view the results for each target.

# This is where you write source(\"R/functions.R\")
tar_source("R/")

# Set target-specific options such as packages that are required:
tar_option_set(packages = c("tidyverse", "ggpubr", "RColorBrewer","qs", 
                            "ggspatial", "sf", "terra", "tidyterra", "rnaturalearth", "mapview",
                            "auk", "vegan", "landscapemetrics", "exactextractr",
                            "targets", "tarchetypes", "geotargets"),
               format = "qs",
               resources = tar_resources(qs = tar_resources_qs(preset = "high"))
)

# List of target objects and the functions to apply to each of them; 
# NOTE this does not RUN the pipeline, only defines it.
list(
  ####read in spatial data
  tar_terra_vect(urbanareas, get_urbanareas(country = "Canada", 
                                            minpopsize = 100000, 
                                            buffersize = 1000)), 
  tar_terra_rast(urban_esarasters, get_urbanESA_grids(urban_polygons = urbanareas, 
                                                      esa_datadir = "F:/DATA/ESA2021LANDCOVER/"),
                 format = "rds"),
  tar_target(urban_ESA_LC, get_urbanESA_LC(rasterlist = urban_esarasters,
                                           outputdir = "/output/ESArasters/",
                                           urban_polygons = urbanareas),
                 format = "rds"),
  tar_target(ghsl_builtLC, get_builtLCvolume(datadir = "/raw_data/ghsl/",
                                             outputdir = "/output/GHSL_builtV/",
                                             urban_polygons = urbanareas),
                 format = "rds"),
  
  ####read in the ebird data
  tar_target(ebd_can, auk_ebd(file = "./raw_data/ebd_can_obsdata.txt", #this is the .txt with the data
                              file_sampling = "./raw_data/ebd_can_sampdata.txt"),
             cue = tar_cue("never")),
  ####create and check list of cities with ebird data
  tar_target(ebird_citydatalist, checkdata_by_city(urbanareas),
             cue = tar_cue("always")), 
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
  #### create a summary table for all the cities combined 
  tar_target(city_ebirdsums, city_summary_table(city_tidyzfdf), 
             pattern = map(city_tidyzfdf), iteration = "vector"),
  #calculate biodiversity indices for every checklist, per city
  tar_target(div_bycity, calc_biodiversity(city_tidyzfdf), pattern = map(city_tidyzfdf), iteration = "vector"),
  
  #calculate landscape metrics for each checklist, per city
  tar_target(lsm_bycity750, calc_landscape_metric(citydiv = div_bycity,
                                               rasterlist = urban_ESA_LC,
                                               bufferradius = 750, #1.5km diameter buffer
                                               crs = "EPSG:3978"), #Canada Atlas Lambert
             pattern = map(div_bycity), iteration = "vector"),
  tar_target(lsm_bycity1500, calc_landscape_metric(citydiv = div_bycity,
                                               rasterlist = urban_ESA_LC,
                                               bufferradius = 1500, #3km diameter buffer
                                               crs = "EPSG:3978"), #Canada Atlas Lambert
             pattern = map(div_bycity), iteration = "vector"),
  tar_target(lsm_bycity2500, calc_landscape_metric(citydiv = div_bycity,
                                               rasterlist = urban_ESA_LC,
                                               bufferradius = 2500, #3km diameter buffer
                                               crs = "EPSG:3978"), #Canada Atlas Lambert
             pattern = map(div_bycity), iteration = "vector"),
  
  #extract the raster values in a buffer around each checklist in a city
  tar_target(builtv_bycity, extract_raster_inpolygon(chkldata = div_bycity,
                                                     rasterlist = ghsl_builtLC,
                                                     bufferradius = 1500),
             pattern = map(div_bycity), iteration = "vector")
)
