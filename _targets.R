library(targets)
library(tarchetypes)

#run tar_make() to run the pipeline
#run tar_manifest() to see the steps in the pipeline
#run tar_visnetwork() to see the visual input-output steps for the pipeline
#tar_read() to view the results for each target.

# This is where you write source(\"R/functions.R\")
tar_source("R/")

# Set target-specific options such as packages that are required:
tar_option_set(packages = c("tidyverse", "dplyr", "ggplot2", "ggpubr","viridis",
                            "rnaturalearth", "terra", "sf", "stars",
                            "auk", "vegan"))

# List of target objects and the functions to apply to each of them; NOTE this does not RUN the pipeline, only defines it.
list(
  tar_target(NA_basemaps, get_naturalearth("ne_basemaps", "North America", "CAN", c("USA", "CAN"))), 
  tar_target(urb_popcentres, get_populationcentres("lpc_000b21a_e.shp", "POPCTR.csv", "EPSG:4326", 100000)),
  tar_target(urb_buffers, st_buffer(urb_popcentres, dist = 5000)), #5km buffer around each urban area

  tar_target(plot_cityarea, plot_urbanarea("Montréal", zoomlevel = 5, basemaplist = NA_basemaps, citylist = urb_popcentres)),
  tar_target(plot_buffarea, plot_urbanarea("Montréal", zoomlevel = 5, basemaplist = NA_basemaps, citylist = urb_buffers)),
  
  tar_target(ebd_can, auk_ebd(file = "./raw_data/ebd_can_obsdata.txt", #this is the .txt with the data
                              file_sampling = "./raw_data/ebd_can_sampdata.txt")),
  tar_target(ebird_citylist, ebird_by_city(urb_buffers)), #create and check list of cities with ebird data
  tar_target(ebird_cityfilters, generate_city_filters(ebird_citylist, ebd_can)), #create ebird filters for cities
  tar_target(data_bycity, ebird_data_bycity(ebird_cityfilters, ebird_citylist)),
  
  tar_target(zf_bycity, auk_zerofill(data_bycity[[1]], collapse = T), pattern = map(data_bycity), iteration = "list"),
  tar_target(zf_checklists, zf_tidy(zf_bycity), pattern = map(zf_bycity), iteration = "list"),
  
  tar_target(table_cityebirdsums, city_summary_table(zf_checklists, urb_popcentres)),
  tar_target(plot_cityebirdsums, city_summary_plot(table_cityebirdsums)),
  
  tar_target(div_bycity, calc_diversity(zf_checklists), pattern = map(zf_checklists), iteration = "list"),
  tar_target(plot_divbycity, biodiv_plot(div_bycity,
                                         "Montréal", zoomlevel = 8, basemaplist = NA_basemaps, citylist = urb_popcentres))
)
