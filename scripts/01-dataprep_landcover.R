## HEADER---------------------------
## Script name: dataprep_landcover.R
##
## Purpose of script: match the ESA landcover data to the selected urban centres for extraction
##
## Author: Andrew Habrich
##
## Date Created: 2024-08-15
## Date last Modified: 2024-08-15
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

## 1. Load relevant packages--------
#install.packages("pacman")
pacman::p_load(tidyverse, ggpubr, RColorBrewer, tictoc,
               ggspatial, sf, terra, tidyterra, rnaturalearth, mapview,
               targets, tarchetypes)

## Import the urban centres ----
fuacan <- read_sf("./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua")
fua_vect <- vect(fuacan)

## 2. landcover data from Copernicus, European space agency ----
# Source: ESA WorldCover 10 m 2021 v200, downloadable here: https://zenodo.org/records/7254221
### 2.1. Identify ESA cells that overlap with urban centres ----
esagrid <- read_sf("./data/esa_worldcover_grid_composites.fgb")

# identify the cells that overlap with the urban centres
esalist_fua <- st_intersects(esagrid, fuacan, sparse = F) %>% 
  #extract all elements that are not 0
  apply(., 1, any) %>% 
  #coerce to dataframe and rename the . column to overlap
  as_tibble() %>%
  rename(urbanoverlap = value) %>%
  #add the tile name to the dataframe
  mutate(tile = esagrid$tile) %>%
  #filter to only the cells that overlap
  filter(urbanoverlap == T)

## A couple of steps required to import the raster files that match the macrotiles
### 2.2. identify macrotiles from esalist_fua ----
## separate the tile column into N and W extracting only the numbers
esalist_fua <- esalist_fua %>% 
  separate(tile, into = c("N", "W"), sep = c(3, 7), remove = F) %>% 
  mutate(N = as.numeric(str_extract(N, "\\d{2}")),
         W = as.numeric(str_extract(W, "\\d{3}"))) %>% 
  # Round the N value up and the W value down to the closest macrotile, which are 3x3
  mutate(Nmacrotile = floor(N/3)*3, 
         Wmacrotile = ceiling(W/3)*3) %>% 
  #create a new column for the macrotile name
  mutate(macrotile = paste0("N", str_pad(Nmacrotile, 2, pad = "0"), "W", 
                            str_pad(Wmacrotile, 3, pad = "0")))

## 3.0 Crop layers to urban areas ----
### 3.1 Identify a list of the relevant macrotile rasters and import ----
datadir <- "E:/DATA/ESA2021LANDCOVER/"
esa_tilelist <- list.files(datadir) %>% 
  str_subset("tif$") %>% #only the raster images
  as_tibble() %>% rename(filename = value) %>% 
  #extract the tile name from the filename, and add it to the dataframe
  mutate(macrotile = str_extract(filename, "(N\\d{2}W\\d{3})")) %>% 
  #filter to only the tiles that overlap with urban centres
  filter(macrotile %in% esalist_fua$macrotile)

# import a list of raster files that match the macrotiles 
esa_rasters <- esa_tilelist %>% 
  mutate(landcover = map(filename, ~rast(paste0(datadir, .x)))) %>% 
  pull(landcover) %>% 
  sprc()

### 4.0 Landsat annual EVI composite ------------------------------------------
outputdir <- "./output/landsat_evi/" #create the output folder if it doesnt exist
if(!dir.exists(outputdir)){
  dir.create(outputdir)
}

# Source: USGS Landsat 8 EVI annual composite, downloaded from Google Earth Engine
evi_rasters <- list.files("./data/LANDSAT8_EVI/", full.names = T) %>%
  as_tibble() %>% rename(path = value) %>%
  filter(str_detect(path, "tif$")) %>%
  mutate(filename = basename(path)) %>%
  mutate(ID = paste0("EVI_", str_extract(filename, "\\d{4}"))) %>%
  mutate(raster = map(path, ~ rast(.x)))

# create a vector of years to assess
years <- evi_rasters %>% pull(ID) %>% str_extract("\\d{4}") %>% unique()
## create a virtual raster dataset for each year of the the EVI data
for (i in 1:length(years)) {
  evi_year <- evi_rasters %>% filter(str_detect(ID, years[i]))
  vrt(
    x = evi_year$path,
    filename = paste0("./output/landsat_evi/evi_", years[i], "_tiles.vrt"),
    overwrite = T
  )
}

## Crop the raster to each of the urban centres for each of the years
for (i in 1:length(years)) {
  evi_year <- evi_rasters %>% filter(str_detect(ID, years[i]))
  vrtfile <- paste0("./output/landsat_evi/evi_", years[i], "_tiles.vrt")
  evi_year_raster <- rast(vrtfile)
  for (j in 1:nrow(fuacan)) {
    fua <- fuacan[j, ]
    fua_name <- fua$eFUA_name
    fua_vect <- vect(fua) %>% buffer(width = 5000)
    evi_fua <- crop(evi_year_raster, fua_vect, filename = paste0(outputdir, fua_name, "_EVI_", years[i], ".tif"),
                    overwrite = T)
  }
}

### 5.0 ECOSTRESS dataset: https://ecostress.jpl.nasa.gov/ ---- NOT USED
#https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/ecostress-overview/#ecostress-naming-conventions
### L3_Meteorological
dir_met <- "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/ECO_L3T_MET"
ecostress_files <- list.files(dir_met, full.names = T) %>%
  as_tibble() %>% rename(path = value) %>% mutate(filename = basename(path)) %>% 
  separate(filename, remove = F,
    into = c("product", "level", "type", "orbitnum", "sceneid",
             "tileid","acqdate","buildid","iterationnum","var"), sep = "_") %>%
  separate(acqdate, into = c("date", "time"), sep = "T") %>% #convert the date column into year, month, day
  mutate(
    year = str_sub(date, 1, 4),
    month = str_sub(date, 5, 6),
    day = str_sub(date, 7, 8)) %>% 
  mutate(raster = map(path, ~ rast(.x)))

## get a dataframe for each variable of interest
eco_ta <- ecostress_files %>% select(path, filename, tileid, year, time, sceneid, var) %>% filter(var == "Ta.tif")
eco_rh <- ecostress_files %>% select(path, filename, tileid, year, time, sceneid, var) %>% filter(var == "RH.tif")

## import the rasters into a spatrasterdatacollection, sprc() object
outcrs <- crs(fua_vect)
eco_ta_rasters <- sprc(eco_ta$path)
## project the rasters to the same crs as the fua_vect
eco_ta_rasters <- project(eco_ta_rasters, outcrs, method = "bilinear", threads = T, partial = T)

## mosaic the rasters
# mosaic(eco_ta_rasters, fun = "mean",
#        filename = "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/ecostress_airtemp_mosaicv2.tif",
#        wopt = list(verbose = T, progress = 5, todisk = T, steps = 50))
# ta <- rast("C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/ecostress_airtemp_mosaicv2.tif")

# eco_rh_rasters <- sprc(eco_rh$path)
# mosaic(eco_rh_rasters, fun = "mean",
#        filename = "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/ecostress_relhumid_mosaic.tif",
#        verbose = F, progress = 10)

### L2_LandSurfaceTemperature
dir_lst <- "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/ECO_L2T_LSTE"
ecolst_files <- list.files(dir_lst, full.names = T) %>%
  as_tibble() %>% rename(path = value) %>% mutate(filename = basename(path)) %>% 
  separate(filename, remove = F,
           into = c("product", "level", "type", "orbitnum", "sceneid",
                    "tileid","acqdate","buildid","iterationnum","var"), sep = "_") %>%
  separate(acqdate, into = c("date", "time"), sep = "T") %>% #convert the date column into year, month, day
  mutate(
    year = str_sub(date, 1, 4),
    month = str_sub(date, 5, 6),
    day = str_sub(date, 7, 8))

# Test LST for some months
lst <- ecolst_files
lst_sprc <- sprc(lst$path)
# mosaic(lst_sprc, fun = "mean", 
#        filename = "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/lst_mosaic.tif",
#        verbose = F, progress = 10)

# lst_vrt <- vrt(x = lst$path,
#       filename = "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/eco_lst2024_march.vrt",
#       overwrite = T)

# merge(lst_vrt, fun = "mean", 
#        filename = "C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/lst_mosaic.tif",
#        verbose = F, progress = 10)

# # convert the fua_vect to the same crs as the vrt_03 using terra package
# fua_v <- project(fua_vect, vrt, partial = T)
# march <- rast("C:/Users/andrewhabrich/OneDrive - Carleton University/GIS storage/ecostress/eco_lst2024_march.vrt")

# # Plot to test
# plot(vrt)
# plot(fua_v)
# plot(fua_v, add = T)

### 6.0 Yale UHI data ----
# find all the files that are tif
uhi_summer <- list.files("./data/YCEO_summeruhi/", full.names = T) %>% #convert to dataframe
  as_tibble() %>% rename(path = value) %>%
  #filter to just the tifs
  filter(str_detect(path, "tif$")) %>%
  #extract the filename
  mutate(filename = basename(path)) %>%
  #create an ID column of just UHI_ and the year number
  mutate(ID = paste0("UHI_", (str_extract(filename, "\\d{4}")))) %>%
  #filter to only the years that match the ebird data (2010-2022)
  filter(str_detect(filename, "201[0-9]")) %>%
  #import the rasters
  mutate(raster = map(path, ~rast(.x))) 
## save to a sprc object and name the layers
uhi_sprc <- sprc(uhi_summer$raster)
names(uhi_sprc) <- uhi_summer$ID

## Plot the first layer
uhi2010 <- uhi_sprc[1]
## crop to the urban centres
uhi2010_fua <- crop(uhi2010, fuacan)
# plot the daytime temperature
plot(fuacan$geom[1])
plot(crop(uhi2010_fua$Daytime, fuacan$geom[1]), add = T)

### Canada HRDEM layers ----
library(biomod2)
library(terra)

####### Check the HRDEM files ## WIP
canvec <- vect(fuacan)
crs(canvec)
dsmtile1 <- rast("./data/HRDEM/Digital_Surface_Model_(VRT).vrt"); mem_info(dsmtile1)
dsmtile1vrt <- vrt("./data/HRDEM/Digital_Surface_Model_(VRT).vrt"); mem_info(dsmtile1vrt)
rprjcanvec <- canvec %>% project(terra::crs(dsmtile1))

## Plot to check
plot(dsmtile1vrt)
plot(rprjcanvec, add = T)

crop1 <- crop(dsmtile1, rprjcanvec, mask = T)

dsmtile2 <- rast("./data/HRDEM/Digital_Surface_Model_(VRT) (2).vrt")
plot(dsmtile1)