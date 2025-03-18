## HEADER---------------------------
## Script name: 03-dataprep_habitatsuitability
##
## Purpose of script: Create habitat suitability layer for bird data
##
## Author: Andrew Habrich
##
## Date Created: 2025-02-18
## Date last Modified: 2025-02-18
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
rm(list = ls())

## 1. Load relevant packages--------
pacman::p_load(sf, terra, tidyterra, rnaturalearth, mapview,
               vegan, landscapemetrics, exactextractr, gdistance,
               tidyverse, ggpubr, RColorBrewer, easystats, tictoc,
               dismo, biomod2, rgbif)

## Import the relevant files
citylist <- read_csv("./cityfiles_summary.csv")
wpgdata <- citylist %>% filter(cityname == "Winnipeg") 
wpg <- readRDS(wpgdata$rds_div) 

### GET CONN METRICS AND RUN MODEL FOR THE CITY
### GET GRAPH BASED, POTENTIAL, FUNCTIONAL, AND STRUCTURAL CONNECTIVITY METRICS

#### GRAPH BASED CONNECTIVITY #############################################
# library(devtools)
# library(remotes)
# remotes::install_github("connectscape/Makurhini", dependencies = TRUE, upgrade = "never")
# library(Makurhini)
library(grainscape)
library(terra)
library(tidyterra)
library(gdistance)

## Landcover code lookup
lcclass_lookup <- tibble(value = c(10, 20, 30, 40, 
                                   50, 60, 70, 80, 
                                   90, 95, 100),
                         category = c("treecover", "shrubland", "grassland", "cropland",
                                      "builtup", "bare", "snow", "water", 
                                      "wetland", "mangrove", "mosslichen"))

lcbinary <- tibble(value = c(10, 20, 30, 40, 
                             50, 60, 70, 80, 
                             90, 95, 100),
                   category = c(1, 1, 1, 0, 
                                0, 0, 0, 0, 
                                0, 0, 0))
lc_veg <- cbind(c(10, 20, 30, 40, 
                  50, 60, 70, 80, 
                  90, 95, 100), 
                c(1, 1, 1, 0, 
                  0, 0, 0, 0, 
                  0, 0, 0))

## Load the relevant targets files
tar_load(cityfiles_summary)
tar_load(urbanareas)

## Load just 1 city to test 
city <- "Winnipeg"
#read in raster LandCover for the city
city_lctile <- cityfiles_summary %>% filter(cityname == city) %>% pull(LCpath) %>% rast()
## read in and project to the raster format
citypoly <- urbanareas %>% 
  filter(eFUA_name == city) %>% 
  project(y = city_lctile)
## read in the city data
citydata <- cityfiles_summary %>% 
  filter(cityname == city) %>% 
  pull(rds_div) %>% readRDS()
city_ploc <- vect(citydata, geom = c("longitude", "latitude"), 
                  crs = "EPSG:4326") %>% 
  project(y = city_lctile)
## Crop the raster to the city polygon
city_lc <- crop(city_lctile, citypoly) %>% mask(citypoly)
is.factor(city_lc)

## get patches of landcover using landscapemetrics
#lsm_patches <- get_patches(city_lc, class = "all", directions = 8, to_disk = T)

## Tidy up the environment
rm(urbanareas, city_lctile, citydata, cityfiles_summary)
#plot(city_lc); plot(city_ploc, add = T)

## coerce to raster for grainscape
rasterlc <- raster(city_lc)
rastercost <- reclassify(rasterlc, rcl = lc_veg)
# rastercostdf <- ggGS(rastercost)
# is.factor(rastercostdf$value)
# rastercostdf$value <- as.factor(rastercostdf$value)
## plot the raster to see how it worked
# ggplot() +
#   geom_raster(
#     data = rastercostdf,
#     aes(x = x, y = y, fill = value)
#   ) +
#   scale_fill_brewer(
#     type = "qual", palette = "Paired",
#     direction = -1, guide = "legend"
#   ) +
#   guides(fill = guide_legend(title = "Resistance")) +
#   theme_grainscape() +
#   theme(legend.position = "right")

## filter the patches to a minimum size to avoid overcomplication
patches <- patchFilter(x = rastercost, area = 2000, directions = 8)
plot(patches)
euclcost <- patches; euclcost[] <- 1
## Get the MPG for the patches of landcover
patchMPG <- MPG(euclcost, patch = (patches == 1))
## Save the MPG to a rds file
#saveRDS(patchMPG, file = "./output/mpg_Winnipeg.rds")
## Test for threshold values in the network (euclidean links)
scalarAnalysis <- threshold(patchMPG, nThresh = 20)
## Quick plot to vis the number of components by the # of thresholds
ggplot(scalarAnalysis$summary, aes(x = maxLink, y = nComponents)) +
  geom_line(colour = "forestgreen") +
  xlab("Link Threshold (resistance units)") +
  ylab("Number of components") +
  theme_light() +
  theme(axis.title = element_text())

# grainscape::export(
#   mpg,
#   dirname = "mpg_Winnipeg",
#   path = "./output/",
#   rasterFormat = "GTiff",
#   overwrite = FALSE,
#   R = FALSE,
#   vorBound = FALSE)

## Grains of connectivity modelling
wpg_goc <- GOC(patchMPG, nThresh = 20) #took 13688 secs
## Save the GOC to a rds file
#saveRDS(wpg_goc, file = "./output/goc_Winnipeg.rds")
plot(grain(wpg_goc, whichThresh = 1), quick = "grainPlot", theme = FALSE)
wpg_goc <- readRDS("./output/goc_Winnipeg.rds")

## check the grains
grain_obj <- grain(wpg_goc, whichThresh = 5)

##### COMMUNITY COMPOSITION ANALYSIS #####################################
### Calculate the geographic distance between all points in the city
city_ploc <- vect(citydata, geom = c("longitude", "latitude"), crs = "EPSG:4326")
dist_mat <- terra::distance(x = city_ploc, unit="m")

### get a matrix for each site to calculate beta-diversity from
sitematrix <- div_bycity %>% filter(cityname == cn) %>% pull(rds_div[1]) %>% readRDS()
## get the site names
site_id <- sitematrix$checklist_id
## remove the columns that are not species observations
sitematrix <- sitematrix %>% select(-c(1:7))
## convert to a matrix
sitematrix <- as.matrix(sitematrix)
# get the site names
rownames(sitematrix) <- site_id
# remove the columns that are not species observations
site_jac <- betadiver(sitematrix, method = "j", order = F)