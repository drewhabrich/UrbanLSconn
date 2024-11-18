## HEADER---------------------------
## Script name: 00-initialization and setup for packages, dependant software, google account
##
## Purpose of script: to faciliate further data extraction and manipulation, this script should be run first if starting on a new computer. Be cautioned though, it can be a pain in the ass to get working properly.
##
## Author: Andrew Habrich
##
## Date Created: 2023-07-26
## Date Updated: 2024-08-30
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

# 1. Load relevant packages--------
#remotes::install_github("r-spatial/rgee") #download the latest version from github
#remotes::install_github("r-earthengine/rgeeExtra")
### INSTALLATION INSTRUCTIONS FROM https://r-spatial.github.io/rgee/articles/rgee01.html ###
### load packages -----------------------------------------------------------
pacman::p_load(tidyverse, ggpubr, RColorBrewer, 
               ggspatial, sf, terra, rnaturalearth, mapview, geojsonio,
               auk, targets, tarchetypes,
               rgee, rgeeExtra, reticulate)

### rgee initialization ----------------------------------------------------------
# Install rgee Python dependencies into a rgee virtual environment
### You need the latest earth engine API, and the numpy package installed. 
### NOTE; if using R studio, select the 'python interpreter' in the Tools > Global options > Python
### Use either anaconda or miniconda to create the 'environment' where the APIs will be installed.

### You also have to install GCLOUD SDK to interface with google credentials. 
#### https://cloud.google.com/sdk/docs/install?hl=en

### You ALSO have to activate your gmail account with earth engine: 
#### https://code.earthengine.google.com/register

Sys.getenv()
## install miniconda if necessary
#reticulate::install_python("3.12:latest")

# Get the username
HOME <- Sys.getenv("HOME")
# 1. Install miniconda
#reticulate::install_miniconda(update = T)

# 2. Install Google Cloud SDK
#https://cloud.google.com/sdk/docs/install#windows

# 3 Set global parameters
#Sys.setenv("RETICULATE_PYTHON" = sprintf("%s/.local/share/r-miniconda/bin/python3", HOME))
#Sys.setenv("EARTHENGINE_GCLOUD" = sprintf("%s/google-cloud-sdk/bin/", HOME))

# Sys.setenv(RETICULATE_PYTHON = "C:/Users/andrh/anaconda3/python.exe")
#rgee::ee_install()
#rgee::ee_install_upgrade(version = "0.1.370")
#rgee::ee_clean_pyenv(Renviron = "global")
rgee::ee_check()
## RUN EE_INSTALL() if you are unfamiliar with python

## Check the python configuration
reticulate::py_available()
reticulate::py_config()
reticulate::py_discover_config()

###  Check installation ----------------------------------------------------
# If everything installed correctly, there should be a python environment with the earth engine and dependent packages installed and a a connection to google earth engine through your google account
# reticulate::py_set_item(name = 'ee', value = "C:/Users/andrh/anaconda3/Lib/site-packages/ee")
#make sure all python dependencies and credentials are valid
rgee::ee_check-tools()
rgee::ee_check_python()
rgee::ee_check_python_packages()

rgee::ee_version() #what version of ee API is being run?
ee_user_info() #authorization token stored here
ee_clean_user_credentials()

# 2. Earth Engine API initialization -----------------------------------------
# session management
rgee::ee_Authenticate(user = 'andrhabr@gmail.com')
rgee::ee_Initialize(user = 'andrhabr@gmail.com', credentials = "persistent")
rgeeExtra:: extra_Initialize()

############# Let's load the GHSL functional urban areas dataset (feature collection)
Map$setCenter(-96.328, 49.382, 4)

# Load the GHSL FUA dataset
fua <- ee$FeatureCollection("projects/ee-ahabrich/assets/ghsl_fua")
# filter FUA to just Canadian cities using the Cntry_name property, and population greater than 100,000
fua_can <- fua$filter(ee$Filter$eq("Cntry_name", "Canada"))$filter(ee$Filter$gt("FUA_p_2015", 100000))
# Buffer each city by 5km and display
fua_can <- fua_can$map(function(f) f$buffer(5000))
Map$addLayer(eeObject = fua_can, name = 'Urban CAN')
# Save the bounds of all cities 
can_bounds <- fua_can$geometry()$bounds()
Map$addLayer(eeObject = can_bounds, name = 'CAN bounds')


############### Let's try to extract the geometry of each city
# What properties can we filter by?
fua_can$first()$propertyNames()$getInfo()
names(fua_can$first()$propertyNames())
## Get a list of cities included
citylist <- fua_can$aggregate_array("eFUA_name")$getInfo()
## write a for loop to extract the geometry bounds of each city and save it in a list with the city name
city_geom <- list()
for (i in 1:length(citylist)){
  city_geom[[i]] <- fua_can[[i]]$geometry()$bounds()
}
# name the list after the city
names(city_geom) <- citylist

## check the single city
Map$addLayer(eeObject = city_geom$Victoria, name = 'Urban CAN') +
  Map$addLayer(eeObject = fua_can[[1]]$geometry(), name = 'Urban CAN')

############# Let's try to download ESA 10m worldcover image set
lc10m <- ee$ImageCollection("ESA/WorldCover/v200")[[1]]
ee_print(lc10m)
names(lc10m)

# Map the land cover map (not specifying vis parameters will use the default of the layer)
Map$addLayer(eeObject = lc10m, name = 'LC')

# Take the first image from the collection and crop it to the extent of the Canadian cities
lc10m_can <- lc10m$clipToCollection(fua_can)$updateMask(lc10m)
Map$addLayer(eeObject = lc10m_can, name = 'LC CAN')
ee_print(lc10m_can)

lc_bounds <- lc10m$clipToBoundsAndScale(geometry = fua_can$geometry(), scale = 10)$updateMask(lc10m)
Map$addLayer(eeObject = lc_bounds, name = 'LC CAN')