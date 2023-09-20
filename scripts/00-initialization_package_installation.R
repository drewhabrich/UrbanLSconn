## HEADER---------------------------
## Script name: 00-initialization and setup for packages, dependant software, google account
##
## Purpose of script: to faciliate furhter data extraction and manipulation, this script should be run first if starting on a new computer.
##
## Author: Andrew Habrich
##
## Date Created: 2023-07-26
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------

## 1. Load relevant packages--------

#install.packages("remotes", "reticulate")
#remotes::install_github("r-spatial/rgee") #download the latest version from github
### INSTALLATION INSTRUCTIONS FROM https://r-spatial.github.io/rgee/articles/rgee01.html ###
# load packages -----------------------------------------------------------
library(rgee)

# initialization ----------------------------------------------------------
# Install rgee Python dependencies into a rgee virtual environment
### You need the latest earth engine API, and the numpy package installed. 
### NOTE; if using R studio, select the 'python interpreter' in the Tools > Global options > Python
### Use either anaconda or miniconda to create the 'environment' where the APIs will be installed.
### You also have to install GCLOUD SDK to interface with google credentials. 
### You ALSO have to activate your gmail account with earth engine: https://code.earthengine.google.com/register

ee_check() #make sure all python dependencies and credentials are valid
ee_check_python()
ee_check_credentials()
ee_check_python_packages()

# session management
ee_get_earthengine_path() #authorization token stored here

ee_Initialize(user = 'andrhabr@gmail.com', credentials = "persistent")
ee_version() #what version of ee API is being run?
ee_user_info()


# Test the installation ----------------------------------------------------
# If everything installed correctly, there should be a python environment with the earth engine and dependent packages installed and a a connection to google earth engine through your google account


## Download elevation data
srtm <- ee$Image("USGS/SRTMGL1_003") #download SRTM elevation data layer from NASA and USGS
## define the colour scheme for mapping viz
viz <- list( 
  max = 6500,
  min = -10,
  palette = c("#000000","#5AAD5A","#A9AD84","#FFFFFF")) #define viz colours
# Visualize data layer interactively (you can scroll and zoom)
Map$addLayer(eeObject = srtm, #object
  visParams =  viz, #viz parameters 
  name = 'SRTM') #name of layer on map

# let's try to download ESA 10m worldcover image set
lc10m <- ee$ImageCollection("ESA/WorldCover/v200")$first()
# Map the land cover map (not specifying vis parameters will use the default of the layer)
Map$addLayer(eeObject = lc10m, 
             name = 'LC')

# What if we estimated cumulative cost based on a LC surface
# A rectangle representing Bangui, Central African Republic.
geometry <- ee$Geometry$Rectangle(list(18.5229, 4.3491, 18.5833, 4.4066))

# Create a source image where the geometry is 1, everything else is 0.
sources <- ee$Image()$toByte()$paint(geometry, 1)

# Mask the sources image with itself.
sources <- sources$updateMask(sources)

# The cost data is generated from classes in ESA/GLOBCOVER.
cover <- ee$Image("ESA/GLOBCOVER_L4_200901_200912_V2_3")$select(0)

# Classes 60, 80, 110, 140 have cost 1.
# Classes 40, 90, 120, 130, 170 have cost 2.
# Classes 50, 70, 150, 160 have cost 3.
cost <- cover$eq(60)$Or(cover$eq(80))$Or(cover$eq(110))$Or(cover$eq(140))$
  multiply(1)$add(
    cover$eq(40)$Or(cover$eq(90))$Or(cover$eq(120))$Or(cover$eq(130))$
      Or(cover$eq(170))$
      multiply(2)$add(
        cover$eq(50)$Or(cover$eq(70))$Or(cover$eq(150))$Or(cover$eq(160))$
          multiply(3)
      )
  )

# Compute the cumulative cost to traverse the lAnd cover.
cumulativeCost <- cost$cumulativeCost(
  source = sources,
  maxDistance = 80 * 1000 # 80 kilometers
)

# Display the results
Map$setCenter(lon = 18.71, lat = 4.2)
Map$setZoom(zoom = 9)

Map$addLayer(
  eeObject = cover,
  visParams = list(),
  name = "Globcover"
) +
  Map$addLayer(
    eeObject = cumulativeCost,
    visParams = list(min = 0, max = 5e4),
    name = "accumulated cost"
  ) +
  Map$addLayer(
    eeObject = geometry,
    visParams = list(color = "FF0000"),
    name = "source geometry"
  )

### get resolution of map
naip <- ee$Image("USDA/NAIP/DOQQ/m_3712213_sw_10_1_20140613")
Map$setCenter(-122.466123, 37.769833, 17)
Map$addLayer(naip, list(bands = c("N", "R", "G")), "NAIP")

naip_resolution <- naip$select("N")$projection()$nominalScale()
cat("NAIP resolution: ", naip_resolution$getInfo())

landsat <- ee$Image("LANDSAT/LC08/C01/T1/LC08_044034_20140318")
landsat_resolution <- landsat$select("B1")$projection()$nominalScale()
cat("Landsat resolution: ", landsat_resolution$getInfo())