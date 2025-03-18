## HEADER---------------------------
## Script name: 03-dataprep_gbif_occurrences
##
## Purpose of script: Get gbif data for bird occurrences in urban areas to parameterize the habitat suitability models
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
pacman::p_load(sf, terra, tidyterra,
               tidyverse, ggpubr,
               rgbif, spocc, CoordinateCleaner)

## Get the fields that can be used for gbif downloads
occ_download_describe("simpleCsv")$fields
occ_download_list() #list the files that have been downloaded
occ_download_list()$results$key #list the download keys
## download the data for all bird species in canada ########
gbif_birds_dl <- occ_download(
  pred("taxonKey", 212), #212 is code for birds in gbif
  pred("gadm", "CAN"), # get data from Canada
  pred_and(pred_gt("year", 2012), 
           pred_lte("year", 2022)
           ), # only get data between 2012-2022
  pred("occurrenceStatus", "PRESENT"),
  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN","LIVING_SPECIMEN",
                                     "MACHINE_OBSERVATION", "PRESERVED_SPECIMEN",
                                     "MATERIAL_CITATION", "MATERIAL_SAMPLE"))),
  pred_not(pred("collectionCode", "EBIRD")), #NOT ebird, to avoid pseudoreplication
  pred_not(pred("datasetKey", "4fa7b334-ce0d-4e88-aaae-2e0c138d049e")), #ebird datasetkey
  pred("hasCoordinate", TRUE), 
  pred("hasGeospatialIssue", FALSE), # remove GBIF default geospatial issues
  format = "SIMPLE_CSV"
) 
occ_download_wait(gbif_birds_dl) # wait for the download to be ready

## Tidy the raw gbif data
gbifdata <- occ_download_get("0002421-250218110819086", path = "./raw_data/") %>% 
            #occ_download_get(gbif_birds_dl, path = "./raw_data/") %>% 
            occ_download_import() 
# What are the columns in the dataset?
glimpse(gbifdata)

# Select only relevant records
gbifmin <- gbifdata %>% select(
  c(gbifID, species, verbatimScientificName, 
    decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, 
    genus, family, order, family, individualCount, countryCode, eventDate, 
    basisOfRecord, institutionCode, collectionCode))

# Tidying
gbifmin <- gbifmin %>% 
  # For the rows with species == "", fill the species column with verbatimScientificName
  mutate(species = ifelse(species == "", verbatimScientificName, species)) %>% 
  # Remove any remaining empty species
  filter(species != "")

## Clean the data using CoordinateCleaner
#flag problems
gbifdat <- data.frame(gbifmin)
glimpse(gbifdat)
names(gbifdat)[4:5] <- c("decimalLongitude", "decimalLatitude")

gbif_clean <- gbifdat %>%
  cc_val() %>%
  cc_equ() %>%
  cc_cen() %>%
  cc_sea() %>%
  cc_zero() %>%
  cc_outl() %>%
  cc_dupl()

## save the tidied dataset to csv
write_csv(gbif_clean, "./data/can_gbif_avesoccurrences.csv")
