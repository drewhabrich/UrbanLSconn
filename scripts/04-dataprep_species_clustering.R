## HEADER---------------------------
## Script name: 04-dataprep_speciesclustering
##
## Purpose of script: Merge gbif data with AVONET2022 traits to cluster into groups for habitat suitability modelling
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
               tidyverse, ggpubr, readxl, DataExplorer, tictoc,
               rgbif, spocc, CoordinateCleaner,
               vegan, FD)

## Read in the data
# Urban area data
fuacan <- read_sf("./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua")

# AVONET trait data
avo_traits <- read_excel("data/AVONET2022/AVONET Supplementary dataset 1.xlsx",
  sheet = "AVONET2_eBird")

# HWI trait data
# hwi_traits <- read_excel("data/AVIAN_HWI2020/Dataset HWI 2020-04-10.xlsx",
#   sheet = "speciesdata") %>% 
#   #remove issues in the column names
#   rename_with(~str_replace_all(., " ", "_")) %>% 
#   rename_with(~str_replace_all(., "-", "_"))

# gbif occurrences
gbif_occ <- read_csv("data/can_gbif_avesoccurrences.csv") %>% st_as_sf(
  coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

## 2. Data preparation----------------
# filter gbif_occ data to occurrences that are within the urban area polygons
gbif_occ_fua <- gbif_occ %>% 
  st_filter(fuacan)

# Join the gbif_occ_fua data with the avo_traits and hiw_traits data
gbif_occ_traits <- gbif_occ_fua %>% 
  left_join(avo_traits, by = c("species" = "Species2")) %>% 
  select(-c(individualCount))
  #left_join(hwi_traits, by = c("species" = "Species_name")) %>% 
  #select(-c(Notes, Synonym))
glimpse(gbif_occ_traits)

## Remove columns that are not needed
gbif_clean <- gbif_occ_traits %>% 
  select(-c(coordinateUncertaintyInMeters, countryCode, basisOfRecord, institutionCode, collectionCode, 
            Avibase.ID2, Order2, Family2,
            Total.individuals, Female, Male, Unknown, Complete.measures, 
            Mass.Source, Mass.Refs.Other, Inference, Traits.inferred, Reference.species 
            #IUCN_name, Tree_name, `Body_mass_(log)`, Order, ,
            #Sample_size, Species_ID, Latitude, Island, Diet,
  )) %>% 
  #coerce to numeric
  #mutate(across(AnnualTemp:PrecipRange, as.numeric)) %>% 
  ## for habitat.density replace 1 with dense, 2 with semiopen, and 3 with open
  mutate(Habitat.Density = case_when(
    Habitat.Density == 1 ~ "dense",
    Habitat.Density == 2 ~ "semi-open",
    Habitat.Density == 3 ~ "open")) %>% 
  ## for migration replace 1 tih sedentary, 2 with partial, and 3 with migratory
  mutate(Migration = case_when(
    Migration == 1 ~ "sedentary",
    Migration == 2 ~ "partial",
    Migration == 3 ~ "migratory")) %>%
  rename(HWI = `Hand-Wing.Index`)
glimpse(gbif_clean)

## Extract all the unique species and their traits in the dataset
unique_species_df <- gbif_clean %>% select(-c(eventDate, geometry)) %>% st_drop_geometry() %>% 
  distinct(species, .keep_all = TRUE)

### 2.1. Var. Correlation check ----------------
## Check for correlations between continuous variables
cor <- correlation::correlation(
  unique_species_df  %>%
    select(c(Beak.Length_Culmen, Beak.Width, Beak.Length_Nares, Beak.Depth, 
             Wing.Length, Secondary1, Kipps.Distance, HWI,
             Tarsus.Length, Tail.Length, Mass)),  na.rm = TRUE)

correlation::cor_sort(as.matrix(cor, distance = "correlation", hclust_method = "complete")) %>%
  correlation::visualisation_recipe() %>% 
  plot() + 
#angle the x-axis labels by 45 degrees
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
#Change the scale fill gradient 
  scale_fill_gradient2(low = "#673AB7", mid = "white", high = "#FF5722", midpoint = 0.5)

## Remove highly correlated morphological traits
spp_traits <- unique_species_df %>%
  select(-c(Beak.Length_Culmen, Beak.Depth,Secondary1, Kipps.Distance, Wing.Length)) %>% 
  column_to_rownames(var = "species") %>% 
  select(-c(gbifID:order))

## Remove rows that are NA
spp_traits <- spp_traits[complete.cases(spp_traits), ]

## 3. Species Trait clustering ----------------
## What are the traits?
colnames(spp_traits)
## Assign weights to each trait, higher values are more important
traitweights <- c(0.05, 0.05, 0.05, 0.15, 0.05, 0.10, 
                  0.125, 0.075, 0.05, 0.05, 0.125, 0.125)
data.frame(trait = colnames(spp_traits), weight = traitweights)
sum(traitweights) #should add up to 1

## Get the combined gower distance matrix with the weights
trait_gowdismat <- sqrt(as.matrix(FD::gowdis(spp_traits, w = traitweights))) 

## Create a distance matrix for each of the traits
# Write a for loop to calculate the distance matrix for each trait column
indtrait_gowdismat <- list()
for (i in 1:ncol(spp_traits)) {
  ind.matrix <- sqrt(as.matrix(FD::gowdis(spp_traits[i])))
  indtrait_gowdismat[[i]] <- ind.matrix
  #name the matrix after the column
  names(indtrait_gowdismat)[i] <- colnames(spp_traits)[i]
}

#Save the distance matrix to rds
saveRDS(trait_gowdismat, "./output/species_clustering/trait_gowdismatrix.rds")
saveRDS(indtrait_gowdismat, "./output/species_clustering/indtrait_gowdismatrix.rds")

### Combine individual matrices and weight them according to importance (weights need to add to 1)
# indtrait_gowdismat <- readRDS("./output/species_clustering/indtrait_gowdismatrix.rds")
## Apply weights to the individual trait matrices
weighted_matrices <- mapply(function(mat, w) mat * w, indtrait_gowdismat, traitweights, SIMPLIFY = FALSE)
comb_weightedmat <- Reduce("+", weighted_matrices)
saveRDS(comb_weightedmat, "./output/species_clustering/indtrait_combined_gowdismatrix.rds")

### 3.1 Cluster species ------------
## Measure of similarity between input matrix distance
## and the one obtained with the clustering (clust.DIST), should be MINIMIZED
## See Mouchet et al. 2008 for rationale
source("./R/04-functions_speciesclustering.R") #import the individual function 

## Get the measures and the dendrogram from the distance matrices
## All traits weighted in one gower distance matrix
fulldendro <- speciesclusters(trait_gowdismat, no_clust_max = 20)
## Individual trait matrices, weighted and combined
combdendro <- speciesclusters(comb_weightedmat, no_clust_max = 20)

### Select a number of clusters and join groups with the original dataframe
k = 5 #selected no.clusters
sppgroups <- cutree(combdendro$clust.dendrograms[[1]], k = k) #get the groups 
sppgroups <- data.frame(species = names(sppgroups), cluster.id = c(sppgroups))

# Plot the tree 
fullcolor = grDevices::colors()[grep('dark', grDevices::colors())]
clustercolors <- sample(fullcolor, length(unique(sppgroups$cluster.id)))
dt.tree <- ape::as.phylo(combdendro$clust.dendrograms[[1]])
plot(dt.tree, type = "cladogram", no.margin = TRUE, 
     tip.color = clustercolors[sppgroups$cluster.id], cex = 0.8)

# Join to the initial dataset and save the table 
# Make sure the rows and columns match
unique_species <- unique_species_df %>% 
  select(-c(Beak.Length_Culmen, Beak.Depth,Secondary1, Kipps.Distance, Wing.Length))
specieslist <- unique_species[complete.cases(unique_species), ]
## Join the groups to the species list
spp_traitclusters <- dplyr::left_join(specieslist, sppgroups, c("species" = "species"))

## Write to csv
# replace the group integer with G1, G2, G3, etc.
spp_traitclusters$cluster.id <- paste0("G", spp_traitclusters$cluster.id)
write_csv(spp_traitclusters, "./output/species_clustering/functionaltraitclusters_k5.csv")

### 3.2 Join with gbif occurrences ------------
# Join the gbif data with the species clusters
sppclustersids <- spp_traitclusters %>% select(species, cluster.id)
gbif_clusters <- left_join(gbif_occ_fua, sppclustersids, c("species" = "species"),
                           keep = FALSE)

## Write to csv
write_csv(gbif_clusters, "./output/species_clustering/gbif_occurrences_sppclusters.csv")

## 4. Overview of traits per cluster ------------
library(gtsummary)
spp_traitclusters <- read_csv("./output/species_clustering/functionaltraitclusters_k5.csv")

# What is the number of species in each cluster.id?
spp_traitclusters %>% 
  group_by(cluster.id) %>% 
  summarise(n = n(), .groups = "drop")

## Count the number of species per group
glimpse(spp_traitclusters)
cat_traits <- c("Habitat", "Habitat.Density", "Migration", 
                "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle")
cont_traits <- c("Beak.Length_Nares", "Beak.Width", "Tarsus.Length", 
                 "Tail.Length", "Mass", "HWI")

## summary table for continuous variables
cont_summarytable <- spp_traitclusters %>%
  select(all_of(cont_traits), cluster.id) %>%
  tbl_summary(
    by = "cluster.id",
    #include median and range in the summary
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{mean} ({sd})", 
                                     "{median} ({p25} - {p75})",
                                     "{range}"),
    missing = c("no")
  ) %>%
  as_gt()  
gt::gtsave(cont_summarytable, file = "./output/species_clustering/cluster_continuousvar_summarytables.docx")

## Create a summary table for categorical variables
cate_summarytable <- spp_traitclusters %>%
    select(all_of(cat_traits), cluster.id) %>%
    tbl_summary(
      by = cluster.id, # No grouping within the cluster
      statistic = list(all_categorical() ~ "{n} ({p}%)"), # Show counts and percentages
      missing = c("no") # Exclude missing data from the summary
    ) %>% 
    as_gt() 
gt::gtsave(cate_summarytable, file = "./output/species_clustering/cluster_categoricalvar_summarytables.docx")
