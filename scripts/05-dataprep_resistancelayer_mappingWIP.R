## HEADER---------------------------
## Script name: 05-dataprep_resistancelayer_mapping
##
## Purpose of script: Use clustered species to generate habitat suitability models 
##              to parameterize the resistance layer for the connectivity analysis
##
## Author: Andrew Habrich
##
## Date Created: 2025-02-18
## Date last Modified: 2025-03-12
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes ---------------------------
rm(list = ls())
pacman::p_load(sf, terra, tidyterra,
               tidyverse, ggpubr, readxl, tictoc,
               biomod2, foreach, doParallel)

## 1.0 Read in the data -----------------------------------------------------
# Urban area data
fuacan <- vect("./data/00-ghs_fua_canada.gpkg", layer = "ghs_fua")

# Raster data locations
landcover <- "./output/ESArasters_factor/"
builtv <- "./output/GHSL_builtv/"
landsatEVI <- "./output/landsat_EVI/"

# gbif occurrences and species clusters
gbif_clusters <- read_csv("./output/species_clustering/gbif_occurrences_sppclusters.csv")
gbif_clusters <- gbif_clusters %>%
  # from the geometry column, remove 'c(' and ')' and split the string into two columns
  mutate(geometry = str_remove_all(as.character(geometry), "c\\(|\\)")) %>%
  separate(geometry, into = c("long", "lat"), sep = ",\\s*") %>% 
  mutate(long = as.numeric(long), lat = as.numeric(lat)) %>%
  # replace the NA in the cluster.id with 0
  mutate(cluster.id = ifelse(is.na(cluster.id), "G0", cluster.id))

# Create subdatasets for the species clusters (cluster.id)
# use a for loop to filter the gbif_clusters data to each cluster.id and save as a csv
# for (i in unique(gbif_clusters$cluster.id)) {
#   gbif_clusters %>%
#     filter(cluster.id == i) %>%
#     select(c(long, lat, cluster.id)) %>% 
#     # replace the values in the cluster.id column with 1 for occurrences
#     mutate(cluster.id = 1) %>%
#     # write the data to a csv
#     write_csv(paste0("./output/habitat_suitability/gbif_occurrences_sppclusters_", i, ".csv"))
# }

### 1.1 Data exploration ----------------------------------------------------
## How many families, order, and genera are there in each cluster.id?
gbif_clusters %>% 
  group_by(cluster.id) %>% 
  summarise(n_obs = n(),
            n_species = n_distinct(species), 
            n_families = n_distinct(family),
            n_order = n_distinct(order),
            n_genus = n_distinct(genus)) 

## Plot a barplot of how many observations are in each cluster.id
gbif_clusters %>% 
  count(cluster.id) %>% 
  ggplot(aes(x = cluster.id, y = n)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  #add the count values over the bars
  geom_text(aes(label = paste0("n = ",n), vjust = -0.3)) +
  labs(title = "Number of observations in each cluster",
       x = "Cluster ID",
       y = "Number of observations")

## 2.0 Filter obs to one city -----------------------
cityname <- "Ottawa"
citypoly <- fuacan %>% 
  filter(eFUA_name == cityname)

## Filter the occurrence data to the city of interest
cityocc <- vect(gbif_clusters, geom = c("long", "lat"), crs = "EPSG:4326", keepgeom = T) %>% 
  #extract the points that are within the city boundary
  .[citypoly] # equivalent to st_intersection(., citypoly)

## 2.1 Get list of cities and observations
# Get the list of cities  
# cities <- unique(fuacan$eFUA_name)
# # create a list of occurrences by the city
# cityocc_list <- list()
# # save a dataframe for each city into the list
# for (city in cities) {
#   cityocc_list[[city]] <- vect(gbif_clusters, geom = c("long", "lat"), 
#                                crs = "EPSG:4326", keepgeom = T) %>% 
#     .[fuacan %>% filter(eFUA_name == city)]
# }

## 3.0 Generate habitat suitability models -----------------------------------
### 3.1 Set parameters for the biomod2 modelling -----------------------
PA.strat <- 'random' #pseudo-absence strategy
PA.n_absences <- 5000
# Params for modelling
lst.mod <- c("RF", "GAM", "ANN") 
nCrossVal <- 5 #number of cross-validation folds
cv.perc <- 0.75 #percentage of data to use for training
nperm.var.imp <- 3 #number of permutations for variable importance
ens.calc <- c('EMca') #ensemble calculation method
nrep.PA <- 5 #number of pseudo-absence replicates
TSS.min <- 0.3 #minimum threshold for TSS

# read in the city specific rater data for the city of interest
city_lc <- rast(paste0(landcover, cityname, ".tif"))
city_bv <- rast(paste0(builtv, cityname,"_builtV.tif")) %>% 
  project(y = city_lc, method = "bilinear")
city_evi <- rast(paste0(landsatEVI, cityname, "_EVI_2020.tif")) %>% 
  project(y = city_lc, method = "bilinear")
## Combine the rasters into 1 object with layers
explrasters <- c(city_bv, city_lc, city_evi)
names(explrasters) <- (c("BuiltVol", "Landcover", "EVI"))
rm(city_lc, city_bv, city_evi, gbif_clusters); gc() #remove them to clear memory

### 3.2 Loop to calculate for each group ########################################
# Create a list of dataframes for each cluster in the cityocc data
clustlist <- list()
for (clust in unique(cityocc$cluster.id)) {
  clustlist[[clust]] <- cityocc %>% 
    project(terra::crs(explrasters)) %>% 
    filter(cluster.id == clust) %>% 
    select(long, lat, cluster.id) %>% 
    mutate(cluster.id = 1)
}
## remove cluster G0
clustlist$G0 <- NULL

# Loop to calculate the habitat suitability models for each cluster
for (clust in names(clustlist)) {
  g <- clustlist[[clust]]
  clustid <- str_extract(clust, "[0-9]+")
  # Create a folder for the sppcluster if it doesn't exist, otherwise continue
  if (!dir.exists(paste0('./output_hs/cl', clustid))) {
    dir.create(paste0('./output_hs/cl', clustid))
  }
  # Format the data for the models
  myBiomodData <- biomod2::BIOMOD_FormatingData(resp.name = paste0("cl", clustid),
                                                resp.var = g[,"cluster.id"], #PRESENCE ONLY
                                                resp.xy = g[, c("long", "lat")],
                                                expl.var = explrasters,
                                                PA.strategy = PA.strat,
                                                PA.nb.rep = nrep.PA,
                                                PA.nb.absences = PA.n_absences,
                                                dir.name = './output_hs',
                                                filter.raster = F)
  saveRDS(myBiomodData, paste0('./output_hs/cl', clustid, '/myBiomodData.rds'))
  
  # Individual model fitting
  bmoptions <- biomod2::bm_ModelingOptions(data.type = "binary", 
                                           models = lst.mod,
                                           bm.format = myBiomodData,
                                           strategy = "bigboss")
  ## Modify some model parameters to avoid model overfitting
  bmoptions@options$RF.binary.randomForest.randomForest@args.values$`_allData_allRun`$nodesize <- 
    round(sum(g$cluster.id)/100)
  
  # Individual model fitting
  myBiomodModelOut <- try(biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                                   modeling.id = 'AllModels',
                                                   models = lst.mod,
                                                   OPT.user = bmoptions,
                                                   OPT.data = 'binary',
                                                   OPT.strategy = 'user.defined',
                                                   OPT.user.base = 'user.defined',
                                                   CV.strategy = 'random',
                                                   CV.nb.rep = nCrossVal,
                                                   CV.perc = cv.perc,
                                                   metric.eval = c('TSS','ROC'), 
                                                   var.import = nperm.var.imp))
  saveRDS(myBiomodModelOut, paste0('./output_hs/cl', clustid, '/myBiomodModelOut.rds'))
  gc()
}

## Check results of the individual models by cluster
# Load the models
clusters <- list.files(paste0('./output_hs/'))

## Create a summary document for the model evaluations for each cluster
for (i in 1:length(clusters)) {
  # Load model
  myBiomodModelOut <- readRDS(paste0('./output_hs/', clusters[i], '/myBiomodModelOut.rds'))
  # Get evaluation scores & variables importance
  Evalplot <- bm_PlotEvalMean(myBiomodModelOut, dataset = "validation")$plot + 
    theme_bw() + 
    scale_color_discrete(name = "Models")
  Evalboxplot <- bm_PlotEvalBoxplot(myBiomodModelOut, group.by = c('algo','run'), dataset = "validation")$plot + 
    theme_bw()
  EvalVarimp <- bm_PlotVarImpBoxplot(myBiomodModelOut)$plot + 
    theme_bw() + ylab("Variable importance") +
    scale_fill_discrete(name = "Variables", labels = c("Landcover", "Vegetation index", "Built Vol. m3"))
  evallist <- list(Evalplot, Evalboxplot, EvalVarimp)
  # Save evaluations to folder
  ggarrange(plotlist = evallist, ncol = 1, labels = clusters[i], hjust = -0.75) %>% 
    ggexport(filename = paste0('./output_hs/cl', i, '_modelevals.png'),
             height = 2400, width = 1800, res = 300)
}

### 3.3 Project single models #################################################
# Project single models based on the best fit model
clustfolders <- list.dirs(path = "./output_hs/", full.names = F, recursive = F)
# Loop to generate the projections for each cluster

for (i in 1:length(clustfolders)) {
  #import model output
  bmout <- readRDS(paste0("./output_hs/", clustfolders[i], "/myBiomodModelOut.rds"))
  #find the best model based on ROC
  bmalgo <- get_evaluations(bmout) %>%
    filter(metric.eval == "ROC") %>%
    filter(run != "allRun") %>%
    group_by(algo) %>%
    summarise(mean = mean(validation, na.rm = T)) %>% 
    #get the algo of the highest mean
    filter(mean == max(mean)) %>% pull(algo)
  #get the name of the best fit model
  bestfitmodel <- get_evaluations(bmout) %>% 
    filter(PA == "allData", run == "allRun", algo == bmalgo, metric.eval == "ROC") %>% 
    pull(full.name)
  #print out which model is being projected
  print(paste0("BIOMOD Projection using", bmalgo))
  #Project single models
  bmProj <- BIOMOD_Projection(bm.mod = bmout,
                               proj.name = 'RF_Proj',
                               new.env = explrasters,
                               models.chosen = bestfitmodel,
                               metric.binary = 'all',
                               metric.filter = 'all',
                               on_0_1000 = T,
                               keep.in.memory = F,
                               build.clamping.mask = F,
                               output.format = ".tif")
  saveRDS(bmProj, paste0("./output_hs/", clustfolders[i], "/Biomod_RF_Proj.rds"))
  #print out the completion of the cluster
  print(paste0("Completed cluster"), clustfolders[i])
  gc()
}

## Ensemble modelling ######################
# bmout1 <- readRDS("./output_hs/cl1/myBiomodModelOut.rds")
# biomod_em1 <- try(biomod2::BIOMOD_EnsembleModeling(bm.mod = bmout1,
#                                                models.chosen = 'all',
#                                                em.by = 'all',
#                                                em.algo = ens.calc ,
#                                                metric.select = c('TSS'),
#                                                metric.select.thresh = TSS.min,
#                                                metric.eval = c('TSS', 'ROC'),
#                                                var.import = nperm.var.imp))
# saveRDS(biomod_em1, "./output_hs/cl1/Biomod_EM.rds")
# 
# # Get evaluation scores & variables importance
# biomod_em1 <- readRDS("./output_hs/cl1/Biomod_EM.rds")
# # Represent evaluation scores
# bm_PlotEvalMean(bm.out = biomod_em1, dataset = 'calibration')
# bm_PlotEvalBoxplot(bm.out = biomod_em1, group.by = c('algo', 'algo'))
# 
# # Represent variables importance
# bm_PlotVarImpBoxplot(bm.out = biomod_em1, group.by = c('expl.var', 'algo', 'algo'))
# bm_PlotVarImpBoxplot(bm.out = biomod_em1, group.by = c('expl.var', 'algo', 'merged.by.PA'))
# bm_PlotVarImpBoxplot(bm.out = biomod_em1, group.by = c('algo', 'expl.var', 'merged.by.PA'))
# 
# # Represent response curves
# bm_PlotResponseCurves(bm.out = biomod_em1,
#                       models.chosen = get_built_models(biomod_em1),
#                       fixed.var = 'median')
# bm_PlotResponseCurves(bm.out = biomod_em1,
#                       models.chosen = get_built_models(biomod_em1),
#                       fixed.var = 'min')
# 
# # Project ensemble models (building single projections)
# get_built_models(biomod_em1)
# EMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = biomod_em1,
#                                               proj.name = 'cl1EM',
#                                               new.env = explrasters,
#                                               models.chosen = 'all',
#                                               metric.binary = 'all',
#                                               metric.filter = 'all',
#                                               on_0_1000 = T,
#                                               keep.in.memory = F,
#                                               build.clamping.mask = F,
#                                               na.rm = T,
#                                               nb.cpu = 1,
#                                               output.format = ".tif")

## 4. Resistance layer parameterization --WIP---------------------------------

# The resistance layers are located in /data/derived-data/ResistanceSurfaces/ folder
# The source layers are located in /data/derived-data/SourceLayers/ folder

p.threshold <- seq(0.5, 0.7, by = 0.1) #threshold to define source from estimated probabilities of suitability
coef.c <- c(2, 4, 8, 16) #habitat suitability to resistance transformation coefficient

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  
  for (c in c(1:k)) {
    
    hab.suit <- terra::rast(here::here(paste0('data/derived-data/outputSDM/', g,'.GroupID.', c, '/proj_CurrentEM/proj_CurrentEM_', g,'.GroupID.', c,'_ensemble.tif')))
    hab.suit <- hab.suit/1000
    res <- 100 - 99*(1-exp(-coef.c*hab.suit))/(1-exp(-coef.c))
    names(res) <- paste0('tranfo.coef.', coef.c)
    terra::writeRaster(res, here::here(paste0('data/derived-data/ResistanceSurfaces/ResistanceSurface_', g, '_GroupID_', c, '_TransfoCoef_', coef.c, '.tif')), overwrite = T)
    
    for (thre in p.threshold) {
      srce <- hab.suit
      srce[srce < thre] <- 0
      terra::writeRaster(srce, here::here(paste0('data/derived-data/SourceLayers/SourceLayer_', g, '_GroupID_', c, '_SuitThreshold_',thre, '.tif')), overwrite = T)
    }
    
  }
}
