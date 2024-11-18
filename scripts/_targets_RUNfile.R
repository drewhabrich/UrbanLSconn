### Load packages
#install.packages("geotargets", repos = c("https://njtierney.r-universe.dev", "https://cran.r-project.org"))
pacman::p_load(tidyverse, ggpubr, RColorBrewer, easystats,
               ggspatial, sf, terra, tidyterra, rnaturalearth, mapview,
               auk, vegan, landscapemetrics, exactextractr,
               targets, tarchetypes, geotargets)

### NOTE ALSO CHECK https://wlandau.github.io/targetopia/ FOR MORE INFO ####

#run tar_make() to run the pipeline
#run tar_manifest() to see the steps in the pipeline
#run tar_visnetwork() to see the visual input-output steps for the pipeline
#tar_read() to view the results for each target.
#tar_load() to load the results for each target.
tar_source("R/")
tar_visnetwork(targets_only = T, label=c("time", "branches", "size"))

## Run pipeline ############################################################
tar_make(reporter = "verbose", as_job = TRUE) #run as_job to run in the background

## Check pipeline PROGRESS ################################################
tar_progress_summary() 
tar_progress_branches()
tar_progress()

tar_poll(
  interval = 30,
  timeout = Inf,
  fields = c("skipped", "dispatched", "completed", "errored", "canceled", "since"),
)

## Check pipeline RESULTS #################################################
## Clear global environment
rm(list=ls())

## Load all target objects (or specific ones)
data <- tar_manifest()
meta <- tar_meta()
#tar_meta_delete(meta = TRUE, progress = TRUE, process = TRUE, verbose = TRUE, delete = "local")

#### Load specific targets ################################################
# SPATIAL DATA -----------------------------------------------------------
### Load the urban areas polygons
targets::tar_load("urbanareas")
glimpse(urbanareas)
plot(urbanareas)
plot(urbanareas %>% filter(eFUA_name == "Ottawa"))
mapview(urbanareas)

### ESA 10m raster data
targets::tar_load("urban_ESA_LC")

### ghsl cropped building volumes raster
targets::tar_load("ghsl_builtLC")

# EBIRD DATA -------------------------------------------------------------
# CITY FILE SUMMARY
targets::tar_load("ebird_citydatalist")

# ebird filters for each city
targets::tar_load("ebird_citydata")

# ebird RDS - this is the raw data for each city
targets::tar_load("ebird_RDS")

# ebird filtered to each city - this is the filtered, combined sampling and checklist data
targets::tar_load("ebird_bycity")
str(ebird_bycity)

# city summary; # of checklists and gamma diversity for the city ---------
targets::tar_load("city_ebirdsums")
city_ebirdsums

# zerofilled data
targets::tar_load("zf_citydata")

# tidy data frame and diversity calculations
targets::tar_load("city_tidyzfdf")

#### Diversity and landscape metrics -------------------------------------
# diversity metrics
targets::tar_load("div_bycity")

# landscape metrics
targets::tar_load("lsm_bycity")

# built environment metrics
targets::tar_load("builtv_bycity")

############################################## TEST ZONE ###############
tar_load("ebird_citydatalist")
# get the city name of the first entry
cn <- ebird_citydatalist$cityname[1]
citysmpl <- ebird_citydatalist %>% filter(cityname == cn) %>% pull(smplRDS[1]) %>% readRDS() %>% 
  #remove useless columns
  select(-c("last_edited_date", "country", "country_code", "state", "state_code", "county", "county_code", 
            "iba_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type",
            "all_species_reported", "project_code", "protocol_code",
            "observer_id", "group_identifier", "trip_comments")) %>% 
  #create a year, month, and ordinal day column
  mutate(year = year(observation_date), 
         month = month(observation_date, label = TRUE, abbr = TRUE, locale = "en"), 
         ordinal_day = yday(observation_date))
  
## checklist data
citydiv <- div_bycity %>% filter(cityname == cn) %>% pull(rds_div[1]) %>% readRDS() %>% 
  select(1:5)
citylsm <- lsm_bycity %>% filter(cityname == cn) %>% pull(rds_lsm[1]) %>% readRDS() %>% 
  select(-starts_with("NA")) %>% #fill 0 for NA values in the columns
  mutate(across(everything(), ~replace_na(., 0)))
citybv <- builtv_bycity %>% filter(cityname == cn) %>% pull(rds_builtv[1]) %>% readRDS() %>% 
  rename(bv_mean = mean,
         bv_stdev = stdev,
         bv_sum = sum)

## Join the dataframes together, by the checklist_id
citydata <- left_join(citydiv, citysmpl, by = "checklist_id") %>% 
            left_join(citylsm, by = "checklist_id") %>% 
            left_join(citybv, by = "checklist_id")

### Check the data ------------------------------------------------------
glimpse(citydata)

# Plot the distribution of specrich, grouped by month 
gghistogram(citydata, x = "specrich", bins = 40,
            add = "mean", rug = F, add_density = T,
            color = "month", fill = "month",
            palette = "Set1", alpha = 0.25)

# Plot the distribution of specrich, grouped by year
citydata %>% mutate(year = as.factor(year)) %>%
gghistogram(x = "specrich", bins = 40, facet.by = "year",
            add = "mean", rug = F, add_density = T,
            color = "year", fill = "year",
            alpha = 0.25) + 
  theme(legend.position = "none")

# Plot the correlation between 'effort' indicators and the species richness
citydata %>% select(specrich, effort_distance_km, duration_minutes, number_observers) %>% 
  correlation() %>% 
  summary(redundant = F) %>% 
  plot() + theme_bw()
plot(cor_test(citydata, "specrich", "duration_minutes"),
     point = list(color = "grey", fill = "black", alpha = 0.5)) + facet_wrap(citydata$year)#checklist duration
plot(cor_test(citydata, "specrich", "effort_distance_km")) #checklist distance travelled

# Plot the correlation between species richness and year
plot(cor_test(citydata, "year", "specrich")) + theme_bw() + 
  #change the x-axis labels to be ech year
  scale_x_continuous(breaks = seq(2010, 2022, 1))

# Plot the correlation between the landscape metrics and raster data
viz <- citydata %>% select(specrich, starts_with("pland"), starts_with("bv")) %>% 
  #remove pland_***_ from the column names, with *** indicating wildcards
  rename_with(~str_remove(., "land_.*_"), starts_with("pland")) %>%
  correlation() %>% summary() %>% 
  visualisation_recipe(labs = list(title = "Correlation bewteen Prop.LC and BuiltEnv. metrics"))

plot(viz) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.70))

# Plot the scatterplot of the pland and bv metrics
citydata %>% select(specrich, starts_with("pland"), starts_with("bv")) %>% 
  rename_with(~str_remove(., "land_.*_"), starts_with("pland")) %>%
  ggscatter(x = "ptreecover", y = "bv_mean", add = "reg.line", 
            alpha = 0.5)

# Plot the distribution of pbuiltup values
citydata %>% rename_with(~str_remove(., "land_.*_"), starts_with("pland")) %>%
  select(pbuiltup, year) %>% 
  gghistogram(x = "pbuiltup", binwidth = 2,
              add = "mean", rug = T, add_density = T,
              alpha = 0.25, color = "black", fill = "pink")

# Plot the built environment metrics
citydata %>% select(bv_mean, year) %>% 
  ggdensity(x = "bv_mean", 
              add = "mean", rug = T, add_density = T,
              alpha = 0.25, color = "black", fill = "grey")

# Plot the distribution of tree cover values
citydata %>% rename_with(~str_remove(., "land_.*_"), starts_with("pland")) %>% 
  select(ptreecover, year) %>% 
  gghistogram(x = "ptreecover", binwidth = 2,
              add = "mean", rug = T, add_density = T,
              alpha = 0.25, color = "black", fill = "green")

## NEED TO SELECT CHECKLISTS WITHIN RANGE OF pLAND VALUES
citydata %>% rename_with(~str_remove(., "land_.*_"), starts_with("pland")) %>% 
  filter(ptreecover >= 25 & ptreecover <= 75) %>%
  #plot scatterplot of tree cover and built volume
  ggscatter(x = "ptreecover", y = "bv_stdev", add = "reg.line", alpha = 0.5, fill = "forestgreen")
  
### NEED TO ACCOUNT FOR SAMPLING EFFORT (distance travelled, duration)
### NEED TO ACCOUNT FOR TEMPORAL SAMPLING