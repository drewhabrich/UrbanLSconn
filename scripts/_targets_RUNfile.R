### Load packages
#install.packages("geotargets", repos = c("https://njtierney.r-universe.dev", "https://cran.r-project.org"))
pacman::p_load(targets, tarchetypes, geotargets, crew, qs2,
               ggspatial, sf, terra, tidyterra, rnaturalearth, mapview,
               auk, vegan, landscapemetrics, exactextractr, gdistance,
               tidyverse, dplyr, ggpubr, RColorBrewer, easystats)

### NOTE ALSO CHECK https://wlandau.github.io/targetopia/ FOR MORE INFO ####
library(targets)
#run tar_make() to run the pipeline
#run tar_manifest() to see the steps in the pipeline
#run tar_visnetwork() to see the visual input-output steps for the pipeline
#tar_read() to view the results for each target.
#tar_load() to load the results for each target.

tar_source("R/")
tar_visnetwork(targets_only = T, label=c("time", "branches", "size"))
tar_visnetwork(targets_only = T)

## Run pipeline ############################################################
targets::tar_make(reporter = "verbose_positives", as_job = T,
         terminate_controller = T) #run as_job to run in the background
targets::tar_make(reporter = "summary", terminate_controller = T, use_crew = F)

## SPECIFIC TARGETS
targets::tar_make("pland500_results", reporter = "verbose_positives", as_job = T, terminate_controller = T)

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
meta <- tar_meta(targets_only = TRUE) 
meta_err <- tar_meta(targets_only = TRUE, fields = error, complete_only = T) 
tar_meta(fields = warnings, complete_only = TRUE)
#tar_meta_delete(meta = TRUE, progress = TRUE, process = TRUE, verbose = TRUE, delete = "local")

#### Load specific targets ################################################
# SPATIAL DATA -----------------------------------------------------------
### Load the urban areas polygons
targets::tar_load("urbanareas")
plot(urbanareas)
plot(urbanareas %>% filter(eFUA_name == "Ottawa"))
mapview(urbanareas)

### ESA 10m raster data
targets::tar_load("urban_ESA_LC")

### ghsl cropped building volumes raster
targets::tar_load("ghsl_builtLC")

# EBIRD DATA -------------------------------------------------------------
targets::tar_load(cityfiles_summary)

############################################## TEST ZONE ###############
tar_load("ebird_citydatalist")
# get the city name of the first entry
cn <- ebird_citydatalist$cityname[23]
citysmpl <- ebird_citydatalist %>% filter(cityname == cn) %>% pull(smplRDS[1]) %>% readRDS() %>% 
  #remove useless columns
  select(-c("last_edited_date", "country", "country_code", 
            "state", "state_code", "county", "county_code", 
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
citylsm <- lsm_bycity2500 %>% filter(cityname == cn) %>% pull(rds_lsm[1]) %>% 
  readRDS() %>% 
  select(-starts_with("NA")) %>% #fill 0 for NA values in the columns
  mutate(across(everything(), ~replace_na(., 0)))
citybv <- builtv_bycity2500 %>% filter(cityname == cn) %>% pull(rds_builtv[1]) %>% readRDS() %>% 
  rename(bv_mean = mean,
         bv_stdev = stdev,
         bv_sum = sum)

## Join the dataframes together, by the checklist_id
citydata <- left_join(citydiv, citysmpl, by = "checklist_id") %>% 
            left_join(citylsm, by = "checklist_id") %>% 
            left_join(citybv, by = "checklist_id") %>% 
  #remove pland_***_ from the column names, with *** indicating wildcards
            rename_with(~str_remove(., "and_.*_"), starts_with("pland")) %>% 
  #replace pl*** with p_***
            rename_with(~str_replace(., "pl", "p_"), starts_with("pl"))

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
     point = list(color = "grey", alpha = 0.25)) + facet_wrap(citydata$year) + theme_bw()#checklist duration
plot(cor_test(citydata, "specrich", "effort_distance_km"),
     point = list(color = "grey", alpha = 0.5)) + facet_wrap(citydata$year) #checklist distance travelled

# Plot the correlation between species richness and year
plot(cor_test(citydata, "year", "specrich")) + theme_bw() + 
  #change the x-axis labels to be ech year
  scale_x_continuous(breaks = seq(2010, 2022, 1))

# Plot the correlation between the landscape metrics and raster data
viz <- citydata %>% select(specrich, starts_with("p_"), starts_with("bv")) %>% 
  correlation() %>% summary() %>% 
  visualisation_recipe(labs = list(title = "Correlation bewteen Prop.LC and BuiltEnv. metrics"))

plot(viz) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.70))

# Plot the scatterplot of the pland and bv metrics
citydata %>% select(specrich, starts_with("p_"), starts_with("bv")) %>% 
  ggscatter(x = "p_treecover", y = "bv_mean", alpha = 0.5)

# Plot the distribution of pbuiltup values
citydata %>%
  gghistogram(x = "p_builtup", binwidth = 2,
              add = "mean", rug = T, add_density = T,
              alpha = 0.25, color = "black", fill = "pink")

# Plot the built environment metrics
citydata %>%
  ggdensity(x = "bv_mean", 
              add = "mean", rug = T, 
              alpha = 0.25, color = "black", fill = "grey")

# Plot the distribution of tree cover values
citydata %>%
  select(p_treecover, year) %>% 
  gghistogram(x = "p_treecover", binwidth = 2,
              add = "mean", rug = T, add_density = T,
              alpha = 0.25, color = "black", fill = "forestgreen")

## NEED TO SELECT CHECKLISTS WITHIN RANGE OF pLAND VALUES
citydata %>% 
  filter(p_treecover >= 25 & p_treecover <= 75) %>%
  #plot scatterplot of tree cover and built volume
  ggscatter(x = "p_treecover", y = "bv_mean", 
            add = "reg.line", add.params = list(color = "black"),
            alpha = 0.5, color = "forestgreen",
            facet.by = "year")
  
# Plot a scatterplot of ptreecover and specrich
citydata %>% 
  filter(p_treecover >= 25 & p_treecover <= 75) %>%
  ggscatter(x = "p_treecover", y = "specrich", alpha = 0.5, color = "forestgreen",
            add = "reg.line", add.params = list(color = "black"), 
            xlab = "Prop. Tree cover in buffer", ylab = "Species richness")

# Plot a scatterplot of bv_mean and specrich
citydata %>% mutate(bv_mean = bv_mean/1000) %>%
  ggscatter(x = "bv_mean", y = "specrich", alpha = 0.5, color = "grey",
            add = "reg.line", add.params = list(color = "black"), 
            xlab = "Mean built heights in buffer", ylab = "Species richness")

##### REGRESSION MODELLING #####
pacman::p_load(mgcv, lme4, MuMIn, easystats, modelsummary)
glimpse(citydata)
### NEED TO ACCOUNT FOR SAMPLING EFFORT (distance travelled, duration)
### NEED TO ACCOUNT FOR TEMPORAL SAMPLING (month, year)

#plot the distribution of the specrich values
citydata %>% gghistogram(x = "specrich", binwidth = 1, alpha = 0.5, fill = "grey")

# Generalized additive model
gam_null <- gam(specrich ~ 1, 
               family = "gaussian", data = citydata) 
gam_mod <- gam(specrich ~ s(p_treecover) + s(bv_mean) + s(duration_minutes) + s(year), 
               family = "gaussian", data = citydata) 
plot(gam_mod, pages = 1)
modelsummary(gam_mod)
coef(gam_mod)

## check the results of the gam
par(mfrow = c(2, 2))
gam.check(gam_mod, k.sample = 5000, k.rep = 200)

## Fit a linear and smoothed term models to compare
gam_l <- gam(specrich ~ p_treecover, 
             family = "gaussian", data = citydata)
gam_s <- gam(specrich ~ s(p_treecover), 
             family = "gaussian", data = citydata)

ggplot(citydata, aes(x = p_treecover, y = specrich)) + geom_point() +
  geom_line(colour = "red", linewidth = 1.2, aes(y = fitted(gam_l))) +
  geom_line(colour = "blue", linewidth = 1.2, aes(y = fitted(gam_s))) +
  theme_bw()
AIC(gam_l, gam_s)
