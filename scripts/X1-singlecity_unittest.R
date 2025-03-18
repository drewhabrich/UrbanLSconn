## HEADER---------------------------
## Script name: X1-singlecity_unittest.R
##
## Purpose of script: This script will tidy the checklist data for Ottawa and merge it with the 
## landscape metrics data to generate some exploratory plots and models - As a single city test.
##
## Author: Andrew Habrich
##
## Date Created: 2025-01-30
## Date last Modified: 2025-03-11
##
## Email: 
## - andrhabr@gmail.com
## - andrewhabrich@cmail.carleton.ca
## 
## Notes --------------------------------------------------------------------
rm(list = ls())
library(chattr)
chattr_use("copilot")

## 1. Load relevant packages ------------------------------------------------
pacman::p_load(sf, terra, tidyterra, rnaturalearth, mapview, 
               mgcv, gam, lme4, mgcv, modelsummary,
               tidyverse, ggpubr, RColorBrewer, easystats, tictoc)

## Import the relevant files
city <- "Ottawa"
citylist <- read_csv("./_unittest/cityfiles_summary.csv") %>% 
  filter(cityname == city) 
cityebird <- readRDS("./_unittest/ebird_smplchkl_data.rds") %>% 
  filter(cityname == city)

## Create the landcover lookup table ----------------------------------------
lcclass_lookup <- tibble(value = c(10, 20, 30, 40, 
                                   50, 60, 70, 80, 
                                   90, 95, 100),
                         category = c("tree", "shrub", "grass", "crop",
                                      "builtup", "bare", "snow", "water", 
                                      "wetland", "mangrove", "mosslichen"))

# Convert landcover codes to the corresponding categories
prefix_map <- c(pland = "p", enn_mn = "ennm", enn_sd = "ennsd")
# Define a function that renames a vector of column names elementwise
rename_fn <- function(x) {
  sapply(x, function(name) {
    # Use regex to extract prefix, numeric code, and suffix
    m <- str_match(name, "^(pland|enn_mn|enn_sd)_(\\d+)_(.*)$")
    if (!is.na(m[1,1])) {
      original_prefix <- m[1,2]      # e.g., "pland"
      code_num <- as.numeric(m[1,3])   # e.g., 10
      suffix <- m[1,4]               # e.g., "500m"
      # Look up the category in lcclass_lookup where code equals code_num
      cat_val <- lcclass_lookup %>%
        filter(value == code_num) %>%
        pull(category)
      # If a match is found, build the new column name; else keep original name
      if (length(cat_val) > 0) {
        paste0(prefix_map[original_prefix], "_", cat_val, "_", suffix)
      } else {
        name
      }
    } else {
      name  # If the pattern doesn't match, return the original name
    }
  })
}

## check the columns in the dataframes
glimpse(citylist)
glimpse(cityebird)

## 2. Read in files ----------------
div <- readRDS(citylist$rds_div)
smpl <- readRDS(cityebird$smplRDS)
chkl <- readRDS(cityebird$chklRDS)

# Read in the LSMetric csv files
pland <- read_csv(paste0("./_unittest/pland_wide_", city, ".csv"))
aggrm <- read_csv(paste0("./_unittest/aggrm_wide_", city, ".csv"))
ennm <- read_csv(paste0("./_unittest/ennm_wide_", city, ".csv"))
dist_bound <- read_csv(paste0("./_unittest/", city, "_chkl_dist.csv")) %>% 
  select(c("checklist_id", "boundarydist_m", "uacentrdist_m"))

## Read in the builtV chkls
## Check what files are in the unittest folder
builtv_files <- list.files("./_unittest") %>% as_tibble() %>% 
  filter(str_detect(value, "builtv")) %>% 
  mutate(buffersize = as.numeric(str_extract(value, "[0-9]+"))) %>% 
  arrange(buffersize)
# import files and remove first column (the "" one)
# and add a suffix to the column names to indicate the buffer size, except for checklist_id
# bv250 <- read_csv(paste0("./_unittest/", builtv_files$value[1]), col_select = -1) %>% 
#   rename_with(~paste0("bv_", ., "_250m"), -checklist_id)
bv500 <- read_csv(paste0("./_unittest/", builtv_files$value[2]), col_select = -1) %>% 
  rename_with(~paste0("bv_", ., "_500m"), -checklist_id)
bv1000 <- read_csv(paste0("./_unittest/", builtv_files$value[3]), col_select = -1) %>% 
  rename_with(~paste0("bv_", ., "_1000m"), -checklist_id)
bv1500 <- read_csv(paste0("./_unittest/", builtv_files$value[4]), col_select = -1) %>% 
  rename_with(~paste0("bv_", ., "_1500m"), -checklist_id)
bv2000 <- read_csv(paste0("./_unittest/", builtv_files$value[5]), col_select = -1) %>% 
  rename_with(~paste0("bv_", ., "_2000m"), -checklist_id)
## Join the bv dataframes together by the checklist_id
bv <- bv500 %>% 
  #left_join(bv250, by = "checklist_id") %>% 
  left_join(bv1000, by = "checklist_id") %>%
  left_join(bv1500, by = "checklist_id") %>%
  left_join(bv2000, by = "checklist_id")

### 2.1 Join ALL the dataframes together using the checklist_id ----
## First join the pland, aggrm, ennm, dist_bound, and bv dataframes to the div dataframe
city_df <- div %>% select(1:7) %>% 
  left_join(pland, by = "checklist_id") %>% 
  left_join(aggrm, by = "checklist_id") %>% 
  left_join(ennm, by = "checklist_id") %>% 
  left_join(dist_bound, by = "checklist_id") %>% 
  left_join(bv, by = "checklist_id") 

## Clean the smpl dataframe 
glimpse(smpl)
smpl_clean <- smpl %>% select(c(checklist_id, observation_date, 
                                time_observations_started, protocol_type, 
                                duration_minutes, effort_distance_km, 
                                effort_area_ha, number_observers))

## Join the smpl_clean dataframe to the city_df dataframe 
#place the smpl dataframe after column 7 in the city_df dataframe
citychkl_lsm <- smpl_clean %>% 
  #inner join to only keep checklists that lsm were extracted for
  inner_join(city_df, by = "checklist_id") %>% 
  #Remove columns for the NA cells
  select(-starts_with("pland_NA")) %>% 
  select(-starts_with("enn_mn_NA")) %>%
  select(-starts_with("enn_sd_NA")) 
glimpse(citychkl_lsm)

## Check the columns in the citychkl_lsm dataframe
citychkl_lsm <- citychkl_lsm %>% 
  #Convert the observation_date to year, month, day columns
  mutate(year = year(observation_date), 
         month = month(observation_date, label = T, abbr = T), 
         day = day(observation_date)) %>% 
  #rearrange the new columns to the front of the dataframe after observation_date
  relocate(year, month, day, .after = observation_date) %>% 
  select(-observation_date) %>% 
  # rename columns with the LC lookup table
  rename_with(rename_fn, .cols = matches("^(pland_|enn_mn_|enn_sd_)\\d+_.*$"))

## Create a subfolder for the derived output if it doesn't exist yet
if (!dir.exists("./_unittest/derived_output")) {
  dir.create("./_unittest/derived_output")
}

#save the city_df to the folder
saveRDS(citychkl_lsm , "./_unittest/derived_output/Ottawa_chkl_lsmdata.rds")
write_csv(citychkl_lsm , "./_unittest/derived_output/Ottawa_chkl_lsmdata.csv")

## 3. Simple exploratory modelling ------------------------------------------
rm(list = ls())

### 3.1 Entire dataset ------------------------------------------------------
citychkl_lsm <- readRDS("./_unittest/derived_output/Ottawa_chkl_lsmdata.rds") 

# Plot the distribution of specrich, grouped by month 
gghistogram(citychkl_lsm, x = "specrich", bins = 40,
            add = "mean", rug = F, add_density = T,
            color = "month", fill = "month",
            palette = "Set1", alpha = 0.25)

# Plot the distribution of specrich, grouped by year
citychkl_lsm %>% mutate(year = as.factor(year)) %>%
  gghistogram(x = "specrich", bins = 40, facet.by = "year",
              add = "mean", rug = F, add_density = T,
              color = "year", fill = "year",
              alpha = 0.25) + 
  theme(legend.position = "none")

# Plot the correlation between 'effort' indicators and the species richness
citychkl_lsm %>% 
  select(specrich, number_observers,
         effort_distance_km, duration_minutes) %>% 
  correlation() %>% 
  summary(redundant = F) %>% 
  plot() + theme_bw()

## checklist duration
plot(cor_test(citychkl_lsm, "duration_minutes", "specrich"),
     point = list(color = "grey", alpha = 0.25),
     smooth = list(color = "blue", alpha = 1)) + 
  facet_wrap(citychkl_lsm$year) + 
  theme_bw() 

## checklist distance travelled
plot(cor_test(citychkl_lsm, "effort_distance_km", "specrich"),
     point = list(color = "grey", alpha = 0.25),
     smooth = list(color = "blue", alpha = 1)) + 
  facet_wrap(citychkl_lsm$year) + 
  theme_bw() 

# Plot the correlation between species richness and year
plot(cor_test(citychkl_lsm, "year", "specrich")) + 
  theme_bw() + 
  #change the x-axis labels to be ech year
  scale_x_continuous(breaks = seq(2010, 2022, 1))

# Plot the correlation between species richness and month
plot(cor_test(citychkl_lsm, "month", "specrich")) + 
  theme_bw() 

### 3.2 Sub-datasets for each bufferradius size -----------------------------
# Green landcover to combine
cols_to_sum <- c("p_tree", "p_shrub", "p_grass")

# individual buffer dataframes
chkl500m <- citychkl_lsm %>% select(-ends_with("_1000m"), -ends_with("_1500m"), -ends_with("_2000m")) %>% 
  rename_with(~str_remove(., "_500m")) %>% #remove the buffer suffix from the column names
  #combine the green pland columns into a new column if they are present in the dataframe
  rowwise() %>%
  mutate(p_greencover = sum(c_across(all_of(intersect(cols_to_sum, names(.)))), na.rm = TRUE)) %>% ungroup()
chkl1000m <- citychkl_lsm %>% select(-ends_with("_500m"), -ends_with("_1500m"), -ends_with("_2000m")) %>% 
  rename_with(~str_remove(., "_1000m")) %>% 
  rowwise() %>%
  mutate(p_greencover = sum(c_across(all_of(intersect(cols_to_sum, names(.)))), na.rm = TRUE)) %>% ungroup()
chkl1500m <- citychkl_lsm %>% select(-ends_with("_500m"), -ends_with("_1000m"), -ends_with("_2000m")) %>% 
  rename_with(~str_remove(., "_1500m")) %>% 
  rowwise() %>%
  mutate(p_greencover = sum(c_across(all_of(intersect(cols_to_sum, names(.)))), na.rm = TRUE)) %>% ungroup()
chkl2000m <- citychkl_lsm %>% select(-ends_with("_500m"), -ends_with("_1000m"), -ends_with("_1500m")) %>% 
  rename_with(~str_remove(., "_2000m")) %>% 
  rowwise() %>%
  mutate(p_greencover = sum(c_across(all_of(intersect(cols_to_sum, names(.)))), na.rm = TRUE)) %>% ungroup()

## Put them into a list
chkl_bufferlist <- list(chkl500m, chkl1000m, chkl1500m, chkl2000m)
names(chkl_bufferlist) <- c("500m", "1000m", "1500m", "2000m")
## Save the list to RDS
saveRDS(chkl_bufferlist, "./_unittest/derived_output/Ottawa_chkl_bufferlist.rds")