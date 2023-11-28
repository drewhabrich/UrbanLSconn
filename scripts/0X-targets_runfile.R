### Load packages
library(targets)
library(tarchetypes)
library(tidyverse)
library(DataExplorer)
library(vegan)
library(patchwork)
library(mapview)

tar_source("R/")
tar_visnetwork(targets_only = T, label=c("time", "branches", "size"))
tar_visnetwork(targets_only = F, label=c("time", "branches", "size"))

## Run pipeline
tar_make()

## Clear global environment
rm(list=ls())

## Load all target objects (or specific ones)
data <- tar_manifest()
meta <- tar_meta()

tar_load(NA_basemaps)
tar_load(urb_popcentres)
tar_load(ebird_citylist)
tar_load(zf_checklists)
tar_load(div_bycity)
tar_load(allcities_divdf)

tar_load(plot_cityebirdsums)
ggsave("./output/figures/tar_fig_urbanebirdchkl_bycity.png", plot = plot_cityebirdsums, scale =2)

tar_load(plot_divbycity)
ggsave("./output/figures/tar_fig_city_diversitymetrics.png", plot = plot_divbycity, scale = 2)

## how many checklists per year
tar_load(allcities_zfdf)
ott_zf <- allcities_zfdf$`Ottawa - Gatineau`
rm(allcities_zfdf)
gc() ##remove and clear memory to remove unnecessary memory usage

ott_zf %>% summarize(count = n_distinct(checklist_id))
str(ott_zf)

ott_zf %>% mutate(year = year(observation_date)) %>% 
  group_by(year) %>% 
  summarize(chkl_count = n_distinct(checklist_id)) %>% 
  ggplot() + 
  geom_bar(aes(y = chkl_count, x = year), fill = "blue", stat = "identity") +
  labs(title = "Complete checklists in Ottawa-Gatineau 2010-23", x = "Year", y = "# of checklists") + 
  theme_bw()
ggsave("./output/figures/tar_fig_ottgat_chklscount.png", plot = last_plot(), scale = 2)

tar_load(div_bycity, branches = 20) #this is ottawa-gatineau
ottdiv <- div_bycity$div_bycity_1d015e62
head(ottdiv, 10)

### How many checklists have NA for diversity metrics?
ottdiv %>% summarize(na_shannon = sum(is.na(shannondiv)),
                     na_specrich = sum(is.na(specrich)),
                     na_evenness = sum(is.na(evenness))) 
## WHY ARE THOSE ROWS EMPTY?

### PLOT CITY WITH DIVERSITY BY YEAR (3 COLS, FOR DIVERSITY, ROWS FOR YEARS)
o_min <- ott_zf %>% select(c("checklist_id","observation_date","effort_hours","effort_distance_km", "effort_speed_kmph"))

o_min <- ott_zf %>% group_by(checklist_id) %>% summarize(obs_date = unique(observation_date),
                                               effort_hours = unique(effort_hours),
                                               effort_distance_km = unique(effort_distance_km),
                                               effort_speed_kmph = unique(effort_speed_kmph),
                                               hour_of_day = unique(hours_of_day))

ott <- left_join(ottdiv, o_min, by = "checklist_id") %>% mutate(year = year(obs_date))

a <- ott %>% ggplot() + 
  geom_violin(aes(x = as.factor(year), y = specrich), fill = "skyblue") + 
  labs(x= "Year", y="Species richness") + theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

b <- ott %>% ggplot() + 
  geom_violin(aes(x = as.factor(year), y = shannondiv), fill = "yellow") + 
  labs(x= "Year", y="Shannon diversity") + theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

c <- ott %>% ggplot() + 
  geom_violin(aes(x = as.factor(year), y = evenness), fill = "coral") + 
  labs(x= "Year", y="Evenness") + theme_bw() +theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
(a | b | c) + plot_annotation(title = "Diversity metrics for complete eBird checklists in Ottawa-Gatineau by year")

ggsave("./output/figures/tar_fig_ottgat_chkldiv_byyear.png", plot = last_plot(), scale = 2)

### CHECK EFFORT VARIABLES
glimpse(ott)
plot_intro(ott)
plot_histogram(ott %>% select(c("effort_distance_km", "effort_hours", "effort_speed_kmph", "hour_of_day")),
               ggtheme = theme_bw())

### BIODIVERSITY MAPS ###
library(gganimate) ##CURRENTLY BROKEN, see 'transformr'
glimpse(ott)
colnames(ott, 10)

city <- urb_popcentres %>% slice(20) %>% sf::st_geometry() %>% sf::st_centroid()
lon_span <- 360 / 2^8
lat_span <- 180 / 2^8
lon_bounds <- c(sf::st_coordinates(city)[1] - lon_span / 2, sf::st_coordinates(city)[1] + lon_span / 2)
lat_bounds <- c(sf::st_coordinates(city)[2] - lat_span / 2, sf::st_coordinates(city)[2] + lat_span / 2)

q <- ott %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% mutate(year = year(obs_date)) 
q %>%   
  ggplot() +
    geom_sf(data = NA_basemaps$ne_land) +
    geom_sf(data = NA_basemaps$ne_state_lines) +
    geom_sf(data = NA_basemaps$ne_country_lines) +
    geom_sf(data = urb_popcentres %>% filter(PCNAME.x == "Ottawa - Gatineau"), fill= "pink") +
    geom_sf(mapping = aes(colour = specrich)) +
    viridis::scale_color_viridis(alpha=0.75, begin = 0.2, end = 1) +
    labs(title = "Species richness - Year: {frame_time}") + transition_time(obs_date) + ease_aes('linear') +
    coord_sf(xlim = lon_bounds, ylim = lat_bounds, crs = sf::st_crs(4326), datum = NA) +
    theme_bw() + theme(legend.position="bottom", legend.title = element_blank()) 

### TRY WITH MAPVIEW ###
## LAYERS BY YEAR?
q <- ott %>% sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% mutate(year = year(obs_date)) 
Y2020 <- q %>% filter(year == 2020)
og_sr <- mapview(Y2020, zcol = "specrich", legend = T, cex = 4, layer.name = "Species richness") %>% 
  removeMapJunk(junk = c("zoomControl", "layersControl", "homeButton", "drawToolbar"))
og_sh <- mapview(Y2020, zcol = "shannondiv", legend = T, cex = 4, layer.name = "Shannon diversity") %>% 
  removeMapJunk(junk = c("zoomControl", "layersControl", "homeButton", "drawToolbar"))
og_ev <- mapview(Y2020, zcol = "evenness", legend = T, cex = 4, layer.name = "Evenness") %>% 
  removeMapJunk(junk = c("zoomControl", "layersControl", "homeButton", "drawToolbar"))

m <- sync(og_sr, og_sh, og_ev, ncol = 3)
print(m)
