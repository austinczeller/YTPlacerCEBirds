library(ggplot2)
library(sf)
library(dplyr)
library(ggspatial)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# Define CRS for the entire analysis
crs_epsg <- 4326

# Load spatial data
surface <- read_sf("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS/YG_SurfaceDisturbance_May2022/ArealFeatures_MostRecent.shp")
linear <- read_sf("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS/YG_SurfaceDisturbance_May2022/LinearFeatures_MostRecent.shp")
locations <- readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/1_CountDataProcessing/out/sites.RDS")

# Filter survey years
locations <- locations %>% filter(!is.na(year))

# Ensure CRS matches
surface <- st_transform(surface, crs = crs_epsg)
linear <- st_transform(linear, crs = crs_epsg)
locations <- st_transform(locations, crs = crs_epsg)

locations <- locations %>% mutate(year = as.numeric(year))%>%filter(!is.na(year))

# Define bounding box
bbox_manual <- 
  st_bbox(c(xmin = -140.8, xmax = -135.5, ymin = 62.8, ymax = 64.5), crs = 4326) %>% 
  st_as_sfc() %>%
  st_transform(crs = crs_epsg)

# Define towns and transform to map CRS
towns <- data.frame(
  name = c("Dawson", "Mayo", "Pelly Crossing"),
  lon = c(-139.432, -135.896, -136.583),   
  lat = c(64.060, 63.595, 62.820)          
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = crs_epsg)

# Load ecoregions data
fgdb <- file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS/Ecoregions_2014_1M.gdb")
ecor <- st_as_sf(vect(fgdb, layer="Ecoregions_2014_1M")) %>%
  st_transform(crs = crs_epsg)

ecor <- ecor %>% filter(ECOREGION_ID %in% c(170, 172, 175, 176, 302)) # Filter ecoregions

# Define a custom green color gradient
custom_greens <- c("#00441b", "#1b7837", "#5aae61", "#a6dba0", "#93C572")

# Plot main map
map <- ggplot() +
  # Ecoregions with a custom green color scale
  geom_sf(data = ecor, aes(fill = ECOREGION)) + 
  scale_fill_manual(name = "Ecoregions", values = custom_greens) +
  
  # Surface and linear features with darker colors for visibility
  geom_sf(data = surface, fill = "#895129", color = "#895129", linewidth = 0.6, aes(group = 1, linetype = "Surface")) +
  geom_sf(data = linear, color = "#895129", linewidth = 0.7, aes(group = 2, linetype = "Linear")) +
  
  # Locations and towns
  geom_sf(data = locations, aes(color = year), size = 2) +
  geom_sf(data = towns, shape = 21, fill = "black", color = "white", size = 3) +  
  
  # Town labels
  geom_text(data = towns, aes(label = name, geometry = geometry), stat = "sf_coordinates", 
            size = 4, fontface = "bold", hjust = 1.2, vjust = 1, color = "black") +
  
  # Scale for survey year colors with a blue gradient
  scale_color_gradient(name = "Survey Year", low = "#A7C7E7", high = "#08306b") +
  
  # Scale for linetypes in the surface and linear features
  scale_linetype_manual(name = "Disturbance", values = c("Surface" = "solid", "Linear" = "solid")) +
  
  # Map details
  annotation_scale(location = "br", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering()) +
  coord_sf(xlim = c(st_bbox(bbox_manual)$xmin, st_bbox(bbox_manual)$xmax), 
           ylim = c(st_bbox(bbox_manual)$ymin, st_bbox(bbox_manual)$ymax)) +
  
  # Customize the theme
  theme_void() +
  theme(legend.position = "right")

# Convert bbox_manual to sf format
study_area <- st_as_sf(bbox_manual)

# Inset Map: Canada with Yukon Highlighted
canada <- ne_states(country = "Canada", returnclass = "sf")
yukon <- canada %>% filter(name == "Yukon")

inset_map <- ggplot() +
  geom_sf(data = canada, fill = "gray80", color = "white") +
  geom_sf(data = yukon, fill = "gray90", color = "black") +
  geom_sf(data = study_area, fill = "red", color = "black", alpha = 0.7) +
  theme_void()

# Adjust inset map size and position
final_map <- map + inset_element(inset_map, 
                                 left = 0.001, bottom = 0.001,  # Move to bottom left
                                 right = 0.30, top = 0.30)  # Increase size more

plot(final_map)
ggsave('C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2023and2024studymap.png', final_map)
