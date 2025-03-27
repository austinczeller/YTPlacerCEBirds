library(tidyverse)
library(climatenaR)
library(sf)
select <- dplyr::select


daysites<-read_rds("D:/CE_birds/2_pipeline/1_CountDataProcessing/out/sites.rds")
daysites <- st_as_sf(daysites)
#sites <- read_sf("CumulativeEffects_GIS/Sites/sites_analysis.shp")
visits <- read_rds("D:/CE_birds/2_pipeline/1_CountDataProcessing/out/visits.RDS")
loc <- daysites %>% group_by(siteID, year) %>% summarise() %>%
  st_centroid()  %>% st_transform(crs = "EPSG:4326")
loc$lon <- st_coordinates(loc)[,1]
loc$lat <- st_coordinates(loc)[,2]
loc$EPSG <- "EPSG:4326"

loc <- loc %>% left_join(visits %>% mutate(year = year(ts)) %>% select(siteID, year) %>% distinct())
###

download.file("https://daymet.ornl.gov/single-pixel/api/data?lat=63.86&lon=-139.22&vars=tmax,prcp&years=2024",
              "D:/CE_birds/0_data/daymet/daymettest.csv")

###  run download first matching unique lat, lon, year.



####weather that year####

#download a csv for the weather for that year
clim <- st_drop_geometry(loc) %>% select(siteID, lat, lon, year) %>% mutate(
                                               latitude = round(lat, digits = 3),
                                               longitude = round(lon, digits = 3),
                                               url=paste0("https://daymet.ornl.gov/single-pixel/api/data?lat=",latitude,
                                                          "&lon=",longitude,"&vars=tmax,prcp&years=",year),
                                               destfile = paste0("D:/CE_birds/0_data/daymet/", siteID,
                                                                 "year", year, ".csv")) %>% distinct() %>% ungroup()

for(s in 1:nrow(clim)){
  download.file(clim$url[s],clim$destfile[s],method="curl")
}


##append to habsum##
#what i want to do here is take the min and max ydays that we surveyed for each location using the visits dataframe 
#and then calculate the average maxtemp and precip during those ydays
visits$yday<-yday(visits$ts)

sitedays<-visits%>%group_by(siteID)%>%
  summarize(minday=min(yday),maxday=max(yday))

hab.sum <- read_rds("D:/CE_birds/2_pipeline/2_HabitatSum/out/hab_sum_site.RDS")


# Initialize a list to store calculated metrics for each site
weather_metrics <- list()

# Loop through downloaded Daymet files
for (file in clim$destfile) {
  # Read Daymet data (skip header rows and ensure proper column names)
  daymet_data <- read_delim(file, delim = ",", skip = 6) %>%
    rename(yday = "yday", tmax = "tmax (deg c)", prcp = "prcp (mm/day)")
  
  # Extract siteID and year from the filename
  file_info <- str_match(basename(file), "([^/]+)year(\\d+)\\.csv")
  siteID <- file_info[2]
  year <- as.numeric(file_info[3])
  
  # Get min and max yday for the site
  site_yday <- sitedays %>% filter(siteID == siteID)
  if (nrow(site_yday) > 0) {
    daymet_filtered <- daymet_data %>%
      filter(yday >= site_yday$minday & yday <= site_yday$maxday)
    
    # Calculate average tmax and prcp
    avg_metrics <- daymet_filtered %>%
      summarise(
        avg_tmax = mean(tmax, na.rm = TRUE),
        avg_prcp = mean(prcp, na.rm = TRUE)
      ) %>%
      mutate(siteID = siteID)
    
    # Append to weather_metrics
    weather_metrics[[siteID]] <- avg_metrics
  }
}

# Combine metrics into a single dataframe
weather_summary <- bind_rows(weather_metrics)

hab.sum <- hab.sum %>%
  left_join(weather_summary, by = "siteID")

write_rds(hab.sum,"D:/CE_birds/2_pipeline/2_HabitatSum/out/hab_sum_site.RDS")
