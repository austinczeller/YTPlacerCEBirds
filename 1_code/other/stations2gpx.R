setwd("E:/FieldSeason2023")
library(tidyverse)
library(sf)
library(readxl)	

opts <- st_read("StationPts/MR2023_FINALclean3/MR2023_FINALclean3.shp")
# st_write(opts %>% mutate(name = StationID) %>% 
#            st_transform(crs = 4326) %>% select(name),
#          "stations2023_original.gpx", driver = "GPX", delete_dsn = TRUE)


fname <- "WCS2023_DataEntry.xlsx"
sheets <- readxl::excel_sheets(fname)
stations <- read_excel(fname, sheet = sheets[1])

pts <- st_as_sf(stations, coords = c("Easting", "Northing"))
st_crs(pts) <- 3155

opts1 <- filter(opts, !SiteID %in% unique(pts$Site))

st_write(pts %>% mutate(stationID = paste(Site, Station, sep = "-")), 
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_shift1.shp", 
         delete_dsn = TRUE)



st_write(pts %>% mutate(name = paste(Site, Station, sep = "-")), "stations2023_shift1.kml", driver = "kml", delete_dsn = TRUE)
st_write(rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
           st_transform(crs = 4326) %>% select(name),
           opts1 %>% mutate(name = StationID) %>% 
             st_transform(crs = 4326) %>% select(name)),
         "stations2023_shift1_all.gpx", driver = "GPX", delete_dsn = TRUE)

heli.pts <- st_read("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/HeliStations_S2_2023.shp")

st_write(heli.pts %>% mutate(name = paste("heli", n(), sep = "-")) %>% 
                 st_transform(crs = 4326) %>% select(name),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_heli.gpx", driver = "GPX", delete_dsn = TRUE)

pts <- st_read("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/MR2023_FINALclean3.shp")

st_write(pts %>% mutate(name = paste(SiteID, Station, sep = "-")) %>% 
           st_transform(crs = 4326) %>% select(name),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_all.gpx", driver = "GPX", delete_dsn = TRUE)


#######shift 2 ----
setwd("C:/Users/jbrow/My Drive/WCSC/MineReclaim/FieldSeason2023")
fname <- "WCS2023_DataEntry_hd2.xlsx"
sheets <- readxl::excel_sheets(fname)
stations <- read_excel(fname, sheet = sheets[1])
stations <- stations %>% rename(UTMzone = "UTM Zone")
unique(stations$UTMzone)

#07V  #3154
#07W #3154
#08W  

pts <- st_as_sf(stations %>% filter(UTMzone == "08W"), coords = c("Easting", "Northing"))
st_crs(pts) <- 3155
pts2 <- st_as_sf(stations %>% filter(UTMzone != "08W"), coords = c("Easting", "Northing"))
st_crs(pts2) <- 3154

st_write(rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
           st_transform(crs = 4326) %>% select(name), 
           pts2 %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
             st_transform(crs = 4326) %>% select(name)),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_shift2.gpx",
         driver = "GPX", delete_dsn = TRUE)

opts2 <- filter(opts, !SiteID %in% unique(pts$Site))

st_write(rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                 st_transform(crs = 4326) %>% select(name), 
               pts2 %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                 st_transform(crs = 4326) %>% select(name)), 
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_shift2.shp", 
         delete_dsn = TRUE)

st_write(rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                 st_transform(crs = 4326) %>% select(name), 
               pts2 %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                 st_transform(crs = 4326) %>% select(name)), 
         "stations2023_shift2.kml", driver = "kml", delete_dsn = TRUE)

pts <- rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
               st_transform(crs = 4326), 
             pts2 %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
               st_transform(crs = 4326))
heli_kml <- pts %>% mutate(Date = ymd(Date)) %>% filter(Date %in% ymd("2023-06-11 UTC", "2023-06-12 UTC"))
st_write(heli_kml %>% mutate(name = paste(Site, Station, sep = "-")), "helipts.kml", driver = "kml", delete_dsn = TRUE)


#########shift 3 ----
fname <- "WCS2023_DataEntry.xlsx"
sheets <- readxl::excel_sheets(fname)
stations <- read_excel(fname, sheet = sheets[1])
stations <- stations %>% rename(UTMzone = "UTM Zone")
unique(stations$UTMzone)

#07V  #3154
#07W #3154
#08W  

pts <- st_as_sf(stations %>% filter(UTMzone %in% c("08W", "08V", "08v")), coords = c("Easting", "Northing"))
st_crs(pts) <- 3155
pts2 <- st_as_sf(stations %>% filter(! UTMzone  %in% c("08W", "08V", "08v")), coords = c("Easting", "Northing"))
st_crs(pts2) <- 3154

aru3 <- rbind(pts %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                st_transform(crs = 4326), 
              pts2 %>% mutate(name = paste(Site, Station, sep = "-")) %>% 
                st_transform(crs = 4326))
aru3 <- filter(aru3, Site %in% c("MR900","MR901","MR902","MR903","MR904","MR905","MR906","MR907","MR908",
                         "MR508","MR17","MR70","MR69","MR509","MR5","MR74","MR502","MR73","MR15"))  


st_write(aru3%>% 
           select(name),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_shift3.gpx",
         driver = "GPX", delete_dsn = TRUE)


st_write(aru3%>% 
           select(name), 
         "stations2023_shift3.kml", driver = "kml", delete_dsn = TRUE)

st_write(aru3 %>% rename(stationID = name), 
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/stations2023_shift3.shp", 
         delete_dsn = TRUE)

HeliDrop <- read_csv("HeliDropPts.csv")

HeliDrop <- st_as_sf(HeliDrop, coords = c("lon", "lat"))
st_crs(HeliDrop) <- 4326
st_write(HeliDrop, 
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/HeliDrop.shp", 
         delete_dsn = TRUE)

HeliDrop <- st_read("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/HeliDrop.shp")
st_write(HeliDrop%>% 
           select(name=Drop),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/HeliDrop.gpx",
         driver = "GPX", delete_dsn = TRUE)

st_write(HeliDrop%>% 
           select(name=Drop),
         "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/HeliDrop.kml",
         driver = "kml", delete_dsn = TRUE)
