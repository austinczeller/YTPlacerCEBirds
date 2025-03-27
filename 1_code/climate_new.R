#get weather data
library(tidyverse)
library(weathercan)

hab.sum <- readRDS("E:/Archive/2_pipeline/2_HabitatSum/out/hab_sum_site.rDS")%>%filter(year!=0)
for(w in 1:length(hab.sum)){
  lat<-hab.sum$lat[w]
  lon<-hab.sum$lon[w]
  coords<-c(lat,lon)
  nearest_station<-stations_search(coords=coords)[1,]
  sdate<-ymd(paste0(hab.sum$year[w],"06-01"))
  edate<-ymd(paste0(hab.sum$year[w],"06-30"))
  weather<-weather_dl(station_ids = nearest_station$station_id,start = sdate, end=edate,interval = "day")
  normals<-normals_dl(nearest_station$climate_id)
  

  }

  weathercan::stations_search(coords=coords)
                            