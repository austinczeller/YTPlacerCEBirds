library(tidyverse)
library(climatenaR)
library(sf)
select <- dplyr::select
# setwd("C:/Users/jbrow/My Drive/WCSC/MineReclaim/ReclaimChrono")
#
#
# Sys.setenv(WT_USERNAME = 'jbrow247@alumni.uwo.ca', WT_PASSWORD = 'Ridpearb3!')
# wt_auth()

# my_projects <- wt_get_download_summary(sensor_id = 'ARU')
# my_projects ##project ID 2078 = SBPPCE 2023, SBPPCE ID 577 = 2021

sites <- read_sf("E:/Archive/2_pipeline/2_HabitatSum/store/sites_analysis.shp")
visits <- read_rds("E:/Archive/2_pipeline/1_CountDataProcessing/out/visits.RDS")
loc <- sites %>% group_by(siteID, year) %>% summarise() %>%
  st_centroid()  %>% st_transform(crs = "EPSG:4326")
loc$lon <- st_coordinates(loc)[,1]
loc$lat <- st_coordinates(loc)[,2]
loc$EPSG <- "EPSG:4326"

loc <- loc %>% left_join(visits %>% mutate(year = year(ts)) %>% select(siteID, year) %>% distinct())
#
#download.file("https://daymet.ornl.gov/single-pixel/api/data?lat=63.86&lon=-139.22&vars=tmax&years=2023",
#             "0_data/daymettest.csv", method = "curl")

###  run download first matching unique lat, lon, year.



clim <- st_drop_geometry(loc) %>% select(siteID, lat, lon, year) %>% mutate(
  latitude = round(lat, digits = 3),
  longitude = round(lon, digits = 3),
  destfile = paste0("0_data/daymet/lat", latitude, "lon", longitude,
                    "year", year, "5yr.csv")) %>%
  
  distinct() %>% ungroup()




##spring
# files <- list.files("0_data/daymet/")
# files <- paste0("0_data/daymet/", files)
# clim %>% select(latitude, longitude, year, destfile) %>%
#   pmap(function(latitude, longitude, year, destfile) {
# if(!destfile %in% files) {
#   download.file(paste0("https://daymet.ornl.gov/single-pixel/api/data?lat=", latitude, "&lon=", longitude,
#                        "&vars=tmax,prcp,swe&start=", year, "-04-01&end=", year,"-06-30"),
#                 destfile = destfile)
# 
#   Sys.sleep(runif(1, 1, 5))
#   }
# 
#   }
#   )

# download.file(paste0("https://daymet.ornl.gov/single-pixel/api/data?lat=", 
#                      filter(clim, destfile == "0_data/daymet/lat64.012lon-135.298year2021spring.csv") %>% pull(latitude), "&lon=", 
#                      filter(clim, destfile == "0_data/daymet/lat64.012lon-135.298year2021spring.csv") %>% pull(longitude),
#                      "&vars=tmax,prcp,swe&start=", 
#                      filter(clim, destfile == "0_data/daymet/lat64.012lon-135.298year2021spring.csv") %>% pull(year), "-04-01&end=", 
#                      filter(clim, destfile == "0_data/daymet/lat64.012lon-135.298year2021spring.csv") %>% pull(year),"-06-30"),
#               destfile = filter(clim, destfile == "0_data/daymet/lat64.012lon-135.298year2021spring.csv") %>% pull(destfile))    

## create summary climate variables for each unique lat, long, year, saves as a single .csv
##summary .csv can be imported into hab.sum and matched to non-unique station lat/long/year

clim <- st_drop_geometry(loc) %>% select(siteID, lat, lon, year) %>% mutate(
  latitude = round(lat, digits = 3),
  longitude = round(lon, digits = 3),
  destfile = paste0("E:/Archive/0_data/daymet/lat", latitude, "lon", longitude,
                    "year", year, "june.csv")) %>%
  distinct() %>% ungroup()

##5 year data
files <- list.files("E:/Archive/0_data/daymet/")
files <- paste0("E:/Archive/0_data/daymet/", files)
clim %>% select(latitude, longitude, year, destfile) %>%
  pmap(function(latitude, longitude, year, destfile) {
    if(!destfile %in% files) { 
      download.file(paste0("https://daymet.ornl.gov/single-pixel/api/data?lat=", latitude, "&lon=", longitude,
                           "&vars=tmax,tmin,prcp,swe&years=", paste(seq(year-4, year, 1), collapse=",")),
                    destfile = destfile)
      
      Sys.sleep(runif(1, 1, 5))
    } 
  })
day

clim <- left_join(clim, map_df(1:length(clim$destfile), function(i) {
  daymet <- read_csv(clim$destfile[i], skip = 7) %>%
    rename(precip = `prcp (mm/day)`, snow = `swe (kg/m^2)`, tmax = `tmax (deg c)`, tmin = `tmin (deg c)`)
  
  mn_june_tmax <- mean(filter(daymet, yday %in% 152:181)$tmax)
  mn_june_tmin <- mean(filter(daymet, yday %in% 152:181)$tmin)
  mn_june_precip <- sum(filter(daymet, yday %in% 152:181)$precip)/5
  mn_total_precip <- sum(daymet$precip)/5
  mn_last_snow <- daymet %>% group_by(year) %>% summarise(last_snow = min(which(snow==0))) %>% 
    summarise(mn_last_snow = mean(last_snow)) %>% pull(mn_last_snow)
  
  
  data.frame(siteID = clim$siteID[i], mn_june_tmax, mn_june_tmin, mn_june_precip, mn_total_precip, mn_last_snow, destfile = clim$destfile[i])
}))
saveRDS(clim, "E:/Archive/0_data/daymet_5yr.rds")


#
# files <- list.files("0_data/daymet")
#
# files = paste0("0_data/daymet/", files)
#
# missing.clim <- clim[!clim$destfile %in% files,]


##spring means. Last snow seems late, might be better to use modis to extract following caribou parturition paper

### variable descriptions can be found here: https://climatebc.ca/help/climateBC/help.htm#_Toc410137601
### Directly calculate annual variables
# MAT = Mean annual temp
# MWMT mean warmest month temp
# MCMT mean coldest month temp
# TD temp diff between MWMT and MCMT
# MAP  mean annual precipitation (mm),
### Derived annual variables:
# NFFD              the number of frost-free days
# FFP                 frost-free period
# bFFP               the day of the year on which FFP begins
# PAS                 precipitation as snow (mm) between August in previous year and July in current year
### Directly calculated seasonal variables:Spring (_sp): March, April and May
# Tave_sp           spring mean temperature (°C)
# PPT_sp            spring precipitation (mm)


