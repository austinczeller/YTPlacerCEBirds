library(stars)
library(sf)
library(tidyverse)

mode <- function(v, na.rm = T) {
  v <- as.vector(v)
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

setwd("C:/Users/jbrow/My Drive/WCSC/MineReclaim/ReclaimChrono")
sites <- read_csv("2_pipeline/dataSummary2023.csv")

sites <- sites %>% mutate(id = paste(site, station, sep = "-"),
                          utm = as.numeric(sub("[a-z]", "", str_to_lower(utm))))
sites <- rbind(st_as_sf(sites %>% filter(utm == 7),
                      coords = c("easting", "northing"),
                      crs = "+proj=utm +zone=7") %>%
               st_transform(crs = "+proj=longlat +datum=NAD83"),
             st_as_sf(sites %>% filter(utm == 8),
                      coords = c("easting", "northing"),
                      crs = "+proj=utm +zone=8") %>%
               st_transform(crs = "+proj=longlat +datum=NAD83"))


reveg <- read_stars("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/RevegetationMapping/Reveg_yod.tif")
crs <- st_crs(reveg)
names(reveg) <- "dyear"
# mlf <- read_stars("CumulativeEffects_GIS/macrolandforms.tif", proxy = F) ##macro landforms. Use this crs
# crs <- st_crs(mlf)

MRstat <- st_transform(sites, crs = crs)

SDisturb <- read_sf("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/CumulativeEffects_GIS/YG_SurfaceDisturbance_May2022/ArealFeatures_May2022_Un_Dis_J.shp")
SDisturb <- st_transform(SDisturb, crs = crs)

##some sites should be grouped into a different hex (either mines span hex boarder, or no available sites of different age within single hex.

##update clustering based on reclaimed site selection ----
# MRstat$ClusterID <- NA
# MRstat[MRstat$SiteID %in% c(59),]$ClusterID <- 1
# MRstat[MRstat$SiteID %in% c(68),]$ClusterID <- 2
# MRstat[MRstat$SiteID %in% c(58,60),]$ClusterID <- 3
# MRstat[MRstat$SiteID %in% c(48, 49),]$ClusterID <- 4
# MRstat[MRstat$SiteID %in% c(47,67,46),]$ClusterID <- 5
# MRstat[MRstat$SiteID %in% c(43,44,45),]$ClusterID <- 6
# MRstat[MRstat$SiteID %in% c(40,41,42),]$ClusterID <- 7
# MRstat[MRstat$SiteID %in% c(37,38),]$ClusterID <- 8
# MRstat[MRstat$SiteID %in% c(39),]$ClusterID <- 9
# MRstat[MRstat$SiteID %in% c(35,34,36,66),]$ClusterID <- 10
# MRstat[MRstat$SiteID %in% c(56,57),]$ClusterID <- 11
# MRstat[MRstat$SiteID %in% c(32,31,62,33),]$ClusterID <- 12
# MRstat[MRstat$SiteID %in% c(28,29,30),]$ClusterID <- 13
# MRstat[MRstat$SiteID %in% c(26,27),]$ClusterID <- 14
# MRstat[MRstat$SiteID %in% c(63,25,24,23),]$ClusterID <- 15
# MRstat[MRstat$SiteID %in% c(65,22,21),]$ClusterID <- 16
# MRstat[MRstat$SiteID %in% c(16,17),]$ClusterID <- 17
# MRstat[MRstat$SiteID %in% c(12,13,15),]$ClusterID <- 18
# MRstat[MRstat$SiteID %in% c(11,10,8),]$ClusterID <- 19
# MRstat[MRstat$SiteID %in% c(7,6),]$ClusterID <- 20
# MRstat[MRstat$SiteID %in% c(18,20,19),]$ClusterID <- 21
# MRstat[MRstat$SiteID %in% c(64,4),]$ClusterID <- 22
# MRstat[MRstat$SiteID %in% c(5,3),]$ClusterID <- 23
# MRstat[MRstat$SiteID %in% c(50, 2, 1),]$ClusterID <- 24
# MRstat[MRstat$SiteID %in% c(61,55,54),]$ClusterID <- 25
# MRstat[MRstat$SiteID %in% c(52,53),]$ClusterID <- 26


reveg.agg <- aggregate(reveg, st_buffer(MRstat, dist = 150), median,na.rm = T)
MRstat <- mutate(MRstat, dyear = reveg.agg$dyear)

##group by site
MRstat.l <- split(MRstat, f = MRstat$site)

MRsite_reveg <- lapply(MRstat.l, function(site) {
  site <- st_buffer(site, dist = 150) %>% st_union
  bbox <- st_bbox(st_buffer(st_union(site), 100))
  reveg.i <-  st_crop(reveg, bbox)
  reveg.i <- aggregate(reveg.i, site, median, na.rm = T)
  reveg.i$dyear
})


#plot reveg with stations and disturbance buffer, and compare to ArcGIS imagery
MRsite_reveg <- lapply(MRstat.l, function(site) {
  site <- st_buffer(site, dist = 150)
  reveg.i <- st_crop(reveg, bbox)

  ggplot() +
    geom_stars(data = reveg.i) +
    geom_sf(data = site, fill = NA, col = "red")
})






##is whole area disturbed, or is there digitization error?
##is the site the same age, or do some stations represent two different disturbance eras?
#Also calculate % active? - may need to use NDVI to distinguish active from old veg

## siteID = 2----
#mix of reveg with small active areas (NA). Reveg mostly covered by mapping (>1985) 
##Assign veg age based on maping of all sites. 
# site <- MRstat.l[["2"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# ##make histograms of reveg age of each site, or groups of stations if of different ages, and assign disturbance year to each statoin. 
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)

##mode likely best estimate
MRstat$veg_age <- NA
MRstat[MRstat$SiteID == 2,]$veg_age <- "1987" #mode

## siteID = 3 ----
## 3 different site ages
# site <- MRstat.l[["3"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# ##make histograms of reveg age of each site, or groups of stations if of different ages, and assign disturbance year to each statoin. 
# 
# ## 162 age 1. Moving to capture different age. 
# reveg.i <- st_crop(reveg, st_buffer(site %>% filter(stationID == 162), 150))
# ##make histograms of reveg age of each site, or groups of stations if of different ages, and assign disturbance year to each statoin. 
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID == 162,]$veg_age <- "2009" #mean, median, mode
#FID 157, 156 = age 2
MRstat[MRstat$stationID %in% c(157,156),]$veg_age  <- "<1985" 
## Remainder = age 3
# reveg.i <- st_crop(reveg, st_buffer(site %>% filter(!stationID %in% c(162, 157,156)), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 3 & is.na(MRstat$veg_age),]$veg_age <- "2000" #mean, median, mode

## siteID = 4 ----
#Recent, with some patches still active and others regrown
# site <- MRstat.l[["4"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# ##make histograms of reveg age of each site, or groups of stations if of different ages, and assign disturbance year to each statoin. 
# hist(reveg.i) ##FID 136 is being mapped as old, but should be same age as others
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 4,]$veg_age <- "2016"

## siteID = 5----
# site <- MRstat.l[["5"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
##make histograms of reveg age of each site, or groups of stations if of different ages, and assign disturbance year to each statoin.
##seems to be a mix of reminant (undisturbed), still disturbed, and 'reclaimed' 
##GIS imagery from 2009 and appears active, sentinel 2021 show some reveg.  2016 might actually be most accurate, but may need some adjustments/manual trend exploration. 
# hist(reveg.i) 
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 5,]$veg_age <- "unk"

## siteID = 6----
# site <- MRstat.l[["6"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# ##definite activity since 2011. some areas still active/bare. Mix of ages... assign based on Station? 
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)

MRstat[MRstat$stationID == 174,]$veg_age <- ">2021"
MRstat[MRstat$stationID %in% c(171:173),]$veg_age <- "unk"

tmp <- sapply(c(167:170, 175), function(x){
  reveg.i <- st_crop(reveg, st_buffer(MRstat.l[["6"]] %>% filter(stationID == x), 60))
  mode(reveg.i$dyear)
})
tmp <- as.character(tmp)
MRstat[MRstat$stationID == 167,]$veg_age <- tmp[1]
MRstat[MRstat$stationID == 168,]$veg_age <- tmp[2]
MRstat[MRstat$stationID == 169,]$veg_age <- tmp[3]
MRstat[MRstat$stationID == 170,]$veg_age <- tmp[4]
MRstat[MRstat$stationID == 175,]$veg_age <- tmp[5]

## siteID = 7----
# site <- MRstat.l[["7"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# Disturbance hugs river closely, with surveys in or near remnant vegetation
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 7,]$veg_age <- "1994"

## siteID = 8----
# site <- MRstat.l[["8"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
##mix of older veg with open areas, reasonably uniform. Some older veg in SW area
## sum per station
# plot(reveg.i)
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat, stationID %in% c(380,  386)), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(380,  386),]$veg_age <- "1990"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat, stationID %in% c(382)), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(382),]$veg_age <- "2007"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat, stationID %in% c(381, 383:385, 387:393)), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(381, 383:385, 387:393),]$veg_age <- "1997"

## siteID = 10----
# site <- MRstat.l[["10"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
MRstat[MRstat$stationID %in% c(395,397,398),]$veg_age <- ">2021"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat, stationID %in% c(394,396)), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(394,396),]$veg_age <- "2012"

## siteID = 11----
# site <- MRstat.l[["11"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 11,]$veg_age <- "1997"

## siteID = 12----
# site <- MRstat.l[["12"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)

tmp <- sapply(c(331, 334:335, 337:345), function(x){
  reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 150))
  mode(reveg.i$dyear)
})
tmp <- as.character(tmp)
MRstat[MRstat$stationID %in% c(331, 334:335, 337:345),]$veg_age <- tmp
MRstat[MRstat$stationID %in% c(332,336),]$veg_age <- ">2021"

## siteID = 13----
# site <- MRstat.l[["13"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i) ## two peaks, plus some old veg = 3 ages. 
MRstat[MRstat$stationID %in% c(348,352,356,358),]$veg_age <- "<1985"
# tmp <- sapply(c(346,347,349,350,351,353,354,355,357), function(x){
#   reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 200))
#   mode(reveg.i$dyear)
# })
# tmp 
MRstat[MRstat$stationID %in% c(346,347,350,351),]$veg_age <- "1990"
MRstat[MRstat$stationID %in% c(349,353,354,355,357),]$veg_age <- "1997"

## siteID = 15----
# site <- MRstat.l[["15"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i, breaks = 20) ## wide range, no clear pattern
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 15,]$veg_age <- "1992"

## siteID = 16----
MRstat[MRstat$SiteID == 16,]$veg_age <- ">2021"

## siteID = 17----
# site <- MRstat.l[["17"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i, breaks = 20)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$SiteID == 17,]$veg_age <- "1986"

# ## siteID = 18----
# site <- MRstat.l[["18"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
MRstat[MRstat$stationID %in% c(57,58,63),]$veg_age <- ">2021"

v <- filter(MRstat, SiteID == 18 & !stationID %in% c(57,58,63)) %>% pull(stationID)
# 
# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% v), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear, na.rm = T)
##probably 2003, but highly variable, may need some tweeking
MRstat[MRstat$stationID %in% v,]$veg_age <- "unk"

## siteID = 19 & 20----
# site <- MRstat.l[["19"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
MRstat[MRstat$stationID %in% c(69,70,71),]$veg_age <- ">2021"
# tmp <- sapply(c(67,68,72), function(x){
#   reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 200))
#   mode(reveg.i$dyear)
# })
# tmp ## 68 was clearly bare in 2017, need to remap
MRstat[MRstat$stationID %in% c(67,68,72),]$veg_age <- "unk"

# site <- MRstat.l[["20"]]
# reveg.20 <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.20)
# v <- c(78,79,80)
# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% v), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear, na.rm = T) ##mode, with some more recent disturbance
MRstat[MRstat$stationID %in% 78:80,]$veg_age <- "1987"
# v <- 73:77
# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% v), 150))
# hist(reveg.i)
# mean(reveg.i$dyear, na.rm = T)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% 73:77,]$veg_age <- "1998"

## siteID = 21----
# site <- MRstat.l[["21"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
#  ## distinct age groups?
# tmp <- sapply(c(17:19, 21), function(x){
#   reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 60))
#   mode(reveg.i$dyear)
# })
# tmp
##17 amd 21 recently disturbed in 2017 (maxar), 18 and 19 likely in forest patches of similar age
##Vegmapping needs more attention
MRstat[MRstat$SiteID == 21,]$veg_age <- "unk"

## siteID = 22----
# site <- MRstat.l[["22"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
MRstat[MRstat$stationID %in% 305:306,]$veg_age <- ">2021"

# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% 301:304), 150))
# hist(reveg.i)
##pre 2000 is wrong based on imagery, most of area recently disturbed in 2017
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
MRstat[MRstat$stationID %in% 301:304,]$veg_age <- "2014" 



## siteID = 23----
MRstat[MRstat$SiteID == 23,]$veg_age <- ">2021" 

## siteID = 24----
# site <- MRstat.l[["24"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)

# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% 311:312), 150))
# hist(reveg.i)
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
# MRstat[MRstat$stationID %in% 311:312,]$veg_age <- "2004" 
# 
# reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID %in% 307:310), 150))
# hist(reveg.i, breaks = 20) ##activity drops off 
# median(reveg.i$dyear, na.rm = T)
# mode(reveg.i$dyear)
# MRstat[MRstat$stationID %in% 307:310,]$veg_age <- "1995" 
MRstat[MRstat$SiteID == 24,]$veg_age <- "unk" 

## siteID = 25----
# site <- MRstat.l[["25"]]
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
# median(reveg.i$dyear, na.rm = T)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 25,]$veg_age <- "1999" ##disturbance peaters out after this date

## siteID = 26----
MRstat[MRstat$SiteID == 26,]$veg_age <- ">2021" 

## siteID = 27----
# site <- MRstat.l[["27"]] ## veg looks uniform, lots of variation in mapping
# reveg.i <- st_crop(reveg, st_buffer(site, 150))
# hist(reveg.i)
# median(reveg.i$dyear, na.rm = T)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 27,]$veg_age <- "unk" 

## siteID = 28----
MRstat[MRstat$SiteID == 28,]$veg_age <- ">2021" 

## siteID = 29----
MRstat[MRstat$SiteID == 29,]$veg_age <- ">2021" 

## siteID = 30----
MRstat[MRstat$stationID %in% c(367, 369, 371, 372, 379),]$veg_age <- "<1985"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(376,378)), 150))
# hist(reveg.i)
# table(reveg.i$dyear) ##2018?
MRstat[MRstat$stationID %in% c(376,378),]$veg_age <- "2018"
# 
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(366,368,370,373,374,375,377)), 150))
# hist(reveg.i)
# table(reveg.i$dyear) 
MRstat[MRstat$stationID %in% c(366,368,370,373,374,375,377),]$veg_age <- "1987"

##siteID = 31---
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(425,426)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(425,426),]$veg_age <- "1987"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(422:424)), 150))
# hist(reveg.i)
# table(reveg.i$dyear) ##edge area = 2002, most of site still open (NA)
MRstat[MRstat$stationID %in% c(422:424),]$veg_age <- ">2021"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(427:429)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(427:429),]$veg_age <- "2000"

##siteID = 32----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 32), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 32,]$veg_age <- "1987" 

##siteID = 33----
MRstat[MRstat$SiteID == 33,]$veg_age <- ">2021"

###siteID = 34----
MRstat[MRstat$stationID %in% c(281,282),]$veg_age <- "<1985"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(283:285)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(283:285),]$veg_age <- "1990" ##some active mining nearby

##siteID = 35----
MRstat[MRstat$SiteID == 35,]$veg_age <- ">2021"

##siteID = 36----
MRstat[MRstat$SiteID == 36,]$veg_age <- ">2021"

##siteID = 37----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID ==37), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$SiteID == 37,]$veg_age <- "1994" ## several consistent areas of 1993-1996

##siteID = 38----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID ==38), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(258, 260:262),]$veg_age <- "1995"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(263:267)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(263:267),]$veg_age <- "1986"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(259)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(259),]$veg_age <- "1998"

## siteID 39
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 39), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
# 
# tmp <- sapply(c(268:272), function(x){
#   reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 150))
#   median(reveg.i$dyear, na.rm = T)
# })
# tmp
# tmp1 <- filter(MRstat, stationID %in% c(268:272))
# tmp1$year <- tmp
# ggplot() + geom_sf(data = tmp1, aes(col = year))
# 
reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(268:270)), 4000))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)


# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(271:272)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$SiteID == 39,]$veg_age <- "2011"

## siteID 40----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 40), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
# 
# tmp <- sapply(c(222:227), function(x){
#   reveg.i <- st_crop(reveg, st_buffer(MRstat %>% filter(stationID == x), 150))
#   median(reveg.i$dyear, na.rm = T)
# })
# tmp 
MRstat[MRstat$SiteID == 40,]$veg_age <- "2002" ##median. Might be different ages, but imagery looks even. 

## site 41----
MRstat[MRstat$stationID %in% c(228:230),]$veg_age <- ">2021"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(231:232)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(231:232),]$veg_age <- "2009"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(233)), 200))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(233),]$veg_age <- "2002" ##could also be 1997

## siteID 42----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(278:279)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(278:279),]$veg_age <- "1994" 
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(280)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(280),]$veg_age <- "1986" 

## siteID 43----
MRstat[MRstat$SiteID == 43,]$veg_age <- "<1985"


## siteID 42----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(278:279)), 150))
# hist(reveg.i)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$stationID %in% c(278:279),]$veg_age <- "1994" 

## siteID 43----
MRstat[MRstat$SiteID == 43,]$veg_age <- "<1985"

## siteID 44----
MRstat[MRstat$SiteID == 44,]$veg_age <- ">2021"

## siteID 45----
MRstat[MRstat$SiteID == 45,]$veg_age <- "<1985"

## siteID 46----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 46), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
# median(reveg.i$dyear, na.rm = T)
MRstat[MRstat$SiteID == 46,]$veg_age <- "1993"

## siteID 47----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 47), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 47,]$veg_age <- "1995"

## siteID 48----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 48), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID == 210), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID == 210,]$veg_age <- "2011"
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(207:209, 212:214)), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(207:209, 212:214),]$veg_age <- "1990"

## siteID 49----
MRstat[MRstat$SiteID == 49,]$veg_age <- ">2021"

## siteID 50----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 50), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(149,150,153)), 100))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 50,]$veg_age <- "unk" ##1998

## siteID 58----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 58), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$SiteID == 58,]$veg_age <- "1987" ##1998

## siteID 59----
MRstat[MRstat$SiteID == 59,]$veg_age <- "unk"

## siteID 60----
MRstat[MRstat$SiteID == 60,]$veg_age <- "<1985"

## siteID 62----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(444:447)), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(444:447),]$veg_age <- "2004"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(448)), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(448),]$veg_age <- "1997" ## based on values north in the same mine stretch.  Should maybe be dropped from site. 


## siteID 63----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(115:117)), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(115:117),]$veg_age <- "<1985"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(118)), 150))
# hist(reveg.i, breaks = 10)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(118),]$veg_age <- "2012"

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(114)), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)
MRstat[MRstat$stationID %in% c(114),]$veg_age <- "2003" ## based on exploring raster, recent disturbance along half of site in 2012, and older veg along other half. c
# more recent

## siteID 64----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID == 64), 150))
# hist(reveg.i, breaks = 20)
# table(reveg.i$dyear)

# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(81:82)), 150))
# hist(reveg.i, breaks = 10)
# table(reveg.i$dyear) ## 95-2000
# 
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(84:85)), 150))
# hist(reveg.i, breaks = 7)
# table(reveg.i$dyear)

MRstat[MRstat$SiteID == 64,]$veg_age <- "1999" 

## siteID 65----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(328,329)), 150))
# hist(reveg.i, breaks = 10)
# table(reveg.i$dyear) ## 95-2000
MRstat[MRstat$stationID %in% c(328:329),]$veg_age <- "2011"
MRstat[MRstat$stationID %in% c(326:327),]$veg_age <- "ref?"
MRstat[MRstat$stationID %in% c(325),]$veg_age <- "<1985"

## siteID 66----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  SiteID %in% 66), 150))
# hist(reveg.i, breaks = 10)
# table(reveg.i$dyear) 
MRstat[MRstat$SiteID == 66,]$veg_age <- "1991"

## siteID 67----
MRstat[MRstat$SiteID == 67,]$veg_age <- ">2021"

## siteID 65----
# reveg.i <- st_crop(reveg, st_buffer(filter(MRstat,  stationID %in% c(328,329)), 150))
# hist(reveg.i, breaks = 10)
# table(reveg.i$dyear) ## 95-2000


#################
as.data.frame(MRstat) %>% select(SiteID, veg_age) %>% distinct() %>% pull(veg_age) %>% table

MRstat$class <- NA
MRstat[MRstat$veg_age == "<1985", "class"] <- 1
MRstat[MRstat$veg_age %in% as.character(1986:1990), "class"] <- 2
MRstat[MRstat$veg_age %in% as.character(1991:1995), "class"] <- 3
MRstat[MRstat$veg_age %in% as.character(1996:2000), "class"] <- 4
MRstat[MRstat$veg_age %in% as.character(2001:2005), "class"] <- 5
MRstat[MRstat$veg_age %in% as.character(2006:2010), "class"] <- 6
MRstat[MRstat$veg_age %in% as.character(2011:2015), "class"] <- 7
MRstat[MRstat$veg_age %in% as.character(2016:2020), "class"] <- 8
MRstat[MRstat$veg_age == ">2021", "class"] <- 9

tmp <- as.data.frame(MRstat) %>% group_by(SiteID) %>% summarise(SiteClass = mode(class))
MRstat <- left_join(MRstat, tmp)

st_write(MRstat, "C:/Users/jbrow/My Drive/WCSC//CumulativeEffects_2022/2_pipeline/MineReclaim/out/MRstations_ageEst.shp", delete_layer=T)

