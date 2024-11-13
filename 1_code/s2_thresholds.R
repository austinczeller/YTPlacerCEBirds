library(tidyverse)
library(terra)
library(sf)
library(irr) ## calculate kappa
library(rgdal) ##to read gdbs

##load in disturbance ----
sites <- readRDS("2_pipeline/1_CountDataProcessing/out/sites.RDS")
fgdb <- "CumulativeEffects_GIS/CumEff_corrected_mjb.gdb"
# List all feature classes in a file geodatabase
# subset(ogrDrivers(), grepl("GDB", name))
# fc_list <- ogrListLayers(fgdb)
# print(fc_list)
dist <- st_as_sf(readOGR(dsn=fgdb, layer="Polygon_disturbance_VegBareWater")) 
## Disturbance separated into bare land, vegetated land, and water based on aerial imagery 
#(spot, sentinel2, Esri basemap, in appropriate year).  
#Only sites selected for analysis, 
#where disturbance was previously identified by YG or patrice, 
#were inspected and digitised. 

dist <- st_intersection(dist %>% st_transform(crs = st_crs(sites)), 
                sites  %>% st_buffer(1500) %>%
                  group_by(year) %>% summarise(geometry = st_union(geometry)))

dist.l <- split(dist, f = dist$year)
v_all <- names(dist.l)
dist.l <- dist.l[1:3]
v <- names(dist.l)

### Watermask B3B8 ----
files <- list.files("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", pattern = "B3B8")
B3B8 <- pmap(list(v, dist.l), function(yr, dist){
  file <- files[str_which(files, as.character(yr))]
  r <- rast(file.path("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", file))
  d <- dist %>% st_transform(crs = crs(r))
  rd <- d %>% 
    split(f = .$simple) %>% ## split by veg, water, bare
    map(function(dtype) {
      rcrop = crop(x = r, y = dtype, mask=T, touches = F)
      rcrop = values(rcrop)
      rcrop = na.omit(rcrop)
      # data.frame(type = unique(dtype$simple), B3B8 = rcrop)
      data.frame(B3B8 = rcrop)
    }) %>% list_rbind(names_to = "type")
  return(rd)
  }) %>% list_rbind(names_to = "year")

B3B8$year <- ifelse(B3B8$year == 1, v[1], ifelse(B3B8$year == 2, v[2], v[3]))
B3B8 <- rename(B3B8, b3b8 = nd)
B3B8$water_dig <- ifelse(B3B8$type == "d_pWater", T, F)

ggplot(B3B8, aes(type, b3b8)) + geom_boxplot()
ggplot(B3B8, aes(b3b8, fill = type)) + geom_histogram()

# B3B8 <- filter(B3B8, type != "d_pVeg")

ggplot(B3B8, aes(b3b8, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.4)+ geom_vline(xintercept = 0.6)
threshold <- seq(-.4, .6, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(B3B8$water_dig, B3B8$b3b8>x), "equal"))$value
})
tmp <- data.frame(threshold, kappa)
ggplot(tmp, aes(threshold, kappa)) + geom_point()
ggplot(B3B8, aes(b3b8, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = tmp[which(tmp$kappa == max(tmp$kappa)),]$threshold)
thresh <- tmp[which(tmp$kappa == max(tmp$kappa)),]$threshold
thresh
max(tmp$kappa) ##kappa =  0.47 at threshold of 0.046
kappa2(cbind(B3B8$water_dig, B3B8$b3b8>thresh))

## NDVI ----
files <- list.files("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", pattern = "ndvi")
ndvi <- pmap(list(v, dist.l), function(yr, dist){
  file <- files[str_which(files, yr)]
  r <- rast(file.path("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", file))
  d <- dist %>% st_transform(crs = crs(r)) 
  rd <- d %>% 
    split(f = .$simple) %>% ## split by veg, water, bare
    map(function(dtype) {
      rcrop = crop(x = r, y = dtype, mask=T, touches = F)
      rcrop = values(rcrop)
      rcrop = na.omit(rcrop)
      # data.frame(type = unique(dtype$simple), ndvi = rcrop)
      data.frame(ndvi_val = rcrop)
    }) %>% list_rbind(names_to = "type")
}) %>% list_rbind(names_to = "year")

ndvi$year <- ifelse(ndvi$year == 1, v[1], ifelse(ndvi$year == 2, v[2], v[3]))
ndvi <- rename(ndvi, ndvi = nd)

bv <- filter(ndvi, type != "d_pWater")

ggplot(bv, aes(ndvi, fill = type)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.2)+ geom_vline(xintercept = .8)
threshold <- seq(-0.2, 0.8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(bv$type== "d_pVeg", bv$ndvi>x)))$value
})
tmp <- data.frame(threshold, kappa)
ggplot(tmp, aes(threshold, kappa)) + geom_point()
ggplot(bv, aes(ndvi, fill = type)) + geom_histogram(alpha = .5) +
  geom_vline(xintercept = tmp[which(tmp$kappa == max(tmp$kappa)),]$threshold)
max(tmp$kappa) ##0.643 at 0.354

set.seed(123)
bv <- rbind(filter(bv, type == "n_dBare"), sample_n(filter(bv, type != "n_dBare"), size=97326))

ggplot(bv, aes(ndvi, fill = type)) + geom_histogram(alpha = .5) +
  geom_vline(xintercept = -0.2)+ geom_vline(xintercept = .8)
threshold <- seq(-0.2, 0.8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(bv$type== "d_pVeg", bv$ndvi>x)))$value
})
tmp2 <- data.frame(threshold, kappa)
ggplot(tmp2, aes(threshold, kappa)) + geom_point()
ggplot(bv, aes(ndvi, fill = type)) + geom_histogram(alpha = .5) +
  geom_vline(xintercept = tmp2[which(tmp2$kappa == max(tmp2$kappa)),]$threshold)
ggplot(bv, aes(ndvi, fill = type)) + geom_histogram(alpha = .5) +
  geom_vline(xintercept = tmp2[which(tmp2$kappa == max(tmp2$kappa)),]$threshold)
max(tmp2$kappa) ##kappa = 0.643 at ndvi threshold = 0.42
tmp2[which(tmp2$kappa == max(tmp2$kappa)),]$threshold
### Thresholds ----
##overall Kappa

threshold_sum <- cbind(ndvi, B3B8 = B3B8[,"b3b8"])
threshold_sum$s2 <- ifelse(threshold_sum$B3B8 > 0.046, "d_pWater",
                           ifelse(threshold_sum$ndvi > 0.42, "d_pVeg", "n_dBare"))
kappa2(cbind(threshold_sum$s2, threshold_sum$type))
# water = B3B8  > 0.046
# bare = ndvi < 0.42
# veg = ndvi > 0.42

threshold_sum$s2 <- paste0("predicted_", threshold_sum$s2)
threshold_sum$type <- paste0("actual_", threshold_sum$type)
table(threshold_sum$s2, threshold_sum$type) ##s2 is y, digitization is vertical
colSums(table(threshold_sum$s2, threshold_sum$type) )
# table(threshold_sum$type

169148/198801##s2 predicts veg correctly 85% of time
5270/13964 ##s2 predicts water correctly 38% of time
75288/97326 ##s2 predicts veg correctly 77% of time

###Mask layers ----
files <- list.files("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", pattern = "B3B8")
watermask <- map(v_all, function(yr){
  file <- files[str_which(files, yr)]
  r <- rast(file.path("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", file))
  r <- r > 0.046
  p <- as.polygons(r)
  names(p) = "water"
  st_as_sf(p)
    }) 
write_rds(watermask, "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/2_pipeline/2_HabitatSum/store/watermask.rds")

files <- list.files("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", pattern = "ndvi")
vegmask <- map(v_all, function(yr){
  file <- files[str_which(files, yr)]
  r <- rast(file.path("C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/0_data/GEE/", file))
  r <- r > 0.42
  p <- as.polygons(r)
  names(p) = "vegetated"
  st_as_sf(p)
}) 
write_rds(vegmask, "C:/Users/jbrow/My Drive/WCSC/CumulativeEffects_2022/2_pipeline/2_HabitatSum/store/vegmask.rds")

