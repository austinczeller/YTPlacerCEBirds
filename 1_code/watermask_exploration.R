### calculate thresholds for NDVI and watermasking
## use digitized disturbance layer to explore thresholds
library(tidyverse)
library(terra)
library(sf)
library(irr) ## calculate kappa

##load in disturbance 
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

dist2021 <- st_intersection(dist %>% st_transform(crs = st_crs(sites)), 
                               sites %>% filter(year == 2021) %>% 
                              st_buffer(1500) %>% st_union())

dist2018 <- st_intersection(dist %>% st_transform(crs = st_crs(sites)), 
                            sites %>% filter(year == 2018) %>% 
                              st_buffer(1500) %>% st_union())
dist2017 <- st_intersection(dist %>% st_transform(crs = st_crs(sites)), 
                            sites %>% filter(year == 2017) %>% 
                              st_buffer(1500) %>% st_union())

## 2021 ----
#### load watermask B3B5 ----
ndwi1 <- rast("0_data/GEE/ndwi_B3B5_2021.tif")

dist2021 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1))

ndwi1 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
    })

ndwi_sum <- pmap(list(ndwi1, names(ndwi1)), 
                  \(v, n) data.frame(type = n, ndwi1 = v))
ndwi_sum <- bind_rows(ndwi_sum)
names(ndwi_sum) <- c("type", "B3B5")


#### load watermask B3B11 ----
ndwi1 <- rast("0_data/GEE/ndwi_B3B11_2021.tif")

ndwi1 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndwiB3B11 <- pmap(list(ndwi1, names(ndwi1)), 
                 \(v, n) data.frame(type = n, ndwi1 = v))
ndwiB3B11 <- bind_rows(ndwiB3B11)
names(ndwiB3B11) <- c("type", "B3B11")
ggplot(ndwiB3B11, aes(type, B3B11)) + geom_boxplot()


#### load watermask B4B11 ----
ndwi1 <- rast("0_data/GEE/ndwi_B4B11_2021.tif")

ndwi1 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndwiB4B11 <- pmap(list(ndwi1, names(ndwi1)), 
                  \(v, n) data.frame(type = n, ndwi1 = v))
ndwiB4B11 <- bind_rows(ndwiB4B11)
names(ndwiB4B11) <- c("type", "B4B11")
ggplot(ndwiB4B11, aes(type, B4B11)) + geom_boxplot()

##load watermask SWM ----

ndwi1 <- rast("0_data/GEE/swm_2021.tif")

ndwi1 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

swm <- pmap(list(ndwi1, names(ndwi1)), 
                  \(v, n) data.frame(type = n, ndwi1 = v))
swm <- bind_rows(swm)
names(swm) <- c("type", "swm")
ggplot(swm, aes(type, swm)) + geom_boxplot()


all.equal(ndwi_sum$type, swm$type)
ndwi_sum$swm <- swm$swm
ndwi_sum$B4B11 <- ndwiB4B11$B4B11
ndwi_sum$B3B11 <- ndwiB3B11$B3B11
ndwi_sum$water_dig <- ifelse(ndwi_sum$type == "d_pWater", T, F)

## B3B5 ----
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.2) + geom_vline(xintercept = 0.4)
threshold <- seq(-0.2, .4, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B3B5>x)))$value
})
b3b5 <- data.frame(threshold, kappa)
ggplot(b3b5, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = b3b5[which(b3b5$kappa == max(b3b5$kappa)),]$threshold)
max(b3b5$kappa) ##0.42 at 0.08

## B3B11 ----
ggplot(ndwi_sum, aes(B3B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -.4)+ geom_vline(xintercept = 0.8)
threshold <- seq(-0.4, .8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B3B11>x)))$value
})
B3B11 <- data.frame(threshold, kappa)
ggplot(B3B11, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B3B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = B3B11[which(B3B11$kappa == max(B3B11$kappa)),]$threshold)
max(B3B11$kappa) ##0.39 at 0.15

## B4B11 ----
ggplot(ndwi_sum, aes(B4B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -.3)+ geom_vline(xintercept = 0.8)
threshold <- seq(-0.3, .8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B4B11>x)))$value
})
B4B11 <- data.frame(threshold, kappa)
ggplot(B4B11, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B4B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = B4B11[which(B4B11$kappa == max(B4B11$kappa)),]$threshold)
max(B4B11$kappa) ##0.31 at 0.13

### SWM ----
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = 0.4)+ geom_vline(xintercept = 2)
threshold <- seq(0.6, 2, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$swm>x)))$value
})
swm <- data.frame(threshold, kappa)
ggplot(swm, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = swm[which(swm$kappa == max(swm$kappa)),]$threshold)
max(swm$kappa) ##0.433 at 1.31



## load NDVI ----

ndwi1 <- rast("0_data/GEE/ndvi_2021.tif")

ndwi1 <- dist2021 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndvi <- pmap(list(ndwi1, names(ndwi1)), 
            \(v, n) data.frame(type = n, ndwi1 = v))
ndvi <- bind_rows(ndvi)
names(ndvi) <- c("type", "ndvi")
ggplot(ndvi, aes(type, ndvi)) + geom_boxplot()

ndwi_sum$ndvi <- ndvi$ndvi

ggplot(ndwi_sum, aes(ndvi, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = 0.4)+ geom_vline(xintercept = -.4)
threshold <- seq(-.4, .4, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$ndvi>x)))$value
})
ndvi <- data.frame(threshold, kappa)
ggplot(ndvi, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(ndvi, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = ndvi[which(ndvi$kappa == max(ndvi$kappa)),]$threshold)
max(ndvi$kappa) 

##NDVI can't seperate water and bare ground.  Need watermask first. 

## 2018 ----

#### load watermask B3B5 ----
ndwi1 <- rast("0_data/GEE/ndwi_B3B5_2018.tif")

dist2018 <- dist2018 %>% 
  st_transform(crs = crs(ndwi1))

ndwi1 <- dist2018 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndwi_sum18 <- pmap(list(ndwi1, names(ndwi1)), 
                 \(v, n) data.frame(type = n, ndwi1 = v))
ndwi_sum18 <- bind_rows(ndwi_sum18)
names(ndwi_sum18) <- c("type", "B3B5")


#### load watermask B3B11 ----
ndwi1 <- rast("0_data/GEE/ndwi_B3B11_2018.tif")

ndwi1 <- dist2018 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndwiB3B11 <- pmap(list(ndwi1, names(ndwi1)), 
                  \(v, n) data.frame(type = n, ndwi1 = v))
ndwiB3B11 <- bind_rows(ndwiB3B11)
names(ndwiB3B11) <- c("type", "B3B11")
ggplot(ndwiB3B11, aes(type, B3B11)) + geom_boxplot()

### load watermask SWM ----

ndwi1 <- rast("0_data/GEE/swm_2018.tif")

ndwi1 <- dist2018 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

swm <- pmap(list(ndwi1, names(ndwi1)), 
            \(v, n) data.frame(type = n, ndwi1 = v))
swm <- bind_rows(swm)
names(swm) <- c("type", "swm")
ggplot(swm, aes(type, swm)) + geom_boxplot()

### NDVI ----
ndwi1 <- rast("0_data/GEE/ndvi_2018.tif")

ndwi1 <- dist2018 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndvi <- pmap(list(ndwi1, names(ndwi1)), 
             \(v, n) data.frame(type = n, ndwi1 = v))
ndvi <- bind_rows(ndvi)
names(ndvi) <- c("type", "ndvi")
ggplot(ndvi, aes(type, ndvi)) + geom_boxplot()


all.equal(ndwi_sum18$type, swm$type)
ndwi_sum18$swm <- swm$swm
ndwi_sum18$B3B11 <- ndwiB3B11$B3B11
ndwi_sum18$water_dig <- ifelse(ndwi_sum18$type == "d_pWater", T, F)
ndwi_sum$year <- 2021
ndwi_sum18$year <- 2018
ndwi_sum18$B4B11 <- NA
ndwi_sum18$ndvi <- ndvi$ndvi
ndwi_sum18 <- select(ndwi_sum18, names(ndwi_sum))

ndwi_sum <- rbind(ndwi_sum, ndwi_sum18)

### B3B5 ----
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.1) + geom_vline(xintercept = 0.4)
threshold <- seq(-0.1, .4, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B3B5>x)))$value
})
b3b5 <- data.frame(threshold, kappa)
ggplot(b3b5, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = b3b5[which(b3b5$kappa == max(b3b5$kappa)),]$threshold)
max(b3b5$kappa) ##0.417 at 0.076

### B3B11 ----
ggplot(ndwi_sum, aes(B3B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -.4)+ geom_vline(xintercept = 0.8)
threshold <- seq(-0.4, .8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B3B11>x)))$value
})
B3B11 <- data.frame(threshold, kappa)
ggplot(B3B11, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B3B11, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = B3B11[which(B3B11$kappa == max(B3B11$kappa)),]$threshold)
max(B3B11$kappa) ##0.407 at 0.14


### SWM ----
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = 0.5)+ geom_vline(xintercept = 2)
threshold <- seq(0.5, 2, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$swm>x)))$value
})
swm <- data.frame(threshold, kappa)
ggplot(swm, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = swm[which(swm$kappa == max(swm$kappa)),]$threshold)
max(swm$kappa) ##0.448 at 1.27


## 2017 ----

#### load watermask B3B5 ----
ndwi1 <- rast("0_data/GEE/ndwi_B3B5_2017.tif")

dist2017 <- dist2017 %>% 
  st_transform(crs = crs(ndwi1))

ndwi1 <- dist2017 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndwi_sum17 <- pmap(list(ndwi1, names(ndwi1)), 
                   \(v, n) data.frame(type = n, ndwi1 = v))
ndwi_sum17 <- bind_rows(ndwi_sum17)
names(ndwi_sum17) <- c("type", "B3B5")

### load watermask SWM ----

ndwi1 <- rast("0_data/GEE/swm_2017.tif")

ndwi1 <- dist2017 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

swm <- pmap(list(ndwi1, names(ndwi1)), 
            \(v, n) data.frame(type = n, ndwi1 = v))
swm <- bind_rows(swm)
names(swm) <- c("type", "swm")
ggplot(swm, aes(type, swm)) + geom_boxplot()

ndwi1 <- rast("0_data/GEE/ndvi_2017.tif")

ndwi1 <- dist2017 %>% 
  st_transform(crs = crs(ndwi1)) %>%
  split(f = .$simple) %>% ## split by veg, water, bare
  map(function(dtype) {
    r = crop(x = ndwi1, y = dtype, mask=T, touches = F)
    r = values(r)
    r = na.omit(r)
    r
  })

ndvi <- pmap(list(ndwi1, names(ndwi1)), 
             \(v, n) data.frame(type = n, ndwi1 = v))
ndvi <- bind_rows(ndvi)
names(ndvi) <- c("type", "ndvi")
ggplot(ndvi, aes(type, ndvi)) + geom_boxplot()

all.equal(ndwi_sum17$type, swm$type)
ndwi_sum17$swm <- swm$swm
ndwi_sum17$B3B11 <- NA
ndwi_sum17$water_dig <- ifelse(ndwi_sum17$type == "d_pWater", T, F)
ndwi_sum$year <- 2021
ndwi_sum17$year <- 2017
ndwi_sum17$B4B11 <- NA
ndwi_sum17$ndvi <- ndvi$ndvi
ndwi_sum17 <- select(ndwi_sum17, names(ndwi_sum))

ndwi_sum <- rbind(ndwi_sum, ndwi_sum17)

### B3B5 ----
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.1) + geom_vline(xintercept = 0.4)
threshold <- seq(-0.1, .4, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$B3B5>x)))$value
})
b3b5 <- data.frame(threshold, kappa)
ggplot(b3b5, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(B3B5, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = b3b5[which(b3b5$kappa == max(b3b5$kappa)),]$threshold)
max(b3b5$kappa) ##0.412 at 0.076

### SWM ----
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = 0.5)+ geom_vline(xintercept = 2)
threshold <- seq(0.5, 2, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(ndwi_sum$water_dig, ndwi_sum$swm>x)))$value
})
swm <- data.frame(threshold, kappa)
ggplot(swm, aes(threshold, kappa)) + geom_point()
ggplot(ndwi_sum, aes(swm, fill = water_dig)) + geom_density(alpha = .5) +
  geom_vline(xintercept = swm[which(swm$kappa == max(swm$kappa)),]$threshold)
max(swm$kappa) ##0.446 at 1.25


# Bare vs Veg ----
bv <- filter(ndwi_sum, !water_dig)
ggplot(bv, aes(ndvi, fill = type)) + geom_density(alpha = .5) +
  geom_vline(xintercept = -0.2)+ geom_vline(xintercept = .8)
threshold <- seq(-0.2, 0.8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(bv$type== "d_pVeg", bv$ndvi>x)))$value
})
ndvi <- data.frame(threshold, kappa)
ggplot(ndvi, aes(threshold, kappa)) + geom_point()
ggplot(bv, aes(ndvi, fill = type)) + geom_density(alpha = .5) +
  geom_vline(xintercept = ndvi[which(ndvi$kappa == max(ndvi$kappa)),]$threshold)
max(ndvi$kappa) ##0.643 at 0.354

bv <- rbind(filter(bv, type == "n_dBare"), sample_n(filter(bv, type != "n_dBare"), size=97326))
            
ggplot(bv, aes(ndvi, fill = type)) + geom_density(alpha = .5) +
geom_vline(xintercept = -0.2)+ geom_vline(xintercept = .8)
threshold <- seq(-0.2, 0.8, length.out = 75)
kappa <- map_dbl(threshold, function(x) {
  (kappa2(cbind(bv$type== "d_pVeg", bv$ndvi>x)))$value
  })
ndvi <- data.frame(threshold, kappa)
ggplot(ndvi, aes(threshold, kappa)) + geom_point()
ggplot(bv, aes(ndvi, fill = type)) + geom_density(alpha = .5) +
  geom_vline(xintercept = ndvi[which(ndvi$kappa == max(ndvi$kappa)),]$threshold)
max(ndvi$kappa) ##0.645 at 0.42


###
# water = swm > 1.25
# bare = ndvi < 0.645
# veg = ndvi > 0.645
            




