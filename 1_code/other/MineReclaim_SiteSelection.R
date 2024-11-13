library(tidyverse)
library(sf)
library(rgdal) ##r geospatial data abstraction library
library(terra)
library(tidyterra)
library(stars)
##############################################################################################################
# Load and Set-up Design Inputs -------------------------------------------
##############################################################################################################

NAME <- 'MineReclaim' ## Name of the R file goes here (without the file extension!)
PROJECT <- 'CumulativeEffects_2022' ## Project folder
# PROJECT_DIR <- 'C:/Users/jbrown1/surfdrive/WCSC' ## Change this to the directory in which your project folder is located, make sure to avoid using single backslashes (i.e. 
PROJECT_DIR <- "C:/Users/jbrow/My Drive/WCSC/"

##Set working directory
setwd(file.path(PROJECT_DIR, PROJECT))
knitr::opts_knit$set(root.dir = file.path(PROJECT_DIR, PROJECT))
# Set  up pipeline folder if missing
### The code below will automatically create a pipeline folder for this code file if it does not exist.

if (dir.exists(file.path('empirical', '2_pipeline'))){
  pipeline <- file.path('empirical', '2_pipeline', NAME)
} else {
  pipeline <- file.path('2_pipeline', NAME)
}

if (!dir.exists(pipeline)) {
  dir.create(pipeline)
  for (folder in c('out', 'store', 'tmp')){
    dir.create(file.path(pipeline, folder))
  }
}


# load ecoregion boundaries
# ogrListLayers(dsn="Habitat/Ecoregions_2014_1M.gdb")
ecoregions <- st_as_sf(readOGR(dsn="CumulativeEffects_GIS/Habitat/Ecoregions_2014_1M.gdb", 
                               layer="Ecoregions_2014_1M"))
ecoregions <- filter(ecoregions, ECOREGION_ID %in% c(172, 302))
##172 = Klondike plateau
##302 = McQuesten Highlands

proj <- st_crs(ecoregions) ##use Yukon Albers projection

##load wetland layer extent

wet.ext <- read_sf("CumulativeEffects_GIS/Arc2R/WetlandLayer_extent.shp")

# Make grid ----
hex.diameter <- 10 #in km
hex <- st_make_grid(
  ecoregions,
  cellsize = sqrt(3)*(hex.diameter*1000/2), ##hexagon with diameter of 10k
  # n = c(round(diff(st_bbox(ecoregions)[c(1,3)])/10000), round(diff(st_bbox(ecoregions)[c(2,4)])/10000)),
  what = "polygons",
  square = F)
hex <- st_as_sf(hex)
hex <-hex %>% mutate(id = 1:length(hex[[1]]))


# Select cells that intersect linear disturbance and Placer claims
##select cells that have a wetland layer
hex <- hex %>%
  filter(st_intersects(x = ., y = st_union(wet.ext), sparse = FALSE)) 


PClaims <- read_sf("C:/Users/jbrow/My Drive/WCSC/MineReclaim/ReclaimChrono/0_data/GeoYukon/Placer_Claims_1M.shp/Placer_Claims_1M.shp")
##prepped in ArcGIS to speedup. Clip using ecorgions (Klondike Plateau and McQuestin)
LDisturb <- read_sf("CumulativeEffects_GIS/Arc2R/YG_LinearDisturbance_May2022_KlondikeMcQuestinClip.shp")

# LDisturb <- read_sf("CumulativeEffects_GIS/YG_SurfaceDisturbance_May2022/LinearFeatures_May2022_NoOverlap.shp")
LDisturb <- st_transform(LDisturb, crs = st_crs(hex))
# LDisturb <- LDisturb %>% filter(st_intersects(x = ., y = st_union(ecoregions, sparse = FALSE))

##select gridcells with placer claims
Phex <- hex %>%
  filter(st_intersects(x = ., y = st_union(PClaims), sparse = FALSE)) 

##select gridclass with linear disturbance (access points?)
Phex <- Phex %>%
  filter(st_intersects(x = ., y = st_union(LDisturb), sparse = FALSE)) 

saveRDS(Phex, file.path(pipeline, "store", "Placer_hexagons.shp")) 
Phex <- readRDS(file.path(pipeline, "store", "Placer_hexagons.shp"))

#annotate Phex with ecoregion
# calculate the geometry for hexagon centroids
cent_hexes <- st_centroid(Phex)
# assign hexagons to their ecoregion by joining the ecoregion to the hexagon centroids
# there may be some NA's where hexagons do not overlap the ecoregion layer
cent_hexes_join <- st_join(cent_hexes, ecoregions[,1:2]) 
# use this to add ecoregion name and number to each hexagon in the jurisdiction
Phex <-
  st_join(Phex, cent_hexes_join[,2:3])

Phex <-
  Phex %>% drop_na() # drop the NA's, hexagons with no overlapping ecoregion
saveRDS(Phex, file.path(pipeline, "store", "final_Placerhex.RDS"))
Phex<- readRDS(file.path(pipeline, "store", "final_Placerhex.RDS"))
#importing as a spatial polygons data frame will allow for the extent to be extracted in the following steps

#Calculate Placer Disturbance per grid---

SDisturb <- read_sf("CumulativeEffects_GIS/YG_SurfaceDisturbance_May2022/ArealFeatures_May2022_Un_Dis_J.shp")
SDisturb <- st_transform(SDisturb, crs = st_crs(Phex))
SDisturb <- st_intersection(SDisturb, Phex) ## add hex ID

Phex <- left_join(Phex, data.frame(id = SDisturb$id, disturbA = SDisturb %>% st_area() %>% as.numeric) %>%
                    group_by(id) %>% summarise(disturbA = sum(disturbA)))
Phex[is.na(Phex$disturbA),]$disturbA <- 0

Phex <- Phex %>% mutate(disturbA = disturbA/64951905*100) ## %cover
ggplot(Phex) + geom_sf(aes(fill = disturbA))
hist(Phex$disturbA, breaks = 6)

##calculate inclusion probability based on total disturbance ----

# first we create a list of ecoregions to be used in the loops
eco.list <- unique(Phex$ECOREGION_ID)

# this function will calculate the habitat-based inclusion probabilities for a given ecoregion
eco.fn <- function(eco) {
  eco.hex <- filter(Phex, ECOREGION_ID == eco) %>% mutate(dbin = round(disturbA))
  d <- density(eco.hex$disturbA)      
  pd <- sapply(eco.hex$disturbA, function(x) d$y[which.min((abs(d$x - x)))])
  pinc <- 1/pd/sum(1/pd) ## probability of inclusion based on total disturbance
  eco.hex$p.inc.disturbA <- pinc
  eco.hex
}

 

# lapply passes each item from our list of ecoregions to the eco.fn
Phex.l <- lapply(eco.list, eco.fn)
Phex <- do.call(rbind, Phex.l) 
Phex %>% ggplot() + geom_sf(aes(fill = p.inc.disturbA))

## reclaim year representation per grid ----
## bring in the Canada Forest Change raster
##pre-clipped in arcGIS for speed. Forest harvest - all change within placer claims. 
## forest fire, all change within placer disturbance (YG). This layer seems to contain a mix of placer clearing not classed as harvest (want to keep), and actual fires (want to remove)

##############################################################################################################
# Calculate Habitat-based Inclusion Probabilities -------------------------
##############################################################################################################

# set up a loop for calculating habitat-based inclusion probabilities
# this is calculated for each hexagon, stratified by ecoregion

# ##canada forest harvest
# CFH <- rast("CumulativeEffects_GIS/Arc2R/CanadaForestChange_Harvest_YA.tif") # read in habitat layer
# CFH <- crop(CFH, Phex) # clip to jurisdiction
# CFH <- mask(CFH, vect(Phex))
# CFH <- mask(CFH, vect(PClaims))
# values(CFH) <- ifelse(values(CFH) == 0, NA, values(CFH)) #0 Value = NA
# 
# ##Canada forest fire
# CFF <- rast("CumulativeEffects_GIS/Arc2R/CanadaForestChange_Fire_YA.tif") # read in habitat layer
# CFF <- crop(CFF, Phex) # clip to jurisdiction
# # ggplot() + geom_spatraster(data = CFF) + geom_sf(data = Phex)
# CFF <- mask(CFF, vect(Phex))
# CFF <- mask(CFF, vect(SDisturb))
# CFF <- mask(CFF, vect(PClaims))
# values(CFF) <- ifelse(values(CFF) == 0, NA, values(CFF)) #0 Value = NA
# 
# #combine. If CFH is not NA, CFH, else, CFF. 
# CFC <- CFH
# values(CFC) <- ifelse(is.na(values(CFC)), values(CFF), values(CFC))
# # ggplot() + geom_spatraster(data = CFC) + geom_sf(data = Phex, fill = NA, col = "red")
# writeRaster(CFC, file.path(pipeline, "store", "CanadaForestChange_combined_placer.tif"), overwrite=T)
# CFC <- rast(file.path(pipeline, "store", "CanadaForestChange_combined_placer.tif"))
# ##NA values should be based neutral probability per raster. 


# first we create a list of ecoregions to be used in the loops
eco.list <- unique(Phex$ECOREGION_ID)

# this function will calculate the habitat-based inclusion probabilities for a given ecoregion
eco.fn <- function(eco) {
  eco.rast <- crop(CFC, vect(subset(Phex, ECOREGION_ID == eco)))
  
  eco.rast <- mask(eco.rast, vect(subset(Phex, ECOREGION_ID == eco))) # clip to the ecoregion
  
  # calculate the habitat-based inclusion probabilities for given pixels
  eco.freq <- as.data.frame(freq(eco.rast)) # create a table with counts of the # of pixels for each habitat class within a given ecoregion
  eco.freq <- subset(eco.freq, !value %in% NA) # remove NA's from the calculations
  eco.freq$habprob <- 1 / (nrow(eco.freq)) / eco.freq$count # calculate probabilities for each class
  
  # eco.freq$count*eco.freq$habprob #each year has combined equal probability of being selected. 
  
  
  # create matrix to reclassify rasters so that pixel value = habitat probability
  eco.old <- eco.freq[["value"]]
  eco.new <- eco.freq[["habprob"]]
  eco.mat <- matrix(c(rbind(eco.old, eco.new)), ncol = 2, byrow = TRUE)
  eco.hab.pix <- classify(eco.rast, eco.mat) # new raster with pixel habitat probability
  Phex.eco <- subset(Phex, ECOREGION_ID == eco)
  
  ##inclusion probability of known reclamation age
  eco.hab.hex <- terra::extract(eco.hab.pix, vect(Phex.eco),
                                 fun = sum, na.rm = TRUE) # new hexagon file with habitat-based inclusion probabilities
  Phex.eco$reclaim.p <- eco.hab.hex[,2]
  
  ##area of known reclamation age
  eco.mat[,2] <- 1
  eco.hab.pix <- classify(eco.rast, eco.mat) # new raster with pixel habitat probabili
  eco.hab.hex <- terra::extract(eco.hab.pix, vect(Phex.eco),
                                fun = sum, na.rm = TRUE) # new hexagon file with habitat-based inclusion probabilities
  Phex.eco$reclaim.A <- eco.hab.hex[,2]
  
  ##total area
  values(eco.hab.pix) <- 1
  eco.hab.hex <- terra::extract(eco.hab.pix, vect(Phex.eco),
                                fun = sum, na.rm = TRUE) # new hexagon file with habitat-based inclusion probabilities
  Phex.eco$total.A <- eco.hab.hex[,2]
  Phex.eco$pc.knownReclaim <- Phex.eco$reclaim.A/Phex.eco$total.A
  
  # sum(Phex.eco$reclaim.p) ##adds to a probability of 1 accross all hexigons
  p.equal <- 1/dim(Phex.eco)[1] ## if all hab the same, each grid would sum to this probability
  
  Phex.eco$p.inc.reclaim <- (Phex.eco$pc.knownReclaim*Phex.eco$reclaim.p) + ##probability is based on reclaim age probability for area where CFC is mapped, 
    (1-Phex.eco$pc.knownReclaim)*(p.equal)##and p.equal for area where not mapped
  # sum(Phex.eco$p.inc.reclaim) ## should be 1
  
  # save
  st_write(obj=Phex.eco, dsn=file.path(pipeline, "store", paste0(eco,"reclaim_hexes.shp")), append=FALSE)
  Phex.eco
}


# lapply passes each item from our list of ecoregions to the eco.fn
Phex.l <- lapply(eco.list, eco.fn)
# hex172<- st_read(file.path(pipeline, "store", "172reclaim_hexes.shp")) 
# hex302<- st_read(file.path(pipeline, "store", "302reclaim_hexes.shp"))
# Phex.l <- list(hex172, hex302)
Phex <- do.call(rbind, Phex.l) 
Phex %>% ggplot() + geom_sf(aes(fill = p.inc.reclaim))

Phex <- Phex %>% group_by(ECOREGION_ID) %>% mutate(p.inc = p.inc.disturbA*p.inc.reclaim/sum(p.inc.disturbA*p.inc.reclaim))
# sum(Phex$p.inc)
         
##aim: 16 grid cells in klondike, 10 in Mayo. (= 52 reclaim sites if 2 per grid)
#overdraw 5 per region in case of accessibility issues. 21 and 15
set.seed(100)
sample.hex <- c(sample(filter(Phex, ECOREGION_ID == 172)%>%pull(id), 21, replace = F, prob = filter(Phex, ECOREGION_ID == 172)%>%pull(p.inc)),
sample(filter(Phex, ECOREGION_ID == 302)%>%pull(id), 25, replace = F, prob = filter(Phex, ECOREGION_ID == 302)%>%pull(p.inc)))

ggplot(filter(Phex, id %in% sample.hex)) + geom_sf(aes(fill = disturbA))
ggplot(filter(Phex, id %in% sample.hex & ECOREGION_ID == 302), aes(disturbA)) +geom_density(col = "red", size = 2) +
  geom_density(data = filter(Phex, ECOREGION_ID == 302))
ggplot(filter(Phex, id %in% sample.hex & ECOREGION_ID == 172), aes(disturbA)) +geom_density(col = "red", size = 2) +
  geom_density(data = filter(Phex, ECOREGION_ID == 172))
##resulting disturbance gradient is more uniform than original selection

ggplot(filter(Phex, id %in% sample.hex & ECOREGION_ID == 302), aes(p.inc.reclaim)) +geom_density(col = "red", size = 2) +
  geom_density(data = filter(Phex, ECOREGION_ID == 302)) + ylim(0,20000)
ggplot(filter(Phex, id %in% sample.hex & ECOREGION_ID == 172), aes(p.inc.reclaim)) +geom_density(col = "red", size = 2) +
  geom_density(data = filter(Phex, ECOREGION_ID == 172)) + ylim(0,20000)

sample.hex <- data.frame(id = sample.hex, order = 1:length(sample.hex))
Shex <- right_join(Phex, sample.hex) %>% select(id, ecoID = ECOREGION_ID, eco=ECOREGION, disturbA, 
                                            reclaim.p, pc.knownReclaim, p.inc, order) %>%
  arrange(order) %>% group_by(eco) %>% mutate(eco.order = row_number()) %>% select(-order)

sample.hex <- sample.hex %>% left_join(Phex %>% select(id, ecoID = ECOREGION_ID, eco=ECOREGION, disturbA, reclaim.p, pc.knownReclaim, p.inc))               
ggplot(ecoregions) + geom_sf(aes(fill = ECOREGION_ID)) + geom_sf(data = wet.ext, alpha = .5, fill = "white") + 
  geom_sf(data = Shex, fill = "red", col = "black") + geom_sf(data = Phex, fill = NA, col = "grey50")

write_sf(Shex, file.path(pipeline, "store", "sampleHex.shp"))
Shex <- st_read(file.path(pipeline, "store", "sampleHex.shp"))

sample.notes <- read.csv("MineReclaim_PolyNotes.csv") %>% rename(id = Hex.ID)

Shex <- left_join(Shex, sample.notes %>% select(id, mapYR = MapYear, notes = Notes, mergeGrid = MergeGrid)) 

Shex %>% filter(mapYR != "x") %>% group_by(ecoID) %>% summarise(n = n()) ## need at least 3 more MAYO sites

Shex[Shex$id == 775,"id"] <- 775729
Shex[Shex$id == 729,"id"] <- 775729
Shex[Shex$id == 567,"id"] <- 567591
Shex[Shex$id == 591,"id"] <- 567591
Shex[Shex$id == 1441,"id"] <- 144113951418
Shex[Shex$id == 1395,"id"] <- 144113951418
Shex[Shex$id == 1418,"id"] <- 144113951418

#Station placement ----
#This was done in ArcGIS.  Site were select within (or near) the drawn hexigons. Hexes without access were discarded. 
#Mine sites around keno were opportunistically added because disturbance data was out of date. 
#Sites were outlined around areas that appeared to have different vegetation regrowth. 
#station points were then placed within the bounds of the site outlines, ensuring points were a minimum of 300m apart. 
#A standardized grid was not used, instead stations were placed centrally within the disturbance outlines while maintaining distance. 
#This was done by placing points (edit>create), creating a 150m buffer around each point (geoanalysis>buffer), then moving
#point + buffer together, ensuring buffers did not over lap with every point. 
#points were then intersected with the station to obtain a siteID and hexID. The 'near' geoanalysis tool was used to 
#calculate nearest distance to water_linear_flow, and the 'extract multi value to point' geoanalysis tool was used to 
#annotate elevation from the gmted raster.
#Annotated station points were exported to CumulativeEffects_GIS/ReclaimStations2023.shp

mlf <- read_stars("CumulativeEffects_GIS/macrolandforms.tif", proxy = F) ##macro landforms. Use this crs
crs <- st_crs(mlf)
MRstat <- read_sf("CumulativeEffects_GIS/ReclaimStations2023.shp")
MRstat <- st_transform(MRstat, crs = crs)

##some sites should be grouped into a different hex (either mines span hex boarder, or no available sites of different age within single hex.

##update clustering based on reclaimed site selection ----
MRstat$ClusterID <- NA
MRstat[MRstat$SiteID %in% c(59),]$ClusterID <- 1
MRstat[MRstat$SiteID %in% c(68),]$ClusterID <- 2
MRstat[MRstat$SiteID %in% c(58,60),]$ClusterID <- 3
MRstat[MRstat$SiteID %in% c(48, 49),]$ClusterID <- 4
MRstat[MRstat$SiteID %in% c(47,67,46),]$ClusterID <- 5
MRstat[MRstat$SiteID %in% c(43,44,45),]$ClusterID <- 6
MRstat[MRstat$SiteID %in% c(40,41,42),]$ClusterID <- 7
MRstat[MRstat$SiteID %in% c(37,38),]$ClusterID <- 8
MRstat[MRstat$SiteID %in% c(39),]$ClusterID <- 9
MRstat[MRstat$SiteID %in% c(35,34,36,66),]$ClusterID <- 10
MRstat[MRstat$SiteID %in% c(56,57),]$ClusterID <- 11
MRstat[MRstat$SiteID %in% c(32,31,62,33),]$ClusterID <- 12
MRstat[MRstat$SiteID %in% c(28,29,30),]$ClusterID <- 13
MRstat[MRstat$SiteID %in% c(26,27),]$ClusterID <- 14
MRstat[MRstat$SiteID %in% c(63,25,24,23),]$ClusterID <- 15
MRstat[MRstat$SiteID %in% c(65,22,21),]$ClusterID <- 16
MRstat[MRstat$SiteID %in% c(16,17),]$ClusterID <- 17
MRstat[MRstat$SiteID %in% c(12,13,15),]$ClusterID <- 18
MRstat[MRstat$SiteID %in% c(11,10,8),]$ClusterID <- 19
MRstat[MRstat$SiteID %in% c(7,6),]$ClusterID <- 20
MRstat[MRstat$SiteID %in% c(18,20,19),]$ClusterID <- 21
MRstat[MRstat$SiteID %in% c(64,4),]$ClusterID <- 22
MRstat[MRstat$SiteID %in% c(5,3),]$ClusterID <- 23
MRstat[MRstat$SiteID %in% c(50, 2, 1),]$ClusterID <- 24
MRstat[MRstat$SiteID %in% c(61,55,54),]$ClusterID <- 25
MRstat[MRstat$SiteID %in% c(52,53),]$ClusterID <- 26

# as.data.frame(MRstat) %>% select(ClusterID, HexID) %>% distinct() %>% 
#   group_by(ClusterID) %>% summarise(nHex = n()) %>%
#   print(n = 30) 
## 8 clusters overlap 2 hexagons
## 14 remain within original hexagons.

#Select within Cluster reference site ----
## Generate new cluster hexagons centered around points. 
# MRClus <- split(MRstat, f = MRstat$ClusterID)
# MRClus <- lapply(MRClus, function(clus){
#   cent <- st_centroid(clus)
#   buff <- st_buffer(cent, dist = 5000)
#   buff
# })
# 


MRClus <- MRstat %>% group_by(ClusterID) %>%
  summarise() ## changes to multipoint by ClusterID
MRClus <- st_centroid(MRClus) %>% st_buffer(dist = 5000)

## generate 300x300m grid of points
MRgrid <- st_make_grid(MRClus, cellsize = 300, what = "corners")
MRgrid <- st_intersection(MRgrid, MRClus)
MRgrid <- st_join(st_as_sf(MRgrid), st_as_sf(MRClus))

## remove points within 600m of disturbance
SDisturb <- read_sf("CumulativeEffects_GIS/YG_SurfaceDisturbance_May2022/ArealFeatures_May2022_Un_Dis_J.shp")
SDisturb <- st_transform(SDisturb, crs = crs)
DistBuff <- st_buffer(SDisturb, dist = 600)
DistBuff <- st_intersection(DistBuff, MRClus)
DistBuff <-st_union(DistBuff)

MRgrid <- st_difference(MRgrid, DistBuff)

## remove points further than 2 SD away from distance to water
fgdb <- "CumulativeEffects_GIS/Habitat/canvec_50k_YT_Hydro.gdb"
waterway <- st_as_sf(readOGR(dsn=fgdb, layer="water_linear_flow_1_ecoregion_union"))
# waterway <- st_transform(waterway, crs = crs)
waterwayClus <- st_intersection(waterway,st_buffer(MRClus, dist = 2500))

waterway.l <- split(waterwayClus, waterwayClus$ClusterID)
MRgrid.l <- split(MRgrid, MRgrid$ClusterID)

MRgrid.l <- mapply(function(grid, water){
  d2w <- st_distance(grid, water)
  grid$d2w <- d2w
  grid
}, MRgrid.l, waterway.l, SIMPLIFY = F)

MRstat.l <- split(MRstat, MRstat$ClusterID)

MRstat.l <- mapply(function(stat, water){

  dist = st_distance(stat, water)
  stat$d2w <- as.numeric(dist[,1])
  stat
}, MRstat.l, waterway.l, SIMPLIFY = F)

wsd <- sapply(MRstat.l, function(x) sd(x$d2w))

MRgrid.l <- mapply(function(grid, sd){
  filter(grid, as.numeric(d2w) <= 2*sd )
}, MRgrid.l, wsd, SIMPLIFY = F)

ggplot(waterway.l[["18"]])+ 
  geom_sf(data = MRClus %>% filter(ClusterID == 18)) + 
  geom_sf(col = "cyan") + 
  geom_sf(data = MRstat %>% filter(ClusterID == 18), col = "red") + 
  geom_sf(data = MRgrid.l[["18"]])

ggplot(waterway.l[["1"]])+ 
  geom_sf(data = MRClus %>% filter(ClusterID == 1)) + 
  geom_sf(col = "cyan") + 
  geom_sf(data = MRstat %>% filter(ClusterID == 1), col = "red") + 
  geom_sf(data = MRgrid.l[["1"]])

#### remove points more with elevation outside 2 SD of mine points?  
#Or based on land form (mostly removing steep slopes?)

##load elevation
gmt <- read_stars("CumulativeEffects_GIS/Elevation/GMTED/50n150w_20101117_gmted_med075_ecoregion.tif", proxy = F)
names(gmt) <- "elev"

MRstat.l <- lapply(MRstat.l, function(pts){
  elev = st_extract(gmt, pts)
  pts$elev <- elev$elev
  pts
})

esd <- sapply(MRstat.l, function(x) sd(x$elev))
emn <- sapply(MRstat.l, function(x) mean(x$elev))
# eup <- emed + 2*esd
# elw <- emed - 2*esd

MRgrid.l <- lapply(MRgrid.l, function(pts) {
  elev = st_extract(gmt, pts)
  pts$elev <- elev$elev
  pts
})

MRgrid.l <- mapply(function(grid, mn, sd){
  f1 = filter(grid, elev <= (mn + 2*sd) & elev >= (mn - 2*sd))
  f2 = filter(grid, elev <= (mn + 3*sd) & elev >= (mn - 3*sd))
  if(length(f1$ClusterID) > 5)  {
    pts = f1
  } else if(length(f2$ClusterID) > 5) {
    pts = f2
  } else {
    pts = grid
  }
  pts
}, MRgrid.l, emn, esd, SIMPLIFY = F)


## check river order to see if similar?
## check landforms to see if similar?
names(mlf) <- "landform"
MRstat.l <- lapply(MRstat.l, function(pts){
  landform = st_extract(mlf, pts)
  pts$landform <- landform$landform
  pts
})

MRgrid.l <- lapply(MRgrid.l, function(pts) {
  landform = st_extract(mlf, pts)
  pts$landform <- landform$landform
  pts
})

mlf_attr <- read_csv("CumulativeEffects_GIS/macrolandforms_attr.csv")

MRstat.l <- lapply(MRstat.l, function(pts){
  left_join(pts, mlf_attr %>% select(VALUE, lf.type = TYPE), by = c("landform" = "VALUE"))
})

MRgrid.l <- lapply(MRgrid.l, function(pts){
  left_join(pts, mlf_attr %>% select(VALUE, lf.type = TYPE), by = c("landform" = "VALUE"))
})

table(unlist(sapply(MRstat.l, function(c) c$lf.type))) ## most common landforms

##landform rank - if > 10 points = 2, between 1 and 10 = 1, otherwise 0

MRgrid.l <- lapply(MRgrid.l, function(pts) {
  pts$lf.rank = ifelse(pts$landform %in% c(42,22,50,41,52,53), 2,
                       ifelse(pts$landform %in% c(31,43,33,34,51,4,23,24), 1, 0))
  pts
})

MRgrid <- do.call(what = sf:::rbind.sf, args = MRgrid.l)
MRstat <- do.call(what = sf:::rbind.sf, args = MRstat.l)
st_write(MRgrid, file.path(pipeline, "out", "inClusterSurveyOptions.shp"), delete_layer=T)
st_write(MRClus, file.path(pipeline, "out", "Clusters.shp"), delete_layer=T)

## select one remaining point probabilistically based on distance to water and elevation distribution of survey points in cluster, 
#this will be central for reference grids (also select adjacent points) - or line of points along stream if possible?


###Disturb-free reference ----

#remove all hex grids overlapping disturbance
Rhex <- hex[!st_intersects(x = hex, y = st_union(SDisturb), sparse = FALSE)[,1],] 

# calculate distance to road & distance to Dawson. 
## 139.3539073°W 64.0392738°N 
Dawson = st_sfc(st_point(c(-139.3521370, 64.0384639)))
st_crs(Dawson) = 4269
Dawson <- st_transform(Dawson, crs = crs)
Rhex$d2Daw <- as.numeric(st_distance(Rhex, Dawson))

ggplot(hex) + geom_sf() + geom_sf(data = Rhex, aes(fill = d2Daw)) + 
  geom_sf(data = Dawson, size = 2, colour = "red")

## select hex cells? ##select randomly based on distance
##Or  select clusters non-randomly based on general area, availability, access, watersheds (less disturbed), river order!


## create grid over selected hex
Rgrid <- st_make_grid(Rhex, cellsize = 300, what = "corners")
Rgrid <- st_intersection(st_as_sf(Rgrid), Rhex)
Rgrid <- st_join(Rgrid, ecoregions %>% select(ecoregion = ECOREGION))
#remove points outside of 1 SD of total distance to water, elevation, landform, river order.
waterwayHex <- st_intersection(waterway,Rhex)
waterway.l <- split(waterwayHex, waterwayHex$id)
Rgrid.l <- split(Rgrid, Rgrid$id)
Rgrid.l <- Rgrid.l[names(waterway.l)]
Rgrid.l <- mapply(function(grid, water){
  d2w <- st_distance(grid, water)
  grid$d2w <- d2w
  grid
}, Rgrid.l, waterway.l, SIMPLIFY = F)



Rgrid.l <- lapply(Rgrid.l, function(grid){
  filter(grid, as.numeric(d2w) < quantile(MRstat$d2w, c(.75)))
})

# gmt <- read_stars("CumulativeEffects_GIS/Elevation/GMTED/50n150w_20101117_gmted_med075_ecoregion.tif", proxy = F)
# names(gmt) <- "elev"
Rgrid.l <- lapply(Rgrid.l, function(pts) {
  elev = st_extract(gmt, pts)
  pts$elev <- elev$elev
  pts
})

ggplot(as.data.frame(MRstat), aes(elev, fill = Ecoregn)) + 
  geom_histogram(binwidth = 25) + xlim(c(250, 1200)) + 
  geom_vline(xintercept = quantile(MRstat$elev, c(.1, .9)))

MRK <- as.data.frame(MRstat) %>% filter(Ecoregn == "Klondike")
MRM <- as.data.frame(MRstat) %>% filter(Ecoregn != "Klondike")

Rgrid.l <-lapply(Rgrid.l, function(grid){
  tmp <- grid %>% mutate(log = ifelse(ecoregion == "Klondike Plateau",
                                      ifelse(elev > quantile(MRK$elev, c(.1)) & 
                                               elev < quantile(MRK$elev, c(.9)), T, F),
                                      ifelse(elev > quantile(MRM$elev, c(.1)) & 
                                               elev < quantile(MRM$elev, c(.9)), T, F)))
  filter(tmp, log) %>% select(-log)
})

Rgrid.l <- lapply(Rgrid.l, function(pts) {
  landform = st_extract(mlf, pts)
  pts$landform <- landform$landform
  pts
})

Rgrid.l <- lapply(Rgrid.l, function(pts){
  left_join(pts, mlf_attr %>% select(VALUE, lf.type = TYPE), by = c("landform" = "VALUE"))
})


Rgrid.l <- lapply(Rgrid.l, function(pts) {
  pts$lf.rank = ifelse(pts$landform %in% c(42,22,50,41,52,53), 2,
                       ifelse(pts$landform %in% c(31,43,33,34,51,4,23,24), 1, 0))
  pts
})

## plot the grids in ArcGIS and select accessible hex cells and sites!
Rgrid <- do.call(what = sf:::rbind.sf, args = Rgrid.l)
st_write(Rgrid, file.path(pipeline, "out", "RefSurveyOptions.shp"), delete_layer=T)
st_write(Rhex, file.path(pipeline, "out", "RefHex.shp"), delete_layer=T)

###make secondary grid ----
# hex774 <- filter(Shex, id == 774) %>% st_intersection(st_union(PClaims))
# stat774 <- st_make_grid(hex774, cellsize = 300, what = "corners")
# stat774 <- st_intersection(stat774, st_union(hex774))
# ggplot(filter(Shex, id == 774)) + geom_sf() + geom_sf(data = hex774) + geom_sf(data = stat774)

# mapID <- Shex %>% filter(mapYR == 2023) %>% as.data.frame() %>% select(id, ecoID) %>% distinct()
# stat.pts <- list()
# for(i in 1:length(mapID$id)){
# hex.id <-  filter(Shex, id == mapID$id[i]) %>% st_intersection(st_union(PClaims)) 
# stat.id<- st_make_grid(hex.id, cellsize = 300, what = "corners")
# stat.id <- st_intersection(stat.id, st_union(hex.id))
# # assign(paste0("stat", mapID$id[i]), stat.id)
# stat.pts[[i]] <- stat.id
# }
# names(stat.pts) <- mapID$id
# mapply(function(x, name) write_sf(x, file.path(pipeline, "store", paste0("stat", name, ".shp"))),
#        stat.pts, mapID$id)
# 
# 
# recSites <- read_sf(file.path(pipeline, "ReclaimedMineSurveySites2324_proposed.shp"))
# recSites$SiteID <- as.numeric(recSites$SiteID)
# recSitesNotes <- read_csv(file.path(pipeline, "ReclaimedMineSurveySites2324_notes.csv"))
# recSites <- left_join(recSites, recSitesNotes)
# recSites %>% group_by(Ecoregion, Proposed_Survey_Year) %>% summarise(n=n())
# ggplot(recSites, aes(fill = Proposed_Survey_Year, col = Proposed_Survey_Year)) + geom_sf(size = 4)
# 
# write_sf(recSites, file.path(pipeline, "ReclaimedMineSurveySites2324_proposedNotes.shp"))



###reveg mapping ----
# reveg.m <- rast("CumulativeEffects_GIS/RevegetationMapping/SiteID_GoldfieldsArea-1985-2022-YOD-Mag-Dur.tif")

reveg.m <- read_stars("CumulativeEffects_GIS/RevegetationMapping/SiteID_GoldfieldsArea-1985-2022-YOD-Mag-Dur.tif",
                      proxy = F)
reveg.m <- split(reveg.m, "band")
reveg.m <- mutate(reveg.m, yod = ifelse(yod <= 0, NA, yod))
reveg.m1 <- read_stars("CumulativeEffects_GIS/RevegetationMapping/SiteID_McQuestenArea-1985-2022-YOD-Mag-Dur.tif",
                      proxy = F)
reveg.m1 <- split(reveg.m1, "band")
reveg.m1 <- mutate(reveg.m1, yod = ifelse(yod <= 0, NA, yod))
reveg <- st_mosaic(reveg.m, reveg.m1)

plot(reveg["yod",,])
hist(reveg$yod) ##big peak in 2007

circle = st_transform(filter(MRClus, ClusterID == 19), crs = st_crs(reveg.m)) 
ggplot()  + geom_stars(data = reveg[circle]) + geom_sf(data = circle, fill = NA) + 
  geom_sf(data = st_crop(st_transform(SDisturb, crs = st_crs(circle)), circle), fill = NA, col = "red", size = 2)

write_stars(reveg, "CumulativeEffects_GIS/RevegetationMapping/Reveg_yod.tif")

