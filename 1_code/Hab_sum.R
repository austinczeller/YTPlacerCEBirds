#New Habitat Summary Extraction#

# Load Packages ----
library(tidyverse)
library(lubridate)
library(sf) ##spatial features
library(readxl)
library(biscale)
library(terra)
library(tidyterra)
library(cowplot)
library(ggbreak)
library(wildrtrax)
library(tictoc)
library(climatenaR)

# Load locations----
crs = "EPSG:3579" # "NAD83 / Yukon Albers"
gis_dir <- "C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS"
pipeline<- "2_pipeline/2_HabitatSum"

## Load sites files compiled from both ARU and PC data, from ECCC and WCS
sites <- readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/1_CountDataProcessing/out/sites.RDS")
sites <- sites%>%
  mutate(stationID = paste(siteID, station, sep = "-"))%>%
  filter(!is.na(year))



ggplot(sites) + geom_sf(aes(color=year))

## Ecoregions and ecozone----
biozone <-  rast(file.path(gis_dir, "Bioclimate_zones_and_subzones/Bioclimate_zones_and_Subzones.tif")) #load raster
cats(biozone)
biozone_cat <- read_csv(file.path(gis_dir, "Bioclimate_zones_and_subzones/Bioclimate_zones_and_Subzone_YA_Att.csv"))

biozone_cat <- biozone_cat %>% mutate(SUBZONE = ifelse(is.na(SUBZONE), ZONE, SUBZONE)) %>% select(-OID_)
biozone_cat <- biozone_cat %>% mutate(zone = str_replace_all(ZONE, " ", "_"))

activeCat(biozone) <- 2 ##switch to biozone, not label
# plot(biozone)

biozone.l <- sites %>% 
  st_buffer(dist = 150) %>%
  st_transform(crs = crs(biozone)) %>%
  split(f = .$siteID) %>%
  map(\(site) crop(x = biozone, y = site, mask=T, touches = F))

biozone <- sites %>% 
  st_buffer(dist = 150) %>%
  st_transform(crs = crs(biozone)) %>%
  crop(x = biozone, y = ., mask=T, touches = F)

## Percent cover of different biozones. boreal_low tends to be surrounded by boreal_high, so using a categorical classification will classify most sites as boreal high... might be better to use % cover of boreal low and/or boreal subalpine 

bz_sum <- map(biozone.l, function(x){
  df <- as.data.frame(x) %>% rename(zone = ZONE)
  df <- df %>% group_by(zone) %>% 
    summarise(n = n()) %>% mutate(prop = n/sum(n))
  df %>% select(-n) %>% 
    pivot_wider(names_from = zone, names_prefix = "P.", values_from = prop)
})
bz_sum <- bind_rows(bz_sum, .id = "siteID")  %>%
  complete(siteID = names(bz_sum))
bz_sum[is.na(bz_sum)] <- 0

names(bz_sum) <- str_replace_all(names(bz_sum), " ", "_")


##classify a site based on the biozone that covers the highest proportion
bz_sum <- bz_sum %>% 
  mutate(ecozone = names(.[,-1])[max.col(.[,-1], 'first')]) %>%
  mutate(ecozone = str_replace(ecozone, 'P.', ""))

bioclim.vars <- c("P.Boreal_Low", "P.Boreal_High", "P.Boreal_Subalpine", "P.Boreal_Alpine_Tundra", "P.Subarctic_Woodland", "P.Subarctic_Alpine_Tundra",    "P.Subarctic_Subalpine")

bz_sum$ecozone <- ordered(bz_sum$ecozone, levels = str_replace(bioclim.vars, 'P.', ''))

##ecoregion 
fgdb <- file.path(gis_dir, "Ecoregions_2014_1M.gdb")
# ogrListLayers(fgdb)

ecor <- st_as_sf(vect(fgdb, layer="Ecoregions_2014_1M"))
ecor <- sites %>% group_by(siteID) %>% summarise() %>% st_centroid() %>% st_transform(crs = st_crs(ecor)) %>%
  st_intersection(ecor, .) %>% mutate(ecoregion=ECOREGION)%>%select(siteID,ecoregion)
table(ecor$ecoregion) 

ecor <- left_join(ecor, bz_sum)

## Only relevent ecoregions
table(ecor$ecoregion) ## mostly Klondike Plateau and McQuesten Highlands
ecor %>% ggplot(aes(col = ecoregion)) + geom_sf() #, Yukon Plateau-North also is likely within range and similar enough to be relevant. 
##Remove ecoregions accept Klondike Plateau, McQuesten Highlands and Yukon Plateau-North
bad.ecor <- filter(ecor, !(ecoregion %in% c("Klondike Plateau", "McQuesten Highlands", "Yukon Plateau-North"))) %>% pull(siteID)

sites <- filter(sites, !(siteID %in% bad.ecor))

## Only relevent Biozones
ecor %>% ggplot(aes(col = ecozone)) + geom_sf()
table(ecor$ecozone) ## mostly boreal low or high
##remove boreal alpine, subarctic woodland, subarctic alpine, and subarctic alpine sites
bad.bz <- bz_sum %>% filter(ecozone %in% c("Boreal_Alpine_Tundra", "Subarctic_Woodland", "Subarctic_Alpine_Tundra", "Subarctic_Subalpine")) %>% pull(siteID)

sites <- sites %>% filter(!(siteID %in% c(bad.bz)))
ecor <- ecor %>% filter(!(siteID %in% c(bad.bz, bad.ecor)))


## Add wetland data----
#DUC
DUC_wet <- rast(file.path(gis_dir, "DUCWetlandInv/Duc Dawson Wetland Inventory_Phase_02/DUC_Dawson_Phase_02_Classification/DUC_DawsonP2L2.tif"))
tictoc::tic()
DUC_wet_sf <- sites %>% 
  st_buffer(dist = 1500) %>%
  st_transform(crs = crs(DUC_wet)) %>%
  crop(x = DUC_wet, y = ., mask=T, touches = F) %>%
  as.polygons(.) %>%
  st_as_sf()
tictoc::toc()

DUC_wet_ext <- DUC_wet > 0 ##extent of raster
DUC_wet_ext <- as.polygons(DUC_wet_ext)
DUC_wet_ext <- st_as_sf(DUC_wet_ext)
DUC_wet_ext <- filter(DUC_wet_ext, Class_Name == 1)

DUC_wet_sf <- DUC_wet_sf %>% filter(Class_Name %in% c("Open Water", "Fen", "Bog", "Swamp")) ## DUC wetland polygons

#MAYO
Mayo_wet <- read_sf(file.path(gis_dir, "MayoWetlandClassification/MayoMcQuesten_L2Wetlands_Mar2022.shp"))

tictoc::tic()
Mayo_wet_sf <- sites %>% 
  st_buffer(dist = 1500) %>%
  st_transform(crs = crs(Mayo_wet)) %>%
  st_intersection(Mayo_wet) %>%
  rename(Class_Name = Class) %>%
  group_by(Class_Name) %>% summarise()
tictoc::toc()

#Beaver River
BR_wet <- rast(file.path(gis_dir, "BeaverRiver_WetlandMap/BRLUP_PEM_Wetland_V1.tif")) 

BR_wet_ext <- BR_wet %>% as.polygons(.) %>% st_as_sf() 
BR_wet_ext <- BR_wet_ext %>% filter(BRLUP_PEM_Wetland_V1 != 128) %>% st_union %>% st_convex_hull()

tictoc::tic()
BR_wet_sf <- sites %>% 
  st_buffer(dist = 1500) %>%
  st_transform(crs = crs(BR_wet)) %>%
  crop(x = BR_wet, y = ., mask=T, touches = F) %>%
  as.polygons(.) %>%
  st_as_sf()
tictoc::toc()

BR_wet_sf  <-  BR_wet_sf %>% left_join(data.frame(BRLUP_PEM_Wetland_V1 = c(1,2, 3,4,5,6,7,8, 128), Class = c("Bog", "Exposed Fluvial", "Fen", "Marsh", "Swamp", "River", "Lake", "Shallow Water", "Other")))

BR_wet_sf <- BR_wet_sf %>% rename(Class_Name = Class) %>%
  filter(!Class_Name %in% c("Other")) 

#Combine wetland layers
## reclassify so wetland classes are the same
DUC_wet_sf <- DUC_wet_sf %>% select(class = Class_Name) %>% 
  mutate(class = if_else(class == "Open Water", "Water", class),
         layer = "DUC") ##Open water
Mayo_wet_sf <- Mayo_wet_sf %>% select(class = Class_Name) %>% 
  mutate(class = if_else(class %in% c("Deep Water", "Shallow Water"), "Water", class), 
         layer = "Mayo")#shallow water
BR_wet_sf <- BR_wet_sf %>% select(class = Class_Name) %>% 
  mutate(class = if_else(class %in% c("River", "Lake", "Shallow Water", "Exposed Fluvial"), "Water", class),
         layer = "BR") ## river, lake, shallow water

## Conver to same CRS
crs <- "EPSG:3579"

DUC_wet_sf <- DUC_wet_sf %>% 
  st_transform(crs = crs)
Mayo_wet_sf <- Mayo_wet_sf %>% 
  st_transform(crs = crs) %>%
  st_difference(st_transform(DUC_wet_ext, crs = crs))  #remove areas where overlaps with DUC
wetland <- bind_rows(list(DUC_wet_sf, Mayo_wet_sf)) # combine into single feature

## Beaver River overlaps with Mayo, so remove overlap area before mergins. 
Mayo_wet_ext <- read_sf(file.path(gis_dir, "MayoWetlandClassification/MayoWetland_extent.shp")) ## Extent of mayo map
wetland_ext <- bind_rows(list(st_transform(DUC_wet_ext, crs = crs),
                              st_transform(Mayo_wet_ext, crs = crs))) # combine into single feature

BR_wet_sf <- BR_wet_sf %>% 
  st_difference(Mayo_wet_ext %>% st_transform(crs = st_crs(BR_wet_sf))) %>%
  st_transform(crs = crs) 

tmp <- BR_wet_ext
BR_wet_ext <- tmp
BR_wet_ext <- st_difference(BR_wet_ext, Mayo_wet_ext%>% st_transform(crs = st_crs(BR_wet_ext)) %>% st_union())

BR_wet_ext <- st_as_sf(BR_wet_ext) %>% rename(geometry = x)
# BR_wet_ext <- st_collection_extract(BR_wet_ext, type = c("POLYGON", "MULTIPOLYGON"))

wetland <- bind_rows(list(wetland, BR_wet_sf))

## same crs as wetland, label the layer
##beaver River only overlaps partially with one site that is half overlapped by Mayo layer, so this site will be assigned 'mayo', and we don't need the BR extend. 
wetland_ext <- bind_rows(
  list(st_transform(DUC_wet_ext, crs = st_crs(wetland)) %>% 
         mutate(layer = "DUC") %>% select(layer), 
       st_transform(Mayo_wet_ext, crs = st_crs(wetland)) %>% 
         mutate(layer = "Mayo")%>% select(layer),
       st_transform(BR_wet_ext, crs = st_crs(wetland)) %>% 
         mutate(layer = "Mayo")%>% select(layer))) #half overlaps one site

ggplot(wetland_ext %>% st_transform(crs = "EPSG:4326")) + 
  geom_sf(aes(fill = layer)) + 
  geom_sf(data = sites, aes(colour = factor(year))) ## not all sites have maps 

#Remove sites without wetland data
no.wet.stat.v <- sites$stationID[!st_intersects(sites, st_union(wetland_ext) %>% 
                                                  st_transform(crs = st_crs(sites)),
                                                sparse = F)]
no.wet.v <- unique(sites$siteID[!st_intersects(sites, st_union(wetland_ext) %>% 
                                                 st_transform(crs = st_crs(sites)),
                                               sparse = F)])
sites <- sites %>% mutate(wetland.layer = ifelse(stationID %in% no.wet.stat.v, F, T))
ggplot(sites) + geom_sf(data = wetland_ext) + geom_sf(aes(col = wetland.layer)) ## a lot of the sites around steward crossing are going to be lost. These also coincide with sites in the Yukon Plateau-North ecoregion

sites<-sites%>%filter(wetland.layer==T)

## Add elevation data----
elev <- rast(file.path(gis_dir, "GMTED/50n150w_20101117_gmted_med075.tif"))
names(elev) <- "elevation"
station_elev <- sites %>% 
  st_buffer(dist = 150) %>%
  st_transform(crs = crs(elev)) %>%
  split(f = .$stationID) %>%
  map(\(stat) crop(x = elev, y = stat, mask=T, touches = F))

#mean elevation at each station
elev_mn <- station_elev %>%
  map(\(stat) global(stat, 'mean', na.rm = T)) %>% bind_rows() %>% mutate(stationID = names(station_elev))

## Add digitized linear and surface disturbance data----
buff_sites<-st_as_sf(sites)%>%
  st_buffer(dist = 1000)%>%
  st_transform(crs=crs)

YG_surfacedisturbance.2021<-read_sf("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS/YG_SurfaceDisturbance_May2022/ArealFeatures_May2022/ArealFeatures_May2022_Un_Dis_J.shp")%>%
  st_transform(crs=crs)%>%
  select(Shape_Leng, Shape_Area, geometry)

Clara_digitizations_2023<-read_sf(file.path(gis_dir,"ClaraDigitization_2023/ArealDist_Added.shp"))%>%
  st_transform(crs = crs)%>%
  select(Shape_Leng, Shape_Area, geometry)


#Add different digitizations here and st_join them to the master

Master_surface_2023<-st_join(YG_surfacedisturbance.2021,Clara_digitizations_2023)

buff_21_23<-buff_sites%>%
  filter(year%in%c(2021:2023))

surface_2023<-Master_surface2023%>%
  st_intersection(buff_21_23)

buff_24<-buff_sites%>%filter(year==2024)


plot(surface_2023,max.plot=1)



# I just manually created 2024 surface and linear files via ARC
lines_2024<-read_sf(file.path(gis_dir,"2024_final_lines.shp")) 
surface_2024<-read_sf(file.path(gis_dir,"Final_2024_surface.shp"))%>%
  st_transform(crs = crs)
plot(surface_2024,max.plot=1)



line_clip_21_23<-line_clip_21_23%>%
  st_difference(surface_2023)


st_write(line_clip_21_23,ds=file.path(gis_dir,"temp.shp"))

site.v <- unique(sites$siteID)
missing.d <- site.v[site.v %in% c(bad.bz, no.wet.v)]

#2023 surface disurbance data is p.dist.23
p.dist.23 <- read_sf(file.path(gis_dir,"ClaraDigitization_2023/Sites_ArealDist_Intersect.shp")) %>%
  st_transform(crs = crs) %>%
  select(Shape_Leng, Shape_Area, geometry) #Clara's digitization

p.dist.23added <- read_sf(file.path(gis_dir,"ClaraDigitization_2023/ArealDist_Added_XY.shp")) %>%
  st_transform(crs = crs)
plot(p.dist.23added%>%select(Shape_Area))
p.dist.23 <- rbind(p.dist.23, p.dist.23added)

### The original 2021 digitization
# Read in disturbance polygons
fgdb <- file.path(gis_dir, "CumEff_corrected_mjb.gdb") ## 2021 digitizations by WCSC
fc_list <- st_layers(fgdb)
print(fc_list)
p.dist <- st_as_sf(vect(fgdb, layer="Polygon_disturbance_VegBareWater"))

fgdb<-file.path(gis_dir,"CumEff_corrected_mjb.gdb")
p_dist_ns <- vect(fgdb, layer = "Polygon_disturbance_NS_NoOverlap")
p_dist_a_ns <- vect(fgdb, layer = "Polygon_disturbance_added_NS")
p_dist_hs <- vect(fgdb, layer = "Polygon_disturbance_HS_cor")
p_dist_a_hs <- vect(fgdb, layer = "Polygon_disturbance_added_HS_cor")

#clean files, and make more uniform
p_dist_ns <- p_dist_ns %>% select(siteID = Id, lat_dist = lat, lon_dist = long, ##lat/lon in
                                  area_m2 = Shape_Area,
                                  type = TYPE_INDUSTRY,
                                  subtype = TYPE_DISTURBANCE,
                                  database = DATABASE,
                                  image = IMAGE_NAME,
                                  image_date = IMAGE_DATE,
                                  image_res = IMAGE_RESOLUTION,
                                  image_sensor = IMAGE_SENSOR,
                                  pc_disturb = disturb, ##ranking from field
                                  pc_intensity = intensity, ##ranking from field
                                  pc_type = type ##classification from field
) %>% mutate(gis_layer = "Polygon_disturbance_NS_NoOverlap")

p_dist_a_ns <- p_dist_a_ns %>% select(siteID = Id, lat_dist = lat, lon_dist = long, ##lat/lon in
                                      area_m2 = Shape_Area,
                                      type = TYPE,
                                      subtype = SUBTYPE,
                                      database = ORIGINAL_DB,
                                      image = IMAGE_NAME,
                                      image_date = IMAGE_DATE,
                                      image_res = IMAGE_RESOLUTION,
                                      image_sensor = IMAGE_SENSOR,
                                      pc_disturb = disturb, ##ranking from field
                                      pc_intensity = intensity, ##ranking from field
                                      pc_type = type_1 ##classification from field
)  %>% mutate(gis_layer = "Polygon_disturbance_added_NS") %>% select(names(p_dist_ns))

p_dist_ns <- rbind(p_dist_ns, p_dist_a_ns)


#clean files, and make more uniform
p_dist_hs <- p_dist_hs %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                  area_m2 = Shape_Area,
                                  type = TYPE_INDUSTRY,
                                  subtype = TYPE_DISTURBANCE,
                                  database = DATABASE,
                                  image = IMAGE_NAME,
                                  image_date = IMAGE_DATE,
                                  image_res = IMAGE_RESOLUTION,
                                  image_sensor = IMAGE_SENSOR
) %>% mutate(gis_layer = "Polygon_disturbance_HS_cor")

p_dist_a_hs <- p_dist_a_hs %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1, ##lat/lon in
                                      area_m2 = Shape_Area,
                                      type = TYPE,
                                      subtype = SUBTYPE,
                                      database = ORIGINAL_DB,
                                      image = IMAGE_NAME,
                                      image_date = IMAGE_DATE,
                                      image_res = IMAGE_RESOLUTION,
                                      image_sensor = IMAGE_SENSOR
)  %>% mutate(gis_layer = "Polygon_disturbance_added_HS_cor") %>% select(names(p_dist_hs))

p_dist_hs <- rbind(p_dist_hs, p_dist_a_hs)

p_dist_hs <- st_as_sf(p_dist_hs)
p_dist_ns<-st_as_sf(p_dist_ns)

p_dist_hs <- p_dist_hs %>% separate(col = siteID, into = c("location_3", "location_4"), sep = "-") %>% ##splits into 6 number code + letter/number
  mutate(siteID = paste(location_3, str_extract(location_4, "[A-Z]+" ), sep = "-")) %>% ##hex code + letter
  select(-starts_with("location_") )

p.dist.2021 <- rbind(p_dist_ns %>% select(type, subtype), p_dist_hs %>% select(type, subtype))


p.dist.2021<- st_transform(p.dist.2021, st_crs(sites))

###this filters out sites already in p.dist
p.dist.2021 <- st_intersection(sites %>% filter(siteID %in% missing.d) %>%
                                 st_buffer(1000),
                               p.dist.2021)


###I AM HERE

## Merge p.dist for sites <= 2021

p.dist <- st_union(p.dist) ## re-digitized in 2022 to identify bare, veg and water
p.dist.2021 <- st_union(p.dist.2021) ## sites from 2021 or earlier that weren't included in p.dist
p.dist <- rbind(st_as_sf(p.dist), st_as_sf(p.dist.2021)) ## all 2021 sites in a single layer

## disturbance is updated annually (e.g. imagery from year of survey or earlier was used)
## This creates a list, where each element corresponds to the survey year. Disturbance data from within this year are cropped within a 2000m buffer of each station surveyed in the given year
tictoc::tic()
p.dist.yr <- sites %>%
  st_buffer(dist = 1000) %>%
  st_transform(crs = st_crs(p.dist)) %>%
  split(sites$year) %>%
  map(\(x) st_intersection(x = p.dist, y = st_union(x)))
tictoc::toc()

##replace 2023 disturbance data with data from Clara's digitization (includes both YG and her added areas)
p.dist.yr[[4]] <- sites %>% 
  filter(year == 2023) %>%
  st_buffer(dist = 1000) %>%
  st_transform(crs = st_crs(p.dist.23)) %>%
  st_union() %>%
  st_intersection(x = p.dist.23, y = .) %>% st_union() %>%
  st_as_sf()

##replace 2024 disturbance data with data from Clara's digitization (includes both YG and her added areas)
p.dist.yr[[5]] <- sites %>% 
  filter(year == 2024) %>%
  st_buffer(dist = 1000) %>%
  st_transform(crs = st_crs(p.dist.23)) %>%
  st_union() %>%
  st_intersection(x = surface_2024, y = .) %>% st_union() %>%
  st_as_sf()

## NDVI: vegmask and Watermask data----
files <- list.files("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GEE", pattern = "B3B8")
files <- files[!str_detect(files, "tif.")]
watermask <- pmap(list(p.dist.yr, c("2017", "2018", "2021","2023","2024")),
                  function(p.dist, yr){
                    file <- files[str_which(files, yr)]
                    r <- rast(file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GEE", file))
                    crs.r <- terra::crs(r)
                    r = crop(x = r, y = p.dist %>% st_transform(crs = crs.r), mask=T, touches = F)
                    r <- r > 0.046
                    r
                  })
write_rds(watermask, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "watermask_dclip_rast.rds"))

files <- list.files("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GEE/", pattern = "ndvi")
files <- files[!str_detect(files, "tif.")] ## removes potential tif.aux files, e.g. if ArcGIS is open
vegmask <- pmap(list(p.dist.yr, c("2017", "2018","2021","2023","2024")), function(p.dist, yr){
  file <- files[str_which(files, yr)]
  r <- rast(file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GEE/", file))
  r = crop(x = r, y = p.dist %>% st_transform(crs = terra::crs(r)), mask=T, touches = F)
  r > 0.42
})
write_rds(vegmask, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "vegmask_dclip_rast.rds"))

p.dist <- pmap(list(vegmask, watermask), function(v, w){
  v <- as.numeric(v)
  mask(x = v, mask = w, maskvalues = TRUE, updatevalue = 3)
})  ## 1 = veg, 0 = bare, 3 = water
write_rds(p.dist, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "p.dist_rast.rds"))

type.df <- data.frame(typeID = c(1, 0, 3), type = c("d_pVeg", "d_pBare", "d_pWater"))
p.dist <- map(p.dist, function(r){
  p <- as.polygons(r)
  names(p) = "typeID"
  p = st_as_sf(p) %>% left_join(type.df)
})
write_rds(p.dist, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/",pipeline, "store", "p.dist_p.rds"))

## add site ID
sites.yr <- sites %>% st_buffer(1000) %>% 
  group_by(year) %>% summarise() %>% 
  st_transform(crs = crs <- "EPSG:3579")%>%
  split(.$year) 

saveRDS(sites.yr,file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",
                  pipeline, "store","sites_yr.RDS"))



## Add Fire data----

fire <- st_read(file.path(gis_dir, "Fire_History.shp/Fire_History.shp"))

names(fire) <- tolower(names(fire))

tictoc::tic()
fire <- sites %>% 
  st_buffer(dist = 1000) %>%
  st_union() %>%
  st_transform(crs = st_crs(fire)) %>%
  st_intersection(x = fire, y = .)
tictoc::toc()

hist(fire$fire_year)
# Attribute values are assigned to sub-geometries; if these are spatially
# constant, as for instance for land use, then this is fine. If they are
# aggregates, such as population count, then this is not fine

## 9000 series seem to correspond with last two numbers of the decade.  Set to mid-decade for an average year?
fire <- mutate(fire, fire_year = ifelse(fire_year > 2030,
                                        decade + 5, fire_year)) %>%
  select(fire_year) %>% 
  group_by(fire_year) %>% summarise()


## NASA ABoVE----
nasa_above <- rast(file.path(gis_dir, "Nasa_Above/Nasa_Above.tif"))
nasa_above <- nasa_above[[31]]# I think that this is taking the last year in the data (2014)
names(nasa_above) <- "lc" ##lc = land classification

landc_cat <- data.frame(
  lc = 1:15,
  class = c("evergreen_forest", "deciduous_forest", "mixed_forest", "woodland", 
            "shrub_low", "shrub_tall", "shrub_open", "herb", "tussock_tundra",
            "sparce_veg", "fen_nasa", "bog_nasa", "shallows", "barren", "water"),
  class10 = c("evergreen_forest", "deciduous_forest", "deciduous_forest",
              "evergreen_forest", "shrub", "shrub", "shrub", "herb", "herb",
              "sparce_veg",  "fen_nasa", "bog_nasa", "shallows", "barren", "water"))

landc_cat[landc_cat$class == "shallows", c("class", "class10")] <- "water"

##clip to 1500m buffer around each station

tic()
landc.ls <- sites %>%
  split(.$siteID) %>%
  map(function(s){
    st_buffer(s, dist = 1000) %>% st_union() %>%
      st_transform(crs = crs(nasa_above)) %>%
      crop(x = nasa_above, y=., mask = T, touches = F) %>%
      as.polygons() %>% st_as_sf() %>%
      left_join(landc_cat)
  })
toc()

landc.ls <- landc.ls %>% bind_rows() %>%
  group_by(class) %>% summarise()  
# reduces each category to a multipolygon. This make computation faster. 

## NASA ABoVE 2020----
tif_files<-dir("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/GIS/Macander",
               pattern=".tif", full.names=T, recursive=F)

landclass<-rast(tif_files)
landclass


## Combined habitat polygons----
site.buff <- sites %>% ## recrop to selected sites, smaller buffer
  st_buffer(dist = 1500) %>%
  group_by(siteID, year) %>% 
  summarise() %>%
  st_transform(crs = crs)

## transform to same crs
p.dist <- map(p.dist, function(site) site %>% st_transform(crs = crs))
sites.yr <- map(sites.yr, function(site) site %>% st_transform(crs = crs))
wetland<-st_as_sf(wetland)
wetland <- wetland %>% st_transform(crs = crs)
fire <- fire %>% st_transform(crs = crs)
landc.ls<- landc.ls%>%st_transform(crs=crs)

## crop all to sites.yr
wetland <- map(sites.yr, function(site){
  st_intersection(site, wetland %>% 
                    select(class, layer)) %>%
    group_by(class, layer) %>% summarize()
})
landc.ls <- map(sites.yr, function(site){
  st_intersection(site, landc.ls %>% select(class)) %>%
    group_by(class) %>% summarize() 
})

## classify fire by age, and remove older burns >= 30 years old
fire <- map(sites.yr, function(site){
  fr = st_intersection(site, fire)
  fr %>% mutate(diff = year - fire_year, 
                class = ifelse(diff %in% 0:9, "Burn-herb",
                               ifelse(diff %in% 10:19,"Burn-shrub", 
                                      ifelse(diff %in% 20:29, "Burn-sappling", "Burn-old"))),
                layer = "fire") %>% 
    filter(class != "Burn-old") %>% ## remove old burns and burns that haven't happened
    arrange(diff) %>% st_difference() %>% #remove overlapping burns, selecting most recent
    group_by(class, layer) %>% summarise()
})

##  separate wetlands into natural and disturbed--
wetland <- pmap(list(p.dist, wetland), function(disturb, wet){
  h_wt <- st_difference(wet, st_union(disturb)) ## remove wetlands in disturbance, leaving natural wetlands
  h_wt$class = paste0("H_", tolower(h_wt$class))
  wt <- st_intersection(wet, st_union(disturb)) ## wetland in disturbance area
  d_water <- st_difference(disturb %>% filter(type == "d_pWater"), st_union(wt)) ## mapped water not in wetland area
  wt$class = paste0("D_", tolower(wt$class))
  d_water %>% select(class = type) %>% mutate(layer = "NDWI") %>% select(names(wt)) %>%
    rbind(wt) %>% rbind(h_wt)
})

### disturb > fire > landclass -
hab.poly <- map(p.dist, function(disturb){
  disturb <- disturb %>% 
    mutate(class = type, layer = "YGfootprint_CWS_update") %>% 
    select(class, layer) 
  
})

## remove disturb from burnt areas
## fire age/class varies by siteID
fire <- pmap(list(fire, hab.poly), function(fr, hp){
  fr = fr %>%  mutate(class = paste0("H_", tolower(class)))
  st_difference(fr, st_union(hp))
})

### merge fire, and disturbance 
hab.poly <- pmap(list(hab.poly, fire), function(x, y){
  y <- y  %>% select(class, layer)
  bind_rows(list(x,y))
})

# erase disturb + fires from nasa above landclass
landc.ls <- pmap(list(hab.poly, landc.ls), function(x, y){
  y <- st_difference(y, st_union(x))
  y$class = paste0("H_", tolower(y$class))
  y
})

hab.poly <- pmap(list(hab.poly, landc.ls), function(x, y){
  y <- y %>% mutate(layer = "nasa_above") %>% select(class, layer)
  bind_rows(list(x,y))
})

## by site
hab.poly.site <- pmap(list(hab.poly, site.buff %>% split(.$year)), 
                      function(hab, sy){
                        site <- split(sy, sy$siteID)
                        map(site, function(s){
                          st_intersection(s, hab) 
                        }) %>% 
                          bind_rows()
                      }) %>% bind_rows()

wetland.site <- pmap(list(wetland, site.buff %>% split(.$year)), 
                     function(hab, sy){
                       site <- split(sy, sy$siteID)
                       map(site, function(s){
                         st_intersection(s, hab) 
                       }) %>% 
                         bind_rows()
                     }) %>% bind_rows()

wetland.site <- wetland.site %>% mutate(class2 = ifelse(str_starts(class,"H_"), "h_wet", "d_wet")) 

##station.  first need to match station to site.  
hab.poly.station <- pmap(list(hab.poly, 
                              sites %>% st_buffer(1500) %>% split(.$year)), 
                         function(hab, sy){
                           site <- split(sy %>% select(siteID, station, stationID), sy$stationID)
                           map(site, function(s){
                             st_intersection(s, hab) 
                           }) %>% 
                             bind_rows()
                         }) %>% bind_rows()

wetland.station <- pmap(list(wetland, 
                             sites %>% st_buffer(1500) %>% split(.$year)), 
                        function(hab, sy){
                          site <- split(sy %>% select(siteID, station, stationID), sy$stationID)
                          map(site, function(s){
                            st_intersection(s, hab) 
                          }) %>% 
                            bind_rows()
                        }) %>% bind_rows()

wetland.station <- wetland.station %>% mutate(class2 = ifelse(str_starts(class,"H_"), "h_wet", "d_wet"))

## Save Combined habitat polys----
saveRDS(hab.poly.station, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "lcmerge_stat.RDS"))
saveRDS(hab.poly.site, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "lcmerge.RDS"))
saveRDS(wetland.site, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "wetland_site.RDS"))
saveRDS(wetland.station, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "wetland_station.RDS"))

## Linear Disturbance----
l_dist_2022 <- read_sf(file.path(gis_dir, "YG_SurfaceDisturbance_May2022/LinearFeatures_May2022_NoOverlap.shp"))

l_dist_2022 <- l_dist_2022 %>% select(length_m = Shape_Leng, 
                                      width_m = WIDTH_M,
                                      width_class = WIDTH_CLAS,
                                      type = TYPE_INDUS,
                                      subtype = TYPE_DISTU,
                                      database = DATABASE,
                                      image = IMAGE_NAME,
                                      image_date = IMAGE_DATE,
                                      image_res = IMAGE_RESO,
                                      image_sensor = IMAGE_SENS) %>%
  st_transform(crs = st_crs(sites))

##2023 digitization, added ontop of 2022 layer 
l_dist_2023 <- read_sf(file.path(gis_dir, "ClaraDigitization_2023/LinDist_Added_XY.shp")) ## can add to l_dist_2022 because it was created with the 2022 layer as the base

l_dist_2023 <- l_dist_2023 %>% rename(length_m= Shape_Leng, 
                                      width_m = Width_M,
                                      width_class = Width_Clas,
                                      type = Type,
                                      subtype = Subtype,
                                      image = Image_Name,
                                      image_date = Image_Date,
                                      image_res = Image_Reso,
                                      image_sensor = Image_Sens) %>%  
  mutate(database = NA)  %>% st_transform(crs = st_crs(sites))

l_dist_2024<- read_sf(file.path(gis_dir,"2024_final_lines.shp"))
head(l_dist_2024)
l_dist_2024<-l_dist_2024%>%rename(length_m= Shape_Leng, 
                                  width_m = Width_M,
                                  width_class = Width_Clas,
                                  type = Type,
                                  subtype = Subtype,
                                  image = Image_Name,
                                  image_date = Image_Date,
                                  image_res = Image_Reso,
                                  image_sensor = Image_Sens)%>%
  mutate(database = NA)%>% st_transform(crs = st_crs(sites))

l_dist_2022 <- bind_rows(list(l_dist_2023, l_dist_2022,l_dist_2024)) ## can add to l_dist_2022 because it was created with the 2022 layer as the base



table(l_dist_2022$width_class) 
l_dist_2022[l_dist_2022$width_class %in% c("High", "High (>8m)"),]$width_class <- "HIGH"
l_dist_2022[l_dist_2022$width_class %in% c("Med (4-8m)"),]$width_class <- "MED"
l_dist_2022[l_dist_2022$width_class %in% c("Low (<4m)"),]$width_class <- "LOW"


### 2021 digitization by WCSC, made over an earlier version of the human footprint layer.
fgdb <- file.path(gis_dir, "CumEff_corrected_mjb.gdb")
l_nrn_ns <- st_as_sf(vect(fgdb, layer="NRNYukonRoadSegments_NS_cor"))
l_nrn_hs <- st_as_sf(vect(fgdb, layer="NRNYukonRoadSegments_HS"))
l_nrn_hs2 <- st_as_sf(vect(fgdb, layer="NRN_Disturbance_368280_D2"))
l_dist_ns <- st_as_sf(vect(fgdb, layer="Line_disturbance_NS_cor"))
l_dist_a_ns <- st_as_sf(vect(fgdb, layer="Line_disturbance_added_NS"))
l_dist_hs <- st_as_sf(vect(fgdb, layer="Line_disturbance_HS_cor"))
l_dist_a_hs <- st_as_sf(vect(fgdb, layer="Line_disturbance_added_HS_cor"))
l_dist_hs2 <-  st_as_sf(vect(fgdb, layer="Line_Disturbance_368280_D2"))

l_nrn_ns <- l_nrn_ns %>% select(siteID = Id, lat_dist = lat, lon_dist = long, ##lat/lon in 
                                length_m = Shape_Leng, 
                                n_lanes = NBRLANES,
                                subtype = ROADCLASS,
                                pc_disturb = disturb, ##ranking from field
                                pc_intensity = intensity, ##ranking from field
                                pc_type = type ##classification from field
) %>% 
  mutate(gis_layer = "NRNYukonRoadSegments_NS_cor", type = "Transportation",
         database = "NRNYukonRoadSegments", 
         image = NA, image_date = NA, image_res = NA,
         image_sensor = NA, 
         width_m = NA,
         width_class = NA)

l_dist_ns <- l_dist_ns %>% select(siteID = Id, lat_dist = lat, lon_dist = long, ##lat/lon in 
                                  length_m = Shape_Length, 
                                  width_m = WIDTH_M,
                                  width_class = WIDTH_CLASS,
                                  type = TYPE_INDUSTRY,
                                  subtype = TYPE_DISTURBANCE,
                                  database = DATABASE,
                                  image = IMAGE_NAME,
                                  image_date = IMAGE_DATE,
                                  image_res = IMAGE_RESOLUTION,
                                  image_sensor = IMAGE_SENSOR,
                                  pc_disturb = disturb, ##ranking from field
                                  pc_intensity = intensity, ##ranking from field
                                  pc_type = type ##classification from field
) %>% mutate(gis_layer = "Line_disturbace_NS_cor", n_lanes = NA) %>%
  select(names(l_nrn_ns))

l_dist_a_ns <- l_dist_a_ns %>% select(siteID = Id, lat_dist = lat, lon_dist = long, ##lat/lon in 
                                      length_m = SHAPE_Length, 
                                      width_m = WIDTH_M,
                                      width_class = WIDTH_CLASS,
                                      type = TYPE,
                                      subtype = SUBTYPE,
                                      database = ORIGINAL_DB,
                                      image = IMAGE_NAME,
                                      image_date = IMAGE_DATE,
                                      image_res = IMAGE_RESOLUTION,
                                      image_sensor = IMAGE_SENSOR,
                                      pc_disturb = disturb, ##ranking from field
                                      pc_intensity = intensity, ##ranking from field
                                      pc_type = type_1 ##classification from field
)  %>% mutate(gis_layer = "Line_disturbance_added_NS", n_lanes = NA) %>%
  select(names(l_nrn_ns))


l_nrn_hs <- l_nrn_hs %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                length_m = Shape_Length, 
                                n_lanes = NBRLANES,
                                subtype = ROADCLASS,
) %>% mutate(gis_layer = "NRNYukonRoadSegments_HS_cor", type = "Transportation",
             database = "NRNYukonRoadSegments", 
             image = NA, image_date = NA, image_res = NA,
             image_sensor = NA, 
             width_m = NA,
             width_class = NA)

l_nrn_hs2 <- l_nrn_hs2 %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                  length_m = Shape_Length, 
                                  n_lanes = NBRLANES,
                                  subtype = ROADCLASS,
) %>% mutate(gis_layer = "NRN_Disturbance_368280_D2", type = "Transportation",
             database = "NRNYukonRoadSegments", 
             image = NA, image_date = NA, image_res = NA,
             image_sensor = NA, 
             width_m = NA,
             width_class = NA)

l_nrn_hs <- rbind(l_nrn_hs, l_nrn_hs2)##buffer missed in first round, so merge clip

l_dist_hs <- l_dist_hs %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                  length_m = Shape_Length, 
                                  width_m = WIDTH_M,
                                  width_class = WIDTH_CLASS,
                                  type = TYPE_INDUSTRY,
                                  subtype = TYPE_DISTURBANCE,
                                  database = DATABASE,
                                  image = IMAGE_NAME,
                                  image_date = IMAGE_DATE,
                                  image_res = IMAGE_RESOLUTION,
                                  image_sensor = IMAGE_SENSOR,
) %>% mutate(gis_layer = "Line_disturbance_HS_cor", n_lanes = NA) %>%
  select(names(l_nrn_hs))

##buffer missed in first round, so merge clip
l_dist_hs2 <- l_dist_hs2 %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                    length_m = Shape_Length, 
                                    width_m = WIDTH_M,
                                    width_class = WIDTH_CLASS,
                                    type = TYPE_INDUSTRY,
                                    subtype = TYPE_DISTURBANCE,
                                    database = DATABASE,
                                    image = IMAGE_NAME,
                                    image_date = IMAGE_DATE,
                                    image_res = IMAGE_RESOLUTION,
                                    image_sensor = IMAGE_SENSOR,
) %>% mutate(gis_layer = "Line_Disturbance_368280_D2", n_lanes = NA) %>%
  select(names(l_nrn_hs))

l_dist_hs <- rbind(l_dist_hs, l_dist_hs2)##buffer missed in first round, so merge clip

l_dist_a_hs <- l_dist_a_hs %>% select(siteID = loc_name, lat_dist = location_l, lon_dist = location_1,
                                      length_m = Shape_Length, 
                                      width_m = WIDTH_M,
                                      width_class = WIDTH_CLASS,
                                      type = TYPE,
                                      subtype = SUBTYPE,
                                      database = ORIGINAL_DB,
                                      image = IMAGE_NAME,
                                      image_date = IMAGE_DATE,
                                      image_res = IMAGE_RESOLUTION,
                                      image_sensor = IMAGE_SENSOR,
)  %>% mutate(gis_layer = "Line_disturbance_added_HS_cor", n_lanes = NA) %>%
  select(names(l_nrn_hs))

## combine the 2021 digitization files
l_dist_ns <- rbind(st_transform(l_nrn_ns, crs = st_crs(l_dist_ns)), l_dist_ns, l_dist_a_ns)

l_dist_hs<-st_transform(l_dist_hs,crs=st_crs(l_dist_a_hs))

#SET SAME CRS?!! 
l_dist_hs <- rbind(st_transform(l_nrn_hs, crs = st_crs(l_dist_hs)), l_dist_hs, l_dist_a_hs)
l_dist_hs <- l_dist_hs %>% separate(col = siteID, into = c("location_3", "location_4"), sep = "-") %>% ##splits into 6 number code + letter/number
  mutate(siteID = paste(location_3, str_extract(location_4, "[A-Z]+" ), sep = "-")) %>% ##hex code + letter
  select(-starts_with("location_"), )

l_dist <- rbind(l_dist_ns %>% select(names(l_dist_hs)), l_dist_hs)

table(l_dist$width_class) 
# table(l_dist$subtype)
# table(is.na(l_dist$subtype))

##### line widths 
## update missing width classes based on measured widths or average width/class for that subtype
table(is.na(l_dist$width_class))
table(is.na(l_dist_2022$width_class))

ggplot(l_dist, aes(width_class, width_m)) + geom_boxplot()

quantile(l_dist_2022[l_dist_2022$width_class == "LOW",]$width_m, c(0, .05, .25, .5, .75, .95, .99, 1), na.rm = T) #1 - 4
quantile(l_dist_2022[l_dist_2022$width_class == "MED",]$width_m, c(0, .05, .25, .5, .75, .95, 1), na.rm = T) #4-8
quantile(l_dist_2022[l_dist_2022$width_class == "HIGH",]$width_m, c(0, .05, .25, .5, .75, .95, 1), na.rm = T) # 8+

table(l_dist_2022[is.na(l_dist_2022$width_class) & l_dist_2022$width_m == 0,]$subtype) ## assume 0 width is NA value
l_dist_2022[is.na(l_dist_2022$width_class) & l_dist_2022$width_m == 0,]$width_m <- NA
l_dist_2022 <- mutate(l_dist_2022, width_class = ifelse(!is.na(width_class), width_class, ## don't replace existing widths
                                                        ifelse(is.na(width_m), NA, ## can't select by NA, so keep NA
                                                               ifelse(width_m <4, "LOW", 
                                                                      ifelse(width_m > 8, "HIGH", "MED")))))

table(l_dist[is.na(l_dist$width_class) & l_dist$width_m == 0,]$subtype) ## assume 0 width is NA value
l_dist <- mutate(l_dist, width_class = ifelse(!is.na(width_class), width_class, ## don't replace existing widths
                                              ifelse(is.na(width_m), NA, ## can't select by NA, so keep NA
                                                     ifelse(width_m <4, "LOW", 
                                                            ifelse(width_m > 8, "HIGH", "MED")))))



## No NRN data have widths or width classes, assume NRN data are all > 4m wide (e.g. med-high class?)?
l_dist[is.na(l_dist$width_class) &
         l_dist$database == "NRNYukonRoadSegments" & 
         !is.na(l_dist$database),]$width_class <- "HIGH"

##remaining unknowns
l_dist_tmp <- rbind(l_dist %>% select(subtype, width_class),
                    l_dist_2022 %>% select(subtype, width_class))

##
# filter(l_dist_tmp, subtype == "Electric Utility Corridor") %>% group_by(width_class) %>% summarise(n = n()) ##high
l_dist[is.na(l_dist$width_class) & 
         l_dist$subtype == "Electric Utility Corridor" & 
         !is.na(l_dist$subtype),]$width_class <- "HIGH"
l_dist_2022[is.na(l_dist_2022$width_class) & 
              l_dist_2022$subtype == "Electric Utility Corridor" & 
              !is.na(l_dist_2022$subtype),]$width_class <- "HIGH"

## All HS local roads unknown, NS local roads mostly medium
# filter(l_dist_tmp, subtype == "Local Road") %>% group_by(width_class) %>% summarise(n = n())##most high
l_dist[is.na(l_dist$width_class) & 
         l_dist$subtype == "Local Road" &
         !is.na(l_dist$subtype),]$width_class <- "HIGH"
l_dist_2022[is.na(l_dist_2022$width_class) & 
              l_dist_2022$subtype == "Local Road" &
              !is.na(l_dist_2022$subtype),]$width_class <- "HIGH"

# filter(l_dist_tmp, subtype == "Arterial Road") %>% group_by(width_class) %>% summarise(n = n())
l_dist_2022[is.na(l_dist_2022$width_class) & 
              l_dist_2022$subtype == "Arterial Road" &
              !is.na(l_dist_2022$subtype),]$width_class <- "HIGH"
try(l_dist[is.na(l_dist$width_class) & 
             l_dist$subtype == "Arterial Road" &
             !is.na(l_dist$subtype),]$width_class <- "HIGH")

# filter(l_dist_tmp, subtype == "Access Road") %>% group_by(width_class) %>% summarise(n = n())
try(l_dist[is.na(l_dist$width_class) & 
             l_dist$subtype == "Access Road" &
             !is.na(l_dist$subtype),]$width_class <- "MED")
l_dist_2022[is.na(l_dist_2022$width_class) & 
              l_dist_2022$subtype == "Access Road" &
              !is.na(l_dist_2022$subtype),]$width_class <- "MED"

# filter(l_dist_tmp, subtype == "Unpaved Road") %>% group_by(width_class) %>% summarise(n = n())
l_dist[is.na(l_dist$width_class) & 
         l_dist$subtype == "Unpaved Road" &
         !is.na(l_dist$subtype),]$width_class <- "MED"
l_dist_2022[is.na(l_dist_2022$width_class) & 
              l_dist_2022$subtype == "Unpaved Road" &
              !is.na(l_dist_2022$subtype),]$width_class <- "MED"

# filter(as.data.frame(l_dist_tmp), str_detect(subtype, "Survey")) %>% group_by(width_class) %>% summarise(n = n()) 
l_dist_2022[is.na(l_dist_2022$width_class) & 
              str_detect(l_dist_2022$subtype, "Survey") &
              !is.na(l_dist_2022$subtype),]$width_class <- "MED"
try(l_dist[is.na(l_dist$width_class) & 
             str_detect(l_dist$subtype, "Survey") &
             !is.na(l_dist$subtype),]$width_class <- "MED")


filter(l_dist_tmp, subtype == "Trail") %>% group_by(width_class) %>% summarise(n = n())
try(l_dist[is.na(l_dist$width_class) & 
             l_dist$subtype == "Trail" &
             !is.na(l_dist$subtype),]$width_class <- "LOW")
try(l_dist_2022[is.na(l_dist_2022$width_class) & 
                  l_dist_2022$subtype == "Trail" &
                  !is.na(l_dist_2022$subtype),]$width_class <- "LOW")

##Most HS unknown features are medium width class
table(filter(l_dist_2022, is.na(width_class)) %>% pull(subtype))
table(filter(l_dist, is.na(width_class)) %>% pull(subtype))


filter(l_dist_tmp, subtype == "Unknown") %>% group_by(width_class) %>% summarise(n = n()) 
filter(l_dist_tmp, is.na(subtype)) %>% group_by(width_class) %>% summarise(n = n())

l_dist[is.na(l_dist$width_class),]$width_class <- "MED"
l_dist_2022[is.na(l_dist_2022$width_class),]$width_class <- "MED"

l_dist_2022$width_class <- factor(l_dist_2022$width_class, levels = c("LOW", "MED", "HIGH"), ordered = T)
l_dist$width_class <- factor(l_dist$width_class, levels = c("LOW", "MED", "HIGH"), ordered = T)

##Select and compare linear distance of YG vs Pat's digitization 

## summarise disturbance within 1000m buffer of each station
l.dist.yg <- st_intersection(l_dist_2022, 
                             sites %>% st_buffer(1000) %>% 
                               group_by(siteID, stationID) %>% summarise()) ## note: includes 2023 WCS digitization ontop of the May 2022 YG footprint data

l.dist.wcs <- st_intersection(l_dist %>% select(-siteID), 
                              sites %>% st_buffer(1000) %>% 
                                group_by(siteID, stationID) %>% summarise()) ## WCSC digitization for sites surveyed in 2017 - 2021, on an older version of the YG layer

##find and compare total length by site
#reduce number of features
l.dist.yg <- l.dist.yg %>% group_by(siteID, stationID, width_class, subtype, width_m) %>% summarise()
l.dist.wcs <- l.dist.wcs %>% group_by(siteID, stationID, width_class, subtype, width_m) %>% summarise()

l.dist.yg$length_yg <- st_length(l.dist.yg)
l.dist.wcs$length_wcs <- st_length(l.dist.wcs)

## calculate total length of all linear features per site from the May 2022 yg layer vs the 2021 wcsc digitizations
l.dist.comp <- full_join(st_drop_geometry(l.dist.yg), 
                         st_drop_geometry(l.dist.wcs)) %>%
  group_by(siteID) %>% summarise(length_yg = sum(length_yg, na.rm = T),
                                 length_wcs = sum(length_wcs, na.rm = T)) %>%
  rowwise() %>% mutate(length = max(c(length_yg, length_wcs), na.rm = T))

## which layer has more disturbance data?
l.dist.comp <- l.dist.comp %>% mutate(layer = ifelse(length == length_yg, "yg", "wcs")) ## layer is the data source with the most linear features for that site, presumably the most up to date

l.dist <- rbind(l.dist.yg %>% rename(length_m = length_yg) %>% mutate(layer = "yg"),
                l.dist.wcs %>% rename(length_m = length_wcs) %>% mutate(layer = "wcs")) %>% 
  right_join(l.dist.comp %>% select(siteID, layer)) ## select rows from right layer

##l.dist: clipped to 1500m buffer, each site has own multilinestring per width_class. 

### reduce subtypes 
local.rd <- c("Access Road", "Transportation - Access Road", "Local Road", "Driveway", "Unpaved Road", "Rural - Driveway", "Transportation - Local Road", "Transportation - Access Assumed", "Local / Street", "Access Assumed", "Service Lane") 
row <- c("Right of Way", "Unknown NRN", "Laneway", "Transportation or Unknown - Right of Way")
trail <- c("Trail", "Transportation - Trail", "Resource / Recreation")
cut <- c("Mining or Unknown - Survey / Cutline", "Survey / Cutline - Placer", "Survey / Cutline - Quartz", "Survey / Cutline", "Survey - Cutline" )
trench <- c("Trench", "Mining - Diversion Channel", "Mining - Trench")
euc <- c("Electric Utility Corridor", "Utility - Electric Utility Corridor")
highway <- c("Highway", "Expressway / Highway")
unk <- c("Unknown", "Mining, Unknown or Utility - Unknown")
airstrip <- c("Transportation - Airstrip")
fuelbreak <- c("Forestry - Fuel Break")

l.dist <- l.dist %>% mutate(subtype2 = ifelse(subtype %in% local.rd, "local road",
                                              ifelse(subtype %in% row, "right of way",
                                                     ifelse(subtype %in% trail, "trail",
                                                            ifelse(subtype %in% cut, "cutline",
                                                                   ifelse(subtype %in% trench, "trench",
                                                                          ifelse(subtype %in% euc, "electric utility corridor",
                                                                                 ifelse(subtype %in% highway, "highway", ifelse(subtype %in% airstrip, "airstrip", ifelse(subtype %in% fuelbreak, "fuel break", "unknown"))))))))))

filter(l.dist, subtype2 == "unknown") %>% pull(subtype) %>% unique()

l.dist.width <- st_drop_geometry(l.dist) %>% group_by(subtype2) %>% 
  summarise(width_mn = mean(width_m, na.rm = T), width_med = median(width_m, na.rm = T))
# what about widths ranges for unknowns? Base this on width class
st_drop_geometry(l.dist) %>% group_by(width_class) %>% summarise(min = min(width_m, na.rm = T), median = median(width_m, na.rm = T), max = max(width_m, na.rm = T))

##For features with unknown widths, use the median value for that subtype.
#*note: might be better to use width_med for all analysis, as width_m is ofter unaccurate and variable throughout the feature, and the imagery resolution and thus measurement resolution is low. 
l.dist <- left_join(l.dist, l.dist.width) %>% 
  mutate(width_m = ifelse(!(is.na(width_m)), width_m, ## if there is a width, keep it
                          ifelse(subtype2 != "unknown", width_med, #if it is a known feature subtype, use the median
                                 ifelse(width_class == "LOW", 3, ifelse(width_class == "MED", 6, 11))))) ## if unknown, use median for the width class

table(l.dist$subtype2, l.dist$width_class) 
table(l.dist$subtype2, l.dist$width_med) 

l.dist$type <- ifelse(l.dist$subtype2 %in% c("highway", "local road", "right of way"), "Road",
                      ifelse(l.dist$width_class == "LOW", "Trail", "Cutline")) ## trail and cutline are generalized terms based on subtype2 median widths, but represent a range of subtypes. Non-road narrow and non-road wide are more accurate terms, but its confusing to change now
table(l.dist$type)
## Linear to polygon 
## if we buffer linear feature by their widths, this can be removed from disturbance surface area which reduce correlation between surface disturbance and linear density. 
## Note: recommended to use the median width vs 
l.dist.p <- st_buffer(l.dist, dist = l.dist$width_med/2, endCapStyle="FLAT") ##width/2 because adding buffer to both sides
l.dist.p <- l.dist.p %>% group_by(siteID, stationID, type) %>% summarise()
## remove overlapping areas, where road > cutline > trail

rm(l_dist, l_dist_2022, l_dist_2023, l_dist_a_hs, l_dist_a_ns, l_dist_hs, l_dist_hs2, l_dist_ns, l_nrn_hs, l_nrn_hs2, l_nrn_ns, l.dist.comp, l.dist.wcs, l.dist.yg)


## remove linear disturbance area from habitat poly

names(l.dist.p)
names(hab.poly.station)

l.dist.p.l <- l.dist.p %>% split(.$stationID)
hab.poly.l <- hab.poly.station %>% split(.$stationID)

all.equal(names(l.dist.p.l), names(hab.poly.l))
##not all sites have linear features? 

hab.poly.l <- hab.poly.l[names(l.dist.p.l)]

## remove linear from current habitat poly (includes surface disturbance)
hab.poly.l <- pmap(list(l.dist.p.l, hab.poly.l), function(disturb, hab.poly){
  st_difference(hab.poly, st_union(disturb))
})

## Merge the two layers into a single layer
hab.poly.l <- pmap(list(l.dist.p.l, hab.poly.l), function(disturb, hp){
  hp <- select(hp, siteID, stationID, class, layer)
  disturb <- disturb %>% 
    mutate(class = paste0("d_", type), layer = "YGfootprint_CWS_update") %>%
    select(names(hp))
  bind_rows(list(disturb, hp))
})

ggplot(hab.poly.l$`MS12-2`, aes(fill = class)) + geom_sf(alpha = .5)

## combine with sites without linear features that were removed from hab.poly.l
hab.poly.station <- bind_rows(bind_rows(hab.poly.l), 
                              hab.poly.station %>% filter(!stationID %in% names(hab.poly.l)))

saveRDS(hab.poly.station, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "lcmerge_stat_linear.RDS"))
## Waterways: Stream Density----
fgdb <- file.path(gis_dir, "canvec_50k_YT_Hydro.gdb")
waterway<-readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/waterway.rds")

##Stream order 

stream_order <- st_read(file.path(gis_dir, "StreamOrder/StreamOrder_CDEMFill1000_shp.shp")) ## created in arcGIS using CDEM elevation raster to calculate flow accumulation and then stream order. 
stream_order <- stream_order %>% rename(order = grid_code) %>% 
  group_by(order) %>% summarize()

stream_order <- st_intersection(stream_order, sites %>% st_buffer(1500) %>% 
                                  st_union() %>% 
                                  st_transform(crs = st_crs(stream_order)))

stream_order <- st_intersection(stream_order, 
                                sites %>% st_buffer(1500) %>% 
                                  group_by(siteID, stationID) %>% summarise() %>% 
                                  st_transform(crs = st_crs(stream_order)))


#### Habitat Summary----
## Decide whether we are using the habitat data with or without polygonized linear features. 
hab.poly.station <- readRDS(file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "lcmerge_stat.RDS")) ## without linear polygons
# hab.poly.station <- readRDS(file.path(pipeline, "store", "lcmerge_stat_linear.RDS")) ## with linear polygons
wetland.station <- readRDS(file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "store", "wetland_station.RDS")) 

## site ID variables
site.sum <- sites %>% st_drop_geometry() %>% 
  select(siteID, organization , year, wetland.layer) %>% 
  distinct() %>% 
  group_by(siteID, organization, year) %>% 
  summarise(wetland.layer = sum(!wetland.layer) == 0) %>% ## if one station is missiing wetland layer, this site will be false 
  ungroup()

# ## surface area
## area per polygon
hab.poly.station$area <- st_area(hab.poly.station)
wetland.station$area <- st_area(wetland.station)

## calculate %cover for each station using a range of buffers. 
hab.poly.site.l <- hab.poly.station %>% split(.$stationID)
wetland.site.l <- wetland.station %>% split(.$stationID)
sites.l <- split(sites, f = sites$stationID)[names(hab.poly.site.l)]

all.equal(names(sites.l), names(hab.poly.site.l))


buff.v <- c(150, 500, 1000)  ## can't go beyon 1000m because 2021 data digitization stoped at 1.5k from site centroid (~1k from the corner sites)
# buff.v <- c(150, 1000)  ## smallest test

##creates a list where each element is a different buffer 
hab.sum.buff <- map(buff.v, function(buff) {
  ##for buffer i, crop habitat polygons around each site
  hab.poly.site.buff.i <- pmap(list(sites.l, hab.poly.site.l),
                               function(site, poly){
                                 site <- st_buffer(site, dist = buff) %>%
                                   st_union() %>% st_transform(crs = st_crs(poly))
                                 st_intersection(poly, site)
                               }) 
  #bind all stations together
  hab.poly.site.buff.i <- bind_rows(hab.poly.site.buff.i)
  hab.poly.site.buff.i$area <- hab.poly.site.buff.i %>% st_area()
  ## calculate the total area per class for that station&buffer
  hab.sum.buff.i <- st_drop_geometry(hab.poly.site.buff.i) %>% group_by(siteID, stationID, class) %>% summarise(area = sum(area))
  ## calculate % cover of each habitat class per station&buffer, join station data together as a df for that buffer
  hab.sum.buff.i <- hab.sum.buff.i %>% group_by(siteID, stationID) %>% 
    summarise(site.area = sum(area)) %>% left_join(hab.sum.buff.i) %>% 
    mutate(prop = as.numeric(area/site.area)) %>%
    select(-area) %>%
    pivot_wider(id_cols = c(siteID, stationID, site.area), names_from = class, 
                values_from = prop, values_fill = 0) %>% 
    mutate(buffer = buff)
  hab.sum.buff.i
})

#join different buffers, and add to site ID info
hab.sum <- right_join(site.sum, list_rbind(hab.sum.buff))
hab.sum$d_surface <- rowSums(hab.sum %>% select(d_pBare, d_pVeg, d_pWater))

## if using hab polygons
#hab.sum$d_linear <- rowSums(hab.sum %>% select(d_Road, d_Trail, d_Cutline))
#hab.sum$d_all <- rowSums(hab.sum %>% select(d_Road, d_Trail, d_Cutline,d_pBare, d_pVeg, d_pWater))

all.equal(names(sites.l), names(wetland.site.l)) ## not all sites have wetlands
sites2.l <- sites.l[names(wetland.site.l)]

### wetland area 
## repeat steps used for hab poly with the wetland polys

gc()


tictoc::tic()
hab.buff <- map(buff.v, function(buff) {
  hab.poly.site.buff.i <- pmap(list(sites2.l, wetland.site.l),
                               function(site, poly){
                                 site <- st_buffer(site, dist = buff) %>%
                                   st_union() %>% st_transform(crs = st_crs(poly))
                                 st_intersection(poly, site)
                               }) 
  hab.poly.site.buff.i <- bind_rows(hab.poly.site.buff.i)
  hab.poly.site.buff.i <- hab.poly.site.buff.i %>% group_by(siteID, stationID, class2) %>% summarise()
  
  hab.poly.site.buff.i$area <- hab.poly.site.buff.i %>% st_area()
  hab.sum.buff.i <- st_drop_geometry(hab.poly.site.buff.i) %>% group_by(siteID, stationID, class2) %>% 
    summarise(area = sum(area))
  
  
  hab.sum.buff.i <- left_join(hab.sum.buff.i, hab.sum %>% filter(buffer == buff) %>% 
                                select(stationID, site.area)) %>%
    mutate(prop = as.numeric(area/site.area)) %>%
    select(-area) %>%distinct()%>%
    pivot_wider(id_cols = c(siteID, stationID, site.area), names_from = class2, 
                values_from = prop, values_fill = 0) %>% 
    mutate(buffer = buff)
  hab.sum.buff.i
})
tictoc::toc()

hab.sum <- left_join(hab.sum, list_rbind(hab.buff))

### add elevation 
hab.sum <- left_join(hab.sum, elev_mn %>% rename(elev = mean))
hab.sum <- left_join(hab.sum, bz_sum %>% select(siteID, hc_borealLow = P.Boreal_Low, hc_borealSubalpine = P.Boreal_Subalpine))
print("done")

##add daymet weather
daymet<-readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/final_weather_summary.RDS")
hab.sum <- left_join(hab.sum, daymet%>%select(!year),by="siteID", relationship = "many-to-one") 

##add linear features
l.dist.l <- l.dist %>% split(.$stationID) ##Not all sites overlap linear features, so those sites will be missing from this data set

sites.l <- sites.l[names(l.dist.l)] 

hab.sum.lf.buff <- map(buff.v, function(buff) {
  l.dist.buff.i <- pmap(list(sites.l, l.dist.l),
                        function(site, poly){
                          site <- st_buffer(site, dist = buff) %>%
                            st_union() %>% st_transform(crs = st_crs(poly))
                          st_intersection(poly, site)
                        }) 
  l.dist.buff.i <- bind_rows(l.dist.buff.i)
  l.dist.buff.i$length <- l.dist.buff.i[st_geometry_type(l.dist.buff.i) %in% c("LINESTRING","MULTILINESTRING"),] %>% st_length()
  
  st_drop_geometry(l.dist.buff.i) %>% 
    group_by(siteID, stationID, type) %>%
    summarise(length = sum(length)) %>% 
    left_join(ungroup(hab.sum) %>% filter(buffer == buff) %>% select(siteID, stationID, site.area)) %>%
    mutate(density = as.numeric(length/site.area),
           type = substring(tolower(type), 1, 1)) %>%distinct()%>%
    pivot_wider(id_cols = c(siteID, stationID),
                names_from = type, 
                names_glue = "d_lden_{type}", 
                values_from = density,
                values_fill = 0) %>% mutate(buffer = buff) 
  
})

hab.sum.lf <- list_rbind(hab.sum.lf.buff)
table(hab.sum.lf$buffer) ## missing rows if there are no linear features that overlap 


#Waterway density
waterway.l <- waterway %>% split(.$stationID)
sites.l <- split(sites, f = sites$stationID)[names(waterway.l)]

valid_indices <- which(!sapply(sites.l, is.null) & !sapply(waterway.l, is.null))
sites.l <- sites.l[valid_indices]
waterway.l <- waterway.l[valid_indices]

water.den <- map(buff.v, function(buff) {
  
  waterway.buff.i <- pmap(list(sites.l, waterway.l),
                          function(site, poly) {
                            if (is.null(site) || is.null(poly)) return(NULL)
                            site <- st_buffer(site, dist = buff) %>%
                              st_union() %>%
                              st_transform(crs = st_crs(poly))
                            st_intersection(poly, site)
                          })
  
  # Remove NULL results from `pmap`
  waterway.buff.i <- waterway.buff.i[!sapply(waterway.buff.i, is.null)]
  
  waterway.buff.i <- bind_rows(waterway.buff.i)
  waterway.buff.i$length <- waterway.buff.i %>% st_length()
  
  st_drop_geometry(waterway.buff.i) %>% 
    group_by(siteID, stationID) %>%
    summarise(length = sum(length)) %>%
    left_join(filter(hab.sum, buffer == buff) %>% select(siteID, stationID, buffer, site.area)) %>%
    mutate(wden = as.numeric(length/site.area)) %>%
    select(siteID, stationID, buffer, wden) %>%
    mutate(buffer = buff)
})


water.den <- list_rbind(water.den)

#stream order
stream_order.l <- stream_order %>% split(.$stationID) 
sites.l <- split(sites, f = sites$stationID)[names(stream_order.l)]

stream_order.sum <- map(buff.v, function(buff) {
  stream_order.buff.i <- pmap(list(sites.l, stream_order.l),
                              function(site, poly){
                                site <- st_buffer(site, dist = buff) %>%
                                  st_union() %>% st_transform(crs = st_crs(poly))
                                st_intersection(poly, site)
                              }) 
  stream_order.buff.i <- bind_rows(stream_order.buff.i)
  stream_order.buff.i$length <- stream_order.buff.i[st_geometry_type(stream_order.buff.i) %in% 
                                                      c("LINESTRING","MULTILINESTRING"),] %>% st_length()
  max.order <- max(unique(stream_order.buff.i$order))
  stream_order.sum.i <- st_drop_geometry(stream_order.buff.i) %>% 
    group_by(siteID, stationID, order) %>%
    summarise(length = sum(length)) %>% 
    left_join(filter(hab.sum, buffer == buff) %>% select(siteID, stationID, buffer, site.area)) %>%
    mutate(density = as.numeric(length/site.area)) %>% distinct()%>%
    pivot_wider(id_cols = c(siteID, stationID, buffer), 
                names_from = order, 
                names_glue = "wden_{order}", 
                values_from = density,
                values_fill = 0) 
  ## the highest order stream within 1000m is 8. Smaller buffers might miss these higher order streams. 
  
  if(max.order < 8){
    
    x <- paste0("wden_", seq(max.order + 1, 8, 1))
    
    for(i in x) {
      stream_order.sum.i <- mutate(stream_order.sum.i, "{i}" := 0)
    }
    
  }
  
  
  stream_order.sum.i %>% mutate(buffer = buff) 
})

stream_order.sum <- list_rbind(stream_order.sum)
table(stream_order.sum$buffer)

water.den <- full_join(water.den, stream_order.sum)

hab.sum <- full_join(hab.sum, water.den) %>% full_join(hab.sum.lf) %>% full_join(hab.sum.lf)

table(hab.sum$buffer) ## all equal now - need to fill missing values with 0
table(hab.sum$buffer, is.na(hab.sum$wden_5))


#### Save Habitat Summary Station-level #####
names(hab.sum)
hab.sum<-hab.sum%>%distinct()
hab.sum$site.area <- as.numeric(hab.sum$site.area)
hab.sum <- hab.sum %>% rename_with(\(x) str_replace(x, pattern = "H_", replacement = "h_"), starts_with('H_'))
hab.sum$wden_high <- rowSums(hab.sum %>% select(wden_6, wden_7, wden_8))
hab.sum <- hab.sum %>% select(-wden_6, -wden_7, -wden_8, - site.area) ## low rep of high order streams, so group into a single class
names(hab.sum)
hab.sum<-hab.sum%>%replace(is.na(.),0)

saveRDS(hab.sum, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "out", "hab_sum.RDS"))


#### Save Habitat Summary Site-level ####
visits <- readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/1_CountDataProcessing/out/visits.RDS") ## time variant survey parameters

visits$stationID <- paste(visits$siteID, visits$station, sep = "-")

hab.sum2 <- left_join(hab.sum, st_drop_geometry(visits), relationship = "many-to-many")
## many to many bc multiple buffers in hab.sum and multiple visits in visits

## Proportional to the number of visits, not the actual area. 
hab.sum2 <- hab.sum2 %>% group_by(organization, siteID, wetland.layer, year, buffer) %>%
  summarise(across(c(starts_with("h_"), starts_with("d_"), elev, starts_with("mn_"), 
                     starts_with("wden")), mean)) 
#add site centroid
sitesC <- sites %>% group_by(siteID) %>% summarise() %>%
  st_centroid() %>% select(siteID)
sitesC$lon <- st_coordinates(sitesC)[,1]
sitesC$lat <- st_coordinates(sitesC)[,2]
sitesC$EPSG <- "EPSG:3579"

sitesC$lon_wgs84 <- st_coordinates(st_transform(sitesC, crs = "EPSG:4326"))[,1]
sitesC$lat_wgs84 <- st_coordinates(st_transform(sitesC, crs = "EPSG:4326"))[,2]


hab.sum2 <- left_join(hab.sum2, st_drop_geometry(sitesC))

saveRDS(hab.sum2, file.path("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds",pipeline, "out", "hab_sum_site.RDS"))

