library(tidyverse)
library(KrigR)
library(sf)
API_User <- 31235
API_Key <- "a487a150-29b4-4aac-ac7b-c4148e4fb1c3"

install.load.package <- function(x){
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c("tidyr", # for turning rasters into ggplot-dataframes
                 "ggplot2", # for plotting
                 "viridis", # colour palettes
                 "cowplot", # gridding multiple plots
                 "ggmap", # obtaining satellite maps
                 "gimms", # to get some pre-existing data to match in our downscaling
                 "rnaturalearth", # for shapefiles
                 "rnaturalearthdata", # for high-resolution shapefiles
                 "mapview" # for generating mapview outputs
)
sapply(package_vec, install.load.package)

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "0_data", "era5") # folder path for data
Dir.Covariates <- file.path(Dir.Base, "2_pipeline", "2_HabitatSum", "tmp", "era5Covariates") # folder path for covariates
Dir.Exports <- file.path(Dir.Base, "2_pipeline", "2_HabitatSum", "tmp", "era5Exports") # folder path for exports
## create directories, if they don't exist yet
Dirs <- sapply(c(Dir.Data, Dir.Covariates, Dir.Exports), 
               function(x) if(!dir.exists(x)) dir.create(x))

source("1_code/FUN_Plotting.R")

sites <- read_sf("CumulativeEffects_GIS/Sites/sites_analysis.shp") 
# loc <- sites %>% group_by(siteID, year) %>% summarise() %>%
#   st_centroid() %>% st_transform(crs = "EPSG:4326")
loc <- sites %>% st_transform(crs = "EPSG:4326")
loc$lon <- st_coordinates(loc)[,1]
loc$lat <- st_coordinates(loc)[,2]
loc$EPSG <- "EPSG:4326"

ggplot(loc) + geom_sf()
ggplot(loc) + geom_histogram(aes(lon, group = year, fill = factor(year)), bins = 30)
ggplot(loc) + geom_histogram(aes(lat, group = year, fill = factor(year)), bins = 10)

ggplot(loc %>% filter(lon > -137.5 & lon < -136.5)) + 
  geom_histogram(aes(lon), bins = 30)

## start with mayo region in 2023 to limit spatial coverage
##

##automate by year

yr <- c("2017, 2018, 2021, 2023")


trial.loc <- filter(loc, year == 2023) 
ext <- as(trial.loc %>% st_buffer(150*100), 'Spatial') ## buffer about the size of a weather system (150km). THis provides more training data vs using just small buffers
extKrig <- as(trial.loc %>% st_buffer(300), 'Spatial')
# Extent_ext <- st_bbox(trial.loc)
# Shape_shp <- st_buffer(trial.loc, 1500)
# bbox <- Extent_ext
# names(bbox) <- c("left", "bottom", "right", "top")
# loc_df <- st_drop_geometry(trial.loc) %>% select()

All23 <- download_ERA(
  Variable = "2m_temperature", # the variable we want to obtain data for
  DataSet = "era5-land", # the data set we want to obtain data from
  DateStart = "2023-04-01", # the starting date of our time-window
  DateStop = "2023-06-30", # the final date of our time-window
  TResolution = "month", ## monthly average
  TStep = 3, ## over 3 months
  Extent = ext, # the spatial preference we are after
  Dir = Dir.Data, # where to store the downloaded data
  FileName = "All23_spring_2mt", # a name for our downloaded file
  API_User = API_User, # your API User Number
  API_Key = API_Key # your API User Key
)

Plot_Raw(All23, Dates = "Spring 2023",  Shp = ext)

##downscaling

Covs_ls <- download_DEM(Train_ras = All23, # the data we want to downscale
                        Target_res = .02, # the resolution we want to downscale to
                        Shape = ext, # extra spatial preferences
                        Dir = Dir.Covariates # where to store the covariate files
)
Plot_Covs(Covs_ls, ext)

## interpolate
All23Krig <- krigR(Data = All23, # data we want to krig as a raster object
                             Covariates_coarse = Covs_ls[[1]], # training covariate as a raster object
                             Covariates_fine = Covs_ls[[2]], # target covariate as a raster object
                             Keep_Temporary = FALSE, # we don't want to retain the individually kriged layers on our hard-drive
                             Cores = 1, # we want to krig on just one core
                             FileName = "All23Krig", # the file name for our full kriging output
                             Dir = Dir.Exports # which directory to save our final input in
)

Plot_Krigs(All23Krig, 
           Shp = ext,
           Dates = c("Spring 2023")
)

##local
AllLoc23Krig <- krigR(Data = All23, # data we want to krig as a raster object
                   Covariates_coarse = Covs_ls[[1]], # training covariate as a raster object
                   Covariates_fine = Covs_ls[[2]], # target covariate as a raster object
                   Keep_Temporary = FALSE, # we don't want to retain the individually kriged layers on our hard-drive
                   Cores = 1, # we want to krig on just one core
                   nmax = 30, ## only use local cells for krig
                   FileName = "AllLoc23Krig", # the file name for our full kriging output
                   Dir = Dir.Exports # which directory to save our final input in
)

Plot_Krigs(AllLoc23Krig, 
           Shp = ext,
           Dates = c("Spring 2023")
) ## more uncertainty (high SD)


trial.loc <- filter(loc, year == 2021) 
ext <- as(trial.loc %>% st_buffer(150*100), 'Spatial') ## buffer about the size of a weather system (150km). THis provides more training data vs using just small buffers
extKrig <- as(trial.loc %>% st_buffer(300), 'Spatial')
# Extent_ext <- st_bbox(trial.loc)
# Shape_shp <- st_buffer(trial.loc, 1500)
# bbox <- Extent_ext
# names(bbox) <- c("left", "bottom", "right", "top")
# loc_df <- st_drop_geometry(trial.loc) %>% select()

All21 <- download_ERA(
  Variable = "2m_temperature", # the variable we want to obtain data for
  DataSet = "era5-land", # the data set we want to obtain data from
  DateStart = "2021-04-01", # the starting date of our time-window
  DateStop = "2021-06-30", # the final date of our time-window
  TResolution = "month", ## monthly average
  TStep = 3, ## over 3 months
  Extent = ext, # the spatial preference we are after
  Dir = Dir.Data, # where to store the downloaded data
  FileName = "All21_spring_2mt", # a name for our downloaded file
  API_User = API_User, # your API User Number
  API_Key = API_Key # your API User Key
)

Plot_Raw(All21, Dates = "Spring 2021",  Shp = ext)

##downscaling

Covs_ls <- download_DEM(Train_ras = All21, # the data we want to downscale
                        Target_res = .02, # the resolution we want to downscale to
                        Shape = ext, # extra spatial preferences
                        Dir = Dir.Covariates # where to store the covariate files
)
Plot_Covs(Covs_ls, ext)

## interpolate
All21Krig <- krigR(Data = All21, # data we want to krig as a raster object
                   Covariates_coarse = Covs_ls[[1]], # training covariate as a raster object
                   Covariates_fine = Covs_ls[[2]], # target covariate as a raster object
                   Keep_Temporary = FALSE, # we don't want to retain the individually kriged layers on our hard-drive
                   Cores = 1, # we want to krig on just one core
                   FileName = "All21Krig", # the file name for our full kriging output
                   Dir = Dir.Exports # which directory to save our final input in
)

Plot_Krigs(All21Krig, 
           Shp = ext,
           Dates = c("Spring 2021")
)

##local
AllLoc21Krig <- krigR(Data = All21, # data we want to krig as a raster object
                      Covariates_coarse = Covs_ls[[1]], # training covariate as a raster object
                      Covariates_fine = Covs_ls[[2]], # target covariate as a raster object
                      Keep_Temporary = FALSE, # we don't want to retain the individually kriged layers on our hard-drive
                      Cores = 1, # we want to krig on just one core
                      nmax = 30, ## only use local cells for krig
                      FileName = "AllLoc21Krig", # the file name for our full kriging output
                      Dir = Dir.Exports # which directory to save our final input in
)

Plot_Krigs(AllLoc21Krig, 
           Shp = ext,
           Dates = c("Spring 2021")
) ## more uncertainty (high SD)


