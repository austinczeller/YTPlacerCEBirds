#new count data processing incorporates the 2024 data
library(tidyverse)
library(lubridate)
library(sf) ##spatial features
library(readxl)
library(biscale)
#library(rgdal) ##r geospatial data abstraction library
# library(raster)
library(terra)
# library(stars)
library(cowplot)
library(ggbreak)
library(wildrtrax)
select <- dplyr::select
theme_set(theme_classic())

# ##Shannon index - all species are represented in a sample, and they are randomly sampled
shannon.index <- function(n) {
  ##n is a vector of abundance per species at a site
  N = sum(n)
  p.i <- n/N
  -sum(p.i * log(p.i))
}

Sys.setenv(WT_USERNAME = 'austinzeller@rocketmail.com', WT_PASSWORD = 'Ke$ha50az') ## enter your own username and password
wt_auth()

my_projects <- wt_get_download_summary(sensor_id = 'ARU')

project_ids <- c(2782)
proj_names <- filter(my_projects, project_id %in% project_ids) %>% pull(project)

aru <- map(project_ids, \(id) wt_download_report(project_id = id, sensor_id = "ARU", report = c("main"), weather_cols = F))
write_rds(aru, file.path("E:/Archive/2_pipeline/1_CountDataProcessing", "tmp", "2024CE_wt.RDS"))

aru <- map(aru, \(df) filter(df, aru_task_status == "Transcribed" & observer != "Patrice Mathieu"))
aru<-aru[[1]]



aru_noise <- read_csv("E:/Archive/0_data/South Beringia Priority Place Cumulative Effects 2024_Tasks_202439.csv")
aru_noise <- aru_noise %>% filter(status %in% c("Transcribed")) %>%
  mutate(recordingDate = as.character(recordingDate), 
         task_duration = paste0(taskLength, "s"),
         industryNoise = factor(industryNoise, c("Unk", "None","Light","Moderate","Heavy"), ordered = T),
         otherNoise = factor(otherNoise, c("Unk","None","Light","Moderate","Heavy"), ordered = T)) %>%
  rowwise %>% mutate(obsNoise = max(industryNoise, otherNoise, na.rm = T)) %>%
  select(location, task_id = internal_task_id, recording_date_time = recordingDate, observer = transcriber, task_method = method, obsNoise) %>%
  ungroup()%>%mutate(recording_date_time=lubridate::ymd_hms(recording_date_time))

aru <- aru%>% left_join(aru_noise)

aru <- aru %>% separate(location, into=c("siteID","station"),sep="-", remove = F) 

aru <- aru %>% select(org = organization, locID = location_id, siteID, station, lat = latitude, lon =  longitude,
                      recID = recording_id, taskID = task_id, recording_date = recording_date_time, 
                      dur = task_duration, task_comments, obs = observer_id, 
                      species_code, species_common_name, tagID = tag_id, 
                      count = individual_count, voc = vocalization, detection_time, obsNoise) %>%
  mutate(recording_date = ymd_hms(recording_date, tz = "Canada/Yukon"), 
         dur = as.numeric(substr(dur,1,nchar(dur)-1)))


tmp <- aru %>% select(siteID, station, recording_date, taskID, obs) %>% distinct() %>% 
  select(-obs, -taskID) %>%
  duplicated() ## is row duplicated?
tmp2 <- aru %>% select(siteID, station, recording_date, taskID, obs) %>% distinct() ## need a matching dataset to filter
table(tmp2[tmp,]$obs) ## approximately random which observer is dropped
aru <- anti_join(aru, tmp2[tmp,]) 


record.sched <- aru %>% select(siteID, station, recording_date, dur) %>% distinct() 
species <- wt_get_species() ## download key to species codes
aru <- anti_join(aru, filter(species, species_class != "AVES" & !(species_common_name %in% c("NONE", "Unidentified"))))


aru %>% filter(count == "TMTT") %>% select(species_code, species_common_name) %>% distinct()
table(aru %>% group_by(locID, species_code) %>% summarise(n.of.ind = n()) %>% filter(species_code %in% c("WWCR", "BANS", "CANG", "VGSW")) %>% pull(species_code),
      aru %>% group_by(locID, species_code) %>% summarise(n.of.ind = n()) %>% filter(species_code %in% c("WWCR", "BANS", "CANG", "VGSW")) %>% pull(n.of.ind))
aru <- aru %>% mutate(count = as.numeric(ifelse(count == "TMTT", 5, count)))


aru.sites <- unique(aru$siteID)
aru$ts <- aru$recording_date
aru$voc <- factor(aru$voc) #call, song, non-vocal

aru <- aru %>%
  select(organization = org, siteID = siteID, station, lat_wt = lat, lon_wt = lon,
         ts, species_code, species_english_name = species_common_name, 
         species_individual_name = tagID, detection_time, abundance = count, dur, observer = obs, obsNoise)

##ceiling to match with PC minute intervals if trying to pair.
##Note recommended, because the PC data with paired recordings isn't very good....
aru$tag_min <- ceiling(aru$detection_time/60)
aru[aru$tag_min == 0 & !is.na(aru$tag_min), "tag_min"] <- 1 ## these tags started at time 0:00, so should be in first minute.

##wildtrax updated GRAJ to CAJA before tag validation, so now mixed
aru[aru$species_code == "GRAJ",]$species_english_name <- "Canada Jay"
aru[aru$species_code == "GRAJ",]$species_code <- "CAJA"

table(aru$station)
aru <- aru %>% mutate(station = as.numeric(gsub("\\D", "", station)))


counts <- aru %>% select(siteID, station, ts, species_code,  species_english_name, individual = species_individual_name, tag_min, abundance) %>% 
  mutate(method = "ARU", detection_type = "A", distance_interval = NA, station = as.numeric(station))


visits <- aru %>% select(siteID, station, ts , dur, obsNoise, observer) %>% distinct() %>% 
  mutate(dur_min = dur/60,  method = "ARU", observer = as.character(observer))

##Unidentified codes
try(counts[counts$species_code == "WARB",]$species_code <- "UNWA") ## Unk warbler
try(counts[counts$species_code == "WOOD",]$species_code <- "UNWO") ## unidentified woodpecker based on comments in eccc
try(counts[counts$species_code == "SPAR",]$species_code <- "UNSP") ## unidentified sparrow
try(counts[counts$species_code == "OWL ",]$species_code <- "UNOW") ## unknown owl 
try(counts[counts$species_code == "OWL",]$species_code <- "UNOW") ## unknown owl 
try(counts[counts$species_code == "LONG",]$species_code <- "UNLO") ##LALO or SMLO
try(counts[counts$species_code == "UTBB",]$species_code <- "UNBT") #Unidentified Black-backed / American Three-toed Woodpecker
try(counts[counts$species_code == "MYRU",]$species_code <- "UNKN") ## no idea Myrtle warbler /yellow rumped warbler?
try(counts[counts$species_code == "RAPT",]$species_code <- "UNRA") ## unknown raptor
try(counts[counts$species_code == "NOOO",]$species_code <- "NONE") ##based on comments
try(counts[counts$species_code == "SWAL",]$species_code <- "UNSW")
try(counts[counts$species_code == "SFAL",]$species_code <- "UNFA") #kestral or merlin
try(counts[counts$species_code == "THRU",]$species_code <- "UNTH")
try(counts[counts$species_code == "PASS",]$species_code <- "UNPA")
try(counts[counts$species_code == "FINC",]$species_code <- "UNFI")
try(counts[counts$species_code == "UNK 'PIK'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK WOODPECKER",]$species_code <- "UNWO")
try(counts[counts$species_code == "UNK 'HAY'",]$species_code <- "UNKN")
try(counts[counts$species_code == "CAJA?",]$species_code <- "CAJA")
try(counts[counts$species_code == "NOWA/GULL?",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'CHIP'",]$species_code <- "UNPA")
try(counts[counts$species_code == "UNK GULL",]$species_code <- "UNGU")
try(counts[counts$species_code == "CORE/PISI?",]$species_code <- "UNFI")
try(counts[counts$species_code == "NOFL?",]$species_code <- "UNWO")
try(counts[counts$species_code == "UNK 'WEEEEE'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'TA TA TA DWI DA DI'",]$species_code <- "UNKN")
try(counts[counts$species_code == "HETH/VATH?",]$species_code <- "UNTH")
try(counts[counts$species_code == "AMRO/NOFL?",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'woo woo waa'" ,]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'WHAAAAA'",]$species_code <- "UNKN")
try(counts[counts$species_code == "CHSP?",]$species_code <- "UNSP")
try(counts[counts$species_code == "YRWA?",]$species_code <- "UNWA")
try(counts[counts$species_code == "UNK 'CHIK CHICK CHICK'",]$species_code <- "UNKN")
try(counts[counts$species_code == "DEJU?" ,]$species_code <- "UNSP")
try(counts[counts$species_code == "UNK DUCK FLIGHT" ,]$species_code <- "UNDU")
try(counts[counts$species_code == "UNK TRILL",]$species_code <- "UNTR")
try(counts[counts$species_code == "UNK 'TI-WHI TI-WHI'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK SPARROW",]$species_code <- "UNSP")
try(counts[counts$species_code == "UNK SHOREBIRD" ,]$species_code <- "UNSH")
try(counts[counts$species_code == "SSHO" ,]$species_code <- "UNSH") ##unidentified small shorebird
try(counts[counts$species_code == "UNK 'WIP WIP WIP'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK RAPTOR",]$species_code <- "UNRA")
try(counts[counts$species_code == "UNK 'I WII I WI IWI IWI'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK SWALLOW",]$species_code <- "UNSW")
try(counts[counts$species_code == "UNK 'PIYO-PIYO-PIPI-PEW'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'HEIIIII'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'TI-OUI  TI-OUI TI-OUI'",]$species_code <- "UNKN")
try(counts[counts$species_code == "RUBL?",]$species_code <- "UNBL")
try(counts[counts$species_code == "UNK 'WEEPA WOO'",]$species_code <- "UNKN")
try(counts[counts$species_code == "CORE?",]$species_code <- "UNFI")
try(counts[counts$species_code == "UNK SPARROW 'TI TI THRIIIIIII IP'",]$species_code <- "UNSP")
try(counts[counts$species_code == "UNK FLYOVER \"CHIT CHIT\"",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'PEER PEER'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'TSSSSST'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'HE-HE HUP HUP HUP'" ,]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'CRT CRT BEET BEET BEET BEET'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK OWL 'WOOOOOOOFFFF'",]$species_code <- "UNOW")
try(counts[counts$species_code == "UNK 'QRI WIIIIII'"  ,]$species_code <- "UNKN")
try(counts[counts$species_code == "PIGR?",]$species_code <- "UNFI")
try(counts[counts$species_code == "UNK 'CHEEEE'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK FLYOVER \"BEER\"",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'ZRRRIIIII'",]$species_code <- "UNKN")
try(counts[counts$species_code == "YRWA/NOWA?",]$species_code <- "UNWA")
try(counts[counts$species_code == "COLO?",]$species_code <- "ULOO")
try(counts[counts$species_code == "TRSW",]$species_code <- "TRES")
try(counts[counts$species_code == "HAFL?",]$species_code <- "UNFL")
try(counts[counts$species_code == "BOCH/CORE?",]$species_code <- "UNPA")
try(counts[counts$species_code == "UNK 'CHIT CHI CHIT AWITA CHU'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK WARBLER 'SIP SIP'",]$species_code <- "UNWA")
try(counts[counts$species_code == "UNKNOWN",]$species_code <- "UNKN")
try(counts[counts$species_code == "GREBE",]$species_code <- "UGRE")
try(counts[counts$species_code == "RUGR/SPGR",]$species_code <- "UGRS")
try(counts[counts$species_code == "UNKNOWN ",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNKNOWN FLYOVER \"CHIT CHIT\"",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'OH WO-OHWO-OHWO'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'CRRT SIP SIP SIP'",]$species_code <- "UNKN")
try(counts[counts$species_code == "AMRO?",]$species_code <- "UNAM")
try(counts[counts$species_code == "RCKI?",]$species_code <- "RCKI")
try(counts[counts$species_code == "UNK 'HEYME HEYME YIP'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'WIT RRR WIT RRR'",]$species_code <- "UNKN")
try(counts[counts$species_code == "WIWA?",]$species_code <- "UNWA")
try(counts[counts$species_code == "SWTH?",]$species_code <- "UNTH")
try(counts[counts$species_code == "UNK 'WIT WIT'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'TIWICK' LOUD",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'WUAA WUAA'",]$species_code <- "UNKN")
try(counts[counts$species_code == "UNK 'SEEP!'",]$species_code <- "UNKN")

alphaCodes <- read_csv("E:/Archive/0_data/species_codes.csv")
counts <- counts %>% left_join(alphaCodes %>% select(species_code = SPEC, scientific_name = SCINAME, COMMONNAME)) %>%
  mutate(species_english_name = ifelse(is.na(species_english_name), COMMONNAME, species_english_name)) %>% 
  select(-COMMONNAME)
counts[counts$species_code == "MEGU",]$scientific_name <- "Larus canus" 
counts[counts$species_code == "MEGU",]$species_english_name <- "Mew Gull"
taxon <- read_csv("E:/Archive/0_data/NACC_list_species.csv")
#seperate scientific_name into genus and species
counts$genus <- str_split_fixed(counts$scientific_name,  " ", 2)[,1]
counts$species <- str_split_fixed(counts$scientific_name,  " ", 2)[,2]

# counts <- mutate(counts, genus = ifelse(species_code %in% c("OCWA", "TEWA"), "Leiothlypis", genus)) ##genus typo
counts <- mutate(counts, genus = ifelse(genus == "", NA, genus),
                 species = ifelse(species == "", NA, species))
counts <- left_join(counts, taxon %>% select(order, family, subfamily, genus)%>% distinct())
unidentified <- filter(counts, (is.na(family)| species_code %in% c("UNSW", "UNRE", "UPCH")) & !(species_code == "NONE")) %>% 
  select(species_code, species_english_name, order, family) %>% distinct()
##UNRE actually removed with removal of non-beringian sites. 

unidentified$genus <- NA


#Unidentified Woodpecker
try(unidentified[unidentified$species_code == "UNWO",]$order <- "Piciformes")
try(unidentified[unidentified$species_code == "UNWO",]$family <- "Picidae")	
#genus either "Picoides" or "Sphyrapicus" or "Colaptes" (woodpecker, sapsucker, flicker)

#Unidentified Blackbird
try(unidentified[unidentified$species_code == "UNWA",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNWA",]$family <- "Icteridae")

#Unidentified Warbler 
try(unidentified[unidentified$species_code == "UNWA",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNWA",]$family <- "Parulidae")

#Unidentified grebe
try(unidentified[unidentified$species_code == "UGRE",]$order <- "Podicipediformes")
try(unidentified[unidentified$species_code == "UGRE",]$family <- "Podicipedidae")

#Unidentified grouse
try(unidentified[unidentified$species_code == "UGRS",]$order <- "Galliformes")
try(unidentified[unidentified$species_code == "UGRS",]$family <- "Phasianidae")

#Unidentified loon
try(unidentified[unidentified$species_code == "ULOO",]$order <- "Gaviiformes")
try(unidentified[unidentified$species_code == "ULOO",]$family <- "Gaviidae")
try(unidentified[unidentified$species_code == "ULOO",]$genus <- "Gavia")


#Unidentified waterfowl (includes ducks geese, swans, etc.)
try(unidentified[unidentified$species_code == "UNWT",]$order <- "Anseriformes")
try(unidentified[unidentified$species_code == "UNWT",]$family <- "Anatidae")

#Unidentified Duck (assuming dabblers, excludes shovelers, wideons, gadwalls, seaducks)
try(unidentified[unidentified$species_code == "UNDU",]$order <- "Anseriformes")
try(unidentified[unidentified$species_code == "UNDU",]$family <- "Anatidae")
try(unidentified[unidentified$species_code == "UNDU",]$genus <- "Anas")

#UNK GOLDENEYE
try(unidentified[unidentified$species_code == "UGOL",]$order <- "Anseriformes")
try(unidentified[unidentified$species_code == "UGOL",]$family <- "Anatidae")
try(unidentified[unidentified$species_code == "UGOL",]$genus <- "Bucephala")

#Unidentified Flycatcher
try(unidentified[unidentified$species_code == "UNFL",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNFL",]$family <- "Tyrannidae")

#Unidentified Gull
try(unidentified[unidentified$species_code == "UNGU",]$order <-  "Charadriiformes")
try(unidentified[unidentified$species_code == "UNGU",]$family <- "Laridae")

#Unidentified Finch
try(unidentified[unidentified$species_code == "UNFI",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNFI",]$family <- "Fringillidae")

#Unidentified Sparrow
try(unidentified[unidentified$species_code == "UNSP",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNSP",]$family <- "Passerellidae")

#unidentified Swallow
try(unidentified[unidentified$species_code == "UNSW",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNSW",]$family <- "Hirundinidae")
try(unidentified[unidentified$species_code == "UNSW",]$genus <- "Hirundo")

#Unidentified Black-backed / American Three-toed Woodpecker
try(unidentified[unidentified$species_code == "UNBT",]$order <- "Piciformes")
try(unidentified[unidentified$species_code == "UNBT",]$family <- "Picidae")
try(unidentified[unidentified$species_code == "UNBT",]$genus <- "Picoides")

#Unidentified Chickadee
try(unidentified[unidentified$species_code == "UPCH",]$order <- "Piciformes")
try(unidentified[unidentified$species_code == "UPCH",]$family <- "Paridae")
try(unidentified[unidentified$species_code == "UPCH",]$genus <- "Poecile")

## Unidentified Thrush
try(unidentified[unidentified$species_code == "UNTH",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNTH",]$family <- "Turdidae")

##unidentified Blackbird
try(unidentified[unidentified$species_code == "UNBL",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNBL",]$family <- "Icteridae")

##unidentified Shorebird
try(unidentified[unidentified$species_code == "UNSH",]$order <-  "Charadriiformes")
try(unidentified[unidentified$species_code == "UNSH",]$family <- "Scolopacidae")

##unidentified Owl
try(unidentified[unidentified$species_code == "UNOW",]$order <-  "Strigiformes")

##unidentified Raptor
try(unidentified[unidentified$species_code == "UNRA",]$order <-  "Accipitriformes")
try(unidentified[unidentified$species_code == "UNRA",]$family <-  "Accipitridae")
try(unidentified[unidentified$species_code == "UNHA",]$order <-  "Accipitriformes")
try(unidentified[unidentified$species_code == "UNHA",]$family <-  "Accipitridae")
##unidentified falcon
try(unidentified[unidentified$species_code == "UNFA",]$order <-  "Falconiformes")
try(unidentified[unidentified$species_code == "UNFA",]$family <-  "Falconidae")
try(unidentified[unidentified$species_code == "UNFA",]$genus <-  "Falco")

##unidentified longspur, either LALO or SMLO
try(unidentified[unidentified$species_code == "UNLO",]$order <-  "Passeriformes")
try(unidentified[unidentified$species_code == "UNLO",]$family <- "Calcariidae")
try(unidentified[unidentified$species_code == "UNLO",]$genus <- "Calcarius")
##Unidentified Trill  ## same as passerine?
##Unidentified passerine in the American Robin Song Complex ## same as passerine?
##Unidentified Passerine
try(unidentified[unidentified$species_code == "UNAM",]$order <- "Passeriformes")
try(unidentified[unidentified$species_code == "UNTR",]$order <-  "Passeriformes")
try(unidentified[unidentified$species_code == "UNPA",]$order <-  "Passeriformes")

try(unidentified[unidentified$species_code == "UNKN",]$species_english_name <- "Unknown")
try(unidentified[unidentified$species_code == "UNPA",]$species_english_name <- "Unidentified Passerine")
try(unidentified[unidentified$species_code == "UNLO",]$species_english_name <- "Unidentified Longspur")
try(unidentified[unidentified$species_code == "UNBT",]$species_english_name <- "Unidentified Black-backed / American Three-toed Woodpecker")
try(unidentified[unidentified$species_code == "UNRA",]$species_english_name <- "Unidentified Raptor")
# unidentified[unidentified$species_code == "UNFA",]$species_english_name <- "Unidentified Falcon"
try(unidentified[unidentified$species_code == "UNFI",]$species_english_name <- "Unidentified Finch")
try(unidentified[unidentified$species_code == "UNGU",]$species_english_name <- "Unidentified Gull")
try(unidentified[unidentified$species_code == "UNTR",]$species_english_name <- "Unidentified Trill")
try(unidentified[unidentified$species_code == "UNSH",]$species_english_name <-  "Unidentified Shorebird")
try(unidentified[unidentified$species_code == "ULOO",]$species_english_name <- "Unidentified Loon")
try(unidentified[unidentified$species_code == "UGRE",]$species_english_name <-  "Unidentified Grebe")
try(unidentified[unidentified$species_code == "UGRS",]$species_english_name <-  "Unidentified Grouse")
try(unidentified[unidentified$species_code == "UNAM",]$species_english_name <-  "Unidentified passerine in the American Robin Song Complex")

unidentified <- distinct(unidentified)
unidentified 

counts <- counts %>% left_join(unidentified %>% rename(order_u = order, family_u = family, genus_u = genus, eng_u = species_english_name))
#remove order/family/genus from UNSW and UNRE, so all unknowns for their latin names
counts[counts$species_code %in% unidentified$species_code, c("order", "family", "genus", "species")] <- NA 
counts <- mutate(counts, species_english_name = ifelse(species_code %in% unidentified$species_code, eng_u, species_english_name)) %>% select(-eng_u)

counts$unk<- ifelse(counts$species_code %in% unidentified$species_code, T, F)


species_code <- c("BANS", "BARS", "CAWA", "CONI", "HUGO", "LEYE", "OSFL", "REKN", "BBSA", "EVGR", "HOGR", "REPH", "RUBL", "SEOW") ##redknot only roselaari type?
status <- c("THR", "SC", "SC", "SC", "THR", "THR", "SC", "THR", "SC", "SC", "SC", "SC", "SC", "THR")
sar <- data.frame(species_code, status)
sar %>% left_join(counts %>% select(species_code, species_english_name) %>% distinct()) ##NAs are SAR that weren't detected
counts$sar <- ifelse(counts$species_code %in% sar$species_code, T, F)

counts %>% filter(sar == T) %>% group_by(species_english_name) %>% summarise(n=n())

##large home ranges
order.HR = c("Accipitriformes", "Strigiformes")
family.HR = c("Laridae", "Gruidae")
sc.HR = c("UNRA", "UNGU", "UNOW", "UNHA", "CORA")
species.HR <- filter(counts, order %in% order.HR | family %in% family.HR | species_code %in% sc.HR) %>% select(species_code, species_english_name) %>% distinct()
species.HR



#N. stations per site
visits$year<-2024
visits$organization<-"WCS"
##wildtrax lat/longs
sites <- aru %>% select(siteID, station, lat = lat_wt, lon = lon_wt, organization) %>%   
  filter(!is.na(lat)) %>% distinct() %>% 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(4269)) %>% ##EPSG for NAD83 datum
  st_transform(crs = st_crs(3579)) ## EPSG for NAD83(CSRS) / Yukon Albers
sites <- left_join(sites, visits %>% select(siteID, year)%>%distinct())
sites$station <- as.numeric(sites$station)
counts_sum <- visits %>% select(siteID, station, year, method, dur_min) %>% distinct() %>% 
  group_by(siteID, year, method) %>% 
  summarise(n.stations = n(), mn.count.dur = mean(dur_min)) ##several stations with just 1 recording? don't have full hex code
hist(counts_sum$n.stations, breaks = seq(.5,max(counts_sum$n.stations) + .5, by = 1)) ## n stations per site. Lots! of sites with just 1 station

counts_sum <- as.data.frame(sites) %>% select(siteID, organization) %>% distinct() %>% left_join(counts_sum)

#recording time per site
counts_sum <- visits %>% group_by(siteID, method) %>% summarise(total.time = sum(dur_min)) %>% left_join(counts_sum)
hist(counts_sum$total.time)


counts <- mutate(counts, site2 = paste(siteID, method)) %>% distinct()
counts.l <- split(counts %>% left_join(visits %>% select(siteID, organization) %>% distinct()), f = counts$site2) ### list per siteID

##when to count unknowns...
counts.l <- lapply(counts.l, function(site){
  #create list of known genuses at that site. Unidentifieds have NA for order/family/genus/species, and *_u is filled in down to the level of identification.  Fully identified species have NAs in *_u columns. 
  genus_known <- unique(site$genus)
  genus_known <- genus_known[!is.na(genus_known)]
  ##if unknown genus (genus_u) is known, and is not in the list of known genuses at that site, set species to 'unk' (will be counted as unique spec)
  if(length(site[!is.na(site$genus_u) & !site$genus_u %in% genus_known,]$species) > 0) {
    site[!is.na(site$genus_u) & !site$genus_u %in% genus_known,]$species <- "unk"
  }
  #Repeat for family
  family_known <- unique(site$family)
  family_known <- family_known[!is.na(family_known)]
  
  if(length(site[!is.na(site$family_u) & !site$family_u %in% family_known,]$species) > 0){
    site[!is.na(site$family_u) & !site$family_u %in% family_known,]$species <- "unk"
    site[!is.na(site$family_u) & !site$family_u %in% family_known,]$genus <- "unk"
  }
  #Repeat for order
  order_known <- unique(site$order)
  order_known <- order_known[!is.na(order_known)]
  
  if(length(site[!is.na(site$order_u) & !site$order_u %in% order_known,]$species)>0){
    site[!is.na(site$order_u) & !site$order_u %in% order_known,]$species <- "unk"
    site[!is.na(site$order_u) & !site$order_u %in% order_known,]$genus <- "unk"
    site[!is.na(site$order_u) & !site$order_u %in% order_known,]$family <- "unk"
  }
  
  ##now copy over the orders, families, genus of unknowns if still missing
  site[is.na(site$order),]$order <- site[is.na(site$order),]$order_u
  site[is.na(site$family),]$family <- site[is.na(site$family),]$family_u
  site[is.na(site$genus),]$genus <- site[is.na(site$genus),]$genus_u
  
  return(site)
})

counts <- data.table::rbindlist(counts.l)

##unknowns are now filled in as species in some sites, and not others, depending on whether they should be counted. 
#filter(aru, species_code == "UNFI") %>% select(order, family, genus, species) %>% distinct()

##abundance
##abundance at sites with no detections is 0
counts[counts$species_code == "NONE", "abundance"] <- 0

##sum of all rows per site/station/time
counts_sum <- counts %>% group_by(siteID,  method) %>% summarise(abundance.all = sum(abundance)) %>%
  left_join(counts_sum) ##cumulative abundance across all stations/visits

##remove large HR birds
counts_sum <- counts %>% mutate(abundance = ifelse(species_code %in% species.HR$species_code, 0, abundance)) %>% #don't count large HR, so set abundance to 0
  group_by(siteID, station, ts, method) %>% summarise(abundance = sum(abundance)) %>% ##count abundance per recording
  group_by(siteID, station, method) %>% summarise(abundance = mean(abundance)) %>% ## mean abundance per station
  group_by(siteID, method) %>% summarise(abundance.mn = mean(abundance)) %>% ## mean abundance across all stations
  left_join(counts_sum)

##number of families
counts_sum <- counts %>% filter(!species_code == "NONE" & !is.na(family) & ##remove unknown families (NA) and counts with no detections 
                                  !(species_code %in% species.HR$species_code)) %>% ## remove large HR birds
  select(siteID, method, order, family) %>% distinct() %>% ##unique families per site
  group_by(siteID, method) %>% summarise(n.family = n()) %>% left_join(counts_sum)

counts_sum <- counts %>% filter(!species_code == "NONE" & !is.na(family)) %>% ##remove unknown families (NA) and counts with no detections
  select(siteID, method, order, family) %>% distinct() %>% ##unique families per site
  group_by(siteID, method) %>% summarise(n.family.all=n()) %>% left_join(counts_sum)

##number of species
counts_sum <- counts %>% filter(!species_code == "NONE" & !is.na(species) & ##remove unknown species (NA) and counts with no detections 
                                  !(species_code %in% species.HR$species_code)) %>% ## remove large HR birds
  select(siteID, method, order, family, genus, species) %>% distinct() %>% ##unique species per site
  group_by(siteID, method) %>% summarise(n.species = n()) %>% left_join(counts_sum)

counts_sum <- counts %>% filter(!species_code == "NONE" & !is.na(species)) %>% ##remove unknown species (NA) and counts with no detections 
  select(siteID, method, order, family, genus, species) %>% distinct() %>% ##unique species per site
  group_by(siteID, method) %>% summarise(n.species.all = n()) %>% left_join(counts_sum)

##abundance per species
species.abundance.l <- lapply(counts.l, function(site){
  site <- filter(site, !species_code == "NONE" & !is.na(species))
  
  ##make a dataframe where for each survey (station/time combination) there is a row for all species detected across all sites/surveys
  site.expand <- expand_grid(site %>% select(station, ts) %>% distinct(), site %>% select(species_code) %>% distinct()) 
  
  ##add abundance s to the surveys. IF the species was not detected in that survey, NA is returned. This should be changed to 0. 
  site.expand <- left_join(site.expand, site %>% select(station, ts, species_code, abundance)) 
  site.expand[is.na(site.expand$abundance),]$abundance <- 0
  
  ##now calculate average abundance across stations 
  site.expand <- site.expand %>% group_by(station, ts, species_code) %>% summarise(abundance = sum(abundance)) %>%##species count per survey
    group_by(station, species_code) %>% summarise(abundance = mean(abundance)) %>% ##highest species count per station
    group_by(species_code) %>% summarise(abundance = sum(abundance), n.stations = n()) #sum across sites
  site.expand$method <- unique(site$method)
  site.expand$organization <- unique(site$organization)
  site.expand$siteID <- unique(site$siteID)
  site.expand
  
})

counts_sum$n.LEYE <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "LEYE",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.OSFL <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "OSFL",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.CONI <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "CONI",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.RUBL <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "RUBL",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.BANS <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "BANS",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.BARS <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "BARS",]$abundance
  ifelse(length(x)==0, 0, x)
})

counts_sum$n.SEOW <- sapply(species.abundance.l, function(site) {
  x <- site[site$species_code == "SEOW",]$abundance
  ifelse(length(x)==0, 0, x)
})

##diversity index
##not accounting for n.station.
counts_sum$div.shan.all <- sapply(species.abundance.l, function(site) shannon.index(site$abundance))
counts_sum$div.shan <- sapply(species.abundance.l, function(site) shannon.index(site %>% 
                                                                                  filter(!(species_code %in% species.HR$species_code))%>%
                                                                                  pull(abundance)))

##top 20 species
counts %>% group_by(species_code, species_english_name) %>% summarise(abundance = n()) %>%
  arrange(desc(abundance)) %>%
  print(n=30)

table(counts_sum$n.LEYE >0)
table(counts_sum$n.OSFL>0)
table(counts_sum$n.CONI>0)
table(counts_sum$n.RUBL>0)
table(counts_sum$n.BANS>0)
table(counts_sum$n.BARS>0)
table(counts_sum$n.SEOW>0)

##number of species where at least 1 individuals was detected at more than 15 sites
dim(data.table::rbindlist(species.abundance.l) %>% group_by(species_code) %>% summarise(n.sites = n()) %>%
      filter(n.sites>15) %>% arrange(desc(n.sites)))
##50 species detected at more than 15 sites
dim(data.table::rbindlist(species.abundance.l) %>% group_by(species_code) %>% summarise(n.sites = n()) %>%
      filter(n.sites>50) %>% arrange(desc(n.sites)))



saveRDS(species.abundance.l,  file.path("E:/Archive/",pipeline, "out", "spec_abun2024data.RDS"))

sites_sf<-st_as_sf(sites)


st_write(sites_sf,'E:/Archive/2_pipeline/1_CountDataProcessing/out/sites24.shp',delete_layer = T)
sites_24<-st_read('E:/Archive/2_pipeline/1_CountDataProcessing/out/sites24.shp')
past_sites<-st_read('E:/Archive/2_pipeline/1_CountDataProcessing/out/sites.shp')
past_sites<-past_sites%>%filter(year!=2024)#%>%mutate(organization=orgnztn)%>%select(!orgnztn)

combined_sites<-past_sites%>%rbind(sites_24)


past_visits<-read.csv("E:/Archive/2_pipeline/1_CountDataProcessing/out/visits.csv")
past_visits<-past_visits%>%filter(year!=2024)
visits<-rbind(past_visits,visits)


st_write(combined_sites, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/", "out", "sites.shp"), delete_layer=T)
saveRDS(combined_sites, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/", "out", "sites.RDS"))
visits <- filter(visits, !(organization == "WCS" & method == "PC"))
write_csv(visits, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/", "out", "visits.csv"))
saveRDS(visits, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/", "out", "visits.RDS"))

past_counts<-read.csv("E:/Archive/2_pipeline/1_CountDataProcessing/out/counts.csv")%>%mutate(ts=ymd_hms(ts))
past_counts<-past_counts%>%mutate(year=lubridate::year(ts))%>%
  filter(year!=2024)%>%select(!year)
counts<-counts%>%mutate(largeHR = ifelse(species_code %in% species.HR$species_code, T, F))

all_counts<-rbind(past_counts,counts)

write_csv(all_counts,"E:/Archive/2_pipeline/1_CountDataProcessing/out/counts.csv")

saveRDS(all_counts,"E:/Archive/2_pipeline/1_CountDataProcessing/out/counts.RDS")


past_counts_sum<-read.csv("E:/Archive/2_pipeline/1_CountDataProcessing/out/counts_sum.csv")%>%filter(year!=2024)
counts_sum<-counts_sum%>%select(!div.shan,!n.species.all,!n.species.all)
all_counts_sum<-rbind(past_counts_sum,counts_sum)


write_csv(all_counts_sum, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/out/counts_sum.csv"))
saveRDS(all_counts_sum, file.path("E:/Archive/2_pipeline/1_CountDataProcessing/out/counts_sum.RDS"))
