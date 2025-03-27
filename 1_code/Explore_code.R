#2024 field data descriptive statistics for reporting

library(wildrtrax)
library(tidyverse)

taxon <- read_csv("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/0_data/NACC_list_species.csv")

Sys.setenv(WT_USERNAME = 'austinzeller@rocketmail.com', WT_PASSWORD = 'Ke$ha50az') ## enter your own username and password
wt_auth()

aru24 <- wt_download_report(project_id = 2782, sensor_id = "ARU", report = c("main"), weather_cols = F)

aru24 <- aru24%>% filter(aru_task_status == "Transcribed")

wt_auth()
aru23 <- wt_download_report(project_id = 2078, sensor_id = "ARU", report = c("main"), weather_cols = F)

aru23 <- aru23%>%
  filter(aru_task_status == "Transcribed")%>%
  filter(observer!="Patrice Mathieu")

#number of speices found 2024
sp_codes24<-aru24%>%
  group_by(species_code)%>%
  summarize(n=n())%>%
  filter(!str_starts(species_code,"UN"))%>%
  filter(!species_code%in%c("BEAV",'WOFR',"UDAB", "UPCH", "NONE","UGOL"))

sp_codes24
##70 bird species in 2024#

sp_codes23<-aru23%>%
  group_by(species_code)%>%
  summarize(n=n())%>%
  filter(!str_starts(species_code,"UN"))%>%
  filter(!species_code%in%c("BEAV",'WOFR',"UDAB","UGOL","UPCH","NONE"))
sp_codes23


sp_site<-aru24%>%select(species_code,location, species_common_name)%>%
  mutate(site = str_extract(location, "^[^-]+"))%>%
  select(-location)%>%
  distinct()%>%
  group_by(species_code)%>%
  summarize(n=n(),species_common_name=first(species_common_name))%>%
  filter(!str_starts(species_code,"UN"))%>%
  filter(!species_code%in%c("BEAV",'WOFR',"UDAB","UPCH", "NONE","RESQ"))

sp_site_23<-aru23%>%select(species_code,location, species_common_name)%>%
  mutate(site = str_extract(location, "^[^-]+"))%>%
  select(-location)%>%
  distinct()%>%
  group_by(species_code)%>%
  summarize(n=n(),species_common_name=first(species_common_name))%>%
  filter(!str_starts(species_code,"UN"))%>%
  filter(!species_code%in%c("BEAV",'WOFR',"UDAB","UPCH", "NONE", "RESQ","UGOL"))


# family
sp_site$common_name<-sp_site$species_common_name
fams<-left_join(sp_site,taxon,by=c("common_name"))
fams%>%select(family)%>%distinct()

#individuals
individual<-aru24%>%select(location,species_code,recording_id,individual_order,individual_count)%>%
  group_by(recording_id)%>%
  summarize(max_count=max(individual_order))%>%
  filter(!is.na(max_count))
sum(individual$max_count)

#SAR INFO

SAR_codes<-c("BANS", "BARS", "CAWA", "CONI", "HUGO", "LEYE", "OSFL", "REKN", "BBSA", "EVGR", "HOGR", "REPH", "RUBL", "SEOW")
SAR_site<-sp_site%>%filter(species_code%in%SAR_codes)
SAR_site


#mine vs unmined richness

aru<-aru24%>%
  mutate(mine=if_else(startsWith(location,"MS"),"Mined","Unmined"))

location_count<-aru%>%
  group_by(location)%>%
  summarize(count = n_distinct(species_code),mine=first(mine))

ggplot(location_count, aes(x = mine, y = count)) +
  geom_boxplot() +
  labs(x="",y = "Species Richness",title="2024 Species Richness") +
  theme_classic()
  

site_richness<-location_count%>%
  mutate(site = str_extract(location, "^[^-]+"))%>%
  group_by(site)%>%
  summarize(richness=mean(count),mine=first(mine))



#SAR
SAR_ARU<-aru%>%
  filter(species_code%in%SAR_codes)%>%
  group_by(species_code)

SARlocation_count<-SAR_ARU%>%
  group_by(location)%>%
  summarize(count = n_distinct(species_code),mine=first(mine))

hab.sum<-readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/2_HabitatSum/out/hab_sum.RDS")
temp<-aru24%>%filter(species_code%in%SAR_codes)%>%mutate(site = str_extract(location, "^[^-]+"))

SARsites<-temp%>%select(species_code,location,site)%>%distinct()
hab.sum<-hab.sum%>%filter(buffer==150)

LEYEsites<-SARsites%>%filter(species_code=="LEYE")
LEYE<-hab.sum%>%filter(siteID%in%LEYEsites$site)

BANSsites<-SARsites%>%filter(species_code=="BANS")
BANS<-hab.sum%>%filter(siteID%in%BANSsites$site)
BANS

OSFLsites<-SARsites%>%filter(species_code=="OSFL")
OSFL<-hab.sum%>%filter(siteID%in%OSFLsites$site)

RUBLsites<-SARsites%>%filter(species_code=="RUBL")
RUBL<-hab.sum%>%filter(siteID%in%RUBLsites$site)

HOGRsites<-SARsites%>%filter(species_code=="HOGR")
HOGR<-hab.sum%>%filter(siteID%in%HOGRsites$site)

## Spp table
sp_site_23<-sp_site_23%>%mutate(Sites_23=n)%>%select(!n)
sp_site_24<-sp_site%>%mutate(Sites_24=n)%>%select(!n)
sp_codes23<-sp_codes23%>%mutate(Detections_23=n)%>%select(!n)
sp_codes24<-sp_codes24%>%mutate(Detections_24=n)%>%select(!n)

sp_table<-full_join(sp_site_23,sp_site_24)%>%
  select(!common_name)%>% #number of sites
  mutate(Prop_site23=Sites_23/41, Prop_site24=Sites_24/38)%>% # % of sites detected sp in
  left_join(sp_codes23)%>%
  left_join(sp_codes24)%>%
  replace(is.na(.),0)%>%
  mutate(Prop_difference= Prop_site24-Prop_site23)%>%
  relocate(Prop_difference, .after= Prop_site24)
  
  
write.csv(file="C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/species_table.csv",
          sp_table)
sp_table
  

sp_table<-sp_site%>%mutate(Sites_24=n)%>%
  mutate(SAR=ifelse(species_code%in%SAR_codes,"Yes","No"))%>%
  select(!n)%>%
  left_join(sp_codes24,by="species_code")%>%
  mutate(Tags_24=n)%>%
  select(!n)%>%
  left_join(sp_codes23, by="species_code")%>%
  mutate(Tags_23=n)%>%
  select(!n)
  
  


sp_table<-select(common_name,Sites_24,SAR)


