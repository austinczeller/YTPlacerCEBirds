
library(tidyverse)
library(sf)
library(vegan)
library(terra)

hab.sum <- readRDS("2_pipeline/4_DataExploration/store/hab.sum_simp.RDS" )
counts <- read_rds("2_pipeline/1_CountDataProcessing/out/counts.RDS")
sites <- readRDS("2_pipeline/2_HabitatSum/store/sites_habpoly.RDS")
sites <- filter(sites, siteID %in% hab.sum$siteID) %>% group_by(siteID, year) %>% summarise()

# write_rds(vegmask, file.path(pipeline, "store", "vegmask_dclip_rast.rds"))

hab23 <- filter(hab.sum, year == 2023 & buffer == 150) 

rip.heli <- c("MR900",  "MR901",  "MR902",  "MR903",  "MR904",  "MR905",  "MR906",  "MR907",  "MR908")
rip.rd <- c("MR501",  "MR502",  "MR508",  "MR509",  "MR520", "MRN501", "MRN502")
hab23$class <- ifelse(hab23$siteID %in% rip.heli, "Riparian - Remote", ifelse(hab23$siteID %in% rip.rd, "Riparian - Near", "Mine"))

## Mn NDVI of mine sites ----
ndvi <- readRDS("2_pipeline/2_HabitatSum/store/disturb_ndvi_median_23.rds")
ndvi <- data.frame(siteID = names(ndvi), med.ndvi = ndvi)
hab23 <- left_join(hab23, ndvi)
hab23$med.ndvi <- ifelse(hab23$class != "Mine", NA, hab23$med.ndvi)
hab23$class <- ifelse(hab23$class != "Mine", hab23$class,
                      ifelse(hab23$med.ndvi >= median(hab23$med.ndvi, na.rm = T), "Mine - Veg", "Mine - Bare"))
hab23  %>% ggplot(aes(d_surface, group=class, fill = class)) + geom_boxplot() ## pretty even %cover of surface between the bare and veg categories
## counts matrix ----
counts <- filter(counts, !largeHR & !unk & siteID %in% hab23$siteID & !species_code == "NONE") %>% 
  group_by(siteID, species_code) %>% summarise(abundance = sum(abundance))

## avoid too many 0s, remove the rarest species
spec.select <- counts %>% 
  select(siteID, species_code) %>% distinct() %>%
  group_by(species_code) %>% summarise(n.sites = n(), p.sites = n.sites/length(unique(hab23$siteID)))
spec.select <- filter(spec.select, n.sites >4) %>% pull(species_code)
counts.m <- pivot_wider(counts %>% filter(species_code %in% spec.select) %>% 
                          select(siteID, species_code, abundance), 
                        id_cols = siteID, names_from=species_code, 
                        values_from = abundance, values_fill = 0)
tmp <- as.matrix(counts.m[,2:length(counts.m[1,])])
rownames(tmp) <- counts.m$siteID
counts.m <- tmp

hab23 <- hab23[order(factor(hab23$siteID, levels=row.names(counts.m))),]

## NMDS ----

set.seed(123)
ord <- metaMDS(counts.m, k = 3, try = 60)

plot(ord)

##1 and 2
plot(ord, type = "n")
text(ord, display = "spec", cex=0.7)
ordihull(ord, hab23$class,   col=1:4, lwd=3, label = T) ## black is mine, red is riparian
# ordihull(ord, hab23$class, col=1:4, lwd=3) ## black is mine, red is riparian
# ordiellipse(ord, hab23$class, col=1:4, kind = "ehull", lwd=3)
# ordispider(ord, hab23$class, col=1:4, label = TRUE, spiders = c("median"))
# ordiellipse(ord, hab23$class, col=1:4, draw="polygon")
# points(ord, display = "sites", cex = 0.8, pch=21, col="red")

## 1 and 3
plot(ord, type = "n", choices = c(1,3))
text(ord, display = "spec", choices = c(1,3), cex=0.7)
ordihull(ord, hab23$class, choices = c(1,3),   col=1:4, lwd=3, label = T) 

## 2 and 3
plot(ord, type = "n", choices = c(2,3))
text(ord, display = "spec", choices = c(2,3), cex=0.7)
ordihull(ord, hab23$class, choices = c(2,3),  col=1:4, lwd=3, label = T) 


# 
min(rowSums(counts.m))
set.seed(1234)
rare.df <- rrarefy(counts.m, 50) 
ord <- metaMDS(rare.df, k = 2, try = 60)

hab23$class_abr <- as.factor(ifelse(hab23$class == "Mine - Bare", "MB", 
                          ifelse(hab23$class == "Mine - Veg", "MV",
                                 ifelse(hab23$class == "Riparian - Remote", "RR", 
                                        "RN"))))
plot(ord, type = "n")
ordiellipse(ord, hab23$class, col=1:4, draw="line", lwd=3,)
# ordihull(ord, hab23$class,   col=1:4, lwd=3, label = F)
text(ord, display="sites", labels = as.character(hab23$class_abr),
                    col=as.numeric(hab23$class_abr))

# rare.df <- rrarefy(counts.m, 20)
rare.rich <- rowSums(rare.df>1)
rare.rich <- data.frame(siteID = names(rare.rich), richness = rare.rich)
rare.rich <- left_join(rare.rich, select(hab23, siteID, class))

ggplot(rare.rich, aes(class, richness, group = class, fill = class)) + geom_boxplot(show.legend = FALSE) + scale_fill_manual(values = 1:4) +
  xlab("Site type") + ylab("Species Richness")
table(rare.rich$class)
