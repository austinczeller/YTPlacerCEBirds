#CCA's

library(tidyverse)
library(vegan)
library(corrplot)



hab.sum<-readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/2_HabitatSum/out/hab_sum.RDS")%>%
  filter(buffer==1000)%>%filter(siteID!="MR40")
counts<-readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/1_CountDataProcessing/out/counts.RDS")
counts<-counts%>%mutate(stationID=paste0(siteID,"-",station))
head(hab.sum)
head(counts)



# Create a species-by-site abundance matrix
species_abundance <- counts %>%
  group_by(stationID, species_code) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = species_code, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("stationID")

species_abundance <- species_abundance[rowSums(species_abundance) > 0, ]



#edit the variables
hab.sum$d_lden_all <- rowSums(hab.sum[,c("d_lden_r", "d_lden_t", "d_lden_c")])

#only keep stuff we will use
hab.sum<-hab.sum%>%select(stationID,elev,d_pVeg,d_pBare,d_lden_all,d_surface,h_evergreen_forest,h_woodland,h_wet)


# Select habitat variables and set stationID as rownames
habitat_data <- hab.sum %>%
  column_to_rownames("stationID")

habitat_data <- habitat_data[rownames(species_abundance), ]

# Ensure stationIDs match between datasets
species_abundance <- species_abundance[rownames(habitat_data), ]

habitat_data <- na.omit(habitat_data)
species_abundance <- species_abundance[rownames(habitat_data), ]  # Ensure species data matches



#check for colinearity
cor_matrix <- cor(habitat_data, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)


# Run the CCA
cca_result <- cca(species_abundance ~ ., data = habitat_data)
disturbance_cca <- cca(species_abundance ~ d_pBare+d_pVeg+d_surface+d_lden_all, data = habitat_data)
hab_cca<-cca(species_abundance ~ h_woodland+h_wet+h_evergreen_forest, data = habitat_data)

# View summary
plot(cca_result, display = c("species","bp"))
species_scores <- scores(cca_result, display = "species")
text(species_scores[, 1], species_scores[, 2], labels = rownames(species_scores), col = "black", cex = 0.7)


plot(disturbance_cca, display = c("species","bp"))
species_scores <- scores(disturbance_cca, display = "species")
text(species_scores[, 1], species_scores[, 2], labels = rownames(species_scores), col = "black", cex = 0.7)


plot(hab_cca, display = c("species","bp"))
species_scores <- scores(hab_cca, display = "species")
text(species_scores[, 1], species_scores[, 2], labels = rownames(species_scores), col = "black", cex = 0.7)

