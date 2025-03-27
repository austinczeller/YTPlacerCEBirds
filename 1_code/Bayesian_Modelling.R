#Bayesian Analysis
library(tidyverse)
library(rjags)
library(coda)


totit <- 4e5    # Number of iterations
thinning <- ceiling(totit / 1e4)  # Thinning rate
nchains <- 4    # Number of chains
K <- 5          # Number of cross-validation folds

# Load habitat covariates
hab.sum <- read_rds("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/2_HabitatSum/out/hab_sum.RDS")

hab.sum$h_open <- rowSums(hab.sum[,c("h_shrub_tall", "h_shrub_open", "h_shrub_low", 
                                     "h_barren", "h_herb", "h_sparce_veg", "h_tussock_tundra",
                                     "h_burn-herb", "h_burn-shrub", "h_burn-sappling",
                                     "h_fen_nasa", "h_bog_nasa")])

hab.sum$d_lden_all <- rowSums(hab.sum[,c("d_lden_r", "d_lden_t", "d_lden_c")])

# Remove unnecessary columns
hab.sum <- hab.sum %>% ungroup()%>%
  select(siteID, buffer, d_surface, h_evergreen_forest, h_woodland,
        h_open, wden, lon_wgs84, lat_wgs84, h_wet, d_wet,avg_prcp,d_lden_all, elev)%>%
  filter(buffer!=500)%>%
  distinct()

head(hab.sum)

hab.sum_wide <- hab.sum %>%
  tidyr::pivot_wider(
    id_cols = siteID,
    names_from = buffer,
    values_from = c(h_wet, h_open, h_evergreen_forest, h_woodland, wden,
                    d_wet,  d_surface, d_lden_all, avg_prcp),
    names_sep = "_"
  )


#hab.sum_wide<-hab.sum_wide%>%mutate(avg_prcp=avg_prcp_150,lat=lat_150,lon=lon_150)%>%
#  select(-lat_1000,-lon_1000,-avg_prcp_1000,-avg_prcp_150,-lat_150,-lon_150)
head(hab.sum_wide)

# Load species counts
counts <- read_rds("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/3_Detectability/out/sitelevel_counts.RDS")

# Join species counts with habitat covariates
# Ensure a common `siteID` column exists in both datasets
model_data <- counts %>%
  inner_join(hab.sum_wide, by = "siteID")

# Check for missing data
model_data <- model_data %>%
  drop_na() # Remove rows with missing data

model_data<-model_data%>%
  select(-perARU,-noise.tag.abun,-noise.linear,-noise.mode,-noise.sq)

head(model_data)

species_list <- unique(model_data$species_code)
results<-list()
set.seed(17)


# Loop over each species
for (sp in species_list) {
  for (buffer in c(150, 1000)) {
    cat("Running model for species:", sp, "Buffer:", buffer, "\n")
  
    # Filter data for the species
    sp_data <- model_data %>%
      filter(species_code == sp)%>%
      select(siteID, species_code, abundance, offset,
             ends_with(paste0(buffer)))
      
  
  # Combine habitat and disturbance variables into a single matrix
  cov_matrix <- sp_data %>%
    select(starts_with("h_"), starts_with("d_"),starts_with("avg")) %>%
    as.matrix()
  
  # Standardize covariates
  cov_matrix <- scale(cov_matrix)

  # Prepare data for JAGS
  jags_data <- list(
    y = sp_data$abundance,
    offset = sp_data$offset,
    cov = cov_matrix,
    nobs = nrow(sp_data),
    ncovs = ncol(cov_matrix)
  )
  
  # Path to model file
  model_path <- "C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/model.txt"
  
  # Initialize and run the JAGS model
  jags_model <- jags.model(file = model_path,
                           data = jags_data,
                           n.chains = 4,
                           n.adapt = totit,
                           quiet = TRUE)
  
  update(jags_model, n.iter=(totit/2))  # Burn-in
  
  samples <- coda.samples(jags_model,
                          variable.names=c("alpha","beta"),
                          thin= thinning,
                          n.iter=totit/thinning)#I need to figure out what I need to do with n.inter here
  
  # Store results
  results[[paste(sp, buffer, sep = "_")]] <- samples
  beepr::beep(4)
  }
}

# Optional: Save results
saveRDS(results, "C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/results.RDS")





