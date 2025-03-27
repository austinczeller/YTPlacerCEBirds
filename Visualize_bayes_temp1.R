library(ggplot2)
library(coda)
results <- readRDS("~/CE_birds/results.RDS")

# Function to visualize results for a single species
visualize_species <- function(samples, species_name, covariate_names) {
  # Combine chains into a single data frame
  combined_samples <- as.data.frame(do.call(rbind, lapply(samples, as.matrix)))
  
  # Extract parameter names and map beta coefficients to covariate names
  param_names <- colnames(combined_samples)
  beta_indices <- grep("^beta", param_names)
  beta_params <- param_names[beta_indices]
  
  # Rename beta parameters with covariate names
  renamed_params <- sapply(beta_params, function(param) {
    idx <- as.numeric(gsub("beta\\[(\\d+)\\]", "\\1", param))
    covariate_names[idx]
  })
  
  # Coefficients with credible intervals
  coef_summary <- summary(samples)$statistics
  credible_intervals <- summary(samples)$quantiles
  
  coef_data <- data.frame(
    Parameter = c("alpha", renamed_params),
    Mean = c(coef_summary["alpha", "Mean"], coef_summary[beta_indices, "Mean"]),
    Lower = c(credible_intervals["alpha", "2.5%"], credible_intervals[beta_indices, "2.5%"]),
    Upper = c(credible_intervals["alpha", "97.5%"], credible_intervals[beta_indices, "97.5%"])
  )
  
  # Highlight consistent variables (95% CI does not overlap zero)
  coef_data$Consistent <- with(coef_data, Lower > 0 | Upper < 0)
  
  # Coefficients plot
  coef_plot <- ggplot(coef_data, aes(x = reorder(Parameter, Mean), y = Mean, color = Consistent)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
    labs(title = paste("Coefficients with 95% Credible Intervals for", species_name),
         x = "Parameter", y = "Mean Estimate") +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = FALSE) +
    theme_minimal() +
    coord_flip()
  
  # Print plot
  print(coef_plot) # Explicit print
}

# Specify species to visualize
selected_species <- c("AMRO", "YRWA", "DEJU", "SWTH", "SAVS")

# Loop through results for selected species and visualize
for (sp in names(results)) {
  if (sp %in% selected_species) {
    cat("Visualizing results for species:", sp, "\n")
    visualize_species(
      samples = results[[sp]], 
      species_name = sp, 
      covariate_names = colnames(cov_matrix) # Pass covariate names
    )
  }
}


ggplot(coef_data, aes(x = reorder(Parameter, Mean), y = Mean, color = Consistent)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  labs(title = paste("Coefficients with 95% Credible Intervals for", species_name),
       x = "Parameter", y = "Mean Estimate") +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = FALSE) +
  theme_minimal() +
  coord_flip()
visualize_species(results$YRWA_150,species_name=YRWA,covariate_names = colnames(cov_matrix))

visualize_species(
  samples = results[[sp]], 
  species_name = sp, 
  covariate_names = colnames(cov_matrix))

  