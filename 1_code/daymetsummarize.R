#make daymet weather good to incorporate into analysis

visits <- readRDS("C:/Users/azeller/OneDrive - Wildlife Conservation Society/Documents/CE_birds/2_pipeline/1_CountDataProcessing/out/visits.RDS") 


library(dplyr)
library(readr)
library(janitor)

# Path to the folder containing the weather data files
weather_data_folder <- "~/CE_birds/0_data/daymet"

# List all CSV files in the folder
csv_files <- list.files(weather_data_folder, pattern = "\\.csv$", full.names = TRUE)

# Summarize the sampling period from the visits dataframe
sampling_period <- visits %>%
  group_by(siteID) %>%
  summarize(
    start = min(yday(ts)),
    end = max(yday(ts)),
    .groups = "drop"
  )

# Initialize an empty list to store dataframes for each site
weather_summary_list <- list()

# Loop over all CSV files
for (file in csv_files) {
  # Extract the siteID and year from the filename
  filename <- basename(file)
  siteID <- sub("^(.*?)year\\d{4}.*$", "\\1", filename)
  year <- as.numeric(sub(".*?year(\\d{4}).*$", "\\1", filename))
  
  # Read the weather data, skipping the first 6 rows and cleaning column names
  weather_data <- read_delim(file, delim = ",", skip = 6) %>%
    clean_names()
  
  # Filter for the site's sampling period
  site_sampling_period <- sampling_period %>% filter(siteID == !!siteID)
  if (nrow(site_sampling_period) == 0) next
  
  start_yday <- site_sampling_period$start
  end_yday <- site_sampling_period$end
  
  weather_filtered <- weather_data %>%
    filter(yday >= start_yday, yday <= end_yday)
  
  # Summarize the weather data
  weather_summary <- weather_filtered %>%
    summarize(
      siteID = siteID,
      year = year,
      avg_prcp = mean(prcp_mm_day, na.rm = TRUE),
      max_prcp = max(prcp_mm_day, na.rm = TRUE),
      avg_tmax = mean(tmax_deg_c, na.rm = TRUE),
      max_tmax = max(tmax_deg_c, na.rm = TRUE)
    )
  
  # Append the summarized data to the list
  weather_summary_list[[length(weather_summary_list) + 1]] <- weather_summary
}

# Combine all site summaries into a single dataframe
final_weather_summary <- bind_rows(weather_summary_list)
final_weather_summary<-final_weather_summary%>%distinct()

saveRDS(final_weather_summary,"~/CE_birds/0_data/final_weather_summary.RDS")
# View the resulting dataframe
print(final_weather_summary)

