task_meta <- read_csv("0_data/South Beringia Priority Place Cumulative Effects 2023_Tasks_202331.csv")
task_meta$recordingDate <- ymd_hms(task_meta$recordingDate, tz = "Canada/Yukon")
recording_meta <- list()
recording_meta[[1]] <- wt_download_report(project_id = project_ids[1], sensor_id = "ARU", report = c("recording"), weather_cols = F)
recording_meta[[2]] <- wt_download_report(project_id = project_ids[2], sensor_id = "ARU", report = c("recording"), weather_cols = F)
recording_meta[[3]] <- wt_download_report(project_id = project_ids[3], sensor_id = "ARU", report = c("recording"), weather_cols = F)

loudness <- recording_meta[[2]] %>% 
  mutate(max_level = (right_max_level+ left_max_level)/2,
         min_level = (right_min_level+ left_min_level)/2,
         peak_level_dbfs = (right_peak_level_dbfs+ left_rms_peak_level_dbfs)/2,
         pk_count = (right_pk_count+ left_pk_count)/2,
         rms_peak_dbfs = (right_rms_peak_dbfs+ left_rms_peak_dbfs)/2,
         rms_trough_dbfs = (right_rms_trough_dbfs+ left_rms_trough_dbfs)/2) 
loudness$recordingDate <- ymd_hms(loudness$recording_date_time, tz = "Canada/Yukon") 
loudness <- loudness %>%  select(location, recordingDate, 
         min_level, max_level, peak_level_dbfs, pk_count, rms_peak_dbfs, 
         rms_peak_dbfs, rms_trough_dbfs)

loudness <- left_join(task_meta, loudness) %>% select(location, recordingDate, industryNoise, otherNoise, audioQuality, taskComments,
                                          min_level, max_level, peak_level_dbfs, pk_count, rms_peak_dbfs, 
                                          rms_peak_dbfs, rms_trough_dbfs)
loudness$industryNoise <- factor(loudness$industryNoise, levels = c("None","Light","Moderate","Heavy"), ordered = T)
loudness$otherNoise <- factor(loudness$otherNoise, levels = c("None","Light","Moderate","Heavy"), ordered = T)
loudness <- loudness %>% rowwise %>% mutate(obsNoise = max(industryNoise, otherNoise, na.rm = T))


ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, min_level)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, max_level)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, peak_level_dbfs)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, pk_count)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, rms_peak_dbfs)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, rms_peak_dbfs)) + geom_boxplot()
ggplot(loudness %>% filter(!is.na(otherNoise)), aes(obsNoise, rms_trough_dbfs)) + geom_boxplot() ## best correlate of observer persived noise


aru <- read_rds(file.path(pipeline, "tmp", "wt_2024_01_10.RDS"))
tags <- bind_rows(aru)
task <- bind_rows(recording_meta)
task <- task %>% 
  mutate(max_level = (right_max_level+ left_max_level)/2,
         min_level = (right_min_level+ left_min_level)/2,
         peak_level_dbfs = (right_peak_level_dbfs+ left_rms_peak_level_dbfs)/2,
         pk_count = (right_pk_count+ left_pk_count)/2,
         rms_peak_dbfs = (right_rms_peak_dbfs+ left_rms_peak_dbfs)/2,
         rms_trough_dbfs = (right_rms_trough_dbfs+ left_rms_trough_dbfs)/2) 

task <- left_join(tags %>% group_by(recording_id) %>% summarise(n.tags = n()), task)

ggplot(task, aes(rms_trough_dbfs, n.tags)) + geom_point() + geom_smooth() ## definitely detect more with lower rms troughdbfs
ggplot(task, aes(rms_peak_dbfs, n.tags)) + geom_point() + geom_smooth()
ggplot(task, aes(min_level, n.tags)) + geom_point() + geom_smooth()
ggplot(task, aes(max_level, n.tags)) + geom_point() + geom_smooth()
ggplot(task, aes(pk_count, n.tags)) + geom_point() + geom_smooth()

loudness <- left_join(loudness, task %>% mutate(recordingDate = ymd_hms(recording_date_time, tz = "Canada/Yukon")) %>%
                        select(location, recordingDate, n.tags))

ggplot(loudness, aes(rms_trough_dbfs, n.tags)) + geom_point() + geom_smooth()
ggplot(loudness %>% filter(!(is.na(obsNoise))), aes(obsNoise, n.tags)) + geom_boxplot()

ggplot(loudness %>% filter(!(is.na(obsNoise))), aes(obsNoise, rms_trough_dbfs)) + geom_boxplot()

ggplot(loudness, aes(rms_trough_dbfs, n.tags, col = obsNoise)) + geom_point() + geom_smooth(se = F)



loudness %>% mutate(Noise = ifelse(obsNoise == "Heavy", "Heavy",
                                   ifelse(obsNoise == "Moderate", "Moderate", "Light")))  %>% 
  group_by(Noise) %>% summarise(mn.tags = mean(n.tags), med.tags = median(n.tags))

