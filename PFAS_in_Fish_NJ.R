library(dataRetrieval)
library(dplyr)
library(tidyverse)
library(factoextra)
library(ggfortify)
library(vegan)
library(RColorBrewer)


# Get PFAS parameters (tissue-based, with units)
pfas_parameters <- parameterCdFile %>%
  filter(grepl("PFAS|perfluoro|polyfluoro", parameter_nm, ignore.case = TRUE),
         parameter_units == "ng/kg",
         casrn != " ",
         srsname != " ")

# Extract unique characteristic names
characteristic <- pcode_to_name(pfas_parameters$parameter_cd)$characteristicname
characteristic <- unique(na.omit(characteristic))

# Download PFAS data for NJ, tissue samples from Morone americana
pfas_data <- readWQPdata(statecode = "New Jersey", characteristicName = characteristic)

# need to fix some unit differences, and estimate BDLs
# using half the detection limit here, which is good enough for small project
# but if this grows into something bigger we should use a better method
pfas_data <- pfas_data %>%
  filter(ActivityMediaName == "Tissue",
         SubjectTaxonomicName == "Morone americana") %>%
  mutate(
    ResultMeasureValue = ifelse(DetectionQuantitationLimitMeasure.MeasureUnitCode == "ng/g",
                                ResultMeasureValue * 1000,
                                ResultMeasureValue),
    DetectionQuantitationLimitMeasure.MeasureValue = ifelse(DetectionQuantitationLimitMeasure.MeasureUnitCode == "ng/g",
                                                            DetectionQuantitationLimitMeasure.MeasureValue * 1000,
                                                            DetectionQuantitationLimitMeasure.MeasureValue),
    ResultMeasureValue = ifelse(ResultDetectionConditionText %in% c("Not Detected", "Present Below Quantification Limit"),
                                DetectionQuantitationLimitMeasure.MeasureValue / 2,
                                ResultMeasureValue)
  )

# Aggregate and reshape data
pfas_grouped <- pfas_data %>%
  group_by(ActivityStartDate, MonitoringLocationIdentifier, CharacteristicName) %>%
  summarise(mean_value = median(ResultMeasureValue), .groups = "drop")

pfas_wide <- pfas_grouped %>%
  pivot_wider(names_from = CharacteristicName, values_from = mean_value)

# Manually drop problematic row if needed
# let's use the 8 most common PFAS compounds
pfas_wide <- pfas_wide[-c(41), 1:10]

# Prepare for clustering
pfas_matrix <- pfas_wide %>%
  select(-ActivityStartDate, -MonitoringLocationIdentifier)

pfas_matrix_imputed <- pfas_matrix %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  mutate(row_id = row_number()) %>%
  column_to_rownames("row_id")

pfas_scaled <- scale(na.omit(pfas_matrix_imputed))

# Cluster and assign cluster groups
dist_matrix <- dist(pfas_scaled)
pfas_hclust <- hclust(dist_matrix, method = "ward.D2")

fviz_nbclust(pfas_scaled, kmeans, method = "wss") +
  labs(title = "Scree Plot Method for Optimal Number of Clusters")

optimal_k <- 2
cluster_assignments <- cutree(pfas_hclust, k = optimal_k)

pfas_wide$cluster <- NA
matching_rows <- as.numeric(rownames(pfas_scaled))
pfas_wide$cluster[matching_rows] <- cluster_assignments

# PCA Visualization
pfas_pca_subset <- pfas_wide[matching_rows, ]

autoplot(prcomp(pfas_scaled), data = pfas_pca_subset, colour = 'cluster') +
  theme_minimal() +
  labs(title = "PCA of PFAS Profiles Colored by Cluster")

# Get site metadata (lat/lon)
site_ids <- unique(pfas_wide$MonitoringLocationIdentifier)
site_info <- whatWQPsites(siteid = site_ids)

site_info_clean <- site_info %>%
  select(MonitoringLocationIdentifier,
         MonitoringLocationName,
         LatitudeMeasure,
         LongitudeMeasure)

# does this shift have to do with time?

pfas_wide_geo <- pfas_wide %>%
  left_join(site_info_clean, by = "MonitoringLocationIdentifier") %>%
  mutate(cluster = as.factor(cluster))

# Temporal cluster shift visualization
ggplot(pfas_wide_geo, aes(x = as.Date(ActivityStartDate), y = cluster, color = cluster)) +
  geom_jitter(width = 0, height = 0.1, alpha = 0.6) +
  labs(title = "Temporal Shift in PFAS Clusters",
       x = "Date",
       y = "Cluster",
       color = "Cluster") +
  theme_minimal()


pfas_comp <- pfas_matrix_imputed  # your original imputed data
adonis_result <- adonis2(pfas_comp ~ cluster, data = pfas_wide_geo[matching_rows, ], method = "bray")

print(adonis_result)

pfas_wide_geo <- pfas_wide_geo %>%
  mutate(
    total_pfas = rowSums(select(., 3:10), na.rm = TRUE)
  )


pfas_relative <- pfas_wide_geo %>%
  mutate(across(3:10, ~ . / total_pfas, .names = "rel_{.col}"))

pfas_relative_long <- pfas_relative %>%
  select(ActivityStartDate, cluster, starts_with("rel_")) %>%
  pivot_longer(cols = starts_with("rel_"), names_to = "compound", values_to = "relative_value")

pfas_acronyms <- c(
  "Perfluorobutanesulfonate" = "PFBS",
  "Perfluorododecanoate" = "PFDoDA",
  "Perfluoroheptanoate" = "PFHpA",
  "Perfluorohexanesulfonate" = "PFHxS",
  "Perfluorohexanoate" = "PFHxA",
  "Perfluorooctanesulfonate" = "PFOS",
  "Perfluoropentanoate" = "PFPeA",
  "Perfluoroundecanoate" = "PFUnDA"
)

pfas_relative_long <- pfas_relative_long %>%
  mutate(compound = sub("^rel_", "", compound))

pfas_relative_long <- pfas_relative_long %>%
  mutate(compound = recode(compound, !!!pfas_acronyms))

ggplot(pfas_relative_long, aes(x = cluster, y = relative_value, fill = compound)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Relative Composition of PFAS through Time",
       y = "Proportion of Total PFAS", x = "Time Period",
       fill = "Compound") +
  scale_x_discrete(labels = c("1" = "Pre 2011", "2" = "Post 2011")) +
  theme_minimal()


pfas_polar <- pfas_relative_long %>%
  group_by(cluster, compound) %>%
  summarise(mean_prop = mean(relative_value, na.rm = TRUE), .groups = "drop")



# Choose an aquatic-friendly palette
fish_nature_palette <- c(
  "#1b4332",  # deep forest green
  "#2a6f4e",  # medium green
  "#40916c",  # softer green
  "#caf0f8",  # pale water blue
  "#00b4d8",  # lighter blue
  "#0077b6",  # ocean blue
  "#7f4f24",  # dark brown
  "#ddb892"   # light tan
)

facet_labels <- c(
  "1" = "Pre 2011\n────────────",
  "2" = "Post 2011\n────────────"
)

ggplot(pfas_polar, aes(x = compound, y = mean_prop, fill = compound)) +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.3) +
  coord_polar() +
  facet_wrap(~ cluster, labeller = as_labeller(facet_labels)) +
  scale_fill_manual(values = fish_nature_palette) +
  scale_y_sqrt(labels = scales::percent_format(accuracy = 1)) +
labs(
    title = "PFAS Fingerprint in White Perch by Time Period",
    subtitle = "Based on mean proportion of PFAS in NJ fish tissue samples",
    y = NULL, x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      angle = 50, vjust = 0.5, hjust = 1, size = 11, face = "bold"
    ),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(color = "gray90", size = 0.1),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    legend.position = "none"
  )
