
install.packages("hdf5r")
install.packages("jsonlite")

library(Seurat)
library(dplyr)
library(hdf5r)
library(jsonlite)

# Define the paths to your .json files
json_paths <- list(
  "dpc2.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc2.5/dpc2.5_summary.json",
  "dpc3.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc3.5/dpc3.5_summary.json",
  "dpc4" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc4/dpc4_summary.json",
  "dpc4.5" = "/Users/lib9aj/git/scRNA-seq/dpc2.5_3.5_4_4.5_Haixiang_Sun/dpc4.5/dpc4.5_summary.json"
)


# Function to extract detailed run information from a JSON file
get_run_info <- function(json_path) {
  json_data <- fromJSON(json_path)
  # Check if 'batches' field exists and is not empty
  if ("batches" %in% names(json_data)) {
    run_info <- data.frame(
      run_id = json_data$batches,
      day = gsub("_summary.json", "", basename(json_path)),
      cell_count = sapply(json_data$batches, function(x) json_data[[paste0(x, "_filtered_bcs")]]),
      feature_count = sapply(json_data$batches, function(x) json_data[[paste0(x, "_pre_normalization_feature_reads")]])
    )
    return(run_info)
  } else {
    return(data.frame(run_id = character(0), day = gsub("_summary.json", "", basename(json_path))))
  }
}

# Iterate over the JSON paths and get detailed run information
run_info_list <- lapply(json_paths, get_run_info)

# Combine all run information into a single data frame
all_run_info <- bind_rows(run_info_list)

# Print detailed run information
print(all_run_info)

# Count the number of runs for each day and summarize other statistics
run_counts <- all_run_info %>%
  group_by(day) %>%
  summarise(
    run_count = n(),
    run_ids = paste(run_id, collapse = ", "),
    total_cells = sum(cell_count, na.rm = TRUE),
    total_features = sum(feature_count, na.rm = TRUE),
    avg_cells_per_run = mean(cell_count, na.rm = TRUE),
    avg_features_per_run = mean(feature_count, na.rm = TRUE)
  )

# Print the summarized information for each day
print(run_counts)






