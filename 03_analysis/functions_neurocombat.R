# Author: Lionel Butry (@lbutry)
# Function to apply neuroComBat to list of connectomes

apply_combat <- function(mat_list, clin_data, batch_var = "scanner", covars = c("Demo_Age", "Demo_Gender")) {
  
  # mat_list:     list of matrices
  # clin_data:    clinical data
  # covars:       character vector of column names to include as covariates
  
  # Reshape to long format (both triangles)
  long_all <- map_dfr(mat_list, ~ {
    rn <- rownames(.x)
    .x %>%
      as.data.frame() %>%
      rownames_to_column("roi1") %>%
      pivot_longer(-roi1, names_to = "roi2", values_to = "value") %>%
      select(roi1, roi2, value)
  }, .id = "participant_id") 
  
  # Extract only upper triangle without diagonal
  long_upper <- {
    roi_order <- rownames(mat_list[[1]])
    long_all |> filter(match(roi1, roi_order) < match(roi2, roi_order))
  }
  
  # Input matrix for combat 
  combat_input_all <- long_upper |> 
    mutate(diagonal = roi1 == roi2) |> 
    filter(diagonal == FALSE) |> 
    mutate(feature = paste(roi1, roi2, sep = "-")) |> 
    select(feature, participant_id, value) |> 
    pivot_wider(names_from = participant_id, values_from = value) |> 
    column_to_rownames("feature") |> 
    as.matrix()
  
  # Filter out features with the same value across all participants
  row_sd <- matrixStats::rowSds(combat_input_all)
  combat_input <- combat_input_all[row_sd > 0, ]
  combat_constant <- combat_input_all[row_sd == 0, ]
  
  # Batch variable
  batch_lookup <- setNames(clin_data[[batch_var]], clin_data$participant_id)
  batch <- colnames(combat_input)
  batch <- batch_lookup[batch]
  
  # Bring clin_data in same order as mat_list
  clin_data <- clin_data |> 
    filter(participant_id %in% names(mat_list)) |> 
    arrange(match(participant_id, names(mat_list)))
  
  # Dynamically build model matrix for covariates
  covars_formula <- paste(c("group", covars), collapse = " + ")
  formula <- as.formula(paste("~", covars_formula))
  mod <- model.matrix(formula, data = clin_data)
  
  # Apply ComBat
  data_harmonized <- neuroCombat::neuroCombat(dat=combat_input, batch=batch, mod = mod)
  
  # Add constant features back to harmonized data
  data_harmonized_all <- rbind(data_harmonized$dat.combat, combat_constant)

  # Restore original feature order
  data_harmonized_all <- data_harmonized_all[rownames(combat_input_all), ]
  
  # Transform ComBat output to list of matrices
  features <- rownames(data_harmonized_all)
  roi_pairs <- str_split_fixed(features, "-", 2)
  
  mat_list <- as.data.frame(data_harmonized_all) |> 
    as.list() |> 
    map2(colnames(data_harmonized_all),
         ~ vec2mat(.x, roi_pairs)) |> 
    set_names(colnames(data_harmonized_all))
  
  # Set NA to 0 (diagonal)
  mat_list <- map(mat_list, ~ {
    .x[is.na(.x)] <- 0
    .x
  })
}
