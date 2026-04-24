# Author: Lionel Butry (@lbutry)
# Purpose: Computes global network properties of the functional connectome
#          for full, DMN, VAN & SMN

source("02_scripts/03_analysis/set_env.R")
source("02_scripts/03_analysis/functions_neurocombat.R")
source("02_scripts/03_analysis/functions_graphtheory.R")

############################################
# ---- Import data ----
############################################

# FC
source("02_scripts/03_analysis/import_fc.R")
mat_list <- fc_mat_list
rm(fc_mat_list)

# Clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") |> 
  filter(participant_id %in% names(mat_list))

# Atlas
bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv")

############################################
# ---- Apply neuroComBat ----
############################################

mat_list_combat <- apply_combat(
  mat_list,
  clin_df,
  batch = "scanner",
  covars = c("Demo_Age", "Demo_Gender"))

###############################################
# Split connectome into subnetworks
###############################################

# DMN (36x36)
rois_dmn <- bn_lut$ROI.Name[bn_lut$network == "DMN"]
mat_dmn_list <- map(mat_list_combat, ~ .x[rois_dmn, rois_dmn])

# SMN (33x33)
rois_smn <- bn_lut$ROI.Name[bn_lut$network == "SMN"]
mat_smn_list <- map(mat_list_combat, ~ .x[rois_smn, rois_smn])

# VAN (22x22)
rois_van <- bn_lut$ROI.Name[bn_lut$network == "VAN"]
mat_van_list <- map(mat_list_combat, ~ .x[rois_van, rois_van])

rm(rois_dmn, rois_smn, rois_van)

############################################
# ---- Apply density threshold ----
############################################

# Wrapper function
thresholding <- function(mat_list, densities) {
  
  mat_list <- map(mat_list, ~ abs(.x))
  
  # Iterate over all subject & across defined density range
  ## Get matrix ranked by connectivity strength + MST as backbone
  ranks_list <- map(mat_list, ~ mst_gthr_order(.x))
  
  ## Apply density thresholds
  thr <- map2(mat_list, ranks_list, function(mat, mat_ranked) {
    map(densities, ~ apply_mst_gthr(mat, mat_ranked, .x)) |>
      set_names(paste0("d", densities)) 
  }) |> list_flatten()
  
  return(thr)
}

# Apply thresholding & store in main_list
densities <- seq(0.1, 0.4, 0.05) # 10-40% in 5% steps

main_list <- list(
  full = thresholding(mat_list_combat, densities),
  dmn = thresholding(mat_dmn_list, densities),
  smn = thresholding(mat_smn_list, densities),
  van = thresholding(mat_van_list, densities)
)

###############################################
# ----Global network properties (GNP) ----
###############################################

# Set up parallel processing & fixed seed
set.seed(42)
plan(multisession, workers = parallel::detectCores() - 1) 

# Iterate over all mat_lists
for (nm in names(main_list)) {
  mat_list <- main_list[[nm]]
  
  # Convert matrices to igraph object
  g_list <- map(mat_list, ~ graph_from_adjacency_matrix(.x, mode = "undirected", weighted = TRUE, diag = FALSE))
  
  # Convert weights to inverse distance (relevant for distance dependend GNPs)
  g_list_invdist <- map(g_list, ~ {E(.x)$weight <- 1 / E(.x)$weight; .x})
  
  # Compute global network properties
  eloc <- future_map_dbl(g_list_invdist, ~ average_local_efficiency(.x)) # Average local efficiency
  eglob <- future_map_dbl(g_list_invdist, ~ global_efficiency(.x)) # Global efficiency
  l <- future_map_dbl(g_list_invdist, ~ mean_distance(.x)) # Characteristic path length
  cc <- future_map_dbl(g_list, ~ avg_cc(.x)) # Clustering coefficient
  str <- future_map_dbl(g_list, ~ avg_strength(.x)) # Strength
  comm <- future_map(g_list, ~ cluster_louvain(.x), .options = furrr_options(seed = 42)) # Community object (get Q, N modules)
  q <- future_map2_dbl(comm, g_list, ~ modularity(.y, membership(.x), E(.y)$weight)) # Modularity
  n_mod <- future_map_dbl(comm, ~ length(.x)) # Number of modules
  comp <- future_map_dbl(g_list, ~ components(.x)$no)
  n_nodes <- future_map_dbl(g_list, ~ length(V(.x)))
  
  # Combine all GNPs
  fc_gnp <- tibble(
    participant_id = names(g_list), 
    cc = cc, 
    eglob = eglob,
    eloc = eloc, 
    l = l, 
    n_mod = n_mod, 
    q = q, 
    str = str,
    comp, n_nodes) |> 
    # Tidy participant_id 
    mutate(
      density = as.numeric(str_extract(participant_id, "(?<=_d)\\d+\\.?\\d*")),
      participant_id = str_extract(participant_id, "PREDICT\\d+"))
  
  # Compute AUC across densities
  fc_gnp_auc <- fc_gnp |> 
    pivot_longer(-c(participant_id, density), names_to = "gnp", values_to = "value") |> 
    summarise(auc = DescTools::AUC(x = density, y = value, method = "trapezoid"),
              .by = c(participant_id, gnp)) |> 
    pivot_wider(names_from = gnp, values_from = auc, names_prefix = "auc_")
  
  # Save results
  saveRDS(fc_gnp, paste0("03_output/gt/fc_", nm, "_w_gnp.rds"))
  saveRDS(fc_gnp_auc, paste0("03_output/gt/fc_", nm, "_w_gnp_auc.rds"))
  cat("Finished GNP calculation for:", nm)
  
  # Clean up
  rm(fc_gnp, fc_gnp_auc, g_list, rand_metrics, cc, eglob, eloc, l, n_mod, q, comm)
}
