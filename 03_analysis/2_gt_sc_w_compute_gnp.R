# Author: Lionel Butry (@lbutry)
# Purpose: Computes global network properties of the structural connectome
#          for full, DMN, VAN & SMN

source("02_scripts/03_analysis/set_env.R")
source("02_scripts/03_analysis/functions_neurocombat.R")
source("02_scripts/03_analysis/functions_graphtheory.R")

############################################
# ---- Import data ----
############################################

# SC
source("02_scripts/03_analysis/import_sc.R")
mat_list <- sc_mat_list
rm(sc_mat_list)

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
  batch = "dwi_protocol",
  covars = c("Demo_Age", "Demo_Gender"))

# Set negative values to zero (cause not biologically meaningful for SIFT2 values)
mat_list_combat <- map(mat_list_combat, ~ pmax(.x, 0))

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

###############################################
# ----Global network properties (GNP) ----
###############################################

# Set up parallel processing & fixed seed
set.seed(42)
plan(multisession, workers = parallel::detectCores() - 1) 

main_list <- list(
  full = mat_list_combat,
  cen = mat_cen_list,
  dan = mat_dan_list,
  dmn = mat_dmn_list,
  lim = mat_lim_list,
  smn = mat_smn_list,
  subc = mat_subc_list,
  van = mat_van_list
)

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
  sc_gnp <- tibble(
    participant_id = names(g_list), 
    eloc = eloc, 
    eglob = eglob,
    l = l, 
    cc = cc, 
    str = str,
    n_mod = n_mod, 
    q = q, 
    comp, 
    n_nodes)
  
  # Save results
  saveRDS(sc_gnp, paste0("03_output/gt/sc_", nm, "_gnp_w.rds"))
}

# Clean up
rm(sc_gnp, g_list, cc, eglob, eloc, l, n_mod, q, comm, str)
