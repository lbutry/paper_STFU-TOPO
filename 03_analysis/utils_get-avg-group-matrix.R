# Author: Lionel Butry (@lbutry)
# Purpose: Helper script; computes accurate average group matrices for visualisation purposes

source("02_scripts/03_analysis/functions_neurocombat.R")
source("02_scripts/03_analysis/functions_graphtheory.R")

############################################
# Average group matrix
############################################

# Load clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") 

# Get IDs of each group
id_all <- clin_df |> pull(participant_id)
id_cg <- clin_df |> filter(group == "CG") |> pull(participant_id)
id_lbpp <- clin_df |> filter(group == "LBP+") |> pull(participant_id)
id_lbpm <- clin_df |> filter(group == "LBP-") |> pull(participant_id)
id_clbp <- clin_df |> filter(group == "LBP+", pain_duration_cat == "chronic") |> pull(participant_id)
id_nclbp <- clin_df |> filter(group == "LBP+", pain_duration_cat %in% c("acute", "subacute")) |> pull(participant_id)

# Import Brainnetome LUT
bn_lut <- read.delim2("01_data/ready4analysis/BN_Atlas_246_LUT.tsv")

########## FC ##########

# Import raw FC
source("02_scripts/03_analysis/import_fc.R")

# Apply Combat
fc_mat_list_combat <- apply_combat(
  fc_mat_list,
  clin_df,
  batch = "scanner",
  covars = c("Demo_Age", "Demo_Gender"))

# Compute average matrix per group
fc_avg_mat_all <- Reduce("+", fc_mat_list_combat) / length(fc_mat_list_combat)

fc_avg_mat_cg <- fc_mat_list_combat[names(fc_mat_list_combat) %in% id_cg]
fc_avg_mat_cg <- Reduce("+", fc_avg_mat_cg) / length(fc_avg_mat_cg)

fc_avg_mat_lbpp <- fc_mat_list_combat[names(fc_mat_list_combat) %in% id_lbpp]
fc_avg_mat_lbpp <- Reduce("+", fc_avg_mat_lbpp) / length(fc_avg_mat_lbpp)

fc_avg_mat_lbpm <- fc_mat_list_combat[names(fc_mat_list_combat) %in% id_lbpm]
fc_avg_mat_lbpm <- Reduce("+", fc_avg_mat_lbpm) / length(fc_avg_mat_lbpm)

fc_avg_mat_clbp <- fc_mat_list_combat[names(fc_mat_list_combat) %in% id_clbp]
fc_avg_mat_clbp <- Reduce("+", fc_avg_mat_clbp) / length(fc_avg_mat_clbp)

fc_avg_mat_nclbp <- fc_mat_list_combat[names(fc_mat_list_combat) %in% id_nclbp]
fc_avg_mat_nclbp <- Reduce("+", fc_avg_mat_nclbp) / length(fc_avg_mat_nclbp)

fc_avg <- list(all = fc_avg_mat_all, CG = fc_avg_mat_cg, `LBP-` =  fc_avg_mat_lbpm, `LBP+` = fc_avg_mat_lbpp, cLBP = fc_avg_mat_clbp, ncLBP = fc_avg_mat_nclbp)
rm(fc_avg_mat_all, fc_avg_mat_cg, fc_avg_mat_lbpp, fc_avg_mat_lbpm, fc_avg_mat_clbp, fc_avg_mat_nclbp)

# Split connectome into subnetworks
rois_van <- bn_lut$ROI.Name[bn_lut$network == "VAN"]
fc_avg_van <- map(fc_avg, ~ .x[rois_van, rois_van])

rois_dmn <- bn_lut$ROI.Name[bn_lut$network == "DMN"]
fc_avg_dmn <- map(fc_avg, ~ .x[rois_dmn, rois_dmn])

rois_smn <- bn_lut$ROI.Name[bn_lut$network == "SMN"]
fc_avg_smn <- map(fc_avg, ~ .x[rois_smn, rois_smn])

rois_3 <- c(rois_van, rois_dmn, rois_smn)
fc_avg_3 <- map(fc_avg, ~ .x[rois_3, rois_3])

# Thresholding at 15%
gthr = 0.15
fc_avg_abs <- map(fc_avg, ~ abs(.x))
ranks_list <- map(fc_avg_abs, ~ mst_gthr_order(.x))
fc_avg_thr <- map2(fc_avg_abs, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

fc_avg_abs_van <- map(fc_avg_van, ~ abs(.x))
ranks_list <- map(fc_avg_abs_van, ~ mst_gthr_order(.x))
fc_avg_thr_van <- map2(fc_avg_abs_van, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

fc_avg_abs_dmn <- map(fc_avg_dmn, ~ abs(.x))
ranks_list <- map(fc_avg_abs_dmn, ~ mst_gthr_order(.x))
fc_avg_thr_dmn <- map2(fc_avg_abs_dmn, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

fc_avg_abs_smn <- map(fc_avg_smn, ~ abs(.x))
ranks_list <- map(fc_avg_abs_smn, ~ mst_gthr_order(.x))
fc_avg_thr_smn <- map2(fc_avg_abs_smn, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

fc_avg_abs_3 <- map(fc_avg_3, ~ abs(.x))
ranks_list <- map(fc_avg_abs_3, ~ mst_gthr_order(.x))
fc_avg_thr_3 <- map2(fc_avg_abs_3, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

fc_avg_thr_list <- list(all = fc_avg_thr, van = fc_avg_thr_van, dmn = fc_avg_thr_dmn, smn = fc_avg_thr_smn, net3 = fc_avg_thr_3)

rm(fc_avg_abs, fc_avg_abs_dmn, fc_avg_abs_smn, fc_avg_abs_van, fc_avg_dmn, fc_avg_smn, fc_avg_thr, fc_avg_thr_dmn,
   fc_avg_thr_smn, fc_avg_thr_van, ranks_list, fc_mat_list_combat, fc_avg_van, fc_avg_abs_3, fc_avg_thr_3)

########## SC ##########

# Import raw SC
source("02_scripts/03_analysis/import_sc.R")

# Apply Combat
sc_mat_list_combat <- apply_combat(
  sc_mat_list,
  clin_df,
  batch = "dwi_protocol",
  covars = c("Demo_Age", "Demo_Gender"))

# Set negative values to zero (cause not biologically meaningful for SIFT2 values)
sc_mat_list_combat <- map(sc_mat_list_combat, ~ pmax(.x, 0))

# Import Brainnetome LUT
bn_lut <- read.delim2("01_data/ready4analysis/BN_Atlas_246_LUT.tsv") |>
  mutate(network = case_when(
    network == "Subcortical" ~ "SUBC",
    network == "Limbic" ~ "LIM",
    .default = network
  ))

# Compute average matrix per group
sc_avg_mat_all <- Reduce("+", sc_mat_list_combat) / length(sc_mat_list_combat)

sc_avg_mat_cg <- sc_mat_list_combat[names(sc_mat_list_combat) %in% id_cg]
sc_avg_mat_cg <- Reduce("+", sc_avg_mat_cg) / length(sc_avg_mat_cg)

sc_avg_mat_lbpp <- sc_mat_list_combat[names(sc_mat_list_combat) %in% id_lbpp]
sc_avg_mat_lbpp <- Reduce("+", sc_avg_mat_lbpp) / length(sc_avg_mat_lbpp)

sc_avg_mat_lbpm <- sc_mat_list_combat[names(sc_mat_list_combat) %in% id_lbpm]
sc_avg_mat_lbpm <- Reduce("+", sc_avg_mat_lbpm) / length(sc_avg_mat_lbpm)

sc_avg_mat_clbp <- sc_mat_list_combat[names(sc_mat_list_combat) %in% id_clbp]
sc_avg_mat_clbp <- Reduce("+", sc_avg_mat_clbp) / length(sc_avg_mat_clbp)

sc_avg_mat_nclbp <- sc_mat_list_combat[names(sc_mat_list_combat) %in% id_nclbp]
sc_avg_mat_nclbp <- Reduce("+", sc_avg_mat_nclbp) / length(sc_avg_mat_nclbp)

sc_avg <- list(all = sc_avg_mat_all, CG = sc_avg_mat_cg, `LBP-` =  sc_avg_mat_lbpm, `LBP+` = sc_avg_mat_lbpp, cLBP = sc_avg_mat_clbp, ncLBP = sc_avg_mat_nclbp)
rm(sc_avg_mat_all, sc_avg_mat_cg, sc_avg_mat_lbpp, sc_avg_mat_lbpm, sc_avg_mat_clbp, sc_avg_mat_nclbp)

# Split connectome into subnetworks
rois_van <- bn_lut$ROI.Name[bn_lut$network == "VAN"]
sc_avg_van <- map(sc_avg, ~ .x[rois_van, rois_van])

rois_dmn <- bn_lut$ROI.Name[bn_lut$network == "DMN"]
sc_avg_dmn <- map(sc_avg, ~ .x[rois_dmn, rois_dmn])

rois_smn <- bn_lut$ROI.Name[bn_lut$network == "SMN"]
sc_avg_smn <- map(sc_avg, ~ .x[rois_smn, rois_smn])

# Thresholding at 15%
gthr = 0.15
sc_avg_abs <- map(sc_avg, ~ .x)
ranks_list <- map(sc_avg_abs, ~ mst_gthr_order(.x))
sc_avg_thr <- map2(sc_avg_abs, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

sc_avg_abs_van <- map(sc_avg_van, ~ .x)
ranks_list <- map(sc_avg_abs_van, ~ mst_gthr_order(.x))
sc_avg_thr_van <- map2(sc_avg_abs_van, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

sc_avg_abs_dmn <- map(sc_avg_dmn, ~ .x)
ranks_list <- map(sc_avg_abs_dmn, ~ mst_gthr_order(.x))
sc_avg_thr_dmn <- map2(sc_avg_abs_dmn, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

sc_avg_abs_smn <- map(sc_avg_smn, ~ .x)
ranks_list <- map(sc_avg_abs_smn, ~ mst_gthr_order(.x))
sc_avg_thr_smn <- map2(sc_avg_abs_smn, ranks_list, ~ apply_mst_gthr(.x, .y, gthr))

sc_avg_thr_list <- list(all = sc_avg_thr, van = sc_avg_thr_van, dmn = sc_avg_thr_dmn, smn = sc_avg_thr_smn)

rm(sc_avg_abs, sc_avg_abs_dmn, sc_avg_abs_smn, sc_avg_abs_van, sc_avg_dmn, sc_avg_smn, sc_avg_thr, sc_avg_thr_dmn,
   sc_avg_thr_smn, sc_avg_thr_van, ranks_list, sc_mat_list_combat, sc_avg_van)

rm(gthr, id_all, id_cg, id_clbp, id_lbpp, id_nclbp, id_lbpm, rois_dmn, rois_smn, rois_van)
