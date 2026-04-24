# Author: Lionel Butry (@lbutry)
# Purpose:
#   (A) Compute structure-function coupling (SFC)
#   (B) Inference stats: Group differences in global, network and regional SFC

source("02_scripts/03_analysis/set_env.R")
source("02_scripts/03_analysis/functions_neurocombat.R")

###############################################
# Import data
###############################################

# Clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds")

# Brainnetome LUT
bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv")

# SC & FC
source("02_scripts/03_analysis/import_fc.R")
source("02_scripts/03_analysis/import_sc.R")

###############################################
# Apply ComBat
###############################################

sc_mat_list <- apply_combat(
  sc_mat_list,
  clin_df |> filter(participant_id %in% names(sc_mat_list)),
  batch = "dwi_protocol",
  covars = c("Demo_Age", "Demo_Gender"))

# Set negative values to zero (cause not biologically meaningful for SIFT2 values)
sc_mat_list <- map(sc_mat_list, ~ pmax(.x, 0))

fc_mat_list <- apply_combat(
  fc_mat_list,
  clin_df |> filter(participant_id %in% names(fc_mat_list)),
  batch = "scanner",
  covars = c("Demo_Age", "Demo_Gender"))

###############################################
# Filter matrices & clinical data
###############################################

# Keep only matrices that have both SC and FC
common_ids <- intersect(names(fc_mat_list), names(sc_mat_list))
fc_mat_list <- fc_mat_list[common_ids]
sc_mat_list <- sc_mat_list[common_ids]

clin_df <- clin_df |> filter(participant_id %in% common_ids)

###############################################
# SFC: Multilinear regression
###############################################
library(brainGraph)

# Function
sfc_mlr <- function(sc_mat, fc_mat, atlas = "brainnetome") {
  
  # Convert matrix to igraph object & set weights from strength to inverse distance
  sc_g <- graph_from_adjacency_matrix(sc_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  E(sc_g)$weight <- 1 / E(sc_g)$weight
  
  # --- PREDICTORS ---
  # Path length
  pl_mat <- distances(sc_g, algorithm = "floyd-warshall")
  diag(pl_mat) <- NA
  
  # Communicability 
  com_mat <- brainGraph::communicability(sc_g)
  dimnames(com_mat) <- dimnames(sc_mat)
  diag(com_mat) <- NA
  
  # Euclidean distance
  ## Fully-connected SC
  full_mat <- sc_mat
  full_mat[full_mat == 0] <- 1
  full_g <- graph_from_adjacency_matrix(full_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  full_g <- set_graph_attr(full_g, "atlas", value = atlas)
  
  ## Get names of all nodes & edges
  nodes <- V(full_g)$name
  edges <- get.edges(full_g, 1:ecount(full_g))

  ## Compute ED & convert to matrix
  ed <- brainGraph::edge_spatial_dist(full_g)
  ed_mat <- matrix(NA, nrow = length(nodes), ncol = length(nodes), dimnames = list(nodes, nodes))
  ed_mat[cbind(edges[,1], edges[,2])] <- ed # upper tri
  ed_mat[cbind(edges[,2], edges[,1])] <- ed # lower tri
  diag(ed_mat) <- NA
  
  # SC
  diag(sc_mat) <- NA
  
  ### SFC per node (= matrix row) ###
  
  r2_vec <- map_dbl(nodes, function(i) {
    summary(lm(fc_mat[i,] ~ sc_mat[i,] + pl_mat[i,] + com_mat[i,] + ed_mat[i,]))$adj.r.squared
  })
  
  names(r2_vec) <- nodes
  return(r2_vec)
}

# Compute SFC
sfc <- map2_dfr(sc_mat_list, fc_mat_list, ~ sfc_mlr(.x, .y), .id = "participant_id") 
saveRDS(sfc, "03_output/sfc/sfc.rds")

###############################################
# Statistics - Region-wise
###############################################
sfc <- readRDS("03_output/sfc/sfc.rds")

# Join clinical & SFC data
sfc_clin <- sfc |>
  left_join(clin_df, by = "participant_id") 

# Get roi names
roi <- names(sfc)[-1]

#### Contrasts ####
# Define contrasts
levels(sfc_clin$group) 

contrasts <- list(
  "CG > LBP-" = c(1, -1, 0),
  "CG > LBP+" = c(1, 0, -1),
  "LBP- > LBP+" = c(0, 1, -1)
)

# Run contrast per roi
results_region_contrasts <- map_dfr(roi, function(region) {
  formula <- as.formula(paste(region, "~ group + Demo_Age + Demo_Gender + scanner"))
  model <- lm(formula, data = sfc_clin)
  emm <- emmeans(model, ~ group)
  emmeans::contrast(emm, method = contrasts) |> 
    broom::tidy() |> 
    mutate(roi = region, .before = term) |> 
    mutate(
      cohen.d = effectsize::t_to_d(
        t = statistic, df_error = df
    ))
})

# Apply FDR per contrast
contrast_cg_lbpm <- results_region_contrasts |> 
  filter(contrast == "CG > LBP-") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_cg_lbpm |> filter(p.value < 0.05) 
contrast_cg_lbpm |> filter(p.fdr < 0.05)
write.csv2(contrast_cg_lbpm, "03_output/sfc/results_roi_contr_cg>lbp-.csv", row.names = FALSE)

contrast_cg_lbpp <- results_region_contrasts |> 
  filter(contrast == "CG > LBP+") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_cg_lbpp |> filter(p.value < 0.05)
contrast_cg_lbpp |> filter(p.fdr < 0.05)
write.csv2(contrast_cg_lbpp, "03_output/sfc/results_roi_contr_cg>lbp+.csv", row.names = FALSE)

contrast_lbpm_lbpp <- results_region_contrasts |> 
  filter(contrast == "LBP- > LBP+") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_lbpm_lbpp |> filter(p.value < 0.05)
contrast_lbpm_lbpp |> filter(p.fdr < 0.05)
write.csv2(contrast_lbpm_lbpp, "03_output/sfc/results_roi_contr_lbp->lbp+.csv", row.names = FALSE)

#### Subgroup analysis ####
sfc_clin_sub <- sfc_clin |> filter(group %in% c("LBP+", "CG")) |>
  mutate(
    group_old = group,
    group = case_when(
      group == "CG" ~ "CG",
      pain_duration_cat == "chronic" ~ "cLBP",
      pain_duration_cat %in% c("subacute", "acute") ~ "ncLBP",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(group)) |>
  mutate(group = factor(group, levels = c("CG", "ncLBP", "cLBP")))

# Define contrasts
levels(sfc_clin_sub$group) # sanity check: factor order

contrasts <- list(
  "CG > ncLBP" = c(1, -1, 0),
  "CG > cLBP" = c(1, 0, -1),
  "ncLBP > cLBP" = c(0, 1, -1)
)

# Run contrast per roi
results_region_contrasts_sub <- map_dfr(roi, function(region) {
  formula <- as.formula(paste(region, "~ group + Demo_Age + Demo_Gender + scanner"))
  model <- lm(formula, data = sfc_clin_sub)
  emm <- emmeans(model, ~ group)
  emmeans::contrast(emm, method = contrasts) |> 
    broom::tidy() |> 
    mutate(roi = region, .before = term) |> 
    mutate(
      cohen.d = effectsize::t_to_d(
        t = statistic, df_error = df
      ))
})

# Apply FDR
contrast_cg_nclbp <- results_region_contrasts_sub |> 
  filter(contrast == "CG > ncLBP") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_cg_nclbp |> filter(p.value < 0.05) 
contrast_cg_nclbp |> filter(p.fdr < 0.05) 
write.csv2(contrast_cg_nclbp, "03_output/sfc/results_roi_contr_cg>nclbp.csv", row.names = FALSE)

contrast_cg_clbp <- results_region_contrasts_sub |> 
  filter(contrast == "CG > cLBP") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_cg_clbp |> filter(p.value < 0.05) 
contrast_cg_clbp |> filter(p.fdr < 0.05)
write.csv2(contrast_cg_clbp, "03_output/sfc/results_roi_contr_cg>clbp.csv", row.names = FALSE)

contrast_nclbp_clbp <- results_region_contrasts_sub |> 
  filter(contrast == "ncLBP > cLBP") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))
contrast_nclbp_clbp |> filter(p.value < 0.05) |> pull(roi)
contrast_nclbp_clbp |> filter(p.fdr < 0.05)
write.csv2(contrast_nclbp_clbp, "03_output/sfc/results_roi_contr_ncLBP>clbp.csv", row.names = FALSE)

###############################################
# Statistics - Global & Network
###############################################

# Compute global & network SFC
sfc_long <- sfc |> 
  pivot_longer(-participant_id, names_to = "ROI.Name", values_to = "sfc") |> 
  left_join(bn_lut)

sfc_network <- sfc_long |> 
  group_by(participant_id, network) |> 
  summarise(sfc = mean(sfc)) |> 
  pivot_wider(names_from = "network", values_from = "sfc")

sfc_global <- sfc_long |> 
  group_by(participant_id) |> 
  summarise(sfc = mean(sfc)) |> 
  mutate(col = "global") |> 
  pivot_wider(names_from = "col", values_from = "sfc")

sfc_agg <- sfc_global |> 
  left_join(sfc_network) |> 
  left_join(clin_df)

saveRDS(sfc_agg, "03_output/sfc/sfc_agg.rds")
rm(sfc_network, sfc_global)

### Main analysis: ANCOVA ###
metrics <- c("global", "DMN", "VAN", "SMN")

sfc_agg_lm <- map(metrics, function(metric) {
  formula <- as.formula(paste(metric, "~ group + Demo_Age + Demo_Gender + scanner"))
  lm(formula, data = sfc_agg)
}); names(sfc_agg_lm) <- metrics

# Tidy results & apply FDR
results_agg_aov <- map_dfr(names(sfc_agg_lm), function(metric) {
  anova(sfc_agg_lm[[metric]]) |> broom::tidy() |> 
    mutate(
      metric = metric,
      eta2 = effectsize::F_to_eta2(
        f = statistic, df = df, df_error = df[term == "Residuals"]
      ))
  }) |> 
  filter(term == "group") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr"))

results_agg_aov |> filter(p.value < 0.05)
results_agg_aov |> filter(p.fdr < 0.05)

write.csv2(results_agg_aov, "03_output/sfc/results_agg_aov.csv", row.names = FALSE)

# Post-hoc tests
emm <- emmeans(sfc_agg_lm$global, spec = "group")
hoc_global <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_global, "03_output/sfc/results_agg_hoc_global.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm$DMN, spec = "group")
hoc_dmn <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_dmn, "03_output/sfc/results_agg_hoc_dmn.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm$VAN, spec = "group")
hoc_van <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_van, "03_output/sfc/results_agg_hoc_van.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm$SMN, spec = "group")
hoc_smn <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_smn, "03_output/sfc/results_agg_hoc_smn.csv", row.names = FALSE)

### Subgroup analysis: ANCOVA ###
sfc_agg_sub <- sfc_agg |> filter(group %in% c("LBP+", "CG")) |>
  mutate(
    group_old = group,
    group = case_when(
      group == "CG" ~ "CG",
      pain_duration_cat == "chronic" ~ "cLBP",
      pain_duration_cat %in% c("subacute", "acute") ~ "ncLBP",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(group)) |>
  mutate(group = factor(group, levels = c("CG", "ncLBP", "cLBP")))

sfc_agg_lm_sub <- map(metrics, function(metric) {
  formula <- as.formula(paste(metric, "~ group + Demo_Age + Demo_Gender + scanner"))
  lm(formula, data = sfc_agg_sub)
}); names(sfc_agg_lm_sub) <- metrics

# Tidy results & apply FDR
results_agg_aov_sub <- map_dfr(names(sfc_agg_lm_sub), function(metric) {
  anova(sfc_agg_lm_sub[[metric]]) |> broom::tidy() |> 
    mutate(
      metric = metric,
      eta2 = effectsize::F_to_eta2(
        f = statistic, df = df, df_error = df[term == "Residuals"]
      ))
}) |> 
  filter(term == "group") |> 
  mutate(p.fdr = p.adjust(p.value, method = "fdr")) 

results_agg_aov_sub |> filter(p.value < 0.05)
results_agg_aov_sub |> filter(p.fdr < 0.05)

write.csv2(results_agg_aov_sub, "03_output/sfc/results_agg_aov_sub.csv", row.names = FALSE)

# Post-hoc tests
emm <- emmeans(sfc_agg_lm_sub$global, spec = "group")
hoc_global <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_global, "03_output/sfc/results_agg_hoc_global_sub.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm_sub$DMN, spec = "group")
hoc_dmn <- pairs(emm, adjust = "none") |>
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_dmn, "03_output/sfc/results_agg_hoc_dmn_sub.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm_sub$VAN, spec = "group")
hoc_van <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_van, "03_output/sfc/results_agg_hoc_van_sub.csv", row.names = FALSE)

emm <- emmeans(sfc_agg_lm_sub$SMN, spec = "group")
hoc_smn <- pairs(emm, adjust = "none") |> 
  broom::tidy() |> 
  mutate(cohen.d = effectsize::t_to_d(t = statistic, df_error = df))
write.csv2(hoc_smn, "03_output/sfc/results_agg_hoc_smn_sub.csv", row.names = FALSE)

###############################################
# Descriptives
###############################################

mean(sfc_agg$global)
sd(sfc_agg$global)

# SFC per group & region (main)
table_network <- sfc_agg |>
  summarise(
    across(
      c(global, DMN, VAN, SMN),
      ~sprintf("%.3f ± %.3f", mean(.x), sd(.x))),
      .by = "group"
  )

table_region <- sfc_clin |>
  summarise(
    across(
      matches("^(SFG|MFG|IFG|OrG|PrG|PCL|STG|MTG|ITG|FuG|PhG|pSTS|SPL|IPL|PCun|PoG|INS|CG|MVO|LOcC|Amyg|Hipp|BG|Tha)"),
      ~sprintf("%.3f ± %.3f", mean(.x), sd(.x))),
    .by = "group"
  )

table_sfc_main <- left_join(table_network, table_region) |> 
  column_to_rownames("group") |>
  t() |>
  as.data.frame() |>
  rownames_to_column("SFC")

# SFC per group & region (sub)
table_network <- sfc_agg_sub |>
  summarise(
    across(
      c(global, DMN, VAN, SMN),
      ~sprintf("%.3f ± %.3f", mean(.x), sd(.x))),
    .by = "group"
  )

table_region <- sfc_clin_sub |>
  summarise(
    across(
      matches("^(SFG|MFG|IFG|OrG|PrG|PCL|STG|MTG|ITG|FuG|PhG|pSTS|SPL|IPL|PCun|PoG|INS|CG|MVO|LOcC|Amyg|Hipp|BG|Tha)"),
      ~sprintf("%.3f ± %.3f", mean(.x), sd(.x))),
    .by = "group"
  )

table_sfc_sub <- left_join(table_network, table_region) |> 
  column_to_rownames("group") |>
  t() |>
  as.data.frame() |>
  rownames_to_column("SFC")

table_sfc_all <- left_join(table_sfc_main, table_sfc_sub) |> 
  select(SFC, CG, `LBP-`, `LBP+`, ncLBP, cLBP)

writexl::write_xlsx(table_sfc_all, "03_output/sfc/suppl_table_sfc.xlsx")
