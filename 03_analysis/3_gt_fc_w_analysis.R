# Author: Lionel Butry (@lbutry)
# Purpose: Group comparisons of FUNCTIONAL global network properties between
#           1) main groups (CG, LBP+, LBP-), and
#           2) subgroups (CG, ncLBP, cLBP)

source("02_scripts/03_analysis/set_env.R")
source("02_scripts/03_analysis/functions_graphtheory.R")
set.seed(42)

###### USER INPUT ######
subgroup <- TRUE # TRUE = CG, ncLBP, cLBP; FALSE = CG, LBP+, LBP-
net <- "full" # Choose connectome: full, dmn, smn, van
########################

############################################
# ---- Import clinical data ----
############################################

clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds")

if (subgroup) {
  clin_df <- clin_df |> filter(group %in% c("LBP+", "CG")) |>
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
}

###############################################
# Explore results
###############################################

# Read global network properties
fc_gnp <- readRDS(paste0("03_output/gt/fc_", net, "_w_gnp_auc.rds"))
skim(fc_gnp)

# Combine gnp data with clinical data
gnp_clin <- fc_gnp |> inner_join(clin_df)
skim(gnp_clin)

gnp_summary <- gnp_clin |>  
  select(group, auc_cc, auc_eglob, auc_eloc, auc_l, auc_n_mod, auc_q, auc_str, auc_btw, auc_close) |>
  mutate(across(c(-group, -auc_l), ~ .x * 1e2)) |> 
  group_by(group) |>
  get_summary_stats() |> 
  select(group, variable, n, mean, sd, ci) |> 
  mutate(value = paste(mean, "±", sd))

# Boxplots (AUC)
plot_auc(gnp_clin, "auc_cc", "CC")
plot_auc(gnp_clin, "auc_eglob", "Eglob")
plot_auc(gnp_clin, "auc_eloc", "Eloc")
plot_auc(gnp_clin, "auc_l", "L")
plot_auc(gnp_clin, "auc_q", "Q")

###############################################
# Group comparison
###############################################

# Set permutation matrix for all analysis
if (subgroup) {
  path_pmat <- "03_output/gt/P_mat_fc_subgroup.rds"
} else {
  path_pmat <- "03_output/gt/P_mat_fc.rds"
}

if (file.exists(path_pmat)) {
  P_mat <- readRDS(path_pmat)
} else {
  P_mat <- permuco::Pmat(np = 50000, n = nrow(gnp_clin))
  saveRDS(P_mat, path_pmat)
}

# Run permutation ANCOVA 
metrics <- c("auc_cc", "auc_eglob", "auc_eloc", "auc_l", "auc_q")
fc_aovperm_results <- map(metrics, function(metric) {
  formula <- as.formula(paste(metric, "~ group + Demo_Age + Demo_Gender + scanner"))
  permuco::aovperm(formula, data = gnp_clin, P = P_mat)
}); names(fc_aovperm_results) <- metrics

# Tidy results
fc_aovperm_results_df <- map_dfr(names(fc_aovperm_results), ~ {
  fc_aovperm_results[[.x]]$table |> 
    rownames_to_column(var = "variable") |> 
    mutate(
      model = .x, .before = 1,
      df.residual = fc_aovperm_results[[.x]]$df.residual
      )})

## apply FDR & compute effect size (partial eta^2)
fdr_controlled <- fc_aovperm_results_df |> 
  filter(variable == "group") |> 
  mutate(
    p.fdr = p.adjust(`resampled P(>F)`, method = "fdr"),
    eta2 = effectsize::F_to_eta2(
      f = `F`,
      df = df,
      df_error = df.residual))

# Save results
if (subgroup) {
  file <- paste0("03_output/gt/fc_w_result_PermANCOVA_subgroup_", net, ".csv")
} else {
  file <- paste0("03_output/gt/fc_w_result_PermANCOVA_", net, ".csv")
}

write.csv2(fdr_controlled, file)

###############################################
# Post-hoc tests; permANCOVA per group-pairs 
###############################################

# Helper function
post_hoc_permaov <- function(pairs, var_metric) {
  
  map_dfr(pairs, \(p) {
    
    dat <- gnp_clin |> 
      filter(group %in% p) |> 
      mutate(group = droplevels(group))
    
    fit <- permuco::aovperm(
      formula = as.formula(paste(var_metric, "~ group + Demo_Age + Demo_Gender + scanner")),
      data = dat,
      P = P_mat
    )
    
    tibble(
      contrast = paste(p, collapse = " vs "),
      F_perm   = fit$table["group", "F"],
      p_perm   = fit$table["group", "resampled P(>F)"],
      cohen.d = effectsize::F_to_d(
        f = fit$table$F[1],
        df = fit$table$df[1],
        df_error = fit$df.residual) 
    )
  })
}

#### Main analysis ####

pairs_list <- list(
  c("CG", "LBP-"),
  c("CG", "LBP+"),
  c("LBP-", "LBP+")
)

# VAN
post_hoc_permaov(pairs_list, "auc_l")
post_hoc_permaov(pairs_list, "auc_eglob")

#### Subgroup analysis ####

pairs_list <- list(
  c("CG", "ncLBP"),
  c("CG", "cLBP"),
  c("ncLBP", "cLBP")
)

# Full FC
post_hoc_permaov(pairs_list, "auc_q")

# VAN
post_hoc_permaov(pairs_list, "auc_eglob")
post_hoc_permaov(pairs_list, "auc_l")
