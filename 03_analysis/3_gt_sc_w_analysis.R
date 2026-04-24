# Author: Lionel Butry (@lbutry)
# Purpose: Group comparisons of STRUCTURAL global network properties between
#           1) main groups (CG, LBP+, LBP-), and
#           2) subgroups (CG, ncLBP, cLBP)

source("02_scripts/03_analysis/set_env.R")
source("02_scripts/03_analysis/functions_graphtheory.R")
set.seed(42)

###### USER INPUT ######
subgroup <- TRUE # TRUE = CG, ncLBP, cLBP; FALSE = CG, LBP+, LBP-
net <- "dmn" # Choose connectome: full, dmn, smn, van
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
sc_gnp <- readRDS(paste0("03_output/gt/sc_w_", net, "_gnp.rds"))

# Combine gnp data with clinical data
gnp_clin <- sc_gnp |> inner_join(clin_df)
skim(gnp_clin)

gnp_summary <- gnp_clin |>  
  select(group, cc, eglob, eloc, l, n_mod, q, str, btw, close) |> 
  mutate(l = l * 1e3) |> 
  group_by(group) |>
  get_summary_stats() |> 
  select(group, variable, n, mean, sd, ci) |> 
  mutate(value = paste(mean, "±", sd))

# Boxplots (AUC)
plot_auc(gnp_clin, "cc", "CC")
plot_auc(gnp_clin, "eglob", "Eglob")
plot_auc(gnp_clin, "eloc", "Eloc")
plot_auc(gnp_clin, "l", "L")
plot_auc(gnp_clin, "q", "Q")

###############################################
# Group comparison
###############################################

# Set permutation matrix for all analysis
if (subgroup) {
  path_pmat <- "03_output/gt/P_mat_sc_subgroup.rds"
} else {
  path_pmat <- "03_output/gt/P_mat_sc.rds"
}

if (file.exists(path_pmat)) {
  P_mat <- readRDS(path_pmat)
} else {
  P_mat <- permuco::Pmat(np = 50000, n = nrow(gnp_clin))
  saveRDS(P_mat, path_pmat)
}

# Run permutation ANCOVA 
metrics <- c("cc", "eglob", "eloc", "l", "q")
sc_aovperm_results <- map(metrics, function(metric) {
  formula <- as.formula(paste(metric, "~ group + Demo_Age + Demo_Gender + dwi_protocol"))
  permuco::aovperm(formula, data = gnp_clin, P = P_mat)
}); names(sc_aovperm_results) <- metrics

# Tidy results
sc_aovperm_results_df <- map_dfr(names(sc_aovperm_results), ~ {
  sc_aovperm_results[[.x]]$table |> 
    rownames_to_column(var = "variable") |> 
    mutate(model = .x, .before = 1)})

## apply FDR
fdr_controlled <- sc_aovperm_results_df |> 
  filter(variable == "group") |> 
  mutate(p.fdr = p.adjust(`resampled P(>F)`, method = "fdr"))

if (subgroup) {
  file <- paste0("03_output/gt/sc_w_result_PermANCOVA_subgroup_", net, ".csv")
} else {
  file <- paste0("03_output/gt/sc_w_result_PermANCOVA_", net, ".csv")
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
      formula = as.formula(paste(var_metric, "~ group + Demo_Age + Demo_Gender + dwi_protocol")),
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

# none significant results

#### Subgroup analysis ####
pairs_list <- list(
  c("CG", "ncLBP"),
  c("CG", "cLBP"),
  c("ncLBP", "cLBP")
)

# DMN
post_hoc_permaov(pairs_list, "l")
post_hoc_permaov(pairs_list, "q")
