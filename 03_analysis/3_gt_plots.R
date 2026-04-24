# Author: Lionel Butry (@lbutry)
# Purpose: Plot the results of the graph theory analysis

######
# Overview of significant results:
# - SC: subgroup: DMN - L
# - SC: subgroup: DMN - Q
# - FC: all: VAN - Eglob
# - FC: all: VAN - L
# - FC: subgroup: VAN - Eglob
# - FC: subgroup: VAN - L
# - FC: subgroup: Full - Q
#
# => FC: VAN-L, VAN-Eglob, Full-Q
# => SC: DMN-L, DMN-Q
######

source("02_scripts/03_analysis/set_env.R")

############################################
# Plot GNP: Compare means (adjusted for age, sex, scanner)
############################################

# Clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") 

# Import significant GNPs
sc_dmn_l_q <- readRDS("03_output/gt/sc_w_dmn_gnp.rds") |> 
  select(participant_id, l, q) |> 
  rename(
    sc_dmn_l = l,
    sc_dmn_q = q,
  )

fc_van_l_eglob <- readRDS("03_output/gt/fc_van_w_gnp_auc.rds") |> 
  select(participant_id, auc_l, auc_eglob) |> 
  rename(
    fc_van_l = auc_l,
    fc_van_eglob = auc_eglob)

fc_full_q <- readRDS("03_output/gt/fc_full_w_gnp_auc.rds") |> 
  select(participant_id, auc_q) |> 
  rename(fc_full_q = auc_q)

all_sig_metrics <- fc_van_l_eglob  |> 
  left_join(fc_full_q) |> 
  left_join(sc_dmn_l_q) |> 
  left_join(clin_df)

rm(sc_dmn_l_q, fc_full_q, fc_van_l_eglob)

# Plotting function
my_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      text = element_text(size = 12, family = "Helvetica", color = "black"),
      line = element_line(linewidth = 1),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 0.75),
      panel.grid = element_blank(),
      axis.ticks.length = unit(1.5, "mm"),
      plot.tag = element_text(face = "bold")
    )
}


plot_adj <- function(data, var_metric, y_label, color = c("#0072B2", "#CC79A7", "#E69F00")) {
  
  # Adjust data by age, sex & scanner
  formula <- as.formula(paste(var_metric, "~ group + Demo_Age + Demo_Gender + scanner"))
  model <- lm(formula, data = data, na.action = na.exclude)
  emm <- emmeans::emmeans(model, "group") |> as_tibble()
  
  data$predicted <- predict(model)
  
  # Define the pairs to connect
  line_segments <- tibble::tibble(
    x = c(1, 2),                  # positions of CG -> LBP-, LBP- -> LBP+
    xend = c(2, 3),
    y = emm$emmean[c(1,2)],
    yend = emm$emmean[c(2,3)]
  )
  
  # Bar plot
  ggplot() +
    geom_jitter(data = data,
                aes(x = group, y = .data[[var_metric]], color = group),
                #aes(x = group, y = predicted, color = group),
                width = 0.15, alpha = 0.15) +
    geom_point(data = emm,
               aes(x = group, y = emmean),
               size = 2) +
    geom_errorbar(data = emm,
                  aes(x = group, ymin = lower.CL, ymax = upper.CL),
                  width = 0.15) +
    geom_segment(data = line_segments,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linetype = 3) +
    scale_color_manual(values = color) +
    labs(y = y_label) +
    my_theme() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank())
}

# Results from main analysis
sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("**", "ns", "ns"), y.position = c(0.01485, 0.0142, 0.0145))
p_fvan_eglob <- plot_adj(all_sig_metrics, "fc_van_eglob", expression(E[glob]~"of fVAN")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("*", "ns", "ns"), y.position = c(13, 12.7, 12.5))
p_fvan_l <- plot_adj(all_sig_metrics, "fc_van_l", "L of fVAN") +
  scale_y_continuous(limits = c(8, 13)) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

# Results from subgroup analysis
all_sig_metrics_subgroup <- all_sig_metrics |> filter(group %in% c("LBP+", "CG")) |>
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

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("**", "*", "ns"), y.position = c(.0359, .03565, .0354))
p_fc_q_sub <- plot_adj(all_sig_metrics_subgroup, "fc_full_q", "Q of FC", color = c("#0072B2", "#66A61E", "#C44500")) +
  scale_y_continuous(limits = c(0.03, 0.0359)) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("*", "ns", "ns"), y.position = c(.015, .0142, .0146))
p_fvan_eglob_sub <- plot_adj(all_sig_metrics_subgroup, "fc_van_eglob", expression(E[glob]~"of fVAN"), color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)
  
sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("*", "*", "ns"), y.position = c(12.8, 12.5, 12.2))
p_fvan_l_sub <- plot_adj(all_sig_metrics_subgroup, "fc_van_l", "L of fVAN", color = c("#0072B2", "#66A61E", "#C44500")) + 
  scale_y_continuous(limits = c(8, 12.8)) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("**", "ns", "*"), y.position = c(0.025, 0.023, 0.024))
p_sdmn_l_sub <- plot_adj(all_sig_metrics_subgroup, "sc_dmn_l", "L of sDMN", color = c("#0072B2", "#66A61E", "#C44500")) +
  scale_y_continuous(limits = c(0.01, 0.026)) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("**", "ns", "ns"), y.position = c(0.7425, 0.735, 0.73))
p_sdmn_q_sub <- plot_adj(all_sig_metrics_subgroup, "sc_dmn_q", "Q of sDMN", color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4)

############################################
# Plot graph on glass-brain
############################################

source("02_scripts/03_analysis/utils_get-avg-group-matrix.R")

# Define color for the networks
color_network <- pal_igv("default")(16)[c(2:8, 16)] |> str_sub(1, -3)
names(color_network) <- c("CEN", "DAN", "DMN", "LIM", "SMN", "SUBC", "VAN", "VIS")

node_size = 1
edge_size = 0.5

# Overall
## All groups: fSMN
brain_all_fsmn <- brainconn(
    atlas = bn_lut |> filter(network == "SMN"),
    conmat = fc_avg_thr_list$smn$all,
    view = c("top"),
    node.color = color_network["SMN"],
    edge.width = edge_size,
    edge.color.weighted = T,
    scale.edge.width = F,
    node.size = node_size,
    show.legend = FALSE,
    background.alpha = 0.8) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## All groups: fDMN
brain_all_fdmn <- brainconn(
    atlas = bn_lut |> filter(network == "DMN"),
    conmat = fc_avg_thr_list$dmn$all,
    view = c("left"),
    node.color = color_network["DMN"],
    edge.width = edge_size,
    edge.color.weighted = T,
    scale.edge.width = F,
    node.size = node_size,
    show.legend = FALSE,
    background.alpha = 0.8) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## All groups: fVAN
brain_all_fvan <- brainconn(
    atlas = bn_lut |> filter(network == "VAN"),
    conmat = fc_avg_thr_list$van$all,
    view = c("left"),
    node.color = color_network["VAN"],
    edge.width = edge_size,
    edge.color.weighted = T,
    scale.edge.width = F,
    node.size = node_size,
    show.legend = FALSE,
    background.alpha = 0.8) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## All groups: Full FC
brain_all_fc <- brainconn(
    atlas = bn_lut,
    conmat = fc_avg_thr_list$all$all,
    view = c("left"),
    node.color = "network",
    edge.width = edge_size,
    edge.color.weighted = T,
    scale.edge.width = FALSE,
    node.size = node_size,
    show.legend = FALSE,
    background.alpha = 0.8) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05)) +
  scale_colour_manual(values = color_network)

# VAN per group
## fVAN - CG
brain_fvan_cg <- brainconn(
  atlas = bn_lut |> filter(network == "VAN"),
  conmat = fc_avg_thr_list$van$CG,
  view = c("left"),
  node.color = "#0072B2",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## fVAN - LBP+
brain_fvan_lbpp <- brainconn(
  atlas = bn_lut |> filter(network == "VAN"),
  conmat = fc_avg_thr_list$van$`LBP+`,
  view = c("left"),
  node.color = "#E69F00",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## fVAN - LBP-
brain_fvan_lbpm <- brainconn(
  atlas = bn_lut |> filter(network == "VAN"),
  conmat = fc_avg_thr_list$van$`LBP-`,
  view = c("left"),
  node.color = "#CC79A7",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## fVAN - cLBP
brain_fvan_clbp <- brainconn(
  atlas = bn_lut |> filter(network == "VAN"),
  conmat = fc_avg_thr_list$van$cLBP,
  view = c("left"),
  node.color = "#C44500",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## fVAN - ncLBP
brain_fvan_nclbp <- brainconn(
  atlas = bn_lut |> filter(network == "VAN"),
  conmat = fc_avg_thr_list$van$ncLBP,
  view = c("left"),
  node.color = "#66A61E",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

# Full FC per group
## Full FC - CG
brain_fc_cg <- brainconn(
  atlas = bn_lut,
  conmat = fc_avg_thr_list$all$CG,
  view = c("left"),
  node.color = "#0072B2",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## Full FC - cLBP
brain_fc_clbp <- brainconn(
  atlas = bn_lut,
  conmat = fc_avg_thr_list$all$cLBP,
  view = c("left"),
  node.color = "#C44500",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

## Full FC - ncLBP
brain_fc_nclbp <- brainconn(
  atlas = bn_lut,
  conmat = fc_avg_thr_list$all$ncLBP,
  view = c("left"),
  node.color = "#66A61E",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 0.05))

# sDMN per group
## sDMN - CG
brain_sdmn_cg <- brainconn(
  atlas = bn_lut |> filter(network == "DMN"),
  conmat = sc_avg_thr_list$dmn$CG,
  view = c("left"),
  node.color = "#0072B2",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 5000))

## sDMN - cLBP
brain_sdmn_clbp <- brainconn(
  atlas = bn_lut |> filter(network == "DMN"),
  conmat = sc_avg_thr_list$dmn$cLBP,
  view = c("left"),
  node.color = "#C44500",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 5000))

## sDMN - ncLBP
brain_sdmn_nclbp <- brainconn(
  atlas = bn_lut |> filter(network == "DMN"),
  conmat = sc_avg_thr_list$dmn$ncLBP,
  view = c("left"),
  node.color = "#66A61E",
  edge.width = edge_size,
  edge.color.weighted = T,
  scale.edge.width = F,
  node.size = node_size,
  show.legend = FALSE,
  background.alpha = 0.5) + 
  scale_edge_color_continuous(low = "white", high = "black", limits = c(0, 5000))

############################################
# Final plot
############################################

# Combine GNP plots
p_gnp_fvan <- p_fvan_l | p_fvan_eglob
p_gnp_fc_q_sub <- p_fc_q_sub
p_gnp_fvan_sub <- p_fvan_l_sub | p_fvan_eglob_sub
p_gnp_sdmn_sub <- p_sdmn_l_sub | p_sdmn_q_sub

# Combine with glass-brain
## Main analysis: fVAN
p_main_fvan <- p_gnp_fvan | (brain_fvan_cg / brain_fvan_lbpm / brain_fvan_lbpp)

## Subgroup: Full FC
p_sub_fc <- p_gnp_fc_q_sub | (brain_fc_cg / brain_fc_nclbp / brain_fc_clbp)

## Subgroup: fVAN
p_sub_fvan <- p_gnp_fvan_sub | (brain_fvan_cg / brain_fvan_nclbp / brain_fvan_clbp) 

## Subgroup: sDMN
p_sub_sdmn <- p_gnp_sdmn_sub | (brain_sdmn_cg / brain_sdmn_nclbp / brain_sdmn_clbp) 

## Correlation
p_cor <- readRDS("03_output/gt/p_cor.rds")

design <- "
AAABBB#
CCDDDEE
"

final_plot <- wrap_elements(p_main_fvan) + 
  wrap_elements(p_sub_fvan) + 
  wrap_elements(p_sub_fc) + 
  wrap_elements(p_sub_sdmn) +
  wrap_elements(p_cor) +
  plot_annotation(tag_levels = "a") +
  plot_layout(design = design) &
  theme(plot.tag = element_text(face = "bold", size = 18))

ggsave("03_output/figures/fig-x_gt_results.svg", final_plot, height = 10, width = 16)
