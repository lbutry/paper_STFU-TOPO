# Author: Lionel Butry (@lbutry)
# Purpose: Create plots for the SFC analysis
#   (A) Surface projected results of regional SFC using "fsbrain"
#   (B) Boxplot-like for group differences

source("02_scripts/03_analysis/set_env.R")
set.seed(42)

###############################################
# Plot global & network results
###############################################

# Import data
sfc_agg <- readRDS("03_output/sfc/sfc_agg.rds")

# Helper functions
my_theme <- function() {
  theme_minimal(base_size = 15) + 
    theme(
      text = element_text(size = 13, family = "Helvetica", color = "black"), 
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

# Plotting - Main analysis
sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("ns", "*", "ns"), y.position = c(0.068, 0.062, 0.065))
p1 <- plot_adj(sfc_agg, "global", "Global SFC") +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.0725, label = expression(ANCOVA:~p[FDR]<0.001), size = 3, hjust = 0)
saveRDS(p1, "03_output/sfc/plots/p_global.rds")

sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("*", "*", "ns"), y.position = c(0.088, 0.081, 0.0845))
p2 <- plot_adj(sfc_agg, "DMN", "DMN-SFC") +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.093, label = expression(ANCOVA:~p[FDR]~"= 0.002"), size = 3, hjust = 0)
saveRDS(p2, "03_output/sfc/plots/p_dmn.rds")

sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("*", "**", "ns"), y.position = c(0.069, 0.062, 0.0655))
p3 <- plot_adj(sfc_agg, "VAN", "VAN-SFC") +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.073, label = expression(ANCOVA:~p[FDR]<0.001), size = 3, hjust = 0) 
saveRDS(p3, "03_output/sfc/plots/p_van.rds")

sig <- tibble(group1 = c("CG", "CG", "LBP-"), group2 = c("LBP+", "LBP-", "LBP+"), p.adj = c("ns", "ns", "ns"), y.position = c(0.091, 0.0865, 0.083))
p4 <- plot_adj(sfc_agg, "SMN", "SMN-SFC") +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.096, label = expression(ANCOVA:~p[FDR]~"= 0.027"), size = 3, hjust = 0)
saveRDS(p4, "03_output/sfc/plots/p_smn.rds")

# Plotting - Subgroup analysis

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

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("*", "ns", "ns"), y.position = c(0.068, 0.062, 0.065))
p5 <- plot_adj(sfc_agg_sub, "global", "Global SFC", color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.0725, label = expression(ANCOVA:~p[FDR]<0.001), size = 3, hjust = 0)
saveRDS(p5, "03_output/sfc/plots/p_global_sub.rds")

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("ns", "ns", "ns"), y.position = c(0.088, 0.081, 0.0845))
p6 <- plot_adj(sfc_agg_sub, "DMN", "DMN-SFC", color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.093, label = expression(ANCOVA:~p[FDR]~"= 0.005"), size = 3, hjust = 0)
saveRDS(p6, "03_output/sfc/plots/p_dmn_sub.rds")

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("**", "*", "ns"), y.position = c(0.069, 0.062, 0.0655))
p7 <- plot_adj(sfc_agg_sub, "VAN", "VAN-SFC", color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.073, label = expression(ANCOVA:~p[FDR]<0.001), size = 3, hjust = 0) 
saveRDS(p7, "03_output/sfc/plots/p_van_sub.rds")

sig <- tibble(group1 = c("CG", "CG", "ncLBP"), group2 = c("cLBP", "ncLBP", "cLBP"), p.adj = c("ns", "ns", "ns"), y.position = c(0.091, 0.0865, 0.083))
p8 <- plot_adj(sfc_agg_sub, "SMN", "SMN-SFC", color = c("#0072B2", "#66A61E", "#C44500")) +
  add_pvalue(sig, tip.length = 0.01, label.size = 4) #+
  #annotate("text", x = 0.5, y = 0.096, label = expression(ANCOVA:~p[FDR]~"= 0.061"), size = 3, hjust = 0)
saveRDS(p8, "03_output/sfc/plots/p_smn_sub.rds")

p_all <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8)

ggsave("03_output/sfc/plots/fig-sfc_means.png", p_all, height = 8, width = 10)

###############################################
# Plot region results on FSAVERAGE
###############################################
 
library(fsbrain)
library(magick)

#### IMPORT DATA ####
# Brainnetome LUT
bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv") |> 
  mutate(network = case_when(
    network == "Limbic" ~ "LIM",
    network == "Subcortical" ~ "SUBC",
    .default = network
  ))

# Results
results_contrast_cg_lbpp <- read.delim("03_output/sfc/results_roi_contr_cg>lbp+.csv", sep = ";")
results_contrast_cg_lbpm <- read.delim("03_output/sfc/results_roi_contr_cg>lbp-.csv", sep = ";")
results_contrast_lbpm_lbpp <- read.delim("03_output/sfc/results_roi_contr_lbp->lbp+.csv", sep = ";")

results_contrast_cg_clbp <- read.delim("03_output/sfc/results_roi_contr_cg>clbp.csv", sep = ";")
results_contrast_cg_nclbp <- read.delim("03_output/sfc/results_roi_contr_cg>nclbp.csv", sep = ";")
results_contrast_nclbp_clbp <- read.delim("03_output/sfc/results_roi_contr_ncLBP>clbp.csv", sep = ";")

#### FUNCTIONS ####
plot_fsbrain <- function(
    statistic_df, 
    color_palette = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu"))),
    color_range = c(-4, 4)
) {
   
   # Join df with Brainnetome LUT
   df <- statistic_df |> 
     rename(ROI.Name = roi) |> 
     left_join(bn_lut)
   
   # Extract statistic in fsbrain readable format
   input_named <- df |> pull(statistic)
   names(input_named) <- df$region
   input_named_lh <- input_named[grepl("_L$", names(input_named))]
   input_named_rh <- input_named[grepl("_R$", names(input_named))]
   
   sjd <- "01_data/ready4analysis/BN_Atlas_freesurfer"
   
   cm <- vis.region.values.on.subject(
     subjects_dir = sjd,
     subject_id = "fsaverage",
     atlas = "brainnetome",
     lh_region_value_list = input_named_lh,
     rh_region_value_list = input_named_rh,
     surface = "pial",
     draw_colorbar = TRUE,
     views = "t4",
     makecmap_options = list(colFn = color_palette, range = color_range))
   
   return(cm)
}
 
save_fsbrain_manually <- function(fsbrain_instance, views, output_file, annot, tmp_dir = "fsbrain_tmp"){
   
  dir.create(tmp_dir)
  files <- character(length(views))
   
  # Snapshot each view
  for (i in seq_along(views)) {
    fsbrain::brainviews(views[i], fsbrain_instance)
    fname <- file.path(tmp_dir, paste0("view_", i, ".png"))
    rgl::snapshot3d(fname, webshot = TRUE, width = 2000, height = 2000)
    files[i] <- fname
  }
  
  # Append all views (2x2 grid)
  imgs <- magick::image_read(files)
  
  ## Trim & keep a bit of padding
  imgs <- magick::image_trim(imgs)
  imgs <- magick::image_border(imgs, color = "white", geometry = "70x70")

  ## Append each pair horizontally
  img_groups <- split(imgs, ceiling(seq_along(imgs) / 2))
  row_images <- lapply(img_groups, function(x) {
    magick::image_append(x, stack = FALSE)
  })
  
  ## Stack rows vertically & add overall padding
  combined <- magick::image_append(magick::image_join(row_images), stack = TRUE)
  combined <- magick::image_border(combined, color = "white", geometry = "120x120") 
  
  # Add annotation
  combined <- image_annotate(
    combined, text = annot, 
    gravity = "north", color = "black", 
    size = 180, location = "+0+10")
  
  # Save img & delete tmp folder
  magick::image_write(combined, output_file)
  unlink(tmp_dir, recursive = TRUE)
}

#### PLOTTING ####
# "sd_dorsal"
views <- c("sd_lateral_lh", "sd_lateral_rh", "sd_medial_lh", "sd_medial_rh")

# Plotting & saving tmp files
### CG > LBP-
p <- plot_fsbrain(results_contrast_cg_lbpm)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_cg>lbp-.png",
  annot = "CG > LBP-"
)
 
### CG > LBP+
p <- plot_fsbrain(results_contrast_cg_lbpp)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_cg>lbp+.png",
  annot = "CG > LBP+"
)
 
### LBP- > LBP+
p <- plot_fsbrain(results_contrast_lbpm_lbpp)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_lbp->lbp+.png",
  annot = "LBP- > LBP+"
)
 
## Subgroup analysis
### CG > cLBP
p <- plot_fsbrain(results_contrast_cg_clbp)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_cg>clbp.png",
  annot = "CG > cLBP"
)

p <- plot_fsbrain(results_contrast_cg_nclbp)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_cg>nclbp.png",
  annot = "CG > ncLBP"
)

p <- plot_fsbrain(results_contrast_nclbp_clbp)
save_fsbrain_manually(
  fsbrain_instance = p, views = views,
  output_file = "03_output/sfc/plots/p_fsbrain_nclbp>clbp.png",
  annot = "ncLBP > cLBP"
)

# Helper functions
# rgl::rglwidget() # view curret fsbrain instance
# view_angles <- get.view.angle.names("t9") # names of view angles

###############################################
# Append all fsbrain plots
###############################################

# Manually annotate ROIs
## Show significant rois
results_contrast_nclbp_clbp |>
  filter(p.value < 0.05) |>
  pull(roi) |> sort()

 testing <- results_contrast_nclbp_clbp |>
  mutate(statistic = case_when(
    str_detect(roi, "PCun_R_4_1") ~ 3,
    str_detect(roi, "PhG_L_6_2") ~ 3,
    str_detect(roi, "PhG_R_6_6") ~ 3,
    .default = 0
  ))

p <- plot_fsbrain(testing)
rgl::rglwidget()

# Get legend
legend <- ggplot(data.frame(x = 0:100, y = 1), aes(x = x, y = y, fill = x)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100), limits = c(-4, 4), name = "t statistic") +
  theme(legend.position = "bottom", 
        legend.title.position = "top", 
        legend.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 13))

legend_grob <- cowplot::get_legend(legend)
legend_plot <- cowplot::ggdraw(legend_grob)
ggsave("03_output/sfc/plots/color_legend_contrast.png", legend_plot)

append_legend <- function(img, path_legend) {
  img_width <- image_info(img)$width
  legend <- image_read(path_legend) |> 
    image_trim() |> 
    image_scale(geometry = paste0("x", round(img_width * 0.15))) |> 
    image_extent(geometry = paste0(img_width, "x"), gravity = "center", color = "white")
  image_append(c(img, legend), stack = TRUE)
}

# Main analysis
sfc_plots1 <- c(
  "03_output/sfc/plots/p_fsbrain_cg>lbp-_annot.png",
  "03_output/sfc/plots/p_fsbrain_cg>lbp+_annot.png",
  "03_output/sfc/plots/p_fsbrain_lbp->lbp+_annot.png"
)

all_plots1 <- image_read(sfc_plots1)
final_plot1 <- image_append(all_plots1, stack = TRUE)
final_plot1 <- append_legend(final_plot1, "03_output/sfc/plots/color_legend_contrast.png")
image_write(final_plot1, "03_output/sfc/plots/sfc_fsbrain_contrasts_main.png")

# Subgroup analysis
sfc_plots2 <- c(
  "03_output/sfc/plots/p_fsbrain_cg>nclbp_annot.png",
  "03_output/sfc/plots/p_fsbrain_cg>clbp_annot.png",
  "03_output/sfc/plots/p_fsbrain_nclbp>clbp_annot.png"
)

all_plots2 <- image_read(sfc_plots2)
final_plot2 <- image_append(all_plots2, stack = TRUE)
final_plot2 <- append_legend(final_plot2, "03_output/sfc/plots/color_legend_contrast.png")
image_write(final_plot2, "03_output/sfc/plots/sfc_fsbrain_contrasts_sub.png")

###############################################
# Append everything to one big beautiful plot
###############################################

path_p_main <- "03_output/sfc/plots/sfc_group_comparison_main.png"
path_p_sub <- "03_output/sfc/plots/sfc_group_comparison_sub.png"

p_main <- p1 + p2 + p3 + p4
ggsave(path_p_main, p_main, height = 9, width = 5, dpi = 1000)

p_sub <- p5 + p6 + p7 + p8
ggsave(path_p_sub, p_sub, height = 9, width = 5, dpi = 1000)

imgs <- magick::image_read(c(
  path_p_main, "03_output/sfc/plots/sfc_fsbrain_contrasts_main.png",
  path_p_sub, "03_output/sfc/plots/sfc_fsbrain_contrasts_sub.png")
)

a <- imgs |> 
  image_scale("x6000") |> 
  image_append() |> 
  image_border(color = "white", geometry = "75x75")

img_width <- image_info(a)$width

img_cor <- image_read("03_output/sfc/plots/p_cor.png") |> 
  image_scale(geometry = paste(img_width, "x"))

final <- image_append(c(a, img_cor), stack = TRUE)

image_write(final, "03_output/sfc/plots/_sfc_group_comparison_.png")
