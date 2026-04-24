# Author: Lionel Butry (@lbutry)
# Purpose: Plot SFC on FSAVERAGE per group (for Figure in Supplements)

source("02_scripts/03_analysis/set_env.R")
library(fsbrain)
library(magick)

###############################################
# Plot SFC on FSAVERAGE per group (Supplements)
###############################################

# Get data
sfc <- readRDS("03_output/sfc/sfc.rds")
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds")

bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv") |> 
  mutate(network = case_when(
    network == "Limbic" ~ "LIM",
    network == "Subcortical" ~ "SUBC",
    .default = network
  ))

# Helper function
get_avg_sfc_per_group <- function(group_name) {
  sfc |>
    pivot_longer(-participant_id, names_to = "roi", values_to = "sfc") |>
    left_join(clin_df, by = "participant_id") |>
    filter(group == group_name) |>
    summarise(statistic = mean(sfc, na.rm=T), .by = roi)
}

get_avg_sfc_per_subgroup <- function(group_name) {
  test <- sfc |> 
    pivot_longer(-participant_id, names_to = "roi", values_to = "sfc") |>
    left_join(clin_df) |> 
    filter(group %in% c("LBP+", "CG")) |>
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
    mutate(group = factor(group, levels = c("CG", "ncLBP", "cLBP"))) |> 
    filter(group == group_name) |>
    summarise(statistic = mean(sfc, na.rm=T), .by = roi)
}

# SFC per group
sfc_cg <- get_avg_sfc_per_group("CG")
sfc_lbpm <- get_avg_sfc_per_group("LBP-")
sfc_lbpp <- get_avg_sfc_per_group("LBP+") 
sfc_nclbp <- get_avg_sfc_per_subgroup("ncLBP")
sfc_clbp <- get_avg_sfc_per_subgroup("cLBP")

sfc_cg$statistic |> max()

# FSBRAIN - Plotting function
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
    views = "sd_dorsal",
    makecmap_options = list(colFn = color_palette, range = color_range))
  
  return(cm)
}

save_fsbrain_manually_row <- function(fsbrain_instance, views, output_file, annot, tmp_dir = "fsbrain_tmp"){
  
  dir.create(tmp_dir)
  files <- character(length(views))
  
  # Snapshot each view
  for (i in seq_along(views)) {
    fsbrain::brainviews(views[i], fsbrain_instance)
    fname <- file.path(tmp_dir, paste0("view_", i, ".png"))
    rgl::snapshot3d(fname, webshot = TRUE, width = 2000, height = 2000)
    files[i] <- fname
  }
  
  # Append all views
  imgs <- magick::image_read(files)
  combined <- magick::image_append(imgs, stack = FALSE)
  
  # Add annotation
  combined <- image_annotate(
    combined, text = annot, 
    gravity = "southwest", color = "black", 
    size = 150, location = "+50+50")
  
  # Save img & delete tmp folder
  magick::image_write(combined, output_file)
  unlink(tmp_dir, recursive = TRUE)
}

# Plotting
views <- c("sd_lateral_lh", "sd_lateral_rh", "sd_dorsal", "sd_medial_lh", "sd_medial_rh")

p <- plot_fsbrain(sfc_cg, color_palette = viridis::magma, color_range = c(0, 0.1))
save_fsbrain_manually_row(
  fsbrain_instance = p, views = views, annot = "CG",
  output_file = "03_output/sfc/plots/sfc_cg.png"
)

p <- plot_fsbrain(sfc_lbpm, color_palette = viridis::magma, color_range = c(0, 0.1))
save_fsbrain_manually_row(
  fsbrain_instance = p, views = views, annot = "LBP-",
  output_file = "03_output/sfc/plots/sfc_lbp-.png"
)

p <- plot_fsbrain(sfc_lbpp, color_palette = viridis::magma, color_range = c(0, 0.1))
save_fsbrain_manually_row(
  fsbrain_instance = p, views = views, annot = "LBP+",
  output_file = "03_output/sfc/plots/sfc_lbp+.png"
)

p <- plot_fsbrain(sfc_nclbp, color_palette = viridis::magma, color_range = c(0, 0.1))
save_fsbrain_manually_row(
  fsbrain_instance = p, views = views, annot = "Subgroup: ncLBP",
  output_file = "03_output/sfc/plots/sfc_nclbp.png"
)

p <- plot_fsbrain(sfc_clbp, color_palette = viridis::magma, color_range = c(0, 0.1))
save_fsbrain_manually_row(
  fsbrain_instance = p, views = views, annot = "Subgroup: cLBP",
  output_file = "03_output/sfc/plots/sfc_clbp.png"
)

# Combine plots & add legend
## Helper function
append_legend_annot <- function(path_plot, path_legend, path_output) {
  a <- image_read(path_plot)
  a_height <- image_info(a)$height
  a_width  <- image_info(a)$width

  b <- image_read(path_legend) |>
    image_trim() |>
    image_scale(geometry = paste0("x", round(a_height * 0.6))) |>
    image_extent(geometry = paste0("x", a_height), gravity = "center", color = "white")

  c <- image_append(c(a, b))

  image_write(c, path_output)
}

## Save legend
legend <- ggplot(data.frame(x = 0:100, y = 1), aes(x = x, y = y, fill = x)) +
  geom_tile() +
  scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, 0.1), name = expression(SFC ~ (R[adj]^2)))

legend_grob <- cowplot::get_legend(legend)
legend_plot <- cowplot::ggdraw(legend_grob); legend_plot
ggsave("03_output/sfc/plots/color_legend_sfc.png", legend_plot)

## Iterate over all groups
input_files <- c(
  "03_output/sfc/plots/sfc_cg.png", 
  "03_output/sfc/plots/sfc_lbp-.png",
  "03_output/sfc/plots/sfc_lbp+.png",
  "03_output/sfc/plots/sfc_nclbp.png",
  "03_output/sfc/plots/sfc_clbp.png")
output_files <- c(
  "03_output/sfc/plots/sfc_legend_cg.png", 
  "03_output/sfc/plots/sfc_legend_lbp-.png", 
  "03_output/sfc/plots/sfc_legend_lbp+.png",
  "03_output/sfc/plots/sfc_legend_nclbp.png",
  "03_output/sfc/plots/sfc_legend_clbp.png")
legend_file <- "03_output/sfc/plots/color_legend_sfc.png"

map2(input_files, output_files, ~ append_legend_annot(.x, legend_file, .y))

# Combine all group-SFC plots
sfc_plots <- output_files

all_plots <- image_read(sfc_plots)
final_plot <- image_append(all_plots, stack = TRUE)
image_write(final_plot, "03_output/sfc/plots/_sfc_final.png")

# rgl::rglwidget() # view curret fsbrain instance
