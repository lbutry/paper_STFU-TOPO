# Requirements: Region-wise results from "4_sfc_analysis.R"
# Author: Lionel Butry (@lbutry)
# Purpose: Alternative result visualisation of regional SFC for the Supplements

source("02_scripts/03_analysis/set_env.R")
library(patchwork)

# Load BN LUT
bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv") |> 
  rename(roi = ROI.Name)

# Load region-wise results
results_paths <- list.files(path = "03_output/sfc", pattern = "results_roi_contr", full.names = T)
results_all <- map_dfr(results_paths, ~ read_csv2(.x)) |> 
  mutate(statistic = as.numeric(statistic)) |> 
  filter(p.value < 0.05)

# Load correlation results & keep only significant roi
sig_cor <- readRDS("03_output/sfc/results_sfc_cor.rds") |> 
  group_by(dv) |> 
  summarise(sig = any(p.value < 0.05)) |> 
  filter(sig == TRUE) |> 
  rename(roi = dv)

# Plot significant results per analysis
main_contrasts <- c("CG > LBP-", "CG > LBP+", "LBP- > LBP+")
sub_contrasts <- c("CG > ncLBP", "CG > cLBP", "ncLBP > cLBP")

my_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      text = element_text(size = 12, family = "Helvetica", color = "black"),
      panel.grid = element_blank(),
      plot.tag = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.x = element_line()
    )
}

plot_compare_sig_results <- function(contr_vec) {
  
  results <- results_all |> 
    dplyr::filter(
      contrast %in% contr_vec,
      p.value < 0.05
    ) |> 
    dplyr::left_join(bn_lut, by = "roi") |> 
    dplyr::arrange(network, roi) |> 
    dplyr::mutate(
      roi = factor(roi, levels = unique(roi)),
      contrast = factor(contrast, levels = contr_vec)
      )
  
  # network colors
  color_network <- ggsci::pal_igv("default")(16)[c(2:8, 16)]
  names(color_network) <- c("CEN", "DAN", "DMN", "LIM", "SMN", "SUBC", "VAN", "VIS")
  
  # bold face vector for y-axis labels
  roi_levels <- levels(results$roi)
  label_face <- ifelse(roi_levels %in% sig_cor$roi, "bold", "plain")
  
  # margin label positions
  margin_labels <- results |>
    dplyr::distinct(roi, network) |>
    dplyr::mutate(roi = factor(roi, levels = levels(results$roi))) |>
    dplyr::arrange(roi) |>
    dplyr::group_by(network) |>
    dplyr::summarise(
      ymin = min(as.numeric(roi)),
      ymax = max(as.numeric(roi)),
      y = mean(as.numeric(roi)),
      .groups = "drop"
    )
  
  ggplot(results, aes(x = contrast, y = roi, fill = statistic)) +
    geom_tile() +
    
    ggh4x::geom_tilemargin(
      aes(x = contrast, species = network),
      sides = "l",
      length = grid::unit(0.15, "npc")
    ) +
    
    ggh4x::scale_listed(
      scalelist = list(
        scale_fill_gradientn(
          name = "t statistic",
          colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
          limits = c(-4, 4)
        ),
        scale_fill_manual(values = color_network, aesthetics = "species", guide = "none")
      ),
      replaces = c("fill", "fill")
    ) +
    
    geom_text(
      data = margin_labels,
      aes(x = 0.4, y = y, label = network),
      angle = 90,
      hjust = 0.5,
      vjust = 1.4,
      size = 3,
      inherit.aes = FALSE
    ) +
    
    labs(y = "SFC") +
    my_theme() +
    theme(axis.text.y = element_text(face = label_face))
}

a <- plot_compare_sig_results(main_contrasts)
b <- plot_compare_sig_results(sub_contrasts)

final <- a + b + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

ggsave("03_output/sfc/plots/fig_sx_sfc_region.svg", final)
