# Author: Lionel Butry (@lbutry)
# Purpose: Explorative partial correlation analysis between SFC and clinical measures
#   + plotting the results in a heatmap

source("02_scripts/03_analysis/set_env.R")

###############################################
# Read data
###############################################

sfc <- readRDS("03_output/sfc/sfc.rds")
sfc_agg <- readRDS("03_output/sfc/sfc_agg.rds")

results_roi <- list.files(path = "03_output/sfc", pattern = "results_roi_contr", full.names = TRUE)

sig_roi <- map(results_roi, ~ read.csv2(.x)) |> 
  bind_rows() |> 
  filter(p.value < 0.05) |> 
  pull(roi) |> 
  unique()

sig_sfc <- sfc_agg |> 
  left_join(sfc) |> 
  select(
    global, DMN, VAN, SMN,
    any_of(sig_roi),
    participant_id,
    group,
    Demo_Age,
    Demo_Gender,
    scanner,
    Pain_Intensity_average_Q1,
    vas_mri,
    Pain_Duration_Q1,
    Quest_ODI_Score,
    Quest_TSK_Score,
    Quest_CSI_Score,
    Quest_PROMIS_Depression,
    Quest_PROMIS_Anxiety
  ) |> 
  mutate(
    Demo_Gender = as.numeric(Demo_Gender == "male"),
    scanner = as.numeric(scanner == "SIEM")
  )

###############################################
# Explorative partial correlations analysis
###############################################

#### On full dataset (PROMIS D, PROMIS A, CSI, TSK, ODI) ####

# Set variables
iv_vars <- c("Quest_TSK_Score", "Quest_PROMIS_Anxiety", "Quest_PROMIS_Depression", "Quest_ODI_Score", "Quest_CSI_Score")
dv_vars <- c("global", "DMN", "VAN", "SMN", sig_roi)
covars <- c("Demo_Age", "Demo_Gender", "scanner")

# Partial correlation analysis (controlling for age, sex, scanner)
a <- map_dfr(dv_vars, ~{
  dv <- .x
  map_dfr(iv_vars, ~{
    iv <- .x
    df <- sig_sfc |>  select(all_of(c(dv, iv, covars))) %>% drop_na()
    res <- pcor.test(df[[dv]], df[[iv]], df |>  select(all_of(covars)), method = "spearman")
    tibble(sample = "all", dv = dv, iv = iv, estimate = res$estimate, statistic = res$statistic, p.value = res$p.value)
  })
})

#### On LBP+ dataset (Avg pain, Cur pain, pain duration) ####
# Set variables
iv_vars <- c("vas_mri", "Pain_Intensity_average_Q1", "Pain_Duration_Q1")
dv_vars <- c("global", "DMN", "VAN", "SMN", sig_roi)
covars <- c("Demo_Age", "Demo_Gender", "scanner")

sig_sfc_lbp <- sig_sfc |> 
  filter(group != "CG")

b <- map_dfr(dv_vars, ~{
  dv <- .x
  map_dfr(iv_vars, ~{
    iv <- .x
    df <- sig_sfc_lbp |>  select(all_of(c(dv, iv, covars))) %>% drop_na()
    res <- pcor.test(df[[dv]], df[[iv]], df |>  select(all_of(covars)), method = "spearman")
    tibble(sample = "lbp", dv = dv, iv = iv, estimate = res$estimate, statistic = res$statistic, p.value = res$p.value)
  })
})

#### Combine & tidy results ####
results_partial <- bind_rows(a, b)
rm(a, b)

# Add significant labels
results_partial <- results_partial |> 
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

max(abs(results_partial$estimate))

# Tidy results
results_partial_neat <- results_partial |> 
  mutate(
    dv = case_when(
      dv == "global" ~ "Global",
      .default = dv),
    iv = case_when(
      iv == "vas_mri" ~ "Cur. pain intensity",
      iv == "Quest_TSK_Score" ~ "TSK",
      iv == "Quest_PROMIS_Depression" ~ "PROMIS D",
      iv == "Quest_PROMIS_Anxiety" ~ "PROMIS A",
      iv == "Quest_ODI_Score" ~ "ODI",
      iv == "Quest_CSI_Score" ~ "CSI",
      iv == "Pain_Intensity_average_Q1" ~ "Avg. pain intensity",
      iv == "Pain_Duration_Q1" ~ "Pain duration",
      .default = NA)
  ) |> 
  mutate(
    iv = factor(iv, level = rev(c("Avg. pain intensity", "Cur. pain intensity", "Pain duration", "CSI", "ODI", "PROMIS A", "PROMIS D", "TSK")))
  ) |> 
  mutate(p.fdr = p.adjust(p.value))

saveRDS(results_partial_neat, "03_output/sfc/results_sfc_cor.rds")

# Explore results manually
explore <- results_partial_neat |> 
  filter(p.value < 0.05) |> 
  filter(iv == "Pain duration") |> 
  filter(estimate < 0)


# Reorder DV by network
bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv") |> 
  rename(dv = ROI.Name)

all_dv <- unique(results_partial_neat$dv)
globals <- c("Global", "DMN", "VAN", "SMN")
others <- setdiff(all_dv, globals)

tmp_order <- results_partial_neat |> 
  left_join(bn_lut) |> 
  filter(dv %in% others) |> 
  arrange(network, dv)

ordered_dv <- c(globals, unique(tmp_order$dv))

results_partial_neat <- results_partial_neat |> 
  left_join(bn_lut) |> 
  mutate(dv = factor(dv, level = ordered_dv)) |> 
  mutate(network = ifelse(dv %in% globals, " ", network))

############################################
# ---- Plot results ----
############################################

my_theme <- function() {
  theme_minimal(base_size = 15) +
    theme(
      text = element_text(size = 15, family = "Helvetica", color = "black"),
      panel.grid = element_blank(),
      plot.tag = element_text(face = "bold"),
      axis.text.y = element_text(size = 13)
    )
}

plot_cor_overview <- function(data) {
  
  ggplot(data, aes(x = dv, y = iv, fill = estimate)) +
    geom_tile() +
    geom_text(aes(label = sig), color = "black", size = 6) +
    scale_fill_continuous(
      name = expression(ρ[partial]),
      palette = RColorBrewer::brewer.pal(11, "PRGn"),
      limits = c(-0.26, 0.26)) +
    my_theme() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.title = element_blank()
    )
}

p_lbp <- plot_cor_overview(results_partial_neat |> filter(sample == "lbp")) +
  theme(axis.text.x = element_blank())

color_network <- c("#FFFFFF00", ggsci::pal_igv("default")(16)[c(2:8, 16)])

margin_labels <- results_partial_neat |>
  filter(sample == "all") |>
  distinct(dv, network) |>
  mutate(dv = factor(dv, levels = levels(results_partial_neat$dv))) |>
  arrange(dv) |>
  group_by(network) |>
  summarise(
    xmin = min(as.numeric(dv)),
    xmax = max(as.numeric(dv)),
    x = mean(as.numeric(dv)),
    .groups = "drop"
  )

p_all <- plot_cor_overview(results_partial_neat |> filter(sample == "all")) +
  labs(x = "SFC") + 
  # Add tilemargins (colored per network)
  ggh4x::geom_tilemargin(aes(x = dv, species = network), sides = "b", length = grid::unit(0.14, "npc")) +
  scale_listed(scalelist = list(
    scale_fill_continuous(name = expression(ρ[partial]),palette = RColorBrewer::brewer.pal(11, "PRGn"),limits = c(-0.25, 0.25)),
    scale_fill_manual(values = color_network, aesthetics = "species", guide = "none")
  ), replaces = c("fill", "fill")) +
  # Network labels
  geom_text(
    data = margin_labels,
    aes(x = x, y = 0, label = network),
    vjust = -0.3, size = 4,
    inherit.aes = FALSE
  ) +
  theme(axis.title.x = element_text())

p_cor <- p_lbp / p_all +
  plot_layout(heights = c(3, 5), guides = "collect")

ggsave("03_output/sfc/plots/p_cor.png", p_cor, height = 4, width = 20, dpi = 1000)
