# Author: Lionel Butry (@lbutry)
# Purpose: Explorative partial correlation analysis between GNPs and clinical measures
#     + plotting the results in a heatmap

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
library(patchwork)

############################################
# ---- Import data ----
############################################

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
  left_join(sc_dmn_l_q)

rm(sc_dmn_l_q, fc_full_q, fc_van_l_eglob)

# Clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") |> 
  select(
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
    ) 

gnp_sig_clin <- all_sig_metrics |> 
  inner_join(clin_df) |> 
  mutate(
    Demo_Gender = as.numeric(Demo_Gender == "male"),
    scanner = as.numeric(scanner == "SIEM")
  )

############################################
# ---- Explorative correlation analysis ----
############################################

#### On full dataset (PROMIS D, PROMIS A, CSI, TSK, ODI) ####

# Set variables
iv_vars <- c("Quest_TSK_Score", "Quest_PROMIS_Anxiety", "Quest_PROMIS_Depression", "Quest_ODI_Score", "Quest_CSI_Score")
dv_vars <- gnp_sig_clin |> select(starts_with("fc_"), starts_with("sc_")) |> names()
covars <- c("Demo_Age", "Demo_Gender", "scanner")

# Partial correlation analysis (controlling for age, sex, scanner)
a <- map_dfr(dv_vars, ~{
  dv <- .x
  map_dfr(iv_vars, ~{
    iv <- .x
    df <- gnp_sig_clin |>  select(all_of(c(dv, iv, covars))) %>% drop_na()
    res <- pcor.test(df[[dv]], df[[iv]], df |>  select(all_of(covars)), method = "spearman")
    tibble(sample = "all", dv = dv, iv = iv, estimate = res$estimate, statistic = res$statistic, p.value = res$p.value)
  })
})

#### On LBP+ dataset (Avg pain, Cur pain, pain duration) ####
# Set variables
iv_vars <- c("vas_mri", "Pain_Intensity_average_Q1", "Pain_Duration_Q1")
dv_vars <- gnp_sig_clin |> select(starts_with("fc_"), starts_with("sc_")) |> names()
covars <- c("Demo_Age", "Demo_Gender", "scanner")

gnp_sig_clin_lbp <- gnp_sig_clin |> 
  filter(group != "CG")

b <- map_dfr(dv_vars, ~{
  dv <- .x
  map_dfr(iv_vars, ~{
    iv <- .x
    df <- gnp_sig_clin_lbp |>  select(all_of(c(dv, iv, covars))) %>% drop_na()
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

# Tidy results
results_partial_neat <- results_partial |> 
  mutate(
    dv = case_when(
      dv == "fc_full_q" ~ "Q of FC",
      dv == "fc_van_eglob" ~ "Eglob of fVAN",
      dv == "fc_van_l" ~ "L of fVAN",
      dv == "sc_dmn_l" ~ "L of sDMN",
      dv == "sc_dmn_q" ~ "Q of sDMN",
      .default = NA),
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
    iv = factor(iv, level = rev(c("Avg. pain intensity", "Cur. pain intensity", "Pain duration", "CSI", "ODI", "PROMIS A", "PROMIS D", "TSK"))),
    dv = factor(dv, level = c("Q of FC", "Eglob of fVAN", "L of fVAN", "Q of sDMN", "L of sDMN"))
  )

write.csv2(results_partial_neat, "03_output/gt/results_partical_cor.csv", row.names = F)

############################################
# ---- Plot results ----
############################################

my_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      text = element_text(size = 12, family = "Helvetica", color = "black"),
      panel.grid = element_blank(),
      plot.tag = element_text(face = "bold")
    )
}

plot_cor_overview <- function(data) {
  
  ggplot(data, aes(x = dv, y = iv, fill = estimate)) +
    geom_tile() +
    geom_text(aes(label = sig), color = "black", size = 6) +
    scale_fill_continuous(
      name = expression(ρ[partial]),
      palette = RColorBrewer::brewer.pal(11, "PRGn"),
      limits = c(-0.2, 0.2)) +
    my_theme() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank()
    )
}

p_all <- plot_cor_overview(results_partial_neat |> filter(sample == "all"))
p_lbp <- plot_cor_overview(results_partial_neat |> filter(sample == "lbp")) +
  theme(axis.text.x = element_blank())

p_cor <- p_lbp / p_all + 
  plot_layout(heights = c(3, 5), guides = "collect")

saveRDS(p_cor, "03_output/gt/p_cor.rds")
