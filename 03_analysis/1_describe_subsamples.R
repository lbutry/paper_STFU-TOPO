# Author: Lionel Butry (@lbutry)

# Describe characteristics of LBP+ subsamples: 
# - pain duration >3 months (cLBP)
# - pain duration <3 months (ncLBP)

############################################
# ---- Set up env ----
############################################

source("02_scripts/03_analysis/set_env.R")

############################################
# ---- Import data ----
############################################

clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") |>
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
  mutate(group = factor(group, levels = c("CG", "ncLBP", "cLBP")))

############################################
# ---- Descriptive statistics ----
############################################

# Overall impression
skim(clin_df)

# Summary: Numeric variables 
clin_summary_num <- clin_df |> 
  select(group,
         Demo_Age,
         Demo_BMI,
         Pain_Duration_Q1,
         Pain_Intensity_average_Q1,
         vas_mri,
         Pain_Past_duration_Q1,
         Quest_ODI_Score,
         Quest_TSK_Score,
         Quest_PROMIS_Anxiety,
         Quest_PROMIS_Depression,
         Quest_CSI_Score
         ) |> 
  group_by(group) |>
  get_summary_stats() |>
  adorn_rounding(digits = 2) |>
  rowwise() |>
  mutate(value = paste0(median, " (", q1, ", ", q3, ")")) |> 
  select(group, variable, n, mean, sd, value)

# Summary: Factorial variables
sex_cg <- clin_df |> filter(group == "CG") |> tabyl(Demo_Gender) |> adorn_rounding(digits = 2) |> mutate(value = paste0(n, " (", percent * 100, " %)"))
sex_clbp <- clin_df |> filter(group == "cLBP") |> tabyl(Demo_Gender) |> adorn_rounding(digits = 2) |> mutate(value = paste0(n, " (", percent * 100, " %)"))
sex_nclbp <- clin_df |> filter(group == "ncLBP") |> tabyl(Demo_Gender) |> adorn_rounding(digits = 2) |> mutate(value = paste0(n, " (", percent * 100, " %)"))

## Pain medication use
meds <- clin_df |> 
  select(participant_id, group, meds) |> 
  mutate(
    meds_nsaid = if_else(str_detect(meds, "NSAID"), TRUE, FALSE),
    meds_meta = if_else(str_detect(meds, "Metamizol"), TRUE, FALSE),
    meds_para = if_else(str_detect(meds, "Paracetamol"), TRUE, FALSE),
    meds_wopi = if_else(str_detect(meds, "Weak opi"), TRUE, FALSE),
    meds_sopi = if_else(str_detect(meds, "Strong opi"), TRUE, FALSE),
    meds_canna = if_else(str_detect(meds, "Cannab"), TRUE, FALSE),
    meds_relax = if_else(str_detect(meds, "Muscle"), TRUE, FALSE),
    meds_other = if_else(str_detect(meds, "Other"), TRUE, FALSE))

meds_per_type <- meds |> 
  summarise(
    nsaid = sum(meds_nsaid),
    meta = sum(meds_meta),
    para = sum(meds_para),
    wopi = sum(meds_wopi),
    sopi = sum(meds_sopi),
    canna = sum(meds_canna),
    relax = sum(meds_relax),
    other = sum(meds_other),
    .by = "group"
  )

############################################
# ---- Inference statistics ----
############################################

# Age
ggpubr::ggqqplot(clin_df, "Demo_Age") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Demo_Age ~ group) # homogeneity of variance
kruskal_test(clin_df, Demo_Age ~ group)
dunn_test(clin_df, Demo_Age ~ group, p.adjust.method = "bonferroni")

# BMI
ggpubr::ggqqplot(clin_df, "Demo_BMI") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Demo_BMI ~ group) # homogeneity of variance
kruskal_test(clin_df, Demo_BMI ~ group)
dunn_test(clin_df, Demo_BMI ~ group, p.adjust.method = "bonferroni")

# Sex
chisq_test(clin_df$Demo_Gender, clin_df$group)

# Avg pain
ggpubr::ggqqplot(clin_df, "Pain_Intensity_average_Q1") + facet_grid(. ~ group) # => no normal distribution
clin_df |>
  filter(group != "CG") |>
  droplevels() |>
  wilcox_test(Pain_Intensity_average_Q1 ~ group)

# Cur pain
ggpubr::ggqqplot(clin_df, "vas_mri") + facet_grid(. ~ group) # => no normal distribution
clin_df |>
  filter(group != "CG") |>
  droplevels() |>
  wilcox_test(vas_mri ~ group)

# Pain duration
ggpubr::ggqqplot(clin_df, "Pain_Duration") + facet_grid(. ~ group) # => no normal distribution
clin_df |>
  filter(group != "CG") |>
  droplevels() |>
  wilcox_test(Pain_Duration ~ group)

# ODI
ggpubr::ggqqplot(clin_df, "Quest_ODI_Score") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Quest_ODI_Score ~ group) # => no homogeneity of variance
kruskal_test(clin_df, Quest_ODI_Score ~ group)
dunn_test(clin_df, Quest_ODI_Score ~ group, p.adjust.method = "bonferroni")

# TSK
ggpubr::ggqqplot(clin_df, "Quest_TSK_Score") + facet_grid(. ~ group) # => approx. normal distribution
levene_test(clin_df, Quest_TSK_Score ~ group) # => homogeneity of variance
anova_test(clin_df, Quest_TSK_Score ~ group)
tukey_hsd(clin_df, Quest_TSK_Score ~ group, p.adjust.method = "bonferroni")

# CSI
ggpubr::ggqqplot(clin_df, "Quest_CSI_Score") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Quest_CSI_Score ~ group) # => no homogeneity of variance
kruskal_test(clin_df, Quest_CSI_Score ~ group)
dunn_test(clin_df, Quest_CSI_Score ~ group, p.adjust.method = "bonferroni")

# PROMIS Anxiety
ggpubr::ggqqplot(clin_df, "Quest_PROMIS_Anxiety") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Quest_PROMIS_Anxiety ~ group) # => homogeneity of variance
kruskal_test(clin_df, Quest_PROMIS_Anxiety ~ group)
dunn_test(clin_df, Quest_PROMIS_Anxiety ~ group, p.adjust.method = "bonferroni")

# PROMIS Depression
ggpubr::ggqqplot(clin_df, "Quest_PROMIS_Depression") + facet_grid(. ~ group) # => no normal distribution
levene_test(clin_df, Quest_PROMIS_Depression ~ group) # => no homogeneity of variance
kruskal_test(clin_df, Quest_PROMIS_Depression ~ group)
dunn_test(clin_df, Quest_PROMIS_Depression ~ group, p.adjust.method = "bonferroni")
