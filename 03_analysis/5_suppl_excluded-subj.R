# Author: Lionel Butry (@lbutry)
# Purpose:
#   (A) Identify participants not included in analysis
#   (B) Test, wheter excluded participants differ from included participants

source("02_scripts/03_analysis/set_env.R")

# Get IDs of participants who have SC, FC & SFC
path_fc <- "01_data/ready4analysis/f04_connectome_partial"
path_sc <- "01_data/ready4analysis/d02_connectome"

id_fc <- list.files(path=path_fc, pattern = "\\.csv$") |> str_extract("PREDICT[0-9]+")
id_sc <- list.files(path=path_sc, pattern = "\\connectome_sift2.csv$") |> str_extract("PREDICT[0-9]+")
id_sfc <- intersect(id_fc, id_sc)

# Read clinical data
clin_df <- readRDS("01_data/ready4analysis/clin_data_353.rds") 

# Subset clin_df for included & excluded participants
clin_df2 <- clin_df |> 
  mutate(
    fc_analysis = ifelse(participant_id %in% id_fc, "included", "excluded"),
    sc_analysis = ifelse(participant_id %in% id_sc, "included", "excluded"),
    sfc_analysis = ifelse(participant_id %in% id_sfc, "included", "excluded")
  ) |> 
  mutate(
    fc_analysis  = factor(fc_analysis, levels = c("included","excluded")),
    sc_analysis  = factor(sc_analysis, levels = c("included","excluded")),
    sfc_analysis = factor(sfc_analysis, levels = c("included","excluded"))
  )

# Differences between included & excluded participants per analysis for
# age, sex, group

## Age: Wilcoxon test
clin_df2 %>% wilcox_test(Demo_Age ~ fc_analysis)
clin_df2 %>% wilcox_test(Demo_Age ~ sc_analysis)
clin_df2 %>% wilcox_test(Demo_Age ~ sfc_analysis)

# Sex: Chi2
chisq_test(clin_df2$Demo_Gender, clin_df2$fc_analysis)
chisq_test(clin_df2$Demo_Gender, clin_df2$sc_analysis)
chisq_test(clin_df2$Demo_Gender, clin_df2$sfc_analysis)

# Group
chisq_test(clin_df2$group, clin_df2$fc_analysis)
chisq_test(clin_df2$group, clin_df2$sc_analysis)
chisq_test(clin_df2$group, clin_df2$sfc_analysis)
