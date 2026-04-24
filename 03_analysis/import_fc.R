# Author: Lionel Butry (@lbutry)
# Purpose: Import functional connectomes in a list

############################################
# ---- Import FC data ----
############################################

# Paths
path_conn <- "01_data/ready4analysis/f04_connectome_partial"

# Load data as named list of matrices
files <- list.files(path=path_conn, pattern = "\\.csv$", full.names = TRUE)

fc_mat_list <- lapply(files, function(file) {
  mat_raw <- read.csv(file)
  mat <- mat_raw[, -1]
  rownames(mat) <- mat_raw$ROI.Name
  mat <- as.matrix(mat)
})

ids <- sub(".*/sub-(.*?)_.*", "\\1", files)
names(fc_mat_list) <- ids

##################################
# Tidy workspace
##################################

rm(files, ids, path_conn)
