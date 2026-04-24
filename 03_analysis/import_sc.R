# Author: Lionel Butry (@lbutry)
# Purpose: Import structural connectomes in a list

############################################
# ---- Import SC data ----
############################################

bn_lut <- read.delim("01_data/ready4analysis/BN_Atlas_246_LUT.tsv")

# Paths
path_conn <- "01_data/ready4analysis/d02_connectome"

# Load data as named list of matrices
files <- list.files(path=path_conn, pattern = "\\connectome_sift2.csv$", full.names = TRUE)

mat_list <- lapply(files, function(file) {
  mat <- as.matrix(read.csv(file, header = FALSE, row.names = NULL))
  dimnames(mat) <- list(bn_lut$ROI.Name, bn_lut$ROI.Name)
  return(mat)
})

ids <- sub(".*/sub-(.*?)_.*", "\\1", files)
names(mat_list) <- ids

##################################
# Tidy workspace
##################################

sc_mat_list <- mat_list
rm(bn_lut, mat_list, ids, files, path_conn)
