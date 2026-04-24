# Author: Lionel Butry (@lbutry)

df2mat <- function(df, input_var) {
  
  #' Converts Long-DF to matrix
  #'
  #' @param df Long-DF. CAVE: First two columns represent regions
  #' @param input_var Name of ..value.. column as string
  #' @return Matrix
  
  rois <- unique(c(df[[1]], df[[2]]))
  mat <- matrix(NA_real_, length(rois), length(rois), dimnames = list(rois, rois))
  with(df, {
    i <- match(df[[1]], rois)
    j <- match(df[[2]], rois)
    values <- df[[input_var]]
    mat[cbind(i, j)] <<- values
    mat[cbind(j, i)] <<- values
  })
  return(mat)
}

# Function to convert one subject vector into matrix
vec2mat <- function(values, roi_pairs) {
  df <- tibble(
    roi1 = roi_pairs[,1],
    roi2 = roi_pairs[,2],
    value = values
  )
  df2mat(df, "value")
}

# Convert matrix to long df
mat2df <- function(mat) {
  mat |> as.data.frame() |> 
    rownames_to_column("roi1") |> 
    pivot_longer(-roi1, names_to = "roi2", values_to = "value") |> 
    filter(roi1 != roi2) # remove diagonal
}
