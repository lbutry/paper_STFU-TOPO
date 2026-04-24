# Author: Lionel Butry (@lbutry)

#######################################
# MST + global density thresholding
#######################################

# Created ranked matrix (MST + gThr)
mst_gthr_order <- function(mat) {
  
  # Build MST using inverse distance (MST tries to minimise cost but higher values in matrix indicate stronger connectivity)
  g <- graph_from_adjacency_matrix(mat, weighted = T, diag = F, mode = "undirected")
  E(g)$weight <- 1 / E(g)$weight # strength -> inverse distance
  
  mask_mst <- mst(g) |> as_adjacency_matrix(sparse = FALSE)
  mask_mst_lower <- mask_mst
  mask_mst_lower[upper.tri(mask_mst_lower)] <- 0 
  num_mst_edges <- sum(mask_mst_lower)
  
  # Get remaining non-MST edges (lower triangle only)
  remaining <- mat
  remaining[upper.tri(remaining) | mask_mst == 1] <- 0
  
  # Rank remaining edges by decreasing weight
  index <- sort(remaining[remaining > 0], decreasing = TRUE, index.return = TRUE)
  
  # Create MST + ranked matrix
  rank_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  rank_mat[remaining > 0][index$ix] <- seq_along(index$ix) # fill matrix with ranks; only lower triangle
  rank_mat[rank_mat > 0] <- rank_mat[rank_mat > 0] + num_mst_edges # add number of MST edges to rank
  rank_mat[mask_mst_lower == 1] <- seq_len(num_mst_edges) # add mst edges (indexed 1 to total MST edges)
  rank_mat[upper.tri(rank_mat)] <- t(rank_mat)[upper.tri(rank_mat)] # make symmetric
  
  return(rank_mat)
}

# Apply rank matrix (MST + global threshold)
apply_mst_gthr <- function(mat, index_mat, density) {
  
  # Get number of edges to keep for given density
  n <- nrow(mat)
  n_edges <- n * (n - 1) / 2
  n_keep <- ceiling(density * n_edges)
  
  # Apply threshold to matrix
  mask <- index_mat <= n_keep
  mat_thr <- mat * mask
  
  return(mat_thr)
}

#######################################
# Global network properties
#######################################

avg_strength <- function(x) {mean(strength(x), na.rm=T)}
avg_cc <- function(x) {mean(transitivity(x, type = "weighted"), na.rm=T)}

#######################################
# Plotting
#######################################

my_theme <- function() {
  theme_minimal(base_size = 14) +
    theme(
      text = element_text(size = 12, face = "bold", family = "Helvetica", color = "black"),
      line = element_line(linewidth = 1),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 0.75),
      panel.grid = element_blank(),
      axis.ticks.length = unit(1.5, "mm") 
    )
}

plot_gnp_line <- function(df, var_metric, y_label) { 
  var_metric <- rlang::sym(var_metric)
  data <- df %>%
    group_by(density, group) %>%
    summarize(mean = mean(!!var_metric),
              se = sqrt(sum((!!var_metric-mean(!!var_metric))^2/(length(!!var_metric)-1)))/sqrt(length(!!var_metric)),
              sd = sd(!!var_metric),
              .groups = "drop")
  
  plot <- ggplot(data, aes(x = density, y = mean, color = group)) +
    scale_color_manual(values = c("#0072B2", "#E69F00", "#CC79A7")) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  width = 0.005) +
    geom_line() +
    labs(x = "Cost", y = y_label) +
    theme(
      axis.text = element_text(size = 12),
    )
  
  return(plot)
}

plot_auc <- function(data, var_metric, y_label) {
  
  var_metric <- rlang::sym(var_metric)
  
  p <- ggplot(data, aes(x = group, y = !!var_metric)) +
    geom_boxplot(alpha = 0.6, aes(fill = group)) +
    geom_jitter(width = 0.15, alpha = 0.4) +
    scale_fill_manual(values = c("#0072B2", "#CC79A7", "#E69F00")) +
    labs(y = y_label) +
    my_theme() +
    theme(
      legend.position = "None",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12)
    )
  return(p)
}

plot_auc_bar <- function(data, var_metric, y_label) {
  var_metric <- rlang::sym(var_metric)
  
  df <- data |> 
    group_by(group) |> 
    summarise(
      mean_auc = mean(!!var_metric),
      se = sqrt(sum((!!var_metric-mean(!!var_metric))^2/(length(!!var_metric)-1)))/sqrt(length(!!var_metric)),
    )
  
  p <- ggplot(df, aes(x = group, y = mean_auc)) +
    geom_col(aes(fill = group)) +
    geom_errorbar(aes(ymin=mean_auc-se, ymax= mean_auc+se), width = 0.2) +
    scale_fill_manual(values = c("#0072B2", "#CC79A7", "#E69F00")) +
    labs(y = y_label) +
    my_theme() +
    theme(
      legend.position = "None",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12)
    )
  return(p)
}
