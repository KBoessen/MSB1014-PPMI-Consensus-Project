################################################################################
# MSB1014 Network Biology
# Final Project Code
################################################################################

# This file contains the main functions used to perform node ranking analysis on
# a precomputed dataset by Xu et al. (https://arxiv.org/abs/2211.12421)

# Used libraries ---------------------------------------------------------------
library(igraph)
library(dplyr)
library(manynet)
library(R.matlab)

#==============================================================================#
# FUNCTIONS
#==============================================================================#

# Subject Networks Function ----------------------------------------------------
threshold_network <- function(corr_matrix, method = "absolute", value = 0.3) {
  if (method == "absolute") {
    adj <- abs(corr_matrix) > value
  } else if (method == "percentile") {
    thresh <- quantile(abs(corr_matrix[upper.tri(corr_matrix)]), value)
    adj <- abs(corr_matrix) >= thresh
  } else if (method == "significance") {
    n <- nrow(corr_matrix)
    pvals <- matrix(NA, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          pvals[i,j] <- cor.test(corr_matrix[i,], corr_matrix[j,])$p.value
        } else {pvals[i,j] <- NA}
      }
    }
    adj <- (pvals < value)
  }
  diag(adj) <- 0
  return(adj)
}

# Jaccard Matrix Function ------------------------------------------------------
compute_jaccard_matrix <- function(networks_list, group_label = "Group") {
  library(igraph)
  
  jaccard_index <- function(g1, g2) {
    A <- as.matrix(as_adjacency_matrix(g1, sparse = FALSE))
    B <- as.matrix(as_adjacency_matrix(g2, sparse = FALSE))
    
    diag(A) <- diag(B) <- 0  # remove self-loops
    A_vec <- A[upper.tri(A)]
    B_vec <- B[upper.tri(B)]
    
    intersection <- sum(A_vec & B_vec)
    union <- sum(A_vec | B_vec)
    
    if (union == 0) return(0)
    return(intersection / union)
  }
  
  n <- length(networks_list)
  J <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in i:n) {
      J[i, j] <- jaccard_index(networks_list[[i]], networks_list[[j]])
      J[j, i] <- J[i, j]  # symmetric
    }
  }
  
  rownames(J) <- colnames(J) <- paste0(group_label, "_", 1:n)
  return(J)
}

# Weighted Consensus Networks Function------------------------------------------
consensus_network <- function(networks_list, Jaccard_mat = NULL, normalize = TRUE,
                              freq_thres = NULL, weight_thres_pct = NULL) {
  library(igraph)
  
  n <- length(networks_list)
  n_nodes <- vcount(networks_list[[1]])
  
  # Compute edge frequency (how many subjects have each edge)
  edge_counts <- matrix(0, n_nodes, n_nodes)
  for (g in networks_list) {
    mat <- as_adjacency_matrix(g, sparse = FALSE)
    edge_counts <- edge_counts + mat
  }
  
  # If no Jaccard provided, use equal weights
  if (is.null(Jaccard_mat)) {
    Jaccard_mat <- matrix(1, n, n)
  }
  
  # Compute average Jaccard weight per edge (subject-pair weighted sum)
  edge_weights <- matrix(0, n_nodes, n_nodes)
  for (i in 1:n) {
    for (j in i:n) {
      mat_i <- as.matrix(as_adjacency_matrix(networks_list[[i]], sparse = FALSE))
      mat_j <- as.matrix(as_adjacency_matrix(networks_list[[j]], sparse = FALSE))
      shared_edges <- mat_i * mat_j
      edge_weights <- edge_weights + Jaccard_mat[i, j] * shared_edges
    }
  }
  
  # Combine edge frequency and Jaccard-weighted overlap
  combined <- edge_counts + edge_weights
  
  # Optional edge frequency threshold
  if (!is.null(freq_thres)) {
    combined[edge_counts < freq_thres] <- 0
  }  
  
  # Optional percentile weight threshold
  if (!is.null(weight_thres_pct)) {
    weight_thres <- quantile(edge_weights, weight_thres_pct)
    combined[edge_weights < weight_thres] <- 0
  }
  
  # Optional normalization
  if (normalize) {
    combined <- combined / max(combined)
  }
  
  # No self-loops
  diag(combined) <- 0
  
  G_weighted <- graph_from_adjacency_matrix(combined, mode = "undirected", weighted = TRUE)
  E(G_weighted)$weight
  
  return(G_weighted)
}

# Centralities Function --------------------------------------------------------
compute_centralities <- function(G) {
  deg <- degree(G)
  btw <- betweenness(G)
  clo <- closeness(G)
  clo[!is.finite(clo)] <- 0
  katz <- node_alpha(G, alpha = 0.85)
  
  centralities <- data.frame(
    node = 1:vcount(G),
    degree = deg,
    betweenness = btw,
    closeness = clo,
    katz = katz
  )
  
  # Rank each centrality (1 = most central)
  ranked <- centralities %>%
    mutate(across(
      c(degree, betweenness, closeness, katz),
      ~ rank(-.x, ties.method = "average")
    ))
  
  # Average the ranks into a composite measure
  ranked$aggregate_rank <- rowMeans(ranked[, c("degree", "betweenness", "closeness", "katz")])
  
  # Re-rank the aggregate so it spans 1 ... N
  ranked$aggregate_rank <- rank(ranked$aggregate_rank, ties.method = "average")
  
  # Make disconnected nodes (degree = 0) N + 1
  ranked$aggregate_rank[centralities$degree == 0] <- vcount(G) + 1
  
  return(ranked)
}

#==============================================================================#
# PIPELINE 
#==============================================================================#

# Import Correlation Matrices --------------------------------------------------

# Path to folder of .mat files
mat_path <- ".../data/corr_mats"
files <- list.files(mat_path, pattern = "_correlation_matrix\\.mat$", full.names = FALSE)
varname <- "data"

# Load all matrices into a list
X_all <- lapply(files, function(f) {
  mat <- readMat(f)
  as.matrix(mat[[varname]])
})

# Extract patient ID numbers (patno)
extract_patno <- function(filename) {
  id <- sub("^sub-(.*?)_correlation_matrix\\.mat$", "\\1", basename(filename))
  return(id)
}
names(X_all) <- sapply(files, extract_patno)

# Identify file indices by class
idx_HC <- grep("control", names(X_all))
idx_PD <- grep("patient", names(X_all))
idx_PR <- grep("prodromal", names(X_all))
idx_SW <- grep("swedd", names(X_all))
idx_per_class <- list(idx_HC,idx_PD,idx_PR,idx_SW)
class_label_list <- list("HC","PD","PR","SW")

# Perform network inference for each class -------------------------------------

G_consensus_list = list()

for (i in 1:4) {
  class_label <- toString(class_label_list[i])
  X_list <- X_all[idx_per_class[[i]]]
  names(X_list) <- gsub("_correlation_matrix\\.mat", "", names(X_list))
  
  # Create networks from correlation matrix (0.95 percentile threshold)
  networks_list <- lapply(X_list, function(mat) graph_from_adjacency_matrix(threshold_network(mat, method="percentile", value=0.95), mode="undirected"))
  
  # Compute Jaccard similarity matrix
  J <- compute_jaccard_matrix(networks_list, group_label = class_label)
  plot_jaccard_heatmap(J, class_label)
  
  # Create consensus network from networks_list and J (0.95 percentile threshold)
  G_consensus <- consensus_network(networks_list, J, weight_thres_pct = 0.95)
  G_consensus_list[[class_label]] <- G_consensus
}

# Compute Node Ranks
centralities_list <- lapply(G_consensus_list, compute_centralities)

# Compute Paired Spearman Rank Coefficient Scores
for (i in 1:(n_groups-1)) {
  for (j in (i+1):n_groups) {
    g1 <- groups[i]; g2 <- groups[j]
    df_compare <- merge(centralities_list[[g1]][,c("node","aggregate_rank")],
                        centralities_list[[g2]][,c("node","aggregate_rank")],
                        by="node", suffixes=c(paste0("_",g1), paste0("_",g2)))
    cor_test <- cor.test(df_compare[[paste0("aggregate_rank_",g1)]],
                         df_compare[[paste0("aggregate_rank_",g2)]],
                         method="spearman",
                         exact = FALSE)
    write.csv(data.frame(group1=g1, group2=g2, rho=cor_test$estimate,
                         pvalue=cor_test$p.value),
              paste0("node_correlation_",g1,"_",g2,".csv"), row.names=FALSE)
  }
}
