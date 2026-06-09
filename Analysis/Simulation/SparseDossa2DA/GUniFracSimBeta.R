
get_script_dir <- function() {
  frames <- sys.frames()
  has_ofile <- vapply(frames, function(x) !is.null(x$ofile), logical(1))
  if (any(has_ofile)) {
    return(dirname(normalizePath(frames[[which(has_ofile)[length(which(has_ofile))]]]$ofile)))
  }
  
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  
  getwd()
}

get_env_value <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

parse_integer_grid <- function(value, default) {
  if (is.na(value) || !nzchar(value)) return(default)
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

SCRIPT_DIR <- get_script_dir()

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(SparseDOSSA2)
  library(MiSPU)   # for throat.otu.tab / throat.tree data in your original workflow
})

required_packages <- c(
  "compositions"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Please install these packages before running the DA simulation: ",
    paste(missing_packages, collapse = ", "),
    "\nFor Bioconductor packages, use BiocManager::install(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))."
  )
}

## -------------------- user / RStudio settings --------------------
FIT_RDS      <- get_env_value("FIT_RDS", file.path(SCRIPT_DIR, "SparseDOSSA2_fit_URT_lambda0.1.rds"))
OUTDIR       <- get_env_value("OUTDIR", file.path(SCRIPT_DIR, "da_sparse12_out"))
BASE_SEED    <- as.integer(get_env_value("BASE_SEED", "12345"))
N_SAMPLES    <- as.integer(get_env_value("N_SAMPLES", "100"))
N_REPS       <- as.integer(get_env_value("N_REPS", "100"))
N_WORKERS    <- as.integer(get_env_value("N_WORKERS", "8"))
PADJUST_METHOD <- get_env_value("PADJUST_METHOD", "BH")
SIGNIFICANCE_LEVEL <- as.numeric(get_env_value("SIGNIFICANCE_LEVEL", "0.05"))
K_COMM_GRID <- parse_integer_grid(get_env_value("K_COMM_GRID", ""), seq(2, 8, 2))

if (!PADJUST_METHOD %in% p.adjust.methods) {
  stop("PADJUST_METHOD must be one of: ", paste(p.adjust.methods, collapse = ", "))
}

if (is.na(N_SAMPLES) || N_SAMPLES < 4L) stop("N_SAMPLES must be at least 4.")
if (is.na(N_REPS) || N_REPS < 1L) stop("N_REPS must be a positive integer.")
if (is.na(N_WORKERS) || N_WORKERS < 1L) stop("N_WORKERS must be a positive integer.")
if (is.na(SIGNIFICANCE_LEVEL) || SIGNIFICANCE_LEVEL <= 0 || SIGNIFICANCE_LEVEL >= 1) {
  stop("SIGNIFICANCE_LEVEL must be between 0 and 1.")
}
if (!length(K_COMM_GRID) || any(is.na(K_COMM_GRID)) || any(K_COMM_GRID < 2L)) {
  stop("K_COMM_GRID must contain integers >= 2.")
}

if (!file.exists(FIT_RDS)) stop("FIT_RDS not found: ", FIT_RDS)

polar_cluster <- function(dist_matrix, k, 
                         max_cluster_size = Inf,
                         first_pole_method = c("median", "max", "random"),
                         subsequent_pole_quantile = 0.75,
                         assignment_method = c("nearest", "balanced"),
                         outlier_filter = TRUE,
                         outlier_threshold = 1.5) {
  #' Polar Ordination Clustering with Robust Pole Selection and Flexible Assignment
  #' 
  #' Performs clustering based on polar ordination with robust pole selection
  #' to avoid outliers and flexible assignment methods.
  #' 
  #' @param dist_matrix A distance matrix (dist object or matrix) with complete dimnames
  #' @param k Number of clusters (poles) to identify
  #' @param max_cluster_size Maximum number of units allowed per cluster (default: Inf)
  #' @param first_pole_method Method for selecting first pole: "median" (median of row sums),
  #'                          "max" (maximum of row sums), or "random" (random selection)
  #' @param subsequent_pole_quantile Quantile to use for subsequent pole selection (default: 0.75).
  #'                                 Use 1.0 for maximum, 0.75 for Q3 (robust), 0.5 for median
  #' @param assignment_method Method for assigning units to poles:
  #'                          "nearest" - assign to nearest pole with size constraints
  #'                          "balanced" - enforce equal cluster sizes (ignores max_cluster_size)
  #' @param outlier_filter If TRUE, exclude outliers when selecting poles (default: TRUE)
  #' @param outlier_threshold IQR multiplier for outlier detection (default: 1.5)
  #' @return A list of class 'polar_cluster' containing clustering results
  #' @examples
  #' dist_mat <- dist(matrix(rnorm(50), ncol = 5))
  #' # Robust method with balanced clusters
  #' result <- polar_cluster(dist_mat, k = 3, 
  #'                         first_pole_method = "median",
  #'                         assignment_method = "balanced",
  #'                         outlier_filter = TRUE)
  
  # Match arguments
  first_pole_method <- match.arg(first_pole_method)
  assignment_method <- match.arg(assignment_method)
  
  # Validate input parameters
  if (!inherits(dist_matrix, "dist") && !is.matrix(dist_matrix)) {
    stop("dist_matrix must be a dist object or matrix")
  }
  
  # Convert to matrix format if necessary
  if (inherits(dist_matrix, "dist")) {
    dist_mat <- as.matrix(dist_matrix)
  } else {
    dist_mat <- dist_matrix
  }
  
  # Check for complete dimnames
  if (is.null(rownames(dist_mat)) || is.null(colnames(dist_mat))) {
    stop("dist_matrix must have complete dimnames")
  }
  
  n <- nrow(dist_mat)
  unit_names <- rownames(dist_mat)
  
  if (k < 1 || k > n) {
    stop("k must be between 1 and the number of units")
  }
  
  if (max_cluster_size < 1) {
    stop("max_cluster_size must be at least 1")
  }
  
  if (subsequent_pole_quantile < 0 || subsequent_pole_quantile > 1) {
    stop("subsequent_pole_quantile must be between 0 and 1")
  }
  
  if (assignment_method == "nearest" && max_cluster_size * k < n) {
    warning("max_cluster_size * k < n: not all units can be assigned. ",
            "Increasing max_cluster_size to ", ceiling(n / k))
    max_cluster_size <- ceiling(n / k)
  }
  
  # ========== Detect and exclude outliers for pole selection ==========
  outlier_units <- character(0)
  pole_candidates <- unit_names
  
  if (outlier_filter) {
    # Calculate average distance for each unit
    row_sums <- rowSums(dist_mat)
    
    # Use IQR method to detect outliers
    Q1 <- quantile(row_sums, 0.25)
    Q3 <- quantile(row_sums, 0.75)
    IQR_val <- Q3 - Q1
    
    # Units with row sum > Q3 + threshold * IQR are considered outliers
    outlier_mask <- row_sums > (Q3 + outlier_threshold * IQR_val)
    outlier_units <- unit_names[outlier_mask]
    
    # Pole candidates exclude outliers
    pole_candidates <- unit_names[!outlier_mask]
    
    if (length(pole_candidates) < k) {
      warning("Too few non-outlier units (", length(pole_candidates), 
              ") for k = ", k, ". Using all units as candidates.")
      pole_candidates <- unit_names
      outlier_units <- character(0)
    }
  }
  
  # ========== Step 1: Select poles from non-outlier candidates ==========
  poles <- character(k)
  
  # Calculate row sums for candidates only
  candidate_row_sums <- rowSums(dist_mat[pole_candidates, , drop = FALSE])
  
  # Select first pole based on specified method
  if (first_pole_method == "median") {
    # Find candidate closest to median of row sums
    median_val <- median(candidate_row_sums)
    first_pole_idx <- which.min(abs(candidate_row_sums - median_val))
    poles[1] <- pole_candidates[first_pole_idx]
    
  } else if (first_pole_method == "max") {
    # Use maximum row sum among candidates
    poles[1] <- pole_candidates[which.max(candidate_row_sums)]
    
  } else if (first_pole_method == "random") {
    # Random selection from candidates
    set.seed(NULL)
    poles[1] <- sample(pole_candidates, 1)
  }
  
  # Select subsequent poles
  for (i in 2:k) {
    if (i == 2) {
      # Second pole: use quantile of distances to first pole (among candidates)
      dist_to_pole1 <- dist_mat[pole_candidates, poles[1]]
      
      if (subsequent_pole_quantile == 1.0) {
        # Use maximum
        poles[2] <- pole_candidates[which.max(dist_to_pole1)]
      } else {
        # Use quantile-based selection
        threshold <- quantile(dist_to_pole1, subsequent_pole_quantile)
        candidates_above_threshold <- pole_candidates[dist_to_pole1 >= threshold]
        
        if (length(candidates_above_threshold) == 1) {
          poles[2] <- candidates_above_threshold
        } else {
          poles[2] <- candidates_above_threshold[which.max(dist_to_pole1[pole_candidates %in% candidates_above_threshold])]
        }
      }
      
    } else {
      # For 3rd+ poles: compute minimum distance to existing poles
      min_dist_to_poles <- apply(dist_mat[pole_candidates, poles[1:(i-1)], drop = FALSE], 1, min)
      
      # Exclude already selected poles
      available_mask <- !pole_candidates %in% poles[1:(i-1)]
      available_dists <- min_dist_to_poles[available_mask]
      available_candidates <- pole_candidates[available_mask]
      
      if (subsequent_pole_quantile == 1.0) {
        # Use maximum
        poles[i] <- available_candidates[which.max(available_dists)]
      } else {
        # Use quantile-based selection
        threshold <- quantile(available_dists, subsequent_pole_quantile)
        candidates_above_threshold <- available_candidates[available_dists >= threshold]
        
        if (length(candidates_above_threshold) == 1) {
          poles[i] <- candidates_above_threshold
        } else {
          poles[i] <- candidates_above_threshold[which.max(available_dists[available_candidates %in% candidates_above_threshold])]
        }
      }
    }
  }
  
  # ========== Step 2: Assign ALL units (including outliers) to poles ==========
  clustering <- integer(n)
  names(clustering) <- unit_names
  reassignments <- 0
  
  # Calculate distance of each unit to each pole
  dist_to_poles_mat <- dist_mat[, poles, drop = FALSE]
  colnames(dist_to_poles_mat) <- 1:k
  
  if (assignment_method == "nearest") {
    # ===== Method 1: Assign to nearest pole with size constraints =====
    
    # Sort units by their minimum distance to any pole
    min_dist_to_any_pole <- apply(dist_to_poles_mat, 1, min)
    unit_order <- order(min_dist_to_any_pole)
    
    # Track current size of each cluster
    current_sizes <- integer(k)
    
    for (idx in unit_order) {
      unit <- unit_names[idx]
      dists <- dist_to_poles_mat[unit, ]
      pole_preference <- order(dists)
      
      assigned <- FALSE
      for (pole_id in pole_preference) {
        if (current_sizes[pole_id] < max_cluster_size) {
          clustering[unit] <- pole_id
          current_sizes[pole_id] <- current_sizes[pole_id] + 1
          assigned <- TRUE
          
          if (pole_id != pole_preference[1]) {
            reassignments <- reassignments + 1
          }
          break
        }
      }
      
      if (!assigned) {
        # Force assign to cluster with minimum current size
        pole_id <- which.min(current_sizes)
        clustering[unit] <- pole_id
        current_sizes[pole_id] <- current_sizes[pole_id] + 1
        reassignments <- reassignments + 1
        warning("Had to exceed max_cluster_size for cluster ", pole_id)
      }
    }
    
  } else if (assignment_method == "balanced") {
    # ===== Method 2: Force balanced cluster sizes =====
    
    # Target size for each cluster
    target_size <- floor(n / k)
    extra_units <- n %% k  # Some clusters will have one extra unit
    
    # Determine target size for each cluster
    cluster_targets <- rep(target_size, k)
    if (extra_units > 0) {
      cluster_targets[1:extra_units] <- target_size + 1
    }
    
    # Track current size of each cluster
    current_sizes <- integer(k)
    
    # Get all pairwise (unit, pole) combinations with distances
    assignment_options <- data.frame(
      unit = rep(unit_names, each = k),
      pole = rep(1:k, times = n),
      distance = as.vector(t(dist_to_poles_mat))
    )
    
    # Sort by distance (prefer closer assignments)
    assignment_options <- assignment_options[order(assignment_options$distance), ]
    
    # Greedy assignment: go through options in order of distance
    for (i in 1:nrow(assignment_options)) {
      unit <- assignment_options$unit[i]
      pole <- assignment_options$pole[i]
      
      # Skip if unit already assigned
      if (clustering[unit] > 0) next
      
      # Skip if this cluster is full
      if (current_sizes[pole] >= cluster_targets[pole]) next
      
      # Assign
      clustering[unit] <- pole
      current_sizes[pole] <- current_sizes[pole] + 1
      
      # Stop if all units assigned
      if (all(clustering > 0)) break
    }
    
    # If any units still unassigned (shouldn't happen), force assign
    unassigned <- which(clustering == 0)
    if (length(unassigned) > 0) {
      for (unit_idx in unassigned) {
        # Assign to cluster with minimum current size
        pole_id <- which.min(current_sizes)
        clustering[unit_names[unit_idx]] <- pole_id
        current_sizes[pole_id] <- current_sizes[pole_id] + 1
      }
    }
  }
  
  # ========== Step 3: Calculate additional statistics ==========
  
  # Distance of each unit to its assigned pole
  dist_to_pole <- numeric(n)
  names(dist_to_pole) <- unit_names
  for (i in 1:n) {
    unit <- unit_names[i]
    pole_id <- clustering[i]
    dist_to_pole[i] <- dist_mat[unit, poles[pole_id]]
  }
  
  # Size of each cluster
  cluster_sizes <- table(factor(clustering, levels = 1:k))
  
  # Within-cluster cohesion: average distance to pole for each cluster
  within_cluster_dist <- sapply(1:k, function(cluster_id) {
    units_in_cluster <- names(clustering)[clustering == cluster_id]
    if (length(units_in_cluster) == 0) return(NA)
    mean(dist_to_pole[units_in_cluster])
  })
  names(within_cluster_dist) <- paste0("Cluster_", 1:k)
  
  # Calculate silhouette coefficient for clustering quality assessment
  silhouette <- numeric(n)
  names(silhouette) <- unit_names
  
  for (i in 1:n) {
    unit <- unit_names[i]
    cluster_id <- clustering[i]
    
    # a(i): average distance to other units in same cluster
    same_cluster <- names(clustering)[clustering == cluster_id]
    if (length(same_cluster) > 1) {
      a_i <- mean(dist_mat[unit, setdiff(same_cluster, unit)])
    } else {
      a_i <- 0
    }
    
    # b(i): minimum average distance to units in other clusters
    other_clusters <- setdiff(1:k, cluster_id)
    if (length(other_clusters) > 0) {
      b_i <- min(sapply(other_clusters, function(other_id) {
        other_cluster <- names(clustering)[clustering == other_id]
        if (length(other_cluster) == 0) return(Inf)
        mean(dist_mat[unit, other_cluster])
      }))
    } else {
      b_i <- 0
    }
    
    # Silhouette coefficient: (b - a) / max(a, b)
    if (max(a_i, b_i) > 0 && is.finite(b_i)) {
      silhouette[i] <- (b_i - a_i) / max(a_i, b_i)
    } else {
      silhouette[i] <- 0
    }
  }
  
  # Average silhouette width across all units
  avg_silhouette <- mean(silhouette)
  
  # ========== Return results ==========
  result <- list(
    clustering = clustering,           # Named vector of cluster assignments
    poles = poles,                     # Names of pole units
    pole_indices = match(poles, unit_names),  # Numeric indices of poles
    dist_to_pole = dist_to_pole,      # Distance to assigned pole for each unit
    cluster_sizes = as.vector(cluster_sizes), # Size of each cluster
    within_cluster_dist = within_cluster_dist, # Average dist to pole per cluster
    silhouette = silhouette,           # Silhouette coefficient per unit
    avg_silhouette = avg_silhouette,   # Overall clustering quality metric
    reassignments = reassignments,     # Number of units reassigned
    max_cluster_size = max_cluster_size, # The size constraint used
    first_pole_method = first_pole_method,  # Method used for first pole
    subsequent_pole_quantile = subsequent_pole_quantile, # Quantile used
    assignment_method = assignment_method,  # Assignment method used
    outlier_filter = outlier_filter,   # Whether outlier filtering was used
    outlier_threshold = if(outlier_filter) outlier_threshold else NULL,
    outliers_detected = outlier_units, # Units identified as outliers
    n_outliers = length(outlier_units), # Number of outliers
    call = match.call(),               # Function call for reference
    k = k,                             # Number of clusters
    n = n                              # Number of units
  )
  
  class(result) <- "polar_cluster"
  return(result)
}

# ========== Print method ==========
print.polar_cluster <- function(x, ...) {
  cat("Polar Ordination Clustering\n")
  cat("===========================\n\n")
  cat("Number of clusters (k):", x$k, "\n")
  cat("Number of units:", x$n, "\n")
  cat("Assignment method:", x$assignment_method, "\n")
  
  if (x$assignment_method == "nearest") {
    cat("Max cluster size:", 
        if (is.infinite(x$max_cluster_size)) "Inf (no limit)" else x$max_cluster_size, 
        "\n")
  }
  
  cat("First pole method:", x$first_pole_method, "\n")
  cat("Subsequent pole quantile:", x$subsequent_pole_quantile, "\n")
  cat("Outlier filter:", x$outlier_filter, "\n")
  
  if (x$outlier_filter) {
    cat("Outliers detected:", x$n_outliers, "units\n")
    if (x$n_outliers > 0 && x$n_outliers <= 5) {
      cat("  (", paste(x$outliers_detected, collapse = ", "), ")\n", sep = "")
    }
  }
  
  cat("Average silhouette width:", round(x$avg_silhouette, 3), "\n")
  
  if (x$reassignments > 0) {
    cat("Reassignments:", x$reassignments, "\n")
  }
  cat("\n")
  
  cat("Poles:\n")
  for (i in 1:x$k) {
    cat(sprintf("  Cluster %d: %s (n=%d, avg.dist=%.3f)\n", 
                i, x$poles[i], x$cluster_sizes[i], x$within_cluster_dist[i]))
  }
  
  cat("\nClustering vector:\n")
  print(head(x$clustering, 10))
  if (x$n > 10) cat("  ... (", x$n - 10, "more units)\n")
  
  invisible(x)
}

# ========== Plot method ==========
plot.polar_cluster <- function(x, ...) {
  par(mfrow = c(1, 2))
  
  # Plot 1: Cluster sizes
  bp <- barplot(x$cluster_sizes, 
                names.arg = paste0("C", 1:x$k),
                main = paste("Cluster Sizes (", x$assignment_method, ")", sep = ""),
                xlab = "Cluster",
                ylab = "Number of Units",
                col = rainbow(x$k),
                ylim = c(0, max(x$cluster_sizes) * 1.1))
  
  # Add max size reference line if using nearest method and limit is finite
  if (x$assignment_method == "nearest" && is.finite(x$max_cluster_size)) {
    abline(h = x$max_cluster_size, col = "red", lty = 2, lwd = 2)
    text(bp[x$k], x$max_cluster_size, 
         labels = paste("Max =", x$max_cluster_size), 
         pos = 3, col = "red", cex = 0.8)
  }
  
  # Add mean line if using balanced method
  if (x$assignment_method == "balanced") {
    abline(h = mean(x$cluster_sizes), col = "blue", lty = 2, lwd = 2)
    text(bp[x$k], mean(x$cluster_sizes), 
         labels = paste("Mean =", round(mean(x$cluster_sizes), 1)), 
         pos = 3, col = "blue", cex = 0.8)
  }
  
  # Plot 2: Silhouette plot
  colors <- rainbow(x$k)[x$clustering]
  barplot(sort(x$silhouette, decreasing = TRUE),
          main = paste0("Silhouette Plot (avg = ", 
                        round(x$avg_silhouette, 3), ")"),
          xlab = "Units",
          ylab = "Silhouette Width",
          col = colors[order(x$silhouette, decreasing = TRUE)],
          border = NA)
  abline(h = 0, lty = 2)
  abline(h = x$avg_silhouette, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

# ========== Summary method ==========
summary.polar_cluster <- function(object, ...) {
  cat("Polar Ordination Clustering Summary\n")
  cat("====================================\n\n")
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Configuration:\n")
  cat("  Number of clusters (k):", object$k, "\n")
  cat("  Assignment method:", object$assignment_method, "\n")
  
  if (object$assignment_method == "nearest") {
    cat("  Max cluster size:", 
        if (is.infinite(object$max_cluster_size)) "Inf" else object$max_cluster_size, 
        "\n")
  }
  
  cat("  First pole method:", object$first_pole_method, "\n")
  cat("  Subsequent pole quantile:", object$subsequent_pole_quantile, "\n")
  cat("  Outlier filter:", object$outlier_filter, "\n")
  
  if (object$outlier_filter) {
    cat("  Outlier threshold (IQR):", object$outlier_threshold, "\n")
    cat("  Outliers detected:", object$n_outliers, "\n")
  }
  
  cat("  Reassignments:", object$reassignments, "\n\n")
  
  cat("Cluster Information:\n")
  cluster_info <- data.frame(
    Pole = object$poles,
    Size = object$cluster_sizes,
    Avg_Dist = round(object$within_cluster_dist, 3)
  )
  rownames(cluster_info) <- paste0("Cluster_", 1:object$k)
  print(cluster_info)
  cat("\n")
  
  cat("Silhouette Statistics:\n")
  cat("  Mean:  ", round(mean(object$silhouette), 3), "\n")
  cat("  Median:", round(median(object$silhouette), 3), "\n")
  cat("  Min:   ", round(min(object$silhouette), 3), "\n")
  cat("  Max:   ", round(max(object$silhouette), 3), "\n")
  
  if (object$n_outliers > 0) {
    cat("\nOutliers detected (", object$n_outliers, "):\n", sep = "")
    cat(" ", paste(object$outliers_detected, collapse = ", "), "\n")
  }
  
  invisible(object)
}


# Comprehensive k selection
select_optimal_k <- function(dist_matrix, k_range = 2:15,
                             first_pole_method = "median",
                             subsequent_pole_quantile = 0.75,
                             assignment_method = "balanced",
                             outlier_filter = TRUE,
                             min_acceptable_size = 5,
                             plot = TRUE) {
  #' Comprehensive selection of optimal k using multiple criteria
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @param min_acceptable_size Minimum acceptable cluster size
  #' @param plot If TRUE, generate diagnostic plots
  #' @return List with recommendations and diagnostic information
  
  cat("Testing k values from", min(k_range), "to", max(k_range), "...\n")
  
  # 1. Silhouette analysis
  cat("1. Running silhouette analysis...\n")
  sil_result <- find_optimal_k(dist_matrix, k_range, 
                               first_pole_method, subsequent_pole_quantile,
                               assignment_method, outlier_filter)
  
  # 2. Elbow method
  cat("2. Calculating within-cluster distances...\n")
  elbow_result <- calculate_wcsd(dist_matrix, k_range,
                                 first_pole_method, subsequent_pole_quantile,
                                 assignment_method, outlier_filter)
  
  # 3. Cluster size diagnostics
  cat("3. Analyzing cluster size distributions...\n")
  diagnostics <- lapply(sil_result$all_results, diagnose_k_quality, min_acceptable_size)
  
  # 4. Find k values with no small clusters
  valid_k <- k_range[sapply(diagnostics, function(d) d$n_small_clusters == 0)]
  
  # 5. Among valid k, find best silhouette
  if (length(valid_k) > 0) {
    valid_indices <- which(k_range %in% valid_k)
    best_valid_idx <- valid_indices[which.max(sil_result$silhouette_scores[valid_indices])]
    recommended_k <- k_range[best_valid_idx]
  } else {
    # If no k has zero small clusters, choose based on quality score
    quality_scores <- sapply(diagnostics, function(d) d$quality_score)
    recommended_k <- k_range[which.max(quality_scores)]
  }
  
  # Generate plots
  if (plot) {
    par(mfrow = c(2, 2))
    
    # Plot 1: Silhouette scores
    plot(k_range, sil_result$silhouette_scores, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Average Silhouette Width",
         main = "Silhouette Analysis", col = "blue", lwd = 2)
    abline(v = sil_result$optimal_k, col = "red", lty = 2)
    abline(h = 0.5, col = "gray", lty = 3)
    text(sil_result$optimal_k, max(sil_result$silhouette_scores), 
         labels = paste("k =", sil_result$optimal_k), pos = 4, col = "red")
    
    # Plot 2: Elbow plot
    plot(k_range, elbow_result$wcsd, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Within-Cluster Sum of Distances",
         main = "Elbow Method", col = "darkgreen", lwd = 2)
    
    # Plot 3: Number of small clusters
    n_small_vec <- sapply(diagnostics, function(d) d$n_small_clusters)
    plot(k_range, n_small_vec, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = paste("# Clusters < ", min_acceptable_size),
         main = "Small Cluster Count", col = "orange", lwd = 2)
    abline(h = 0, col = "red", lty = 2)
    
    # Plot 4: Quality score
    quality_scores <- sapply(diagnostics, function(d) d$quality_score)
    plot(k_range, quality_scores, type = "b", pch = 19,
         xlab = "Number of Clusters (k)", ylab = "Composite Quality Score",
         main = "Overall Quality", col = "purple", lwd = 2)
    abline(v = recommended_k, col = "red", lty = 2)
    text(recommended_k, max(quality_scores), 
         labels = paste("Recommended k =", recommended_k), pos = 4, col = "red")
    
    par(mfrow = c(1, 1))
  }
  
  # Print summary
  cat("\n=== K SELECTION SUMMARY ===\n")
  cat("Silhouette-optimal k:", sil_result$optimal_k, 
      "(silhouette =", round(sil_result$optimal_silhouette, 3), ")\n")
  cat("Valid k values (no small clusters):", 
      if(length(valid_k) > 0) paste(valid_k, collapse = ", ") else "None", "\n")
  cat("RECOMMENDED k:", recommended_k, "\n")
  
  cat("\nDetailed diagnostics by k:\n")
  for (i in seq_along(k_range)) {
    d <- diagnostics[[i]]
    cat(sprintf("k=%2d: Sil=%.3f, Small=%d/%d, CV=%.2f, Status=%s\n",
                d$k, d$avg_silhouette, d$n_small_clusters, d$k, 
                d$size_cv, d$recommendation))
  }
  
  return(list(
    recommended_k = recommended_k,
    silhouette_optimal_k = sil_result$optimal_k,
    valid_k_values = valid_k,
    silhouette_analysis = sil_result,
    elbow_analysis = elbow_result,
    diagnostics = diagnostics,
    all_results = sil_result$all_results
  ))
}

# Function to find optimal k using silhouette analysis
find_optimal_k <- function(dist_matrix, k_range = 2:10, 
                           first_pole_method = "median",
                           subsequent_pole_quantile = 0.75,
                           assignment_method = "balanced",
                           outlier_filter = TRUE) {
  #' Find optimal number of clusters using silhouette analysis
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test (default: 2:10)
  #' @param ... Other parameters passed to polar_cluster
  #' @return List with silhouette scores and recommended k
  
  results <- list()
  silhouette_scores <- numeric(length(k_range))
  cluster_size_stats <- list()
  
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    
    # Run clustering
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    
    results[[i]] <- result
    silhouette_scores[i] <- result$avg_silhouette
    
    # Track cluster size statistics
    cluster_size_stats[[i]] <- list(
      k = k_val,
      sizes = result$cluster_sizes,
      min_size = min(result$cluster_sizes),
      max_size = max(result$cluster_sizes),
      sd_size = sd(result$cluster_sizes),
      n_small = sum(result$cluster_sizes < 5)  # Number of small clusters
    )
  }
  
  # Find optimal k (maximum silhouette)
  optimal_idx <- which.max(silhouette_scores)
  optimal_k <- k_range[optimal_idx]
  
  return(list(
    k_range = k_range,
    silhouette_scores = silhouette_scores,
    optimal_k = optimal_k,
    optimal_silhouette = silhouette_scores[optimal_idx],
    all_results = results,
    cluster_size_stats = cluster_size_stats
  ))
}

# Function to calculate within-cluster sum of distances
calculate_wcsd <- function(dist_matrix, k_range = 2:10,
                           first_pole_method = "median",
                           subsequent_pole_quantile = 0.75,
                           assignment_method = "balanced",
                           outlier_filter = TRUE) {
  #' Calculate Within-Cluster Sum of Distances for different k
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @return List with WCSD values
  
  wcsd_values <- numeric(length(k_range))
  
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    
    # Calculate total within-cluster sum of distances
    wcsd_values[i] <- sum(result$dist_to_pole)
  }
  
  return(list(
    k_range = k_range,
    wcsd = wcsd_values
  ))
}

# Calculate gap statistic
calculate_gap_statistic <- function(dist_matrix, k_range = 2:10,
                                    B = 50,  # Number of bootstrap samples
                                    first_pole_method = "median",
                                    subsequent_pole_quantile = 0.75,
                                    assignment_method = "balanced",
                                    outlier_filter = TRUE) {
  #' Calculate Gap Statistic for optimal k selection
  #' 
  #' @param dist_matrix Distance matrix
  #' @param k_range Range of k values to test
  #' @param B Number of bootstrap reference datasets
  #' @return List with gap statistics
  
  n <- nrow(as.matrix(dist_matrix))
  dist_mat <- as.matrix(dist_matrix)
  
  observed_wcsd <- numeric(length(k_range))
  expected_wcsd <- matrix(0, nrow = length(k_range), ncol = B)
  
  # Calculate observed WCSD
  for (i in seq_along(k_range)) {
    k_val <- k_range[i]
    result <- polar_cluster(
      dist_matrix, 
      k = k_val,
      first_pole_method = first_pole_method,
      subsequent_pole_quantile = subsequent_pole_quantile,
      assignment_method = assignment_method,
      outlier_filter = outlier_filter
    )
    observed_wcsd[i] <- sum(result$dist_to_pole)
  }
  
  # Generate B reference datasets and calculate expected WCSD
  for (b in 1:B) {
    # Generate random distance matrix with similar properties
    # Using uniform distribution in the range of original distances
    rand_dist <- matrix(runif(n * n, min(dist_mat), max(dist_mat)), n, n)
    rand_dist <- (rand_dist + t(rand_dist)) / 2  # Make symmetric
    diag(rand_dist) <- 0
    rownames(rand_dist) <- colnames(rand_dist) <- rownames(dist_mat)
    
    for (i in seq_along(k_range)) {
      k_val <- k_range[i]
      result <- polar_cluster(
        rand_dist, 
        k = k_val,
        first_pole_method = first_pole_method,
        subsequent_pole_quantile = subsequent_pole_quantile,
        assignment_method = assignment_method,
        outlier_filter = FALSE  # Don't filter outliers in random data
      )
      expected_wcsd[i, b] <- sum(result$dist_to_pole)
    }
  }
  
  # Calculate gap statistic
  gap <- log(rowMeans(expected_wcsd)) - log(observed_wcsd)
  
  # Calculate standard deviation
  sdk <- apply(log(expected_wcsd), 1, sd)
  sk <- sdk * sqrt(1 + 1/B)
  
  # Find optimal k using 1-SE rule
  # Choose smallest k such that Gap(k) >= Gap(k+1) - s_{k+1}
  optimal_k <- k_range[1]
  for (i in 1:(length(k_range) - 1)) {
    if (gap[i] >= gap[i + 1] - sk[i + 1]) {
      optimal_k <- k_range[i]
      break
    }
  }
  
  return(list(
    k_range = k_range,
    gap = gap,
    gap_se = sk,
    optimal_k = optimal_k
  ))
}

# Diagnose cluster quality
diagnose_k_quality <- function(result, min_acceptable_size = 5) {
  #' Diagnose the quality of a clustering result
  #' 
  #' @param result A polar_cluster result object
  #' @param min_acceptable_size Minimum acceptable cluster size
  #' @return List with diagnostic information
  
  sizes <- result$cluster_sizes
  n_small <- sum(sizes < min_acceptable_size)
  prop_small <- n_small / result$k
  
  # Calculate size balance (coefficient of variation)
  cv_size <- sd(sizes) / mean(sizes)
  
  # Identify problematic clusters
  problematic <- which(sizes < min_acceptable_size)
  
  quality_score <- result$avg_silhouette * (1 - prop_small) * (1 / (1 + cv_size))
  
  diagnosis <- list(
    k = result$k,
    n_units = result$n,
    avg_silhouette = result$avg_silhouette,
    n_small_clusters = n_small,
    prop_small_clusters = prop_small,
    size_cv = cv_size,
    problematic_clusters = problematic,
    quality_score = quality_score,
    recommendation = ifelse(n_small > result$k * 0.2, 
                            "TOO MANY SMALL CLUSTERS - Consider reducing k",
                            ifelse(result$avg_silhouette < 0.25,
                                   "WEAK CLUSTERING - Try different k or parameters",
                                   "ACCEPTABLE"))
  )
  
  return(diagnosis)
}

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## beta = 0 is the null / type-I-error setting.
beta_grid <- seq(0, 1, 0.05)
signal_mode_grid <- c("abundance", "prevalence", "both")

## -------------------- differential-abundance helpers --------------------
row_tss <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  rs <- rowSums(x)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  keep <- is.finite(rs) & rs > 0
  if (any(keep)) out[keep, ] <- x[keep, , drop = FALSE] / rs[keep]
  out
}

replace_zeros_matrix <- function(X, constant = 0.5) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  X[X == 0] <- constant
  X
}

normalize_rows <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  rs <- rowSums(X)
  if (any(rs <= 0)) stop("Some rows have non-positive row sums.")
  sweep(X, 1, rs, FUN = "/")
}

normalize_rows_allow_all_zero <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  n <- nrow(X)
  p <- ncol(X)
  rs <- rowSums(X)
  Pi <- matrix(NA_real_, nrow = n, ncol = p, dimnames = dimnames(X))
  positive_rows <- rs > 0
  zero_rows <- rs == 0
  
  if (any(positive_rows)) {
    Pi[positive_rows, ] <- sweep(
      X[positive_rows, , drop = FALSE],
      1,
      rs[positive_rows],
      FUN = "/"
    )
  }
  
  if (any(zero_rows)) {
    Pi[zero_rows, ] <- matrix(1 / p, nrow = sum(zero_rows), ncol = p)
  }
  
  Pi
}

get_component_variance_after_normalization <- function(data,
                                                       group_factor = 1,
                                                       zero_constant = 0.5) {
  X <- data[, -group_factor, drop = FALSE]
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  var_vec <- apply(Pi, 2, var)
  names(var_vec) <- colnames(Pi)
  var_vec
}

ahc_transformation <- function(data, group_factor = 1) {
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  Pi <- normalize_rows_allow_all_zero(X)
  k <- ncol(Pi)
  
  sqrt_Pi <- sqrt(Pi)
  sqrt_pi0 <- rep(1 / sqrt(k), k)
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  s_pi <- sqrt(1 - c_pi^2)
  scale_factor <- asin(s_pi) / s_pi
  scale_factor[s_pi == 0] <- 1
  
  ahc_data <- (sqrt_Pi - outer(c_pi, sqrt_pi0)) * scale_factor
  colnames(ahc_data) <- colnames(Pi)
  
  data.frame(Group = groups, ahc_data, check.names = FALSE)
}

reference_ahc_transformation <- function(data,
                                         group_factor = 1,
                                         ref_component) {
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  Pi <- normalize_rows_allow_all_zero(X)
  k <- ncol(Pi)
  
  if (missing(ref_component) || is.null(ref_component)) ref_component <- k
  if (ref_component < 1 || ref_component > k) stop("ref_component is out of range.")
  
  sqrt_Pi <- sqrt(Pi)
  sqrt_pi0 <- rep(1 / sqrt(k), k)
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  s_pi <- sqrt(1 - c_pi^2)
  scale_factor <- asin(s_pi) / s_pi
  scale_factor[s_pi == 0] <- 1
  
  ref_sqrt <- sqrt_Pi[, ref_component]
  rahc_data <- sweep(
    sqrt_Pi[, -ref_component, drop = FALSE],
    1,
    ref_sqrt,
    FUN = "-"
  )
  rahc_data <- rahc_data * scale_factor
  colnames(rahc_data) <- colnames(Pi)[-ref_component]
  
  data.frame(Group = groups, rahc_data, check.names = FALSE)
}

alr_transformation <- function(data,
                               group_factor = 1,
                               ref_component,
                               zero_constant = 0.5) {
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  p <- ncol(Pi)
  
  if (missing(ref_component) || is.null(ref_component)) ref_component <- p
  if (ref_component < 1 || ref_component > p) stop("ref_component is out of range.")
  
  alr_data <- log(Pi[, -ref_component, drop = FALSE] / Pi[, ref_component])
  colnames(alr_data) <- colnames(Pi)[-ref_component]
  
  data.frame(Group = groups, alr_data, check.names = FALSE)
}

clr_transformation <- function(data,
                               group_factor = 1,
                               zero_constant = 0.5) {
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  gm <- exp(rowMeans(log(Pi)))
  clr_data <- log(Pi / gm)
  colnames(clr_data) <- colnames(Pi)
  
  data.frame(Group = groups, clr_data, check.names = FALSE)
}

ilr_transformation <- function(data,
                               group_factor = 1,
                               zero_constant = 0.5) {
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  ilr_data <- as.matrix(compositions::ilr(compositions::acomp(Pi)))
  colnames(ilr_data) <- paste0("ILR", seq_len(ncol(ilr_data)))
  
  data.frame(Group = groups, ilr_data, check.names = FALSE)
}

choose_nondiff_reference <- function(reference_effects) {
  nondiff_index <- which(!reference_effects)
  if (length(nondiff_index) == 0L) stop("No non-differential feature is available as reference.")
  tail(nondiff_index, 1)
}

choose_random_reference <- function(reference_effects) {
  sample(seq_along(reference_effects), size = 1)
}

choose_signal_reference <- function(reference_effects) {
  signal_index <- which(reference_effects)
  if (length(signal_index) == 0L) stop("No intended signal feature is available as reference.")
  tail(signal_index, 1)
}

choose_last_reference <- function(reference_effects) {
  length(reference_effects)
}

choose_minvar_reference <- function(data,
                                    group_factor = 1,
                                    zero_constant = 0.5) {
  var_vec <- get_component_variance_after_normalization(
    data = data,
    group_factor = group_factor,
    zero_constant = zero_constant
  )
  which.min(var_vec)
}

get_reference_table <- function(data,
                                reference_effects,
                                true_effects,
                                zero_constant = 0.5) {
  var_vec <- get_component_variance_after_normalization(
    data = data,
    group_factor = 1,
    zero_constant = zero_constant
  )
  
  ref <- c(
    choose_nondiff_reference(reference_effects),
    choose_random_reference(reference_effects),
    choose_signal_reference(reference_effects),
    choose_last_reference(reference_effects),
    choose_minvar_reference(data, group_factor = 1, zero_constant = zero_constant)
  )
  
  data.frame(
    ref_strategy = c("nondiff_ref", "random_ref", "signal_ref", "last_ref", "minvar_ref"),
    ref = ref,
    ref_is_signal = reference_effects[ref],
    ref_is_true_effect = true_effects[ref],
    ref_variance = as.numeric(var_vec[ref]),
    stringsAsFactors = FALSE
  )
}

get_alr_reference_table <- function(data,
                                    reference_effects,
                                    true_effects,
                                    zero_constant = 0.5) {
  ref_tbl <- get_reference_table(data, reference_effects, true_effects, zero_constant)
  ref_tbl$transformation <- c(
    "ALR_nondiff_ref",
    "ALR_random_ref",
    "ALR_signal_ref",
    "ALR_last_ref",
    "ALR_minvar_ref"
  )
  ref_tbl[, c("transformation", "ref_strategy", "ref", "ref_is_signal",
              "ref_is_true_effect", "ref_variance")]
}

get_rahc_reference_table <- function(data,
                                     reference_effects,
                                     true_effects,
                                     zero_constant = 0.5) {
  ref_tbl <- get_reference_table(data, reference_effects, true_effects, zero_constant)
  ref_tbl$transformation <- c(
    "RAHC_nondiff_ref",
    "RAHC_random_ref",
    "RAHC_signal_ref",
    "RAHC_last_ref",
    "RAHC_minvar_ref"
  )
  ref_tbl[, c("transformation", "ref_strategy", "ref", "ref_is_signal",
              "ref_is_true_effect", "ref_variance")]
}

run_feature_ttests <- function(transformed_data,
                               p_adjust_method = PADJUST_METHOD) {
  groups <- factor(transformed_data$Group)
  X <- as.matrix(transformed_data[, -1, drop = FALSE])
  storage.mode(X) <- "numeric"
  group_tab <- table(groups)
  
  if (length(group_tab) != 2L) stop("t-test requires exactly two groups.")
  if (any(group_tab < 2L)) stop("Each group must have at least two samples.")
  
  raw_p <- vapply(seq_len(ncol(X)), function(j) {
    y <- X[, j]
    if (all(!is.finite(y)) || length(unique(y[is.finite(y)])) < 2L) return(NA_real_)
    tryCatch(stats::t.test(y ~ groups)$p.value, error = function(e) NA_real_)
  }, numeric(1))
  
  adj_p <- p.adjust(raw_p, method = p_adjust_method)
  
  data.frame(
    feature = colnames(X),
    raw_p = raw_p,
    adj_p = adj_p,
    rejected = !is.na(adj_p) & adj_p < SIGNIFICANCE_LEVEL,
    score = NA_real_,
    stringsAsFactors = FALSE
  )
}

align_truth <- function(p_table, true_effects) {
  if (!is.null(names(true_effects)) && all(p_table$feature %in% names(true_effects))) {
    return(as.logical(true_effects[p_table$feature]))
  }
  
  if (length(true_effects) != nrow(p_table)) {
    stop("Truth vector and p-value table have incompatible lengths.")
  }
  
  as.logical(true_effects)
}

safe_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) <= 1L || all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

calculate_detection_metrics <- function(p_table,
                                        true_effects,
                                        significance_level = SIGNIFICANCE_LEVEL) {
  truth <- align_truth(p_table, true_effects)
  rejected <- as.logical(p_table$rejected)
  rejected[is.na(rejected)] <- FALSE
  
  TP <- sum(rejected & truth)
  FP <- sum(rejected & !truth)
  FN <- sum(!rejected & truth)
  n_signal <- sum(truth)
  n_null <- sum(!truth)
  
  power <- if (n_signal == 0L) NA_real_ else TP / n_signal
  fdr <- if ((TP + FP) == 0L) 0 else FP / (TP + FP)
  type1_error <- if (n_signal == 0L) as.integer(FP > 0L) else NA_real_
  false_positive_rate <- if (n_null == 0L) NA_real_ else FP / n_null
  
  list(
    Power = power,
    FDR = fdr,
    Type1Error = type1_error,
    FalsePositiveRate = false_positive_rate,
    TP = TP,
    FP = FP,
    FN = FN,
    Rejections = sum(rejected)
  )
}

make_method_result <- function(transformation,
                               p_table,
                               true_effects,
                               runtime_sec,
                               alr_ref = NA_integer_,
                               alr_ref_strategy = NA_character_,
                               alr_ref_is_signal = NA,
                               alr_ref_is_true_effect = NA,
                               alr_ref_variance = NA_real_,
                               method_error = NA_character_) {
  metrics <- calculate_detection_metrics(p_table, true_effects)
  truth <- align_truth(p_table, true_effects)
  
  signal_raw <- p_table$raw_p[truth]
  null_raw <- p_table$raw_p[!truth]
  signal_adj <- p_table$adj_p[truth]
  null_adj <- p_table$adj_p[!truth]
  
  data.frame(
    transformation = transformation,
    power = metrics$Power,
    fdr = metrics$FDR,
    type1_error = metrics$Type1Error,
    false_positive_rate = metrics$FalsePositiveRate,
    TP = metrics$TP,
    FP = metrics$FP,
    FN = metrics$FN,
    rejections = metrics$Rejections,
    n_tested = nrow(p_table),
    n_signal_tested = sum(truth),
    n_null_tested = sum(!truth),
    n_na_p = sum(is.na(p_table$raw_p)),
    min_raw_p_signal = safe_min(signal_raw),
    min_raw_p_null = safe_min(null_raw),
    min_adj_p_signal = safe_min(signal_adj),
    min_adj_p_null = safe_min(null_adj),
    n_raw_p_signal_less_005 = sum(signal_raw < 0.05, na.rm = TRUE),
    n_adj_p_signal_less_005 = sum(signal_adj < SIGNIFICANCE_LEVEL, na.rm = TRUE),
    alr_ref = alr_ref,
    alr_ref_strategy = alr_ref_strategy,
    alr_ref_is_signal = alr_ref_is_signal,
    alr_ref_is_true_effect = alr_ref_is_true_effect,
    alr_ref_variance = alr_ref_variance,
    runtime_sec = runtime_sec,
    method_error = method_error,
    stringsAsFactors = FALSE
  )
}

make_failed_method_result <- function(transformation,
                                      true_effects,
                                      runtime_sec,
                                      method_error,
                                      alr_ref = NA_integer_,
                                      alr_ref_strategy = NA_character_,
                                      alr_ref_is_signal = NA,
                                      alr_ref_is_true_effect = NA,
                                      alr_ref_variance = NA_real_) {
  n_tested <- length(true_effects)
  
  data.frame(
    transformation = transformation,
    power = NA_real_,
    fdr = NA_real_,
    type1_error = NA_real_,
    false_positive_rate = NA_real_,
    TP = NA_integer_,
    FP = NA_integer_,
    FN = NA_integer_,
    rejections = NA_integer_,
    n_tested = n_tested,
    n_signal_tested = sum(true_effects),
    n_null_tested = sum(!true_effects),
    n_na_p = n_tested,
    min_raw_p_signal = NA_real_,
    min_raw_p_null = NA_real_,
    min_adj_p_signal = NA_real_,
    min_adj_p_null = NA_real_,
    n_raw_p_signal_less_005 = NA_integer_,
    n_adj_p_signal_less_005 = NA_integer_,
    alr_ref = alr_ref,
    alr_ref_strategy = alr_ref_strategy,
    alr_ref_is_signal = alr_ref_is_signal,
    alr_ref_is_true_effect = alr_ref_is_true_effect,
    alr_ref_variance = alr_ref_variance,
    runtime_sec = runtime_sec,
    method_error = method_error,
    stringsAsFactors = FALSE
  )
}

run_method_safely <- function(transformation,
                              p_fun,
                              true_effects,
                              alr_ref = NA_integer_,
                              alr_ref_strategy = NA_character_,
                              alr_ref_is_signal = NA,
                              alr_ref_is_true_effect = NA,
                              alr_ref_variance = NA_real_) {
  t0 <- proc.time()[["elapsed"]]
  
  out <- tryCatch({
    p_table <- p_fun()
    t1 <- proc.time()[["elapsed"]]
    make_method_result(
      transformation = transformation,
      p_table = p_table,
      true_effects = true_effects,
      runtime_sec = t1 - t0,
      alr_ref = alr_ref,
      alr_ref_strategy = alr_ref_strategy,
      alr_ref_is_signal = alr_ref_is_signal,
      alr_ref_is_true_effect = alr_ref_is_true_effect,
      alr_ref_variance = alr_ref_variance
    )
  }, error = function(e) {
    t1 <- proc.time()[["elapsed"]]
    make_failed_method_result(
      transformation = transformation,
      true_effects = true_effects,
      runtime_sec = t1 - t0,
      method_error = conditionMessage(e),
      alr_ref = alr_ref,
      alr_ref_strategy = alr_ref_strategy,
      alr_ref_is_signal = alr_ref_is_signal,
      alr_ref_is_true_effect = alr_ref_is_true_effect,
      alr_ref_variance = alr_ref_variance
    )
  })
  
  out
}

run_one_dataset <- function(data,
                            true_effects,
                            reference_effects,
                            seed) {
  set.seed(seed)
  alr_ref_original <- choose_nondiff_reference(reference_effects)
  rahc_ref_tbl <- get_rahc_reference_table(data, reference_effects, true_effects, zero_constant = 0.5)
  alr_ref_tbl <- get_alr_reference_table(data, reference_effects, true_effects, zero_constant = 0.5)
  method_results <- list()
  idx <- 0L
  
  idx <- idx + 1L
  method_results[[idx]] <- run_method_safely(
    transformation = "AHC",
    p_fun = function() run_feature_ttests(ahc_transformation(data)),
    true_effects = true_effects
  )
  
  for (i in seq_len(nrow(rahc_ref_tbl))) {
    this_ref <- rahc_ref_tbl$ref[i]
    truth_RAHC <- true_effects[-this_ref]
    
    idx <- idx + 1L
    method_results[[idx]] <- run_method_safely(
      transformation = rahc_ref_tbl$transformation[i],
      p_fun = function(ref_i = this_ref) {
        run_feature_ttests(reference_ahc_transformation(data, ref_component = ref_i))
      },
      true_effects = truth_RAHC,
      alr_ref = this_ref,
      alr_ref_strategy = rahc_ref_tbl$ref_strategy[i],
      alr_ref_is_signal = rahc_ref_tbl$ref_is_signal[i],
      alr_ref_is_true_effect = rahc_ref_tbl$ref_is_true_effect[i],
      alr_ref_variance = rahc_ref_tbl$ref_variance[i]
    )
  }
  
  for (i in seq_len(nrow(alr_ref_tbl))) {
    this_ref <- alr_ref_tbl$ref[i]
    truth_ALR <- true_effects[-this_ref]
    
    idx <- idx + 1L
    method_results[[idx]] <- run_method_safely(
      transformation = alr_ref_tbl$transformation[i],
      p_fun = function(ref_i = this_ref) {
        run_feature_ttests(alr_transformation(data, ref_component = ref_i))
      },
      true_effects = truth_ALR,
      alr_ref = this_ref,
      alr_ref_strategy = alr_ref_tbl$ref_strategy[i],
      alr_ref_is_signal = alr_ref_tbl$ref_is_signal[i],
      alr_ref_is_true_effect = alr_ref_tbl$ref_is_true_effect[i],
      alr_ref_variance = alr_ref_tbl$ref_variance[i]
    )
  }
  
  idx <- idx + 1L
  method_results[[idx]] <- run_method_safely(
    transformation = "CLR",
    p_fun = function() run_feature_ttests(clr_transformation(data)),
    true_effects = true_effects
  )
  
  idx <- idx + 1L
  method_results[[idx]] <- run_method_safely(
    transformation = "ILR",
    p_fun = function() run_feature_ttests(ilr_transformation(data)),
    true_effects = true_effects[-alr_ref_original],
    alr_ref = alr_ref_original,
    alr_ref_strategy = "nondiff_ref_dimension_match",
    alr_ref_is_signal = reference_effects[alr_ref_original],
    alr_ref_is_true_effect = true_effects[alr_ref_original],
    alr_ref_variance = NA_real_
  )
  
  out <- do.call(rbind, method_results)
  rownames(out) <- NULL
  out
}

make_metadata_matrix <- function(y_vec) {
  y_vec <- as.integer(y_vec)
  matrix(
    y_vec,
    ncol = 1,
    dimnames = list(paste0("sample_", seq_along(y_vec)), "Y")
  )
}

make_spike_config <- function(feature_ids, signal_mode, effect_size) {
  props <- switch(
    signal_mode,
    abundance = "abundance",
    prevalence = "prevalence",
    both = c("abundance", "prevalence"),
    stop("Unknown signal_mode: ", signal_mode)
  )
  
  do.call(rbind, lapply(props, function(prop) {
    data.frame(
      metadata_datum = 1,
      feature_spiked = feature_ids,
      associated_property = prop,
      effect_size = effect_size,
      stringsAsFactors = FALSE
    )
  }))
}

upper_tri_values <- function(mat) {
  if (nrow(mat) <= 1L) return(numeric(0))
  mat[upper.tri(mat, diag = FALSE)]
}

get_boxplot_summary <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(data.frame(
      within_cop_n_pair = 0L,
      within_cop_mean = NA_real_,
      within_cop_sd = NA_real_,
      within_cop_min = NA_real_,
      within_cop_q1 = NA_real_,
      within_cop_median = NA_real_,
      within_cop_q3 = NA_real_,
      within_cop_max = NA_real_,
      within_cop_iqr = NA_real_,
      within_cop_lower_whisker = NA_real_,
      within_cop_upper_whisker = NA_real_,
      within_cop_n_outlier = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  bp <- boxplot.stats(x, coef = 1.5, do.conf = FALSE, do.out = TRUE)
  
  data.frame(
    within_cop_n_pair = length(x),
    within_cop_mean = mean(x),
    within_cop_sd = if (length(x) > 1L) sd(x) else NA_real_,
    within_cop_min = min(x),
    within_cop_q1 = safe_quantile(x, 0.25),
    within_cop_median = safe_quantile(x, 0.50),
    within_cop_q3 = safe_quantile(x, 0.75),
    within_cop_max = max(x),
    within_cop_iqr = IQR(x),
    within_cop_lower_whisker = bp$stats[1],
    within_cop_upper_whisker = bp$stats[5],
    within_cop_n_outlier = length(bp$out),
    stringsAsFactors = FALSE
  )
}

get_est_mean_abs <- function(fit_obj) {
  f <- fit_obj$EM_fit$fit
  if (is.null(f)) f <- fit_obj$fit
  if (is.null(f)) stop("Cannot find EM_fit$fit or fit in FIT_RDS object.")
  
  pi0 <- f$pi0
  mu <- f$mu
  
  if (!is.null(f$sigma)) {
    sigma <- f$sigma
  } else if (!is.null(f$sigma2)) {
    sigma <- sqrt(f$sigma2)
  } else {
    stop("Cannot find sigma or sigma2 in fit object.")
  }
  
  if (is.null(names(mu)) || is.null(names(pi0)) || is.null(names(sigma))) {
    stop("mu/pi0/sigma must be named by OTU IDs.")
  }
  
  common_ids <- Reduce(intersect, list(names(mu), names(pi0), names(sigma)))
  if (!length(common_ids)) stop("No overlapping OTU IDs across mu/pi0/sigma.")
  
  mu <- mu[common_ids]
  pi0 <- pi0[common_ids]
  sigma <- sigma[common_ids]
  
  (1 - pi0) * exp(mu + 0.5 * sigma^2)
}

make_analysis_frame <- function(counts_mat, y_vec) {
  counts_mat <- as.matrix(counts_mat)
  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("sample_", seq_len(nrow(counts_mat)))
  }
  
  groups <- factor(ifelse(as.integer(y_vec) == 0L, "A", "B"), levels = c("A", "B"))
  group_tab <- table(groups)
  
  if (length(group_tab) != 2L || any(group_tab < 2L)) {
    stop("Both groups must be present with at least two samples each.")
  }
  
  data.frame(Group = groups, counts_mat, check.names = FALSE)
}

make_seed <- function(rep_id, setting_id, offset = 0L) {
  seed <- (as.numeric(BASE_SEED) +
             as.numeric(offset) +
             as.numeric(rep_id) * 1000003 +
             as.numeric(setting_id) * 1009) %% .Machine$integer.max
  seed <- as.integer(seed)
  if (is.na(seed) || seed <= 0L) seed <- 1L
  seed
}

add_setting_metadata <- function(method_results,
                                 rep_id,
                                 setting_id,
                                 K_COMM,
                                 job,
                                 cluster_info,
                                 counts_sim,
                                 Y,
                                 signal_abs,
                                 seed_sim,
                                 seed_test,
                                 sim_runtime_sec) {
  meta <- data.frame(
    rep_id = rep_id,
    setting_id = setting_id,
    k = K_COMM,
    signal_cluster = job$cluster_idx,
    cluster_ordered = cluster_info$cluster_ordered,
    cluster_abund_sum_template = cluster_info$cluster_abund_sum_template,
    cluster_n_otus_template = cluster_info$cluster_n_otus_template,
    signal_mode = as.character(job$signal_mode),
    beta = job$beta,
    is_null = as.integer(job$beta == 0),
    n_samples = nrow(counts_sim),
    y_mean = mean(Y),
    signal_mean_abs = mean(signal_abs, na.rm = TRUE),
    signal_sd_abs = safe_sd(signal_abs),
    seed_sim = seed_sim,
    seed_test = seed_test,
    sim_runtime_sec = sim_runtime_sec,
    stringsAsFactors = FALSE
  )
  
  cluster_cols <- c(
    "cluster_avg_silhouette",
    "cluster_mean_relative_abundance",
    "cluster_mean_prevalence",
    "within_cop_n_pair",
    "within_cop_mean",
    "within_cop_sd",
    "within_cop_min",
    "within_cop_q1",
    "within_cop_median",
    "within_cop_q3",
    "within_cop_max",
    "within_cop_iqr",
    "within_cop_lower_whisker",
    "within_cop_upper_whisker",
    "within_cop_n_outlier"
  )
  
  out <- cbind(
    meta[rep(1L, nrow(method_results)), , drop = FALSE],
    method_results,
    cluster_info[rep(1L, nrow(method_results)), cluster_cols, drop = FALSE]
  )
  
  rownames(out) <- NULL
  out
}

run_one_setting <- function(rep_id,
                            setting_id,
                            K_COMM,
                            job,
                            pam_k,
                            cluster_summary_k) {
  cluster_idx_chr <- as.character(job$cluster_idx)
  C_otus <- names(pam_k)[pam_k == job$cluster_idx]
  if (!length(C_otus)) stop("Signal cluster has no OTUs.")
  
  cluster_info <- cluster_summary_k[
    cluster_summary_k$signal_cluster == job$cluster_idx,
    ,
    drop = FALSE
  ]
  
  seed_sim <- make_seed(rep_id, setting_id, offset = 0L)
  seed_test <- make_seed(rep_id, setting_id, offset = 2000000L)
  
  set.seed(make_seed(rep_id, setting_id, offset = 500000L))
  Y <- rbinom(N_SAMPLES, 1, 0.5)
  metadata_matrix <- make_metadata_matrix(Y)
  spike_config <- make_spike_config(
    feature_ids = C_otus,
    signal_mode = as.character(job$signal_mode),
    effect_size = job$beta
  )
  
  t0_sim <- proc.time()[["elapsed"]]
  set.seed(seed_sim)
  sim <- SparseDOSSA2::SparseDOSSA2(
    template = fit,
    n_sample = N_SAMPLES,
    new_features = FALSE,
    spike_metadata = spike_config,
    metadata_matrix = metadata_matrix,
    verbose = FALSE
  )
  t1_sim <- proc.time()[["elapsed"]]
  
  A_abs <- t(sim$simulated_matrices$a_spiked)
  counts_sim <- t(sim$simulated_data)
  
  tree <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(A_abs)))
  A_abs <- A_abs[, tree$tip.label, drop = FALSE]
  counts_sim <- counts_sim[, tree$tip.label, drop = FALSE]
  
  if (is.null(rownames(counts_sim))) {
    rownames(counts_sim) <- paste0("sample_", seq_len(nrow(counts_sim)))
  }
  rownames(A_abs) <- rownames(counts_sim)
  
  libsize <- rowSums(counts_sim)
  keep_samp <- is.finite(libsize) & libsize > 0
  if (!all(keep_samp)) {
    A_abs <- A_abs[keep_samp, , drop = FALSE]
    counts_sim <- counts_sim[keep_samp, , drop = FALSE]
    Y <- Y[keep_samp]
  }
  
  C_curr <- intersect(C_otus, colnames(counts_sim))
  if (!length(C_curr)) stop("No signal OTUs remain after simulation alignment.")
  
  analysis_data <- make_analysis_frame(counts_sim, Y)
  groups <- analysis_data$Group
  reference_effects <- colnames(counts_sim) %in% C_curr
  true_effects <- if (job$beta == 0) {
    rep(FALSE, ncol(counts_sim))
  } else {
    reference_effects
  }
  names(reference_effects) <- colnames(counts_sim)
  names(true_effects) <- colnames(counts_sim)
  
  signal_abs <- rowSums(A_abs[, C_curr, drop = FALSE])
  method_results <- run_one_dataset(
    data = analysis_data,
    true_effects = true_effects,
    reference_effects = reference_effects,
    seed = seed_test
  )
  
  add_setting_metadata(
    method_results = method_results,
    rep_id = rep_id,
    setting_id = setting_id,
    K_COMM = K_COMM,
    job = job,
    cluster_info = cluster_info,
    counts_sim = counts_sim,
    Y = Y,
    signal_abs = signal_abs,
    seed_sim = seed_sim,
    seed_test = seed_test,
    sim_runtime_sec = t1_sim - t0_sim
  )
}

run_one_replicate <- function(rep_id) {
  rows_out <- vector("list", total_settings)
  row_idx <- 0L
  
  for (K_COMM in K_COMM_GRID) {
    pam_k <- pam_list[[as.character(K_COMM)]]
    settings <- settings_list[[as.character(K_COMM)]]
    cluster_summary_k <- cluster_summary_list[[as.character(K_COMM)]]
    
    for (s in seq_len(nrow(settings))) {
      row_idx <- row_idx + 1L
      job <- settings[s, , drop = FALSE]
      rows_out[[row_idx]] <- run_one_setting(
        rep_id = rep_id,
        setting_id = row_idx,
        K_COMM = K_COMM,
        job = job,
        pam_k = pam_k,
        cluster_summary_k = cluster_summary_k
      )
    }
  }
  
  res_df <- do.call(rbind, rows_out[seq_len(row_idx)])
  rownames(res_df) <- NULL
  res_df
}

make_replicate_status <- function(rep_id,
                                  status,
                                  n_rows = NA_integer_,
                                  message = NA_character_) {
  data.frame(
    rep_id = rep_id,
    status = status,
    n_rows = n_rows,
    message = message,
    stringsAsFactors = FALSE
  )
}

run_one_replicate_safe <- function(rep_id) {
  tryCatch(
    {
      res_df <- run_one_replicate(rep_id)
      list(
        status = make_replicate_status(rep_id, "ok", nrow(res_df)),
        results = res_df
      )
    },
    error = function(e) {
      list(
        status = make_replicate_status(rep_id, "error", NA_integer_, conditionMessage(e)),
        results = NULL
      )
    }
  )
}

summarize_results <- function(all_replicate_results) {
  group_cols <- c(
    "k",
    "signal_cluster",
    "cluster_ordered",
    "signal_mode",
    "beta",
    "is_null",
    "transformation",
    "alr_ref_strategy"
  )
  
  cluster_cols <- c(
    "cluster_abund_sum_template",
    "cluster_n_otus_template",
    "cluster_avg_silhouette",
    "cluster_mean_relative_abundance",
    "cluster_mean_prevalence",
    "within_cop_n_pair",
    "within_cop_mean",
    "within_cop_sd",
    "within_cop_min",
    "within_cop_q1",
    "within_cop_median",
    "within_cop_q3",
    "within_cop_max",
    "within_cop_iqr",
    "within_cop_lower_whisker",
    "within_cop_upper_whisker",
    "within_cop_n_outlier"
  )
  
  metric_cols <- c(
    "power",
    "fdr",
    "type1_error",
    "false_positive_rate",
    "TP",
    "FP",
    "FN",
    "rejections",
    "n_tested",
    "n_signal_tested",
    "n_null_tested",
    "n_na_p",
    "min_raw_p_signal",
    "min_raw_p_null",
    "min_adj_p_signal",
    "min_adj_p_null",
    "n_raw_p_signal_less_005",
    "n_adj_p_signal_less_005",
    "alr_ref",
    "alr_ref_variance",
    "runtime_sec",
    "sim_runtime_sec"
  )
  
  key_parts <- lapply(all_replicate_results[group_cols], function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "<NA>"
    x
  })
  key <- do.call(paste, c(key_parts, sep = "\r"))
  idx_list <- split(seq_len(nrow(all_replicate_results)), key, drop = TRUE)
  
  summary_list <- lapply(idx_list, function(idx) {
    g <- all_replicate_results[idx, , drop = FALSE]
    out <- g[1L, group_cols, drop = FALSE]
    out <- cbind(out, g[1L, cluster_cols, drop = FALSE])
    
    out$n_reps <- length(unique(g$rep_id))
    out$n_rows <- nrow(g)
    out$n_method_error <- sum(!is.na(g$method_error) & nzchar(g$method_error))
    out$prop_method_error <- out$n_method_error / nrow(g)
    out$prop_alr_ref_signal <- safe_mean(as.numeric(g$alr_ref_is_signal))
    out$prop_alr_ref_true_effect <- safe_mean(as.numeric(g$alr_ref_is_true_effect))
    
    for (col in metric_cols) {
      out[[paste0("mean_", col)]] <- safe_mean(g[[col]])
      out[[paste0("sd_", col)]] <- safe_sd(g[[col]])
    }
    
    out
  })
  
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  list(
    summary = summary_df
  )
}

## -------------------- fixed manuscript-style setup --------------------
data("throat.otu.tab")
data("throat.tree")

counts <- as.matrix(throat.otu.tab)
common <- intersect(colnames(counts), throat.tree$tip.label)
counts <- counts[, common, drop = FALSE]
tree0  <- ape::drop.tip(throat.tree, setdiff(throat.tree$tip.label, colnames(counts)))

keep <- colSums(counts > 0) > 1
counts_f <- counts[, keep, drop = FALSE]
tree0    <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(counts_f)))

fit <- readRDS(FIT_RDS)

est_abs_mean <- get_est_mean_abs(fit)
missing_otus <- setdiff(colnames(counts_f), names(est_abs_mean))
if (length(missing_otus) > 0) {
  stop(
    "These OTUs are in counts_f but missing from est_abs_mean: ",
    paste(head(missing_otus, 20), collapse = ", ")
  )
}
est_abs_mean <- est_abs_mean[colnames(counts_f)]

## template-level relative abundance and prevalence
counts_f_tss <- row_tss(counts_f)
otu_mean_rel_template <- colMeans(counts_f_tss, na.rm = TRUE)
otu_mean_prev_template <- colMeans(counts_f > 0, na.rm = TRUE)

## -------------------- community distance matrix reused across K --------------------
cop <- ape::cophenetic.phylo(tree0)
cop <- cop[colnames(counts_f), colnames(counts_f)]

total_settings <- sum(K_COMM_GRID) * length(beta_grid) * length(signal_mode_grid)

pam_list              <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
rank_map_list         <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
abund_sum_list        <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
n_otus_list           <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
settings_list         <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
cluster_summary_list  <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))

cat(sprintf(
  "[START] DA simulation | reps=%d | settings/rep=%d | K grid=%s | N=%d | workers=%d\n",
  N_REPS, total_settings, paste(K_COMM_GRID, collapse = ","), N_SAMPLES, N_WORKERS
))
cat(sprintf("[P-ADJUST] method=%s | alpha=%.3f\n",
            PADJUST_METHOD, SIGNIFICANCE_LEVEL))

for (K_COMM in K_COMM_GRID) {
  
  ## -------------------- community definition --------------------
  pam_k <- polar_cluster(
    dist_matrix = as.matrix(cop),
    k = K_COMM,
    max_cluster_size = Inf,
    first_pole_method = "median",
    subsequent_pole_quantile = 1,
    assignment_method = "balanced",
    outlier_filter = TRUE,
    outlier_threshold = 1.5
  )$clustering
  
  abund_sum <- tapply(est_abs_mean[names(pam_k)], pam_k, sum, na.rm = TRUE)
  n_otus    <- tapply(pam_k, pam_k, length)
  
  ord_labels <- names(sort(abund_sum, decreasing = TRUE))
  rank_map   <- setNames(seq_along(ord_labels), ord_labels)
  
  ## silhouette
  sil_width_by_otu <- rep(NA_real_, length(pam_k))
  names(sil_width_by_otu) <- names(pam_k)
  sil_width_by_otu[] <- tryCatch({
    sil_obj <- cluster::silhouette(as.integer(pam_k), dist = as.dist(cop))
    as.numeric(sil_obj[, "sil_width"])
  }, error = function(e) {
    rep(NA_real_, length(pam_k))
  })
  names(sil_width_by_otu) <- names(pam_k)
  
  ## cluster summary table for this K
  cluster_summary_k <- vector("list", K_COMM)
  
  for (g in seq_len(K_COMM)) {
    otus_g <- names(pam_k)[pam_k == g]
    cop_g  <- cop[otus_g, otus_g, drop = FALSE]
    within_vals <- upper_tri_values(cop_g)
    bp_sum <- get_boxplot_summary(within_vals)
    
    cluster_summary_k[[g]] <- cbind(
      data.frame(
        k = K_COMM,
        signal_cluster = g,
        cluster_ordered = as.integer(rank_map[as.character(g)]),
        cluster_abund_sum_template = as.numeric(abund_sum[as.character(g)]),
        cluster_n_otus_template = as.integer(n_otus[as.character(g)]),
        cluster_avg_silhouette = safe_mean(sil_width_by_otu[otus_g]),
        cluster_mean_relative_abundance = safe_mean(otu_mean_rel_template[otus_g]),
        cluster_mean_prevalence = safe_mean(otu_mean_prev_template[otus_g]),
        stringsAsFactors = FALSE
      ),
      bp_sum
    )
  }
  
  cluster_summary_k <- do.call(rbind, cluster_summary_k)
  
  ## settings: one signal community at a time, binary Y, no confounder
  settings <- expand.grid(
    cluster_idx = seq_len(K_COMM),
    beta = beta_grid,
    signal_mode = signal_mode_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  pam_list[[as.character(K_COMM)]]             <- pam_k
  rank_map_list[[as.character(K_COMM)]]        <- rank_map
  abund_sum_list[[as.character(K_COMM)]]       <- abund_sum
  n_otus_list[[as.character(K_COMM)]]          <- n_otus
  settings_list[[as.character(K_COMM)]]        <- settings
  cluster_summary_list[[as.character(K_COMM)]] <- cluster_summary_k
  
  cat(sprintf("[K=%d] settings=%d\n", K_COMM, nrow(settings)))
}

rep_ids <- seq_len(N_REPS)
available_cores <- parallel::detectCores(logical = TRUE)
if (is.na(available_cores) || available_cores < 1L) available_cores <- 1L
actual_workers <- max(1L, min(N_WORKERS, length(rep_ids), available_cores))

start_time <- Sys.time()
cat(sprintf("[RUN] Launching %d replicate(s) on %d worker(s).\n", N_REPS, actual_workers))

combined_file <- file.path(OUTDIR, "replicate_results_da_sparse12_all.csv")
status_file <- file.path(OUTDIR, "da_sparse12_replicate_status.csv")
if (file.exists(combined_file)) file.remove(combined_file)
if (file.exists(status_file)) file.remove(status_file)

if (actual_workers > 1L) {
  cl <- parallel::makeCluster(actual_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(ape)
      library(cluster)
      library(SparseDOSSA2)
      library(MiSPU)
      library(compositions)
    })
    NULL
  })
  
  export_names <- setdiff(ls(), "cl")
  parallel::clusterExport(cl, varlist = export_names, envir = environment())
}

rep_batches <- split(rep_ids, ceiling(seq_along(rep_ids) / actual_workers))
status_accum <- vector("list", length(rep_batches))

for (b in seq_along(rep_batches)) {
  batch_ids <- rep_batches[[b]]
  cat(sprintf(
    "[BATCH %d/%d] reps=%s\n",
    b, length(rep_batches), paste(batch_ids, collapse = ",")
  ))
  
  batch_list <- if (actual_workers > 1L) {
    parallel::parLapplyLB(cl, batch_ids, run_one_replicate_safe)
  } else {
    lapply(batch_ids, run_one_replicate_safe)
  }
  
  batch_status <- do.call(rbind, lapply(batch_list, `[[`, "status"))
  status_accum[[b]] <- batch_status
  
  status_df_so_far <- do.call(rbind, status_accum[seq_len(b)])
  write.csv(status_df_so_far, status_file, row.names = FALSE)
  
  if (any(batch_status$status != "ok")) {
    print(batch_status[batch_status$status != "ok", , drop = FALSE])
    stop("At least one replicate failed. See status file: ", status_file)
  }
  
  batch_results <- do.call(rbind, lapply(batch_list, `[[`, "results"))
  rownames(batch_results) <- NULL
  append_to_combined <- file.exists(combined_file)
  
  write.table(
    batch_results,
    file = combined_file,
    append = append_to_combined,
    col.names = !append_to_combined,
    row.names = FALSE,
    sep = ",",
    quote = TRUE,
    qmethod = "double"
  )
  
  rm(batch_list, batch_results)
  invisible(gc())
}

status_df <- do.call(rbind, status_accum)
write.csv(status_df, status_file, row.names = FALSE)

all_replicate_results <- read.csv(combined_file, stringsAsFactors = FALSE, check.names = FALSE)
summary_obj <- summarize_results(all_replicate_results)
rm(all_replicate_results)
invisible(gc())
summary_file <- file.path(OUTDIR, "summary_results_da_sparse12.csv")
write.csv(summary_obj$summary, summary_file, row.names = FALSE)

summary_rds <- file.path(OUTDIR, "summary_results_da_sparse12.rds")
saveRDS(
  list(
    summary = summary_obj$summary,
    replicate_status = status_df,
    settings_by_k = settings_list,
    pam_by_k = pam_list,
    rank_map_by_k = rank_map_list,
    abund_sum_by_k = abund_sum_list,
    n_otus_by_k = n_otus_list,
    cluster_summary_by_k = cluster_summary_list,
    K_COMM_grid = K_COMM_GRID,
    beta_grid = beta_grid,
    signal_mode_grid = signal_mode_grid,
    n_reps = N_REPS,
    n_samples = N_SAMPLES,
    n_workers = actual_workers,
    p_adjust_method = PADJUST_METHOD,
    significance_level = SIGNIFICANCE_LEVEL,
    combined_replicate_file = combined_file
  ),
  summary_rds
)

end_time <- Sys.time()

cat("\n[DONE] Differential-abundance simulation finished.\n")
cat("Start time:", as.character(start_time), "\n")
cat("End time:", as.character(end_time), "\n")
cat("Total runtime:", as.character(end_time - start_time), "\n\n")
cat("Replicate status:\n")
cat(status_file, "\n\n")
cat("Combined replicate CSV:\n")
cat(combined_file, "\n\n")
cat("Summary CSV:\n")
cat(summary_file, "\n\n")
cat("Summary RDS:\n")
cat(summary_rds, "\n\n")

cat("Quick check: overall method summary\n")
overall_cols <- c(
  "transformation",
  "mean_power",
  "mean_fdr",
  "mean_type1_error",
  "mean_rejections",
  "n_method_error",
  "prop_method_error"
)
overall_check <- summary_obj$summary[, intersect(overall_cols, colnames(summary_obj$summary)), drop = FALSE]
overall_check <- aggregate(
  overall_check[, setdiff(colnames(overall_check), "transformation"), drop = FALSE],
  by = list(transformation = overall_check$transformation),
  FUN = safe_mean
)
print(overall_check)
