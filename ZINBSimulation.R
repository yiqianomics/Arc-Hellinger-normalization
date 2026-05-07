############################################################
## Simulation: New AHC normalization vs ALR, CLR, ILR
## 100 replications, same setting as original script
############################################################

## =========================
## 1. Library path and packages
## =========================

.libPaths(c("/home/zhang.16383/RPackage_r451", .libPaths()))

user_lib <- "/home/zhang.16383/RPackage_r451"

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

is_installed <- function(pkg, lib) {
  pkg %in% rownames(installed.packages(lib.loc = lib))
}

install_and_load <- function(pkg, lib) {
  if (!is_installed(pkg, lib)) {
    install.packages(pkg, lib = lib)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE, lib.loc = lib)
  )
}

packages <- c("compositions", "dplyr", "tibble")

for (pkg in packages) {
  install_and_load(pkg, user_lib)
}


## =========================
## 2. Data generation
## =========================

generate_data <- function(n = 100, p = 50, beta0 = 5, beta = 1,
                          alpha = 3, q = 0.3) {
  
  if (n %% 2 != 0) stop("n must be even.")
  if (p %% 2 != 0) stop("p must be even.")
  
  n1 <- n / 2
  n2 <- n / 2
  
  x <- c(rep(1, n1), rep(0, n2))
  
  ## Significant columns: first p / 2 features
  mu_nonzero <- exp(beta0 + x * beta)
  size_nonzero <- 1 / alpha
  prob_nonzero <- 1 / (alpha * mu_nonzero + 1)
  
  xn <- replicate(
    p / 2,
    rnbinom(n, size = size_nonzero, prob = prob_nonzero)
  )
  
  ## Non-significant columns: last p / 2 features
  beta_zero <- 0
  mu_zero <- exp(beta0 + x * beta_zero)
  size_zero <- 1 / alpha
  prob_zero <- 1 / (alpha * mu_zero + 1)
  
  xz <- replicate(
    p / 2,
    rnbinom(n, size = size_zero, prob = prob_zero)
  )
  
  X <- cbind(xn, xz)
  
  ## Zero inflation
  pi_zero <- replicate(p, rbinom(n, 1, 1 - q))
  Xc <- X * pi_zero
  
  ## Library-size effect
  s <- rlnorm(n, meanlog = 1)
  Xc <- sweep(Xc, 1, s, FUN = "*")
  
  group <- ifelse(x == 1, "A", "B")
  
  colnames(Xc) <- paste0("X", seq_len(p))
  
  out <- data.frame(Group = group, Xc, check.names = FALSE)
  return(out)
}


## =========================
## 3. Helper functions
## =========================

replace_zeros_matrix <- function(X, constant = 0.5) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  X[X == 0] <- constant
  return(X)
}


normalize_rows <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  rs <- rowSums(X)
  
  if (any(rs <= 0)) {
    stop("Some rows have non-positive row sums.")
  }
  
  sweep(X, 1, rs, FUN = "/")
}


normalize_rows_allow_all_zero <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  n <- nrow(X)
  p <- ncol(X)
  
  rs <- rowSums(X)
  Pi <- matrix(NA_real_, nrow = n, ncol = p)
  
  positive_rows <- rs > 0
  zero_rows <- rs == 0
  
  Pi[positive_rows, ] <- sweep(X[positive_rows, , drop = FALSE],
                               1, rs[positive_rows], FUN = "/")
  
  ## If a row is all zero, define its composition as the uniform composition.
  ## This matches the behavior that zero replacement would give for an all-zero row.
  if (any(zero_rows)) {
    Pi[zero_rows, ] <- matrix(1 / p, nrow = sum(zero_rows), ncol = p)
  }
  
  colnames(Pi) <- colnames(X)
  rownames(Pi) <- rownames(X)
  
  return(Pi)
}


## =========================
## 4. Transformation functions
## =========================

## ---- ALR ----
alr_transformation <- function(data, group_factor = 1,
                               ref_component = NULL,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  p <- ncol(Pi)
  
  if (is.null(ref_component)) {
    ref_component <- p
  }
  
  alr_data <- log(Pi[, -ref_component, drop = FALSE] /
                    Pi[, ref_component])
  
  colnames(alr_data) <- colnames(Pi)[-ref_component]
  
  out <- data.frame(Group = groups, alr_data, check.names = FALSE)
  return(out)
}


## ---- CLR ----
clr_transformation <- function(data, group_factor = 1,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  gm <- exp(rowMeans(log(Pi)))
  clr_data <- log(Pi / gm)
  
  colnames(clr_data) <- colnames(Pi)
  
  out <- data.frame(Group = groups, clr_data, check.names = FALSE)
  return(out)
}


## ---- ILR ----
ilr_transformation <- function(data, group_factor = 1,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  ilr_data <- compositions::ilr(compositions::acomp(Pi))
  ilr_data <- as.matrix(ilr_data)
  
  colnames(ilr_data) <- paste0("ILR", seq_len(ncol(ilr_data)))
  
  out <- data.frame(Group = groups, ilr_data, check.names = FALSE)
  return(out)
}


## ---- New method: Angular Hellinger Compositional Normalization ----
ahc_transformation <- function(data, group_factor = 1) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  ## No zero replacement for ordinary zeros.
  ## Only all-zero rows are mapped to the uniform composition.
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
  
  out <- data.frame(Group = groups, ahc_data, check.names = FALSE)
  return(out)
}


## =========================
## 5. Testing functions
## =========================

run_feature_ttests <- function(transformed_data) {
  
  groups <- factor(transformed_data$Group)
  X <- transformed_data[, -1, drop = FALSE]
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  p_values <- numeric(ncol(X))
  
  for (j in seq_len(ncol(X))) {
    
    y <- X[, j]
    
    p_values[j] <- tryCatch({
      t.test(y ~ groups)$p.value
    }, error = function(e) {
      NA_real_
    })
  }
  
  adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
  
  return(adjusted_p_values)
}


get_truth_vector <- function(method_name, p, ref_component = p) {
  
  truth_full <- c(rep(TRUE, p / 2), rep(FALSE, p / 2))
  
  if (method_name == "ALR") {
    return(truth_full[-ref_component])
  }
  
  if (method_name == "ILR") {
    ## ILR returns p - 1 coordinates.
    ## This keeps the dimension consistent with the original simulation framework.
    return(truth_full[-ref_component])
  }
  
  return(truth_full)
}


calculate_power_fdr <- function(p_values, true_effects,
                                significance_level = 0.05) {
  
  if (length(p_values) != length(true_effects)) {
    stop("Length mismatch between p_values and true_effects.")
  }
  
  rejected <- !is.na(p_values) & p_values < significance_level
  
  TP <- sum(rejected & true_effects)
  FP <- sum(rejected & !true_effects)
  FN <- sum(!rejected & true_effects)
  
  total_signal <- sum(true_effects)
  
  power <- TP / total_signal
  
  fdr <- if ((TP + FP) == 0) {
    NA_real_
  } else {
    FP / (TP + FP)
  }
  
  out <- list(
    Power = power,
    FDR = fdr,
    TP = TP,
    FP = FP,
    FN = FN,
    Rejections = sum(rejected)
  )
  
  return(out)
}


run_one_replicate <- function(rep_id, n, p, beta0, beta, alpha, q,
                              significance_level = 0.05) {
  
  data <- generate_data(
    n = n,
    p = p,
    beta0 = beta0,
    beta = beta,
    alpha = alpha,
    q = q
  )
  
  transformation_functions <- list(
    AHC = ahc_transformation,
    ALR = alr_transformation,
    CLR = clr_transformation,
    ILR = ilr_transformation
  )
  
  replicate_results <- tibble()
  
  for (method_name in names(transformation_functions)) {
    
    trans_func <- transformation_functions[[method_name]]
    
    transformed_data <- trans_func(data)
    
    p_values <- run_feature_ttests(transformed_data)
    
    true_effects <- get_truth_vector(
      method_name = method_name,
      p = p,
      ref_component = p
    )
    
    metrics <- calculate_power_fdr(
      p_values = p_values,
      true_effects = true_effects,
      significance_level = significance_level
    )
    
    replicate_results <- bind_rows(
      replicate_results,
      tibble(
        rep_id = rep_id,
        transformation = method_name,
        power = metrics$Power,
        fdr = metrics$FDR,
        TP = metrics$TP,
        FP = metrics$FP,
        FN = metrics$FN,
        rejections = metrics$Rejections
      )
    )
  }
  
  return(replicate_results)
}


## =========================
## 6. Main simulation
## =========================

run_multiple_simulations <- function(num_reps = 100,
                                     n = 100,
                                     p = 50,
                                     significance_level = 0.05,
                                     seed = 12345) {
  
  set.seed(seed)
  
  alpha_grid <- seq(1, 10, by = 2)
  beta0_grid <- seq(1, 10, by = 2)
  beta_grid <- seq(1, 10, by = 2)
  q_grid <- seq(0, 0.9, by = 0.2)
  
  total_settings <- length(alpha_grid) *
    length(beta0_grid) *
    length(beta_grid) *
    length(q_grid)
  
  all_results <- tibble()
  setting_id <- 0
  
  for (alpha in alpha_grid) {
    for (beta0 in beta0_grid) {
      for (beta in beta_grid) {
        for (q in q_grid) {
          
          setting_id <- setting_id + 1
          
          cat(
            "Setting", setting_id, "/", total_settings,
            "| alpha =", alpha,
            "| beta0 =", beta0,
            "| beta =", beta,
            "| q =", q, "\n"
          )
          
          setting_results <- tibble()
          
          for (rep_id in seq_len(num_reps)) {
            
            one_result <- run_one_replicate(
              rep_id = rep_id,
              n = n,
              p = p,
              beta0 = beta0,
              beta = beta,
              alpha = alpha,
              q = q,
              significance_level = significance_level
            )
            
            one_result <- one_result %>%
              mutate(
                alpha = alpha,
                beta0 = beta0,
                beta = beta,
                q = q
              )
            
            setting_results <- bind_rows(setting_results, one_result)
          }
          
          all_results <- bind_rows(all_results, setting_results)
        }
      }
    }
  }
  
  return(all_results)
}


## =========================
## 7. Run simulation
## =========================

num_reps <- 100
n <- 100
p <- 50

all_simulation_results <- run_multiple_simulations(
  num_reps = num_reps,
  n = n,
  p = p,
  significance_level = 0.05,
  seed = 12345
)


## =========================
## 8. Summarize results
## =========================

summary_results <- all_simulation_results %>%
  group_by(alpha, beta0, beta, q, transformation) %>%
  summarize(
    mean_power = mean(power, na.rm = TRUE),
    sd_power = sd(power, na.rm = TRUE),
    mean_fdr = mean(fdr, na.rm = TRUE),
    sd_fdr = sd(fdr, na.rm = TRUE),
    mean_TP = mean(TP, na.rm = TRUE),
    mean_FP = mean(FP, na.rm = TRUE),
    mean_rejections = mean(rejections, na.rm = TRUE),
    .groups = "drop"
  )


## =========================
## 9. Save outputs
## =========================

write.csv(
  all_simulation_results,
  file = file.path(getwd(), "replicate_results_AHC_ALR_CLR_ILR_100rep.csv"),
  row.names = FALSE
)

write.csv(
  summary_results,
  file = file.path(getwd(), "summary_results_AHC_ALR_CLR_ILR_100rep.csv"),
  row.names = FALSE
)

cat("\nSimulation finished.\n")
cat("Replicate-level results saved to:\n")
cat(file.path(getwd(), "replicate_results_AHC_ALR_CLR_ILR_100rep.csv"), "\n")
cat("Summary results saved to:\n")
cat(file.path(getwd(), "summary_results_AHC_ALR_CLR_ILR_100rep.csv"), "\n")