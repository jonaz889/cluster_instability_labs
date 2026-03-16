standardize_data_unit_variance <- function(X) {
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- X / sqrt(sum(X^2) / (nrow(X) - 1))
  X
}

run_rim_once_over_d_nested_scaled <- function(
    n = 3000,
    ds = c(2, 4, 8, 16, 32, 64),
    ks = 1:10,
    nstart = 50,
    seed = 1
) {
  set.seed(seed)
  
  ds <- sort(unique(ds))
  d_max <- max(ds)
  
  X_max <- matrix(rnorm(n * d_max), nrow = n, ncol = d_max)
  
  runs <- vector("list", length(ds))
  
  for (i in seq_along(ds)) {
    d <- ds[i]

    X <- X_max[, 1:d, drop = FALSE]
    X <- standardize_data_unit_variance(X)
    #if (seed == 2) {
    cat("d =", d, " total variance =", sum(X^2) / (nrow(X) - 1), "\n")
    #}
    kms <- vector("list", max(ks))
    for (k in ks) {
      kms[[k]] <- kmeans(X, centers = k, nstart = nstart)
    }
    
    rim_rates <- numeric(length(ks) - 1)
    for (j in seq_len(length(ks) - 1)) {
      k <- ks[j]
      k1 <- ks[j + 1]
      rim_rates[j] <- rim_rate_between_k_and_k1(kms[[k]], kms[[k1]])
    }
    
    runs[[i]] <- list(
      d = d,
      n = n,
      ks = ks,
      rim_rates = rim_rates
    )
  }
  
  runs
}

run_rim_reps_over_d_nested_scaled <- function(
    n = 3000,
    ds = c(2, 4, 8, 16, 32, 64),
    ks = 1:10,
    nstart = 50,
    reps = 50,
    seed = 1
) {
  out <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    out[[r]] <- run_rim_once_over_d_nested_scaled(
      n = n,
      ds = ds,
      ks = ks,
      nstart = nstart,
      seed = seed + r
    )
  }
  
  out
}

rim_reps_over_d_nested_scaled_sum <- function(runs_nested, conf = 0.95) {
  ds <- sapply(runs_nested[[1]], function(x) x$d)
  
  summary_list <- vector("list", length(ds))
  
  for (i in seq_along(ds)) {
    runs_for_one_d <- lapply(runs_nested, function(rep_run) rep_run[[i]])
    tbl <- rim_reps_sum(runs_for_one_d, conf = conf)
    tbl$d <- ds[i]
    summary_list[[i]] <- tbl
  }
  
  do.call(rbind, summary_list)
}