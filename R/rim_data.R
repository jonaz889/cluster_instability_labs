# R/rim_resampling.R

# Euclidian dist
dist_euclid <- function(a, b) sqrt(sum((a - b)^2))

# Alignment with dist mat
align_k_to_k1_leg <- function(km_k, km_k1) {
  k <- nrow(km_k$centers)
  
  dist_mat <- matrix(NA_real_, nrow = k, ncol = k + 1)
  
  for (i in 1:k) {
    for (j in 1:(k + 1)) {
      dist_mat[i, j] <- dist_euclid(km_k$centers[i, ], km_k1$centers[j, ])
    }
  }
  
  perms <- gtools::permutations(k, k)
  
  best_sum_distances <- numeric(k + 1)
  best_perm_for_i <- matrix(NA_integer_, nrow = k + 1, ncol = k)
  
  for (i in 1:(k + 1)) {
    dist_mat_wo_i <- as.matrix(dist_mat[, -i])
    
    best_sum <- Inf
    best_perm <- rep(NA_integer_, k)
    
    for (perm in 1:nrow(perms)) {
      total <- 0
      for (j in 1:k) {
        total <- total + dist_mat_wo_i[j, perms[perm, j]]
      }
      
      if (total < best_sum) {
        best_sum <- total
        best_perm <- perms[perm, ]
      }
    }
    
    best_sum_distances[i] <- best_sum
    best_perm_for_i[i, ] <- best_perm
  }
  
  best_removed_i <- which.min(best_sum_distances)
  
  pair_old_to_new <- numeric(k)
  perm <- best_perm_for_i[best_removed_i, ]
  
  for (i in 1:k) {
    pair_old_to_new[i] <- (1:(k + 1))[-best_removed_i][perm[i]]
  }
  
  new_to_old_map <- numeric(k + 1) + (k + 1)
  for (i in 1:k) {
    new_to_old_map[pair_old_to_new[i]] <- i
  }
  
  new_to_old_map[km_k1$cluster]
}

align_k_to_k1 <- function(km_k, km_k1) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required. Install it with install.packages('clue').")
  }
  
  old_centers <- as.matrix(km_k$centers)
  new_centers <- as.matrix(km_k1$centers)
  
  k <- nrow(old_centers)
  
  if (nrow(new_centers) != k + 1) {
    stop("km_k1 must have exactly one more center than km_k.")
  }
  
  if (ncol(old_centers) != ncol(new_centers)) {
    stop("Old and new centers must have the same number of columns.")
  }
  
  # Build k x (k+1) cost matrix of Euclidean distances
  cost_mat <- matrix(0, nrow = k, ncol = k + 1)
  for (i in 1:k) {
    for (j in 1:(k + 1)) {
      cost_mat[i, j] <- sqrt(sum((old_centers[i, ] - new_centers[j, ])^2))
    }
  }
  
  # Add dummy row with zero cost to every new cluster
  cost_aug <- rbind(cost_mat, rep(0, k + 1))
  
  # Solve assignment: each row gets one column
  assignment <- clue::solve_LSAP(cost_aug)
  
  # assignment[r] = column assigned to row r
  assigned_cols <- as.integer(assignment)
  
  # Map new cluster index -> old cluster label
  # Default unmatched cluster gets label k+1
  new_to_old_map <- rep(k + 1, k + 1)
  
  # Real rows 1:k correspond to old clusters
  for (old_idx in 1:k) {
    new_idx <- assigned_cols[old_idx]
    new_to_old_map[new_idx] <- old_idx
  }
  
  # Relabel each observation in km_k1$cluster
  new_to_old_map[km_k1$cluster]
}

rim_rate_between_k_and_k1 <- function(km_k, km_k1) {
  k <- nrow(km_k$centers)
  km_k1_aligned <- align_k_to_k1(km_k, km_k1)
  
  is_rim <- (km_k1_aligned != km_k$cluster) & (km_k1_aligned != k + 1)
  mean(is_rim)
}

run_rim_once <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  kms <- vector("list", max(ks))
  for (k in ks) {
    kms[[k]] <- kmeans(X, centers = k, nstart = nstart)
  }
  
  rim_rates <- numeric(max(ks) - 1)
  for (k in 1:(max(ks) - 1)) {
    rim_rates[k] <- rim_rate_between_k_and_k1(kms[[k]], kms[[k + 1]])
  }
  
  list(
    d = d,
    n = n,
    ks = ks,
    rim_rates = rim_rates
  )
}

run_rim_reps <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, reps = 50, seed = 1) {
  runs <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    runs[[r]] <- run_rim_once(
      n = n,
      d = d,
      ks = ks,
      nstart = nstart,
      seed = seed + r
    )
  }
  
  runs
}

run_rim_reps_over_d <- function(ds, n = 3000, ks = 1:10, nstart = 50, reps = 50, seed = 1) {
  runs <- vector("list", length(ds))
  
  for (d in seq_along(ds)) {
    runs[[d]] <- run_rim_reps(
      n = n,
      d = ds[d],
      ks = ks,
      nstart = nstart,
      reps=reps,
      seed = seed + 1000*d
    )
  }
      
  runs
}


standardize_data_unit_variance <- function(X)
{
  # total variance = 1, so coordinate variance = 1/d
    X <- scale(X, center = TRUE, scale = FALSE)
    X <- X / sqrt(sum(X^2) / (nrow(X) - 1))
}

run_rim_once_scaled <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, seed = 1) {
  set.seed(seed)
  
  X <- matrix(rnorm(n * d), nrow = n)
  
  X <- standardize_data_unit_variance(X)
  
  kms <- vector("list", max(ks))
  for (k in ks) {
    kms[[k]] <- kmeans(X, centers = k, nstart = nstart)
  }
  
  rim_rates <- numeric(length(ks) - 1)
  for (i in seq_len(length(ks) - 1)) {
    k <- ks[i]
    k1 <- ks[i + 1]
    rim_rates[i] <- rim_rate_between_k_and_k1(kms[[k]], kms[[k1]])
  }
  
  list(
    d = d,
    n = n,
    ks = ks,
    rim_rates = rim_rates
  )
}

run_rim_reps_scaled <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, reps = 50, seed = 1) {
  runs <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    runs[[r]] <- run_rim_once_scaled(
      n = n,
      d = d,
      ks = ks,
      nstart = nstart,
      seed = seed + r
    )
  }
  
  runs
}

run_rim_reps_over_d_scaled <- function(ds, n = 3000, ks = 1:10, nstart = 50, reps = 50, seed = 1) {
  runs <- vector("list", length(ds))
  
  for (i in seq_along(ds)) {
    runs[[i]] <- run_rim_reps_scaled(
      n = n,
      d = ds[i],
      ks = ks,
      nstart = nstart,
      reps = reps,
      seed = seed + 1000 * i
    )
  }
  
  runs
}


rim_reps_sum <- function(runs, conf = 0.95) {
  reps <- length(runs)
  ks <- runs[[1]]$ks
  
  rim_mat <- sapply(runs, function(run) run$rim_rates)
  # rows = k values, cols = replications
  
  mean_rim <- rowMeans(rim_mat)
  sd_rim <- apply(rim_mat, 1, sd)
  se_rim <- sd_rim / sqrt(reps)
  tcrit <- qt(1 - (1 - conf) / 2, df = reps - 1)
  
  data.frame(
    k = ks[-length(ks)],
    mean = mean_rim,
    sd = sd_rim,
    lower = mean_rim - tcrit * se_rim,
    upper = mean_rim + tcrit * se_rim
  )
}


rim_reps_over_d_sum <- function(rim_reps_over_d, conf = 0.95) {
  summaries <- vector("list", length(rim_reps_over_d))
  
  for (i in seq_along(rim_reps_over_d)) {
    s <- rim_reps_sum(rim_reps_over_d[[i]], conf = conf)
    d_val <- rim_reps_over_d[[i]][[1]]$d
    s$d <- d_val
    summaries[[i]] <- s
  }
  
  do.call(rbind, summaries)
}
across_dimension_means_sum <- function(summary_all, conf = 0.95) {
  ks <- sort(unique(summary_all$k))
  ds <- sort(unique(summary_all$d))
  
  out <- data.frame(
    k = ks,
    mean = NA_real_,
    sd = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )
  
  tcrit <- qt(1 - (1 - conf) / 2, df = length(ds) - 1)
  
  for (i in seq_along(ks)) {
    vals <- summary_all$mean[summary_all$k == ks[i]]
    m <- mean(vals)
    s <- sd(vals)
    se <- s / sqrt(length(vals))
    
    out$mean[i] <- m
    out$sd[i] <- s
    out$lower[i] <- m - tcrit * se
    out$upper[i] <- m + tcrit * se
  }
  
  out
}

plot_overall_dimension_mean_ci <- function(summary_dim_mean) {
  plot(summary_dim_mean$k, summary_dim_mean$mean,
       type = "l", lwd = 2,
       ylim = range(c(summary_dim_mean$lower, summary_dim_mean$upper)),
       xlab = "k",
       ylab = "Mean rim data frequency",
       main = "Mean rim curve across dimensions",
       cex.main = 1.4, cex.lab = 1.2)
  
  lines(summary_dim_mean$k, summary_dim_mean$lower, lty = 2)
  lines(summary_dim_mean$k, summary_dim_mean$upper, lty = 2)
  points(summary_dim_mean$k, summary_dim_mean$mean, pch = 16)
}
plot_rim_curves_by_d <- function(summary_all) {
  ds <- sort(unique(summary_all$d))
  ks <- sort(unique(summary_all$k))
  
  ylim <- range(c(summary_all$lower, summary_all$upper))
  
  plot(NULL,
       xlim = range(ks),
       ylim = ylim,
       xlab = "k",
       ylab = "Mean rim data frequency",
       main = "Mean rim curves for different dimensions",
       cex.main = 1.5, cex.lab = 1.3)
  
  for (i in seq_along(ds)) {
    sub <- summary_all[summary_all$d == ds[i], ]
    lines(sub$k, sub$mean, lwd = 2, col = i)
    points(sub$k, sub$mean, pch = 16, col = i)
  }
  
  legend("topleft",
         legend = paste0("d = ", ds),
         col = seq_along(ds),
         lwd = 2,
         pch = 16,
         bty = "n")
}
plot_all_rim_curves <- function(runs) {
  ks <- runs[[1]]$ks[-length(runs[[1]]$ks)]
  
  plot(NULL,
       xlim = range(ks),
       ylim = c(0, max(sapply(runs, function(run) max(run$rim_rates)))),
       xlab = "k",
       ylab = "Rim data frequency",
       main = paste0("Rim data frequency across different data, d = ", runs[[1]]$d),
       cex.main=1.5, cex.lab=1.3)
  
  for (r in seq_along(runs)) {
    lines(ks, runs[[r]]$rim_rates)
  }
}

plot_rim_ci <- function(summary_tbl, d) {
  plot(summary_tbl$k, summary_tbl$mean,
       type = "l", lwd = 2,
       ylim = range(c(summary_tbl$lower, summary_tbl$upper)),
       xlab = "k",
       ylab = "Rim data frequency",
       main = paste0("Mean and CI of rim data frequency, d = ", d),
       cex.main = 1.5, cex.lab = 1.3)
  
  lines(summary_tbl$k, summary_tbl$lower, lty = 2)
  lines(summary_tbl$k, summary_tbl$upper, lty = 2)
  points(summary_tbl$k, summary_tbl$mean, pch = 16)
}