# R/rim_dimension_variability.R

dist_euclid <- function(a, b) sqrt(sum((a - b)^2))

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
  old_centers <- as.matrix(km_k$centers)
  new_centers <- as.matrix(km_k1$centers)
  
  k <- nrow(old_centers)
  
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

run_rim_once <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, seed = 1,
                         scale_mode = c("identity", "constant_total_variance")) {
  scale_mode <- match.arg(scale_mode)
  
  set.seed(seed)
  
  if (scale_mode == "identity") {
    X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  } else {
    X <- matrix(rnorm(n * d, sd = 1 / sqrt(d)), nrow = n, ncol = d)
  }
  
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
    rim_rates = rim_rates,
    scale_mode = scale_mode
  )
}

run_rim_reps <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, reps = 50, seed = 1,
                         scale_mode = c("identity", "constant_total_variance")) {
  scale_mode <- match.arg(scale_mode)
  
  runs <- vector("list", reps)
  for (r in seq_len(reps)) {
    runs[[r]] <- run_rim_once(
      n = n,
      d = d,
      ks = ks,
      nstart = nstart,
      seed = seed + r,
      scale_mode = scale_mode
    )
  }
  
  runs
}

run_rim_dimension_grid <- function(ds, n = 3000, ks = 1:10, nstart = 50, reps = 50, seed = 1,
                                   scale_mode = c("identity", "constant_total_variance")) {
  scale_mode <- match.arg(scale_mode)
  
  out <- vector("list", length(ds))
  names(out) <- paste0("d_", ds)
  
  for (i in seq_along(ds)) {
    out[[i]] <- run_rim_reps(
      n = n,
      d = ds[i],
      ks = ks,
      nstart = nstart,
      reps = reps,
      seed = seed + 1000 * i,
      scale_mode = scale_mode
    )
  }
  
  out
}

summarise_rim_reps <- function(runs, conf = 0.95) {
  reps <- length(runs)
  ks <- runs[[1]]$ks
  d <- runs[[1]]$d
  scale_mode <- runs[[1]]$scale_mode
  
  rim_mat <- sapply(runs, function(run) run$rim_rates)
  
  mean_rim <- rowMeans(rim_mat)
  sd_rim <- apply(rim_mat, 1, sd)
  se_rim <- sd_rim / sqrt(reps)
  tcrit <- qt(1 - (1 - conf) / 2, df = reps - 1)
  
  data.frame(
    d = d,
    scale_mode = scale_mode,
    k = ks[-length(ks)],
    mean = mean_rim,
    sd = sd_rim,
    lower = mean_rim - tcrit * se_rim,
    upper = mean_rim + tcrit * se_rim
  )
}

summarise_rim_grid <- function(grid_runs, conf = 0.95) {
  do.call(rbind, lapply(grid_runs, summarise_rim_reps, conf = conf))
}

plot_rim_mean_curves_by_dimension <- function(summary_tbl) {
  ds <- sort(unique(summary_tbl$d))
  ks <- sort(unique(summary_tbl$k))
  
  plot(NULL,
       xlim = range(ks),
       ylim = range(summary_tbl$mean),
       xlab = "k",
       ylab = "Mean rim data frequency",
       main = "Mean rim data frequency by dimension")
  
  for (i in seq_along(ds)) {
    df <- summary_tbl[summary_tbl$d == ds[i], ]
    lines(df$k, df$mean, type = "b", pch = 16, col = i, lwd = 2)
  }
  
  legend("topright",
         legend = paste0("d = ", ds),
         col = seq_along(ds),
         lwd = 2,
         pch = 16,
         bty = "n")
}

plot_rim_sd_curves_by_dimension <- function(summary_tbl) {
  ds <- sort(unique(summary_tbl$d))
  ks <- sort(unique(summary_tbl$k))
  
  plot(NULL,
       xlim = range(ks),
       ylim = range(summary_tbl$sd),
       xlab = "k",
       ylab = "SD of rim data frequency",
       main = "SD of rim data frequency by dimension")
  
  for (i in seq_along(ds)) {
    df <- summary_tbl[summary_tbl$d == ds[i], ]
    lines(df$k, df$sd, type = "b", pch = 16, col = i, lwd = 2)
  }
  
  legend("topright",
         legend = paste0("d = ", ds),
         col = seq_along(ds),
         lwd = 2,
         pch = 16,
         bty = "n")
}

plot_rim_ci_panels_by_dimension <- function(summary_tbl) {
  ds <- sort(unique(summary_tbl$d))
  n_panels <- length(ds)
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)
  
  op <- par(mfrow = c(nrow, ncol), mar = c(4, 4, 3, 1))
  on.exit(par(op), add = TRUE)
  
  for (d in ds) {
    df <- summary_tbl[summary_tbl$d == d, ]
    
    plot(df$k, df$mean,
         type = "l",
         ylim = range(c(df$lower, df$upper)),
         xlab = "k",
         ylab = "Rim frequency",
         main = paste0("d = ", d),
         lwd = 2)
    
    lines(df$k, df$lower, lty = 2)
    lines(df$k, df$upper, lty = 2)
    points(df$k, df$mean, pch = 16)
  }
}