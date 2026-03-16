# R/e1_gauss_distance.R
source("R/lall")

# two Gaussians, same noise across distance, only dist changes
e1_run_sigma <- function(
    sigma,
    seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
    kmeans_nstart = 40,
    iter.max = 50,
    algorithm = "Hartigan-Wong"
) {
  set.seed(seed)
  
  k <- 2
  N <- n * k
  g <- rep(1:k, each = n)
  
  # fixed random component
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  dist <- seq(dist_min, dist_max, length.out = exper_amount)
  
  mis_rate <- numeric(exper_amount)
  tot_within <- numeric(exper_amount)
  
  for (t in seq_len(exper_amount)) {
    mx <- c(dist[t] / 2, -dist[t] / 2)
    X <- cbind(
      x_r * sigma + mx[g],
      y_r * sigma
    )
    
    km <- kmeans(
      X,
      centers = 2,
      nstart = kmeans_nstart,
      iter.max = iter.max,
      algorithm = algorithm
    )
    
    tot_within[t] <- km$tot.withinss
    mis_rate[t] <- mis_rate_aligned_k2(g, km$cluster)
  }
  
  list(
    sigma = sigma,
    seed = seed,
    n = n,
    dist = dist,
    mis_rate = mis_rate,
    tot_within = tot_within,
    kmeans_nstart = kmeans_nstart,
    iter.max = iter.max,
    algorithm = algorithm
  )
}

# Run for each sigma
e1_sigma_sweep <- function(
    sigmas,
    seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
    kmeans_nstart = 40,
    iter.max = 50,
    algorithm = "Hartigan-Wong"
) {
  runs <- lapply(sigmas, function(s) {
    e1_run_sigma(
      sigma = s,
      seed = seed + round(1000 * s),
      n = n,
      dist_min = dist_min,
      dist_max = dist_max,
      exper_amount = exper_amount,
      kmeans_nstart = kmeans_nstart,
      iter.max = iter.max,
      algorithm = algorithm
    )
  })
  names(runs) <- paste0("sigma_", sigmas)
  runs
}