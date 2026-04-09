# This file contains functions to simulate data from gaussian distribution with
# equal variance and run kmeans. Various types and wrappers are represented 
# used for different experiments.
run_kmeans_separation <- function(
    sigma,
    x_r,
    y_r,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    n_dist = 1000,
    kmeans_nstart = 40,
    iter.max = 50,
    algorithm = "Hartigan-Wong",
    keep_all = TRUE) {
  k <- 2
  N <- n * k
  true_lab <- rep(1:k, each = n)
  
  dist <- seq(dist_min, dist_max, length.out = n_dist)
  
  mis_rate <- numeric(n_dist)
  tot_within <- numeric(n_dist)
  
  if (keep_all) {
    kms <- vector("list", n_dist)
    Xs  <- vector("list", n_dist)
  }
  
  for (t in seq_len(n_dist)) {
    mx <- c(dist[t] / 2, -dist[t] / 2)
    
    X <- cbind(
      x_r * sigma + mx[true_lab],
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
    mis_rate[t] <- mis_rate_aligned_k2(true_lab, km$cluster)
    
    if (keep_all) {
      Xs[[t]] <- X
      kms[[t]] <- km
    }
  }
  
  out <- list(
    sigma = sigma,
    n = n,
    dist = dist,
    mis_rate = mis_rate,
    tot_within = tot_within,
    kmeans_nstart = kmeans_nstart,
    iter.max = iter.max,
    algorithm = algorithm
  )
  
  if (keep_all) {
    out$Xs <- Xs
    out$kms <- kms
  }
  
  out
}
# Same as above but save seed
run_kmeans_separation_seeded <- function(
    sigma,
    rep_seed,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    n_dist = 1000,
    kmeans_nstart = 25,
    keep_all = FALSE
) {
  N <- 2 * n
  
  set.seed(rep_seed)
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  run_kmeans_separation(
    sigma = sigma,
    x_r = x_r,
    y_r = y_r,
    n = n,
    dist_min = dist_min,
    dist_max = dist_max,
    n_dist = n_dist,
    kmeans_nstart = kmeans_nstart,
    keep_all = keep_all
  )
}

# Run above for each sigma
run_kmeans_separation_sigmas <- function(
    sigmas,
    seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    n_dist = 1000,
    kmeans_nstart = 40,
    reps = 1,
    iter.max = 50,
    algorithm = "Hartigan-Wong",
    keep_all = TRUE
) {
  k <- 2
  N <- n * k
  
  reps_list <- vector("list", reps)
  
  for (rep in seq_len(reps)) {
    set.seed(seed + rep)
    x_r <- rnorm(N)
    y_r <- rnorm(N)
    
    reps_list[[rep]] <- lapply(sigmas, function(s) {
      run_kmeans_separation(
        sigma = s,
        x_r = x_r,
        y_r = y_r,
        n = n,
        dist_min = dist_min,
        dist_max = dist_max,
        n_dist = n_dist,
        kmeans_nstart = kmeans_nstart,
        iter.max = iter.max,
        algorithm = algorithm,
        keep_all = keep_all
      )
    })
  }
  
  reps_list
}

run_kmeans_separation_sigmas_rand_across_sigmas <- function(
    sigmas,
    seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    n_dist = 1000,
    kmeans_nstart = 40,
    reps = 1,
    iter.max = 50,
    algorithm = "Hartigan-Wong",
    keep_all = TRUE
) {
  k <- 2
  N <- n * k
  
  reps_list <- vector("list", reps)
  
  for (rep in seq_len(reps)) {
    reps_list[[rep]] <- lapply(seq_along(sigmas), function(j) {
      s <- sigmas[j]
      
      set.seed(seed + 1000 * rep + j)
      x_r <- rnorm(N)
      y_r <- rnorm(N)
      
      run_kmeans_separation(
        sigma = s,
        x_r = x_r,
        y_r = y_r,
        n = n,
        dist_min = dist_min,
        dist_max = dist_max,
        n_dist = n_dist,
        kmeans_nstart = kmeans_nstart,
        iter.max = iter.max,
        algorithm = algorithm,
        keep_all = keep_all
      )
    })
  }
  
  reps_list
}

# Same as above but set both fixed random component
# AND seed
run_kmeans_separation_with_reps <- function(
    sigma,
    x_r,
    y_r,
    reps = 100,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    n_dist = 1000,
    kmeans_nstart = 25,
    seed = 1,
    keep_all = FALSE
) {
  rep_runs <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    set.seed(seed + r)
    
    rep_runs[[r]] <- run_kmeans_separation(
      sigma = sigma,
      x_r = x_r,
      y_r = y_r,
      n = n,
      dist_min = dist_min,
      dist_max = dist_max,
      n_dist = n_dist,
      kmeans_nstart = kmeans_nstart,
      keep_all = keep_all    )  }
  
  rep_runs
}