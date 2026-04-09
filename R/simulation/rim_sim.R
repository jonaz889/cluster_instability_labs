# Simulate rim data
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

# Generate normalized data.... ------

# Identity
std_data_1 <- function(n,d)
{
  X <- matrix(rnorm(n * d, sd = 1), nrow = n)
  X  
}

# Subtract mean and divide by sd
std_data_2 <- function(n,d)
{
  X <- matrix(rnorm(n * d, sd = 1), nrow = n)
  
  X <- scale(X)
  X
}

# Divide by norm
std_data_3 <- function(n,d)
{
  X <- matrix(rnorm(n * d), nrow = n)
  X / sqrt(rowSums(X^2))  
}



run_rim_once_scaled <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, seed = 1, scale_func=std_data_1) {
  set.seed(seed)
  X <- scale_func(n,d)
  
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

run_rim_reps_scaled <- function(n = 3000, d = 3, ks = 1:10, nstart = 50, reps = 50, seed = 1
                                , scale_func=std_data_1) {
  runs <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    runs[[r]] <- run_rim_once_scaled(
      n = n,
      d = d,
      ks = ks,
      nstart = nstart,
      seed = seed + r,
      scale_func = scale_func
    )
  }
  
  runs
}

run_rim_reps_over_d_scaled <- function(ds, n = 3000, ks = 1:10, nstart = 50, 
                                       reps = 50, seed = 1, scale_func=std_data_1) {
  runs <- vector("list", length(ds))
  
  for (i in seq_along(ds)) {
    runs[[i]] <- run_rim_reps_scaled(
      n = n,
      d = ds[i],
      ks = ks,
      nstart = nstart,
      reps = reps,
      seed = seed + 1000 * i,
      scale_func = scale_func
    )
  }
  
  runs
}
