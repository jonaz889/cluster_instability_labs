#R/e3_func.R
e3_func <- function(sigma){
  gc()
  n <- 300
  k_true <- 3
  k <- 2
  N <- n * k_true
  exper_amount <- 200
  # Base noise
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  # True labels
  g <- rep(1:k_true, each = n)
  
  #dist <- seq(0, 2000, length.out=exper_amount)
  side_length <- 20
  pull <- seq(0,5, length.out = exper_amount)
  # Cluster parameters
  sx <- lapply(1:exper_amount, function(x){c(sigma,sigma, sigma)}) # x scale
  sy <- lapply(1:exper_amount, function(x){c(sigma,sigma, sigma)})     # y scale
  mx <- lapply(1:exper_amount, function(x){c(side_length/2, -side_length/2, 0)}) # x shift
  my <- lapply(1:exper_amount, function(x) {c(0, 0, (side_length * sqrt(3))/2 + pull[x])})  
  # Construct dataset
  X <- lapply(1:exper_amount, function(x){
    cbind(
      x_r * sx[[x]][g] + mx[[x]][g],
      y_r * sy[[x]][g] + my[[x]][g]
    )
  }
  )
  ## cluster that dataK-means
  kms <- lapply(1:exper_amount, function(x){kmeans(X[[x]], centers = k, nstart = 20)})
  
  ari <- unlist(lapply(1:exper_amount, function(t)
    mclust::adjustedRandIndex(g, kms[[t]]$cluster)
  ))
  
  return(list(
    X = X, side_length = side_length, pull = pull, g = g,
    kms = kms,
    clusters = lapply(1:exper_amount, function(t) kms[[t]]$cluster),
    ari = ari
  ))
}


bin_clusterings_by_ari <- function(sigma, R = 50, efunc) {
  run0 <- efunc(sigma)   
  
  pull <- run0$pull
  exper_amount <- length(pull)
  
  bins_per_pull <- vector("list", exper_amount)
  
  for (t in 1:exper_amount) {
    
    X_t <- run0$X[[t]]                         
    k <- nrow(run0$kms[[t]]$centers)    
    
    reps_clusters <- list()
    reps_X <- list()
    bin_counts <- integer(0)
    
    for (r in 1:R) {
      km <- kmeans(X_t, centers = k, nstart = 1)
      cl <- km$cluster
      
      if (length(reps_clusters) == 0) {
        reps_clusters[[1]] <- cl
        reps_X[[1]] <- X_t
        bin_counts[1] <- 1
      } else {
        ari_to_reps <- sapply(reps_clusters,function(repcl) mclust::adjustedRandIndex(repcl, cl))
        
        j <- which(ari_to_reps == 1)
        
        if (length(j) > 0) {
          bin_counts[j[1]] <- bin_counts[j[1]] + 1
        } else {
          reps_clusters[[length(reps_clusters) + 1]] <- cl
          reps_X[[length(reps_X) + 1]] <- X_t
          bin_counts[length(bin_counts) + 1] <- 1
        }
      }
    }
    
    bins_per_pull[[t]] <- list(
      pull = pull[t],
      n_bins = length(reps_clusters),
      counts = bin_counts,
      reps_clusters = reps_clusters,
      reps_X = reps_X
    )
  }
  
  list(pull = pull,bins_per_pull = bins_per_pull)
}


plot_reps_for_pull <- function(func_dat, t, cols, pch = 16) {
  bp <- func_dat$bins_per_pull[[t]]
  reps_X <- bp$reps_X
  reps_cl <- bp$reps_clusters
  counts <- bp$counts
  pull_val <- bp$pull
  
  # Square
  ncol <- ceiling(sqrt(bp$n_bins))
  nrow <- ceiling(bp$n_bins / ncol)
  
  op <- par(mfrow = c(ceiling(bp$n_bins / ceiling(sqrt(bp$n_bins))), ceiling(sqrt(bp$n_bins))), mar = c(3.5, 3.5, 3, 1))
  on.exit(par(op), add = TRUE)
  
  for (i in 1:bp$n_bins) {
    X <- reps_X[[i]]
    cl <- reps_cl[[i]]
    plot(X,
         col = cols[cl], pch = pch, asp = 1,
         xlab = expression(x[1]),ylab = expression(x[2]),
         cex.main = 1.5, cex.lab = 1.3,
         main = sprintf("pull = %.3f   freq = %d", pull_val, counts[i]))
  }
  
}


