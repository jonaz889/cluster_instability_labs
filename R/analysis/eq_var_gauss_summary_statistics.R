
summarise_mis_ci_over_reps <- function(runs, conf = 0.95) {
  n_reps <- length(runs)
  n_sigmas <- length(runs[[1]])
  
  out <- vector("list", n_sigmas)
  
  for (j in seq_len(n_sigmas)) {
    mis_mat <- sapply(runs, function(rep_run) rep_run[[j]]$mis_rate)
    
    mean_mis <- rowMeans(mis_mat)
    sd_mis <- apply(mis_mat, 1, sd)
    se_mis <- sd_mis / sqrt(n_reps)
    tcrit <- qt(1 - (1 - conf) / 2, df = n_reps - 1)
    
    dist <- runs[[1]][[j]]$dist
    sigma <- runs[[1]][[j]]$sigma
    
    out[[j]] <- data.frame(
      sigma = sigma,
      dist = dist,
      dist_over_sigma = dist / sigma,
      mean = mean_mis,
      sd = sd_mis,
      se = se_mis,
      lower = mean_mis - tcrit * se_mis,
      upper = mean_mis + tcrit * se_mis,
      theory = theory_misclassification(dist, sigma)
    )
  }
  
  names(out) <- paste0("sigma_", sapply(out, function(x) x$sigma[1]))
  out
}

summarise_mis_ci_single_sigma <- function(rep_runs, conf = 0.95) {
  n_reps <- length(rep_runs)
  
  mis_mat <- sapply(rep_runs, function(run) run$mis_rate)# rows:distance points, cols:reps
  
  mean_mis <- rowMeans(mis_mat)
  sd_mis <- apply(mis_mat, 1, sd)
  se_mis <- sd_mis / sqrt(n_reps)
  
  tcrit <- qt(1 - (1 - conf) / 2, df = n_reps - 1)
  
  dist <- rep_runs[[1]]$dist
  sigma <- rep_runs[[1]]$sigma
  
  data.frame(
    sigma = sigma,
    dist = dist,
    dist_over_sigma = dist / sigma,
    mean = mean_mis,
    sd = sd_mis,
    se = se_mis,
    lower = mean_mis - tcrit * se_mis,
    upper = mean_mis + tcrit * se_mis
  )
}


summarise_selected_points <- function(sum_tab, x_points) {
  # Pick closest to selected points
  i <- sapply(x_points, function(x) which.min(abs(sum_tab$dist_over_sigma - x)))
  
  out <- sum_tab[i, c("sigma", "dist", "dist_over_sigma", "mean", "sd", "lower", "upper")]
  rownames(out) <- NULL
  out
}


theory_misclassification <- function(dist, sigma) {
  1 - pnorm(dist / (2 * sigma))
}


