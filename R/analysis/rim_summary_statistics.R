
rim_reps_sum <- function(runs, conf = 0.95) {
  reps <- length(runs)
  ks <- runs[[1]]$ks
  d <- runs[[1]]$d
  
  rim_mat <- sapply(runs, function(run) run$rim_rates)# rows:k , cols:replications
  
  mean_rim <- rowMeans(rim_mat)
  sd_rim <- apply(rim_mat, 1, sd)
  se_rim <- sd_rim / sqrt(reps)
  tcrit <- qt(1 - (1 - conf) / 2, df = reps - 1)
  
  data.frame(
    d=d,
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
    s$d <- rim_reps_over_d[[i]][[1]]$d
    summaries[[i]] <- s
  }
  
  do.call(rbind, summaries)
}
across_dimension_means_sum <- function(rim_reps_over_d_summary, conf = 0.95) {
  ks <- sort(unique(rim_reps_over_d_summary$k))
  ds <- sort(unique(rim_reps_over_d_summary$d))
  
  out <- data.frame(
    k = ks,
    mean = NA_real_,
    sd = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )
  
  tcrit <- qt(1 - (1 - conf) / 2, df = length(ds) - 1)
  
  for (i in seq_along(ks)) {
    vals <- rim_reps_over_d_summary$mean[rim_reps_over_d_summary$k == ks[i]]
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
