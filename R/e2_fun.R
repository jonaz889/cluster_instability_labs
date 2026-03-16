# R/e1_bin_story.R
source("R/lall.R")

# Fixed data one run
e1_run_sigma <- function(
    sigma,
    x_r,
    y_r,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
    kmeans_nstart = 40,
    iter.max = 50,
    algorithm = "Hartigan-Wong",
    keep_all = TRUE
) {
  k <- 2
  N <- n * k
  g <- rep(1:k, each = n)
  
  dist <- seq(dist_min, dist_max, length.out = exper_amount)
  
  mis_rate <- numeric(exper_amount)
  tot_within <- numeric(exper_amount)
  
  if (keep_all) {
    kms <- vector("list", exper_amount)
    Xs  <- vector("list", exper_amount)
  }
  
  for (t in seq_len(exper_amount)) {
    mx <- c(dist[t] / 2, -dist[t] / 2)
    
    X <- cbind(
      x_r * sigma + mx[g],
      y_r * sigma
    )
    
    km <- kmeans(X,centers = 2,nstart = kmeans_nstart,iter.max = iter.max,algorithm = algorithm)
    
    tot_within[t] <- km$tot.withinss
    mis_rate[t] <- mis_rate_aligned_k2(g, km$cluster)
    
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
run_single_sigma_rep <- function(
    sigma,
    rep_seed,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
    kmeans_nstart = 25,
    keep_all = FALSE
) {
  N <- 2 * n
  
  set.seed(rep_seed)
  x_r <- rnorm(N)
  y_r <- rnorm(N)
  
  e1_run_sigma(
    sigma = sigma,
    x_r = x_r,
    y_r = y_r,
    n = n,
    dist_min = dist_min,
    dist_max = dist_max,
    exper_amount = exper_amount,
    kmeans_nstart = kmeans_nstart,
    keep_all = keep_all
  )
}


# Run above for each sigma
e1_sigma_sweep <- function(
    sigmas,
    seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
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
      e1_run_sigma(
        sigma = s,
        x_r = x_r,
        y_r = y_r,
        n = n,
        dist_min = dist_min,
        dist_max = dist_max,
        exper_amount = exper_amount,
        kmeans_nstart = kmeans_nstart,
        iter.max = iter.max,
        algorithm = algorithm,
        keep_all = keep_all
      )
    })
    
    names(reps_list[[rep]]) <- paste0("sigma_", sigmas)
  }
  
  reps_list
}

# Same as above but set both fixed random component
# AND seed
e1_run_sigma_kmeans_reps <- function(
    sigma,
    x_r,
    y_r,
    reps = 100,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 1000,
    kmeans_nstart = 25,
    seed = 1,
    keep_all = FALSE
) {
  rep_runs <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    set.seed(seed + r)
    
    rep_runs[[r]] <- e1_run_sigma(
      sigma = sigma,
      x_r = x_r,
      y_r = y_r,
      n = n,
      dist_min = dist_min,
      dist_max = dist_max,
      exper_amount = exper_amount,
      kmeans_nstart = kmeans_nstart,
      keep_all = keep_all
    )
  }
  
  rep_runs
}

# Identifies drops over a certain threshold
# Computes dist to center from clusterings
identify_drops <- function(runs, drop_size = 0.2) {
  
  dist_to_center <- vector("list", length(runs))
  drop_indices <- vector("list", length(runs))
  
  for (run in seq_along(runs)) {
    
    dist_to_center[[run]] <- t(vapply(seq_along(runs[[run]]$kms),
      function(i) {
        center <- runs[[run]]$kms[[i]]$centers
        c(abs(center[1,2]), abs(center[2,2])) # get only |y-values|
      },
      numeric(2)
    ))
    
    # Do +1 since diff
    drop_indices[[run]] <- which(abs(diff(runs[[run]]$mis_rate)) >= drop_size) + 1
  }
  
  list(drop_indices = drop_indices,dist_to_center = dist_to_center)
}

bin_story <- function(
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
  runs <- e1_sigma_sweep(
    sigmas = sigmas,
    seed = seed,
    n = n,
    dist_min = dist_min,
    dist_max = dist_max,
    exper_amount = exper_amount,
    kmeans_nstart = kmeans_nstart,
    reps = 1,
    iter.max = iter.max,
    algorithm = algorithm,
    keep_all = TRUE
  )
  
  drop_info <- identify_drops(runs, exper_amount, drop_size = 0.05)
  
  list(runs = runs, drop_info = drop_info)
}


screen_drop_runs <- function(
    sigma,
    reps = 1000,
    base_seed = 1,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 400,
    kmeans_nstart = 25,
    drop_threshold = 0.05
) {
  out <- vector("list", reps)
  
  for (r in seq_len(reps)) {
    rep_seed <- base_seed + r
    
    run <- run_single_sigma_rep(
      sigma = sigma,
      rep_seed = rep_seed,
      n = n,
      dist_min = dist_min,
      dist_max = dist_max,
      exper_amount = exper_amount,
      kmeans_nstart = kmeans_nstart,
      keep_all = FALSE
    )
    
    drop_sizes <- -diff(run$mis_rate)
    drop_indice <- which(drop_sizes >= drop_threshold)
    
    if (length(drop_indice) == 0) {
      out[[r]] <- data.frame(
        rep = r,
        rep_seed = rep_seed,
        has_drop = FALSE,
        drop_from_index = NA_integer_,
        drop_to_index = NA_integer_,
        drop_size = NA_real_,
        dist_before = NA_real_,
        dist_after = NA_real_
      )
    } else {
      out[[r]] <- data.frame(
        rep = r,
        rep_seed = rep_seed,
        has_drop = TRUE,
        drop_from_index = drop_indice[1],
        drop_to_index = drop_indice[1] + 1L,
        drop_size = drop_indice[i],
        dist_before = run$dist[i],
        dist_after = run$dist[i + 1L]
      )
    }
  }
  
  # Very smart way to turn into dataframe!!!
  do.call(rbind, out)
}

# RUN WITH SPECIFIC SEED KEEP ALL IF NOT USING THIS!!!! 
# OR PROBLEMS AGAIN......
rerun_drop_rep <- function(
    sigma,
    rep_seed,
    n = 1000,
    dist_min = 0,
    dist_max = 3,
    exper_amount = 400,
    kmeans_nstart = 25
) {
  run_single_sigma_rep(
    sigma = sigma,
    rep_seed = rep_seed,
    n = n,
    dist_min = dist_min,
    dist_max = dist_max,
    exper_amount = exper_amount,
    kmeans_nstart = kmeans_nstart,
    keep_all = TRUE
  )
}


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
  
  mis_mat <- sapply(rep_runs, function(run) run$mis_rate)
  # rows = distance points, cols = repetitions
  
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


summarise_selected_points <- function(df, x_points) {
  idx <- sapply(x_points, function(x) which.min(abs(df$dist_over_sigma - x)))
  
  out <- df[idx, c("sigma", "dist", "dist_over_sigma", "mean", "sd", "lower", "upper")]
  rownames(out) <- NULL
  out
}


find_first_drop <- function(run, drop_threshold = 0.05) {
  drop_sizes <- -diff(run$mis_rate)
  drop_index <- which( drop_sizes>= drop_threshold)

  if (length(drop_index) == 0) {
    return(NULL)
  }
  
  list(
    drop_from_index = drop_index[1],
    drop_to_index = drop_index[1] + 1L,
    drop_size = drop_sizes[drop_index[1]],
    dist_before = run$dist[drop_index[1]],
    dist_after = run$dist[drop_index[1] + 1],
    dist_over_sigma_before = run$dist[drop_index[1]] / run$sigma,
    dist_over_sigma_after = run$dist[drop_index[1] + 1] / run$sigma
  )
}


make_drop_table <- function(rep_runs, drop_threshold = 0.05) {
  out <- lapply(seq_along(rep_runs), function(r) {
    info <- find_first_drop(rep_runs[[r]], drop_threshold = drop_threshold)
    
    if (is.null(info)) {
      data.frame(
        rep = r,
        has_drop = FALSE,
        drop_from_index = NA_integer_,
        drop_to_index = NA_integer_,
        drop_size = NA_real_,
        dist_before = NA_real_,
        dist_after = NA_real_,
        dist_over_sigma_before = NA_real_,
        dist_over_sigma_after = NA_real_
      )
    } else {
      data.frame(
        rep = r,
        has_drop = TRUE,
        drop_from_index = info$drop_from_index,
        drop_to_index = info$drop_to_index,
        drop_size = info$drop_size,
        dist_before = info$dist_before,
        dist_after = info$dist_after,
        dist_over_sigma_before = info$dist_over_sigma_before,
        dist_over_sigma_after = info$dist_over_sigma_after
      )
    }
  })
  
  do.call(rbind, out)
}

theory_misclassification <- function(dist, sigma) {
  1 - pnorm(dist / (2 * sigma))
}

add_theory_to_drop_table <- function(drop_tbl, sigma) {
  drop_tbl$theory_before <- theory_misclassification(drop_tbl$dist_before, sigma)
  drop_tbl$theory_after  <- theory_misclassification(drop_tbl$dist_after, sigma)
  drop_tbl
}

center_axis_displacement <- function(km) {
  mean(abs(km$centers[, 2]))
}