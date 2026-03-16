# R/e1_drop_theory.R

theory_misclassification <- function(dist, sigma) {
  1 - pnorm(dist / (2 * sigma))
}

find_first_drop <- function(run, drop_threshold = 0.05) {
  drop_sizes <- -diff(run$mis_rate)
  drop_index <- which(drop_sizes >= drop_threshold)
  
  if (length(drop_index) == 0) {
    return(NULL)
  }
  
  i <- drop_index[1]
  
  list(
    drop_from_index = i,
    drop_to_index = i + 1L,
    drop_size = drop_sizes[i],
    dist_before = run$dist[i],
    dist_after = run$dist[i + 1L],
    dist_over_sigma_before = run$dist[i] / run$sigma,
    dist_over_sigma_after = run$dist[i + 1L] / run$sigma,
    theory_before = theory_misclassification(run$dist[i], run$sigma),
    theory_after = theory_misclassification(run$dist[i + 1L], run$sigma)
  )
}

make_drop_table_from_sweep <- function(runs, drop_threshold = 0.05) {
  # runs[[rep]][[sigma_index]]
  out <- list()
  idx <- 1L
  
  for (r in seq_along(runs)) {
    for (j in seq_along(runs[[r]])) {
      run <- runs[[r]][[j]]
      info <- find_first_drop(run, drop_threshold = drop_threshold)
      
      if (is.null(info)) {
        out[[idx]] <- data.frame(
          rep = r,
          sigma = run$sigma,
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
        out[[idx]] <- data.frame(
          rep = r,
          sigma = run$sigma,
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
      
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, out)
}

summarise_drop_table <- function(drop_tbl) {
  keep <- drop_tbl[drop_tbl$has_drop, ]
  if (nrow(keep) == 0) return(data.frame())
  
  split_tbl <- split(keep, keep$sigma)
  
  out <- lapply(split_tbl, function(df) {
    data.frame(
      sigma = df$sigma[1],
      n_drops = nrow(df),
      mean_dist_before = mean(df$dist_before, na.rm = TRUE),
      sd_dist_before = sd(df$dist_before, na.rm = TRUE),
      mean_r_before = mean(df$dist_over_sigma_before, na.rm = TRUE),
      sd_r_before = sd(df$dist_over_sigma_before, na.rm = TRUE)
      )
  })
  
  do.call(rbind, out)
}

summarise_mis_ci_over_reps_with_theory <- function(runs, conf = 0.95) {
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

extract_center_y_displacement <- function(runs) {
  # Only works if keep_all = TRUE
  out <- list()
  idx <- 1L
  
  for (r in seq_along(runs)) {
    for (j in seq_along(runs[[r]])) {
      run <- runs[[r]][[j]]
      
      if (is.null(run$kms)) next
      
      disp <- vapply(run$kms, function(km) {
        mean(abs(km$centers[, 2]))
      }, numeric(1))
      
      out[[idx]] <- data.frame(
        rep = r,
        sigma = run$sigma,
        dist = run$dist,
        dist_over_sigma = run$dist / run$sigma,
        center_y_disp = disp
      )
      
      idx <- idx + 1L
    }
  }
  
  do.call(rbind, out)
}

summarise_center_y_displacement <- function(center_tbl) {
  split_tbl <- split(center_tbl, center_tbl$sigma)
  
  out <- lapply(split_tbl, function(df) {
    agg_mean <- aggregate(center_y_disp ~ dist + dist_over_sigma + sigma, data = df, FUN = mean)
    agg_sd   <- aggregate(center_y_disp ~ dist + dist_over_sigma + sigma, data = df, FUN = sd)
    
    names(agg_mean)[names(agg_mean) == "center_y_disp"] <- "mean_disp"
    names(agg_sd)[names(agg_sd) == "center_y_disp"] <- "sd_disp"
    
    merge(agg_mean, agg_sd, by = c("dist", "dist_over_sigma", "sigma"))
  })
  
  do.call(rbind, out)
}

# R/plots_drop_theory.R

plot_mis_with_theory <- function(summary_list) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n <- length(summary_list)
  nr <- ceiling(n / 2)
  nc <- min(2, n)
  par(mfrow = c(nr, nc))
  
  for (j in seq_along(summary_list)) {
    df <- summary_list[[j]]
    
    ylim <- range(c(df$lower, df$upper, df$theory), na.rm = TRUE)
    
    plot(
      df$dist_over_sigma, df$mean,
      type = "l", lwd = 2,
      ylim = ylim,
      xlab = "d / sigma",
      ylab = "Misclassification",
      main = paste("sigma =", df$sigma[1])
    )
    
    lines(df$dist_over_sigma, df$lower, lty = 2)
    lines(df$dist_over_sigma, df$upper, lty = 2)
    lines(df$dist_over_sigma, df$theory, lwd = 2, lty = 3)
    
    legend(
      "topright",
      legend = c("Empirical mean", "95% CI", "Theory"),
      lty = c(1, 2, 3),
      lwd = c(2, 1, 2),
      bty = "n"
    )
  }
}

plot_drop_hist_r <- function(drop_tbl) {
  keep <- drop_tbl[drop_tbl$has_drop, ]
  sigmas <- sort(unique(keep$sigma))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n <- length(sigmas)
  nr <- ceiling(n / 2)
  nc <- min(2, n)
  par(mfrow = c(nr, nc))
  
  for (s in sigmas) {
    x <- keep$dist_over_sigma_before[keep$sigma == s]
    hist(
      x,
      breaks = 20,
      main = paste("Drop locations, sigma =", s),
      xlab = "d / sigma at first drop"
    )
    abline(v = mean(x, na.rm = TRUE), lwd = 2, lty = 2)
  }
}

plot_drop_hist_theory <- function(drop_tbl) {
  keep <- drop_tbl[drop_tbl$has_drop, ]
  sigmas <- sort(unique(keep$sigma))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n <- length(sigmas)
  nr <- ceiling(n / 2)
  nc <- min(2, n)
  par(mfrow = c(nr, nc))
  
  for (s in sigmas) {
    x <- keep$theory_before[keep$sigma == s]
    hist(
      x,
      breaks = 20,
      main = paste("Theory at drop, sigma =", s),
      xlab = "1 - Phi(d / (2 sigma)) at first drop"
    )
    abline(v = mean(x, na.rm = TRUE), lwd = 2, lty = 2)
  }
}

plot_center_disp_summary <- function(center_summary) {
  sigmas <- sort(unique(center_summary$sigma))
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  n <- length(sigmas)
  nr <- ceiling(n / 2)
  nc <- min(2, n)
  par(mfrow = c(nr, nc))
  
  for (s in sigmas) {
    df <- center_summary[center_summary$sigma == s, ]
    plot(
      df$dist_over_sigma, df$mean_disp,
      type = "l", lwd = 2,
      xlab = "d / sigma",
      ylab = "Mean vertical center displacement",
      main = paste("sigma =", s)
    )
  }
}


rm(list = ls())

source("R/utils-runlog.R")
source("R/lall.R")
source("R/load_all.R")

params <- list(
  sigmas = c(1),
  n = 500,
  dist_min = 0,
  dist_max = 8,
  exper_amount = 150,
  reps = 200,
  kmeans_nstart = 20,
  drop_threshold = 0.05
)

seed <- 1

meta <- start_run("drop_vs_theory", params = params, seed = seed)

runs <- e1_sigma_sweep(
  sigmas = params$sigmas,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  exper_amount = params$exper_amount,
  kmeans_nstart = params$kmeans_nstart,
  reps = params$reps,
  keep_all = FALSE
)

mis_summary <- summarise_mis_ci_over_reps_with_theory(runs)
drop_tbl <- make_drop_table_from_sweep(runs, drop_threshold = params$drop_threshold)
drop_summary <- summarise_drop_table(drop_tbl)

res <- list(
  runs = runs,
  mis_summary = mis_summary,
  drop_tbl = drop_tbl,
  drop_summary = drop_summary
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")

png(draft_fig_path(meta, "mis_vs_theory.png"), width = 1200, height = 900)
plot_mis_with_theory(mis_summary)
dev.off()

png(draft_fig_path(meta, "drop_hist_r.png"), width = 1200, height = 900)
plot_drop_hist_r(drop_tbl)
dev.off()

png(draft_fig_path(meta, "drop_hist_theory.png"), width = 1200, height = 900)
plot_drop_hist_theory(drop_tbl)
dev.off()

write.csv(drop_tbl, file = result_path(meta, "drop_table.csv"), row.names = FALSE)
write.csv(drop_summary, file = result_path(meta, "drop_summary.csv"), row.names = FALSE)