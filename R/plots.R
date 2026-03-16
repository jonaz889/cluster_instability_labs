# R/plots.R


## Experiment 1
e1_plot_mis_vs_dist <- function(runs, main = "Distance vs misclassification rate") {
  sigmas <- vapply(runs, function(x) x$sigma, numeric(1))
  
  plot(runs[[1]]$dist, runs[[1]]$mis_rate,
       type = "l", lwd = 2,
       xlab = expression(paste("Dist (", Delta, ")")),
       ylab = "Misclassification rate",
       main = main,
       cex.lab = 1.3,
       cex.main = 1.5, 
       col=1)
  
  if (length(runs) > 1) {
    for (i in 2:length(runs)) lines(runs[[i]]$dist, runs[[i]]$mis_rate, lwd = 2, col=i)
  }
  legend("topright",
         legend = paste0("\u03C3 = ", sigmas),
         lty = 1, lwd = 2, bty = "n", col=1:length(sigmas))
}

e1_plot_mis_vs_dist_over_sigma <- function(runs,
                                           xlim = c(0, 6),
                                           ylim = c(0, 0.5),
                                           main = expression(paste("Mis-rate vs ", Delta, "/", sigma))) {
  sigmas <- vapply(runs, function(x) x$sigma, numeric(1))
  
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = expression(Delta / sigma),
       ylab = "Misclassification rate",
       main = main,
       cex.lab = 1.3,
       cex.main = 1.5)
  
  for (i in seq_along(runs)) {
    lines(runs[[i]]$dist / runs[[i]]$sigma, runs[[i]]$mis_rate, lwd = 2, col=i)
  }
  legend("topright",
         legend = paste0("\u03C3 = ", sigmas),
         lty = 1, lwd = 2, bty = "n")
}





##  Experiment 2A

# THIS IS REALLY BAD DONT USE IT POLYGONS AREUGLY
e1_add_ci_band <- function(x, lower, upper, col) {
  polygon(x = c(x, rev(x)),y = c(lower, rev(upper)),col = col,border = NA
  )
}

## BAD
e1_plot_mis_ci_raw <- function(ci_tbls, main = "Misclassification vs raw separation") {
  sigmas <- vapply(ci_tbls, function(tab) unique(tab$sigma), numeric(1))
  n_sigmas <- length(ci_tbls)
  
  cols <- seq_len(n_sigmas)
  band_cols <- vapply(cols, function(i) adjustcolor(i, alpha.f = 0.20), character(1))
  
  xlim <- range(unlist(lapply(ci_tbls, function(tab) tab$dist)))
  ylim <- c(0, max(unlist(lapply(ci_tbls, function(tab) tab$upper)), na.rm = TRUE))
  
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = expression(Delta),
       ylab = "Mean misclassification rate",
       main = main,
       cex.lab = 1.3,
       cex.main = 1.5)
  
  for (i in seq_along(ci_tbls)) {
    tab <- ci_tbls[[i]]
    e1_add_ci_band(tab$dist, tab$lower, tab$upper, band_cols[i])
    lines(tab$dist, tab$mean, col = cols[i], lwd = 2)
  }
  
  legend("topright",
         legend = paste0("\u03C3 = ", sigmas),
         col = cols, lwd = 2, bty = "n")
}

#BAD
e1_plot_mis_ci_normalized <- function(ci_tbls,main = "Misclassification vs normalized separation") {
  sigmas <- vapply(ci_tbls, function(tab) unique(tab$sigma), numeric(1))
  n_sigmas <- length(ci_tbls)
  
  cols <- seq_len(n_sigmas)
  band_cols <- vapply(cols, function(i) adjustcolor(i, alpha.f = 0.20), character(1))
  
  xlim <- range(unlist(lapply(ci_tbls, function(tab) tab$dist_over_sigma)))
  ylim <- c(0, max(unlist(lapply(ci_tbls, function(tab) tab$upper)), na.rm = TRUE))
  
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = expression(Delta / sigma),
       ylab = "Mean misclassification rate",
       main = main,
       cex.lab = 1.3,
       cex.main = 1.5)
  
  for (i in seq_along(ci_tbls)) {
    tab <- ci_tbls[[i]]
    e1_add_ci_band(tab$dist_over_sigma, tab$lower, tab$upper, band_cols[i])
    lines(tab$dist_over_sigma, tab$mean, col = cols[i], lwd = 2)
  }
  
  legend("topright",
         legend = paste0("\u03C3 = ", sigmas),
         col = cols, lwd = 2, bty = "n")
}

#BAD
e1_plot_collapse_residuals <- function(ci_tbls,main = "Deviation from pooled normalized curve") {
  sigmas <- vapply(ci_tbls, function(tab) unique(tab$sigma), numeric(1))
  n_sigmas <- length(ci_tbls)
  cols <- seq_len(n_sigmas)
  
  x_min <- max(vapply(ci_tbls, function(tab) min(tab$dist_over_sigma), numeric(1)))
  x_max <- min(vapply(ci_tbls, function(tab) max(tab$dist_over_sigma), numeric(1)))
  x_grid <- seq(x_min, x_max, length.out = 300)
  
  mean_mat <- sapply(ci_tbls, function(tab) {
    approx(tab$dist_over_sigma, tab$mean, xout = x_grid)$y
  })
  
  pooled <- rowMeans(mean_mat)
  resid_mat <- sweep(mean_mat, 1, pooled, "-")
  
  ylim <- range(resid_mat, na.rm = TRUE)
  
  plot(NULL, xlim = range(x_grid), ylim = ylim,
       xlab = expression(Delta / sigma),
       ylab = "Mean curve - pooled mean curve",
       main = main,
       cex.lab = 1.3,
       cex.main = 1.5)
  
  abline(h = 0, lty = 2)
  
  for (i in seq_len(n_sigmas)) {
    lines(x_grid, resid_mat[, i], col = cols[i], lwd = 2)
  }
  
  legend("topright",
         legend = paste0("\u03C3 = ", sigmas),
         col = cols, lwd = 2, bty = "n")
}


## Experiment 2B

e1_plot_drop_mechanism <- function(res,
                                   sigma_index = 1,
                                   main_left = "Center orientation vs separation",
                                   main_right = "Misclassification with drop location") {
  run_obj <- res$runs[[sigma_index]]
  drop_idx <- res$drop_info$drop_indices[[sigma_index]]
  
  if (length(drop_idx) == 0) {
    stop("No drops found for this sigma.")
  }
  
  d <- drop_idx[1]
  
  par(mfrow = c(1, 2))
  
  plot(
    run_obj$dist,
    apply(res$drop_info$dist_to_center[[sigma_index]], 1, mean),
    type = "l", lwd = 2,
    xlab = expression(Delta),
    ylab = "Mean |y-coordinate| of centers",
    main = main_left,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  abline(v = run_obj$dist[d], lty = 2)
  
  plot(
    run_obj$dist,
    run_obj$mis_rate,
    type = "l", lwd = 2,
    xlab = expression(Delta),
    ylab = "Misclassification rate",
    main = main_right,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  abline(v = run_obj$dist[d], lty = 2)
}


e1_plot_before_after_drop <- function(res,
                                      sigma_index = 1,
                                      main_before = "Before drop",
                                      main_after = "After drop") {
  run_obj <- res$runs[[sigma_index]]
  drop_idx <- res$drop_info$drop_indices[[sigma_index]]
  
  if (length(drop_idx) == 0) {
    stop("No drops found for this sigma.")
  }
  
  d <- drop_idx[1]
  
  if (d <= 1) {
    stop("Drop occurs at first grid point; cannot plot before/after.")
  }
  
  par(mfrow = c(1, 2))
  
  plot(
    run_obj$Xs[[d - 1]],
    col = run_obj$kms[[d - 1]]$cluster,
    pch = 16, cex = 0.5,
    xlab = "", ylab = "", axes = FALSE,
    main = main_before,
    cex.main = 1.5
  )
  points(run_obj$kms[[d - 1]]$centers, pch = 19, cex = 1.5, col = 5)
  box()
  
  plot(
    run_obj$Xs[[d]],
    col = run_obj$kms[[d]]$cluster,
    pch = 16, cex = 0.5,
    xlab = "", ylab = "", axes = FALSE,
    main = main_after,
    cex.main = 1.5
  )
  points(run_obj$kms[[d]]$centers, pch = 19, cex = 1.5, col = 5)
  box()
}


## Experiment 2C

e1_plot_drop_distribution <- function(drop_tbl,
                                      main = "Distribution of normalized drop location") {
  sigmas <- sort(unique(drop_tbl$sigma))
  
  boxplot(
    split(drop_tbl$r_star, drop_tbl$sigma),
    names = sigmas,
    xlab = expression(sigma),
    ylab = expression(r^"*" == Delta^"*" / sigma),
    main = main,
    cex.lab = 1.3,
    cex.main = 1.5
  )
}


e1_plot_mean_rstar <- function(drop_summary,
                               main = "Mean normalized drop location") {
  ylim <- range(drop_summary$lower_r_star, drop_summary$upper_r_star, na.rm = TRUE)
  
  plot(
    drop_summary$sigma,
    drop_summary$mean_r_star,
    ylim = ylim,
    pch = 19,
    type = "b",
    xlab = expression(sigma),
    ylab = expression("Mean " * r^"*"),
    main = main,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  
  arrows(
    x0 = drop_summary$sigma,
    y0 = drop_summary$lower_r_star,
    x1 = drop_summary$sigma,
    y1 = drop_summary$upper_r_star,
    angle = 90, code = 3, length = 0.05
  )
}


e1_plot_mean_diststar <- function(drop_summary,
                                  main = "Mean raw drop location") {
  ylim <- range(drop_summary$lower_dist_star, drop_summary$upper_dist_star, na.rm = TRUE)
  
  plot(
    drop_summary$sigma,
    drop_summary$mean_dist_star,
    ylim = ylim,
    pch = 19,
    type = "b",
    xlab = expression(sigma),
    ylab = expression("Mean " * Delta^"*"),
    main = main,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  
  arrows(
    x0 = drop_summary$sigma,
    y0 = drop_summary$lower_dist_star,
    x1 = drop_summary$sigma,
    y1 = drop_summary$upper_dist_star,
    angle = 90, code = 3, length = 0.05
  )
}