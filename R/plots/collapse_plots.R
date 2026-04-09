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

