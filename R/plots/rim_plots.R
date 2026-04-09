
plot_overall_dimension_mean_ci <- function(summary_dim_mean) {
  plot(summary_dim_mean$k, summary_dim_mean$mean,
       type = "l", lwd = 2,
       ylim = range(c(summary_dim_mean$lower, summary_dim_mean$upper)),
       xlab = "k",
       ylab = "Mean rim data frequency",
       main = "Mean rim curve across dimensions",
       cex.main = 1.4, cex.lab = 1.2)
  
  lines(summary_dim_mean$k, summary_dim_mean$lower, lty = 2)
  lines(summary_dim_mean$k, summary_dim_mean$upper, lty = 2)
  points(summary_dim_mean$k, summary_dim_mean$mean, pch = 16)
}

plot_rim_curves_by_d <- function(summary_across_d) {
  ds <- sort(unique(summary_across_d$d))
  ks <- sort(unique(summary_across_d$k))
  
  ylim <- range(c(summary_across_d$lower, summary_across_d$upper))
  
  plot(NULL,
       xlim = range(ks),
       ylim = ylim,
       xlab = "k",
       ylab = "Mean rim data frequency",
       main = "Mean rim curves for different dimensions",
       cex.main = 1.5, cex.lab = 1.3)
  
  for (i in seq_along(ds)) {
    sub <- summary_across_d[summary_across_d$d == ds[i], ]
    lines(sub$k, sub$mean, lwd = 2, col = i)
    points(sub$k, sub$mean, pch = 16, col = i)
  }
  
  legend("topleft",
         legend = paste0("d = ", ds),
         col = seq_along(ds),
         lwd = 2,
         pch = 16,
         bty = "n")
}

plot_all_rim_curves <- function(runs) {
  ks <- runs[[1]]$ks[-length(runs[[1]]$ks)]
  
  plot(NULL,
       xlim = range(ks),
       ylim = c(0, max(sapply(runs, function(run) max(run$rim_rates)))),
       xlab = "k",
       ylab = "Rim data frequency",
       main = paste0("Rim data frequency across different data, d = ", runs[[1]]$d),
       cex.main=1.5, cex.lab=1.3)
  
  for (r in seq_along(runs)) {
    lines(ks, runs[[r]]$rim_rates)
  }
}

plot_rim_ci <- function(rim_sum, col=1, ylim=c(0,0.5)) {
  plot(rim_sum$k, rim_sum$mean,
       type = "b", lwd = 2,
       ylim = ylim,
       xlab = "k",
       ylab = "Rim data frequency",
       main = paste0("Mean and CI of rim data frequency"),
       cex.main = 1.5, cex.lab = 1.3, col=col)
  
  lines(rim_sum$k, rim_sum$lower, lty = 2, col=col)
  lines(rim_sum$k, rim_sum$upper, lty = 2,col=col)
}

plot_rim_ci_add_line <- function(rim_sum, col=1) {
  lines(rim_sum$k, rim_sum$mean, type="b", lwd = 2, col = col)  
  lines(rim_sum$k, rim_sum$lower, lty = 2, col=col)
  lines(rim_sum$k, rim_sum$upper, lty = 2,col=col)
}

plot_rim_sd_vs_k <- function(rim_sum,col=1, ylim=c(0,0.2)) {
  plot(rim_sum$k, rim_sum$sd,type="b",lwd=2,
       ylim = ylim,
       xlab = "k",
       ylab = "SD of rim-data frequency",
       main = "SD of rim-data frequency vs. k",
       cex.main = 1.3,
       cex.lab = 1.2, col=col)
}
plot_rim_sd_vs_k_add_line <- function(rim_sum,col=1) {
  lines(rim_sum$k, rim_sum$sd, type="b",col=col,lwd=2)
}
