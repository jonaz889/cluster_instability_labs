rm(list = ls())

source("R/load_all.R")

sigmas <- c(1,2,3)

params <- list(
  sigmas = sigmas,
  n = 500,
  dist_min = 0,
  dist_max = 7*max(sigmas),
  n_dist = 500,
  reps = 200,#200,
  kmeans_nstart = 20
)

seed <- 1

meta <- start_run("exp_02_03_data_variance", params = params, seed = seed)

runs <- run_kmeans_separation_sigmas_rand_across_sigmas(
  sigmas = params$sigmas,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  n_dist = params$n_dist,
  kmeans_nstart = params$kmeans_nstart,
  reps = params$reps,
  keep_all = FALSE
)

# Summaries across replications
sigma_runs <- summarise_mis_ci_over_reps(runs)
names(sigma_runs) <- paste0("sigma_", sigmas)


mean_runs_for_plot <- lapply(sigma_runs, function(x) {
  list(
    sigma = x$sigma[1],
    dist = x$dist,
    mis_rate = x$mean
  )
})


## Plot misclassification vs distance
pdf(
  make_fig_path(meta, "misclass_vs_dist_theo.pdf"),
  width = 10,
  height = 6
)
par(mfrow=c(1,1))
plot_mis_dist(mean_runs_for_plot)
x_theory <- seq(0, 20, length.out = 5000)
lines(x_theory, theory_misclassification(x_theory, 1),
      lty = 2, lwd = 3, col = "black")
lines(x_theory, theory_misclassification(x_theory, 2),
      lty = 2, lwd = 3, col = "black")
lines(x_theory, theory_misclassification(x_theory, 3),
      lty = 2, lwd = 3, col = "black")
legend("right",
       legend = c("Theory"),
       col = c("black"),
       lwd = c(3),
       lty = c(2),
       cex = 0.9,
       bty = "n")
dev.off()
## Plot misclassification vs distance/variance
pdf(
  make_fig_path(meta, "misclass_vs_dist_over_sigma_theo.pdf"),
  width = 10,
  height = 6
)

plot_mis_dist_over_sigma(mean_runs_for_plot, xlim=c(0,7))
x_theory <- seq(0, 7, length.out = 500)
lines(x_theory, theory_misclassification(x_theory, 1),
      lty = 2, lwd = 3, col = "black")
legend("right",
       legend = c("Theory"),
       col = c("black"),
       lwd = c(3),
       lty = c(2),
       cex = 0.9,
       bty = "n")

dev.off()


mean_minus_theory_for_plot <- lapply(sigma_runs, function(x) {
  list(
    sigma = x$sigma[1],
    dist = x$dist,
    mis_rate = x$mean - x$theory
  )
})

pdf(
  make_fig_path(meta, "misclass_vs_dist_over_sigma_error_theo.pdf"),
  width = 10,
  height = 6
)

plot_mis_dist_over_sigma(mean_minus_theory_for_plot, xlim = c(0, 7))
abline(h = 0, lty = 2, lwd = 2)
which.max(mean_minus_theory_for_plot$sigma_1$mis_rate)
dev.off()

#  CI plot for each sigma
for (nm in names(sigma_runs)) {
  pdf(
    make_fig_path(meta, paste0(nm, "_ci_curve.pdf")),
    width = 10,
    height = 6
  )
  
  plot_single_sigma_ci(sigma_runs[[nm]])
  
  dev.off()
}

#zoomed CI plot for each sigma
for (nm in names(sigma_runs)) {
  pdf(
    make_fig_path(meta, paste0(nm, "_ci_curve_zoomed.pdf")),
    width = 12,
    height = 4.5
  )
  
  par(mfrow = c(1, 3))
  plot_single_sigma_ci(sigma_runs[[nm]], xlim = c(0, 2))
  plot_single_sigma_ci(sigma_runs[[nm]], xlim = c(2, 4))
  plot_single_sigma_ci(sigma_runs[[nm]], xlim = c(4, 6))
  
  dev.off()
}

#SD all sigmas together
pdf(
  make_fig_path(meta, "sigma_sd_vs_dist_over_sigma.pdf"),
  width = 10,
  height = 6
)

par(mfrow = c(1, 1))
plot(NULL,
     xlim = c(0, 7),
     ylim = c(0, 0.03),
     xlab = expression(d(mu[1], mu[2]) / sigma),
     ylab = "SD of Misclassification rate",
     main = "SD of misclassification rate",
     cex.lab = 1.3, cex.main = 1.5)

for (i in seq_along(sigma_runs)) {
  lines(sigma_runs[[i]]$dist_over_sigma,
        sigma_runs[[i]]$sd,
        col = i,
        lwd = 2)
}

peak_x <- sigma_runs[[1]]$dist_over_sigma[which.max(sigma_runs[[1]]$sd)]

abline(v = peak_x, lty = 2, lwd = 2)

legend("topright",
       legend = c("σ = 1", "σ = 2", "σ = 3",
                  paste0("SD peak = ", round(peak_x, 2))),
       col = c(seq_along(sigma_runs), "black"),
       lty = c(rep(1, length(sigma_runs)), 2),
       lwd = 2,
       bty = "n")

dev.off()

# Summary at below selected points 
x_points <- c(0.1, 0.5, 1, 4, 6.1)

summary_tbl <- do.call(
  rbind,
  lapply(sigma_runs, function(df) summarise_selected_points(df, x_points))
)

write.csv(
  summary_tbl,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_summary_points.csv")
  ),
  row.names = FALSE
)

# Save results and params
res <- list(
  sigma_runs = sigma_runs,
  summary_tbl = summary_tbl,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
