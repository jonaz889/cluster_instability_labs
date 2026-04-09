rm(list = ls())

source("R/load_all.R")

sigmas <- c(1, 2, 3)

n <- 1000
N <- n * 2

seed <- 1

set.seed(seed)
x_r <- rnorm(N)
y_r <- rnorm(N)

params <- list(
  sigmas = sigmas,
  n = n,
  dist_min = 0,
  dist_max = 7 * max(sigmas),
  n_dist = 500,
  reps = 500,
  kmeans_nstart = 1
)

meta <- start_run("exp_02_02_kmeans_variance", params = params, seed = seed)

sigma_runs <- lapply(sigmas, function(s) {
  
  rep_runs <- run_kmeans_separation_with_reps(
    sigma = s,
    x_r = x_r,
    y_r = y_r,
    reps = params$reps,
    n = params$n,
    dist_min = params$dist_min,
    dist_max = params$dist_max,
    n_dist = params$n_dist,
    kmeans_nstart = params$kmeans_nstart,
    seed = seed + 1 * s
  )
  
  summarise_mis_ci_single_sigma(rep_runs)
})

names(sigma_runs) <- paste0("sigma_", sigmas)

# CI plots for each sigma
pdf(
  make_fig_path(meta, "sigma_ci_curves.pdf"),
  width = 12,
  height = 4.5
)

par(mfrow = c(1, 3))

for (i in seq_along(sigma_runs)) {
  plot_single_sigma_ci(sigma_runs[[i]])
}

dev.off()

# Summary table at points specified below
x_points <- c(0.1, 0.5, 6.19)

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

pdf(make_fig_path(meta, "sd_miss_rate.pdf"), width=7, height=5)

plot(NULL,
     xlim = c(0,7),
     ylim = c(0,0.01),
     xlab = expression(d(mu[1],mu[2]) / sigma),
     ylab = "SD of misclassification rate",
     main = "SD across k-means initializations")

for (i in seq_along(sigma_runs)) {
  lines(sigma_runs[[i]]$dist_over_sigma,
        sigma_runs[[i]]$sd,
        col = i,
        lwd = 2)
}

legend("topright",
       legend = paste0("σ = ", c(1,2,3)),
       col = 1:3,
       lwd = 2,
       bty = "n")

dev.off()

# Save results and params
res <- list(
  sigma_runs = sigma_runs,
  summary_tbl = summary_tbl,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")