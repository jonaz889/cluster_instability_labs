rm(list = ls())

source("R/load_all.R")

params <- list(
  sigmas = c(1/2, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5),
  n =  1000,
  dist_min = 0,
  dist_max = 35,
  n_dist = 200,#4000,
  reps = 1,
  kmeans_nstart = 25
)

# Seed 2 er drop seed
# Seed 1 er clean
seed <- 2

meta <- start_run("exp_02_01", params = params, seed = seed)

runs <- run_kmeans_separation_sigmas(
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

res <- list(
  runs = runs,
  params = params
)


# misclass. vs dist and misclass, vs dist/sigma
pdf(
  make_fig_path(meta, "02_01_exploratory_misclas_vs_dist.pdf"),
  width = 10,
  height = 4.5
)

par(mfrow = c(1, 2))
plot_mis_dist(runs[[1]])
plot_mis_dist_over_sigma(runs[[1]], xlim = c(0, 7))

# Where misclassification first hits 0 for each sigma
zero_mis_points <- sapply(runs[[1]], function(run) {
  i <- which(run$mis_rate == 0)[1]
  if (is.na(i)) NA_real_ else run$dist[i] / run$sigma
})

mean_zero_mis_points <- mean(zero_mis_points, na.rm = TRUE)
sd_zero_mis_points <- sd(zero_mis_points, na.rm = TRUE)
se_zero_mis_points <- sd_zero_mis_points / sqrt(sum(!is.na(zero_mis_points)))

abline(v = mean_zero_mis_points, lty = 3)

dev.off()

# summary table for the zero-mis points
zero_tbl <- data.frame(
  sigma = params$sigmas,
  first_zero_dist_over_sigma = zero_mis_points
)

write.csv(
  zero_tbl,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_zero_points.csv")
  ),
  row.names = FALSE
)


pdf(
  make_fig_path(meta, "02_01_exploratory_misclas_vs_dist.pdf"),
  width = 10,
  height = 4.5
)

par(mfrow = c(1, 2))

plot(NULL,
     xlim = c(params$dist_min, params$dist_max),
     ylim = c(0, 0.5),
     xlab = expression(d(mu[1], mu[2])),
     ylab = "Misclassification rate",
     main = "Misclassification vs distance")

for (run in runs[[1]]) {
  lines(run$dist, run$mis_rate, lwd = 2, col=1)
  lines(run$dist,
        theoretic_misclass(run$dist, run$sigma),
        lwd = 2, col=2)
}

legend("topright",
       legend = c("Empirical", "Theoretical"),
       col = c(1,2),
       lwd = 2,
       bty = "n")
rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
