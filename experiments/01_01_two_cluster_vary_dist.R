# experiments/01_e1_gauss_distance_sigma_sweep.R
rm(list = ls())

source("R/utils-runlog.R")
source("R/plots.R")
source("R/e1_bin_story.R")

params <- list(
  sigmas = c(0.1, 0.5, 1, 2, 3),
  n = 1000,
  dist_min = 0,
  dist_max = 3,
  exper_amount = 1000,
  kmeans_nstart = 40,
  iter.max = 50,
  algorithm = "Hartigan-Wong"
)
seed <- 1

meta <- start_run("e1_gauss_sigma_sweep", params = params, seed = seed)

runs <- e1_sigma_sweep(
  sigmas = params$sigmas,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  exper_amount = params$exper_amount,
  kmeans_nstart = params$kmeans_nstart,
  iter.max = params$iter.max,
  algorithm = params$algorithm
)

# Save results
save_result(meta, list(params = params, runs = runs))

# Misclass. rate vs distance
png(draft_fig_path(meta, "mis_vs_dist.png"), width = 1600, height = 1000)
e1_plot_mis_vs_dist(runs)
dev.off()

# Collapse by dist/sigma
png(draft_fig_path(meta, "mis_vs_dist_over_sigma.png"), width = 1600, height = 1000)
e1_plot_mis_vs_dist_over_sigma(runs)
dev.off()

cat("Done:", meta$experiment, meta$run_id, "\n")