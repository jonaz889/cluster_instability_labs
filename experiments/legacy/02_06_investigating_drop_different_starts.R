rm(list = ls())

source("R/load_all.R")

params <- list(
  sigmas = c(1,2,4,8),
  n = 500,
  dist_min = 0,
  dist_max = 8,
  n_dist = 150,
  reps = 200,
  kmeans_nstart = 20
)

seed <- 1

meta <- start_run("exp2C_drop_variability", params = params, seed = seed)

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

drop_tbl <- extract_drop_table(runs)
drop_summary <- summarise_drop_table(drop_tbl)

res <- list(
  runs = runs,
  drop_tbl = drop_tbl,
  drop_summary = drop_summary
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")

png(draft_fig_path(meta, "drop_distribution.png"), width = 900, height = 700)
e1_plot_drop_distribution(drop_tbl)
dev.off()

png(draft_fig_path(meta, "mean_drop_locations.png"), width = 900, height = 700)
par(mfrow = c(1,2))
e1_plot_mean_rstar(drop_summary)
e1_plot_mean_diststar(drop_summary)
dev.off()