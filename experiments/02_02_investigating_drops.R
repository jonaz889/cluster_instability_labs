rm(list = ls())

source("R/utils-runlog.R")
source("R/plots.R")
source("R/e1_bin_story.R")

params <- list(
  sigma = 5,
  n = 1000,
  dist_min = 0,
  dist_max = 8,
  exper_amount = 200,
  kmeans_nstart = 40
)

seed <- 1

meta <- start_run("exp2B_drop_mechanism", params = params, seed = seed)

res <- bin_story(
  sigmas = params$sigma,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  exper_amount = params$exper_amount,
  kmeans_nstart = params$kmeans_nstart
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")

png(draft_fig_path(meta, "drop_mechanism.png"), width = 1400, height = 600)
e1_plot_drop_mechanism(res)
dev.off()

png(draft_fig_path(meta, "before_after_drop.png"), width = 1200, height = 600)
e1_plot_before_after_drop(res)
dev.off()