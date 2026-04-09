rm(list = ls())
source("R/load_all.R")
params <- list(
  n = 500,
  ds = c(2, 4, 8),# 16, 32, 64),
  ks = 1:10,
  reps = 100,
  kmeans_nstart = 10
)

seed <- 1

meta <- start_run(
  "exp_04_03_rim_same_data_nested_d",
  params = params,
  seed = seed
)

runs_nested <- run_rim_reps_over_d_nested_scaled(
  n = params$n,
  ds = params$ds,
  ks = params$ks,
  nstart = params$kmeans_nstart,
  reps = params$reps,
  seed = seed
)

summary_by_d <- rim_reps_over_d_nested_scaled_sum(runs_nested)

write.csv(
  summary_by_d,
  file.path("outputs", "tables",
            paste0(meta$experiment, "_", meta$run_id, "_summary_by_d.csv")),
  row.names = FALSE
)

pdf(make_fig_path(meta, "rim_mean_curves_by_d.pdf"), width = 8, height = 6)
par(mfrow=c(1,1))
plot_rim_curves_by_d(summary_by_d)
dev.off()

res <- list(
  runs_nested = runs_nested,
  summary_by_d = summary_by_d,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")