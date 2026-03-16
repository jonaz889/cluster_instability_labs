rm(list = ls())

source("R/load_all.R")

params <- list(
  n = 3000,
  d = 3,
  ks = 1:10,
  reps = 50,
  kmeans_nstart = 50
)

seed <- 1

meta <- start_run("exp_04_01_rim_resampling_fixed_d", params = params, seed = seed)

runs <- run_rim_reps(
  n = params$n,
  d = params$d,
  ks = params$ks,
  nstart = params$kmeans_nstart,
  reps = params$reps,
  seed = seed
)

summary_tbl <- summarise_rim_reps(runs)

write.csv(
  summary_tbl,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_summary.csv")
  ),
  row.names = FALSE
)

# all curves
pdf(draft_fig_path(meta, "all_rim_curves.pdf"), width = 8, height = 6)
plot_all_rim_curves(runs)
dev.off()

# mean && CI
pdf(draft_fig_path(meta, "rim_ci.pdf"), width = 8, height = 6)
plot_rim_ci(summary_tbl, d = params$d)
dev.off()

res <- list(
  runs = runs,
  summary_tbl = summary_tbl,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
