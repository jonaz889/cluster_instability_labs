rm(list = ls())

source("R/load_all.R")
library("bench")

params <- list(
  n = 1000,
  ds = c(2,4,8,16,32,64),
  ks = 1:10,
  reps = 100,
  kmeans_nstart = 25
)

seed <- 1
set.seed(1)

meta <- start_run(
  "exp_04_02_rim_resampling_many_d",
  params = params,
  seed = seed
)
res_bench <- bench::mark(
  run_rim_reps_over_d(
    ds = params$ds,
    n = params$n,
    ks = params$ks,
    nstart = params$kmeans_nstart,
    reps = params$reps,
    seed = seed
  ),
  iterations = 5,
  check = FALSE
)

res_bench

runs_by_d <- run_rim_reps_over_d(
  ds = params$ds,
  n = params$n,
  ks = params$ks,
  nstart = params$kmeans_nstart,
  reps = params$reps,
  seed = seed
)

# one summary per dimension
summary_list <- lapply(runs_by_d, rim_reps_sum)

# attach dimension column and bind together
summary_by_d <- do.call(
  rbind,
  lapply(seq_along(summary_list), function(i) {
    tbl <- summary_list[[i]]
    tbl$d <- params$ds[i]
    tbl
  })
)

# CI across dimension-specific means
summary_across_d <- across_dimension_means_sum(summary_by_d)

# save combined table
write.csv(
  summary_by_d,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_summary_by_d.csv")
  ),
  row.names = FALSE
)

# save across-dimension summary table
write.csv(
  summary_across_d,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_summary_across_d.csv")
  ),
  row.names = FALSE
)

# plot all mean curves in one figure
pdf(draft_fig_path(meta, "rim_mean_curves_by_d.pdf"), width = 8, height = 6)
plot_rim_curves_by_d(summary_by_d)
dev.off()

# plot CI across dimension-specific means
pdf(draft_fig_path(meta, "rim_ci_across_d.pdf"), width = 8, height = 6)
plot_overall_dimension_mean_ci(summary_across_d)
dev.off()

res <- list(
  runs_by_d = runs_by_d,
  summary_by_d = summary_by_d,
  summary_across_d = summary_across_d,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
