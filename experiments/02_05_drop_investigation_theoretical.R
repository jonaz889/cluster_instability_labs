rm(list = ls())
source("R/load_all.R")

params <- list(
  sigma = 1,
  n = 400,
  dist_min = 0,
  dist_max = 5,
  n_dist = 200,
  reps = 10,
  kmeans_nstart = 25,
  drop_threshold = 0.05,
  seed = 7
)

meta <- start_run(
  "exp_drop_switch_analysis",
  params = params,
  seed = params$seed
)

N <- 2 * params$n

set.seed(params$seed)
x_r <- rnorm(N)
y_r <- rnorm(N)

rep_runs <- run_kmeans_separation_with_reps(
  sigma = params$sigma,
  x_r = x_r,
  y_r = y_r,
  reps = params$reps,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  n_dist = params$n_dist,
  kmeans_nstart = params$kmeans_nstart,
  seed = params$seed,
  keep_all = TRUE
)

analysis <- drops_summary(
  rep_runs = rep_runs,
  drop_threshold = params$drop_threshold)

first_drop_rep_i <- which(analysis$has_drop)[1]
first_drop_rep <- rep_runs[[first_drop_rep_i]]
analysis_first_drop_rep <- analysis[first_drop_rep_i,]

## Plot mis vs dist
pdf(make_fig_path(meta, "mis_vs_dist_first_drop_rep"), width = 10, height = 5)

plot_mis_dist_over_sigma(list(first_drop_rep))
abline(v=(first_drop_rep$dist/first_drop_rep$sigma)[analysis_first_drop_rep$drop_from_index], lty=3)
abline(v=(first_drop_rep$dist/first_drop_rep$sigma)[analysis_first_drop_rep$drop_to_index], lty=3)

dev.off()

## Plot drop
first_drop_ana <- drops_summary(list(first_drop_rep), params$drop_threshold)

pdf(make_fig_path(meta, "drop_switch_example.pdf"), width = 10, height = 5)

par(mfrow = c(1, 2))

true_lab <- rep(1:nrow(first_drop_rep$kms[[1]]$centers), each = first_drop_rep$n)

plot_kmeans(
  X = first_drop_rep$Xs[[first_drop_ana$drop_from_index]],
  pred_lab = first_drop_rep$kms[[first_drop_ana$drop_from_index]]$cluster,
  true_lab = true_lab,
  dist = first_drop_rep$dist[first_drop_ana$drop_from_index],
  mis_rate = first_drop_rep$mis_rate[first_drop_ana$drop_from_index]
)
plot_kmeans(
  X = first_drop_rep$Xs[[first_drop_ana$drop_to_index]],
  pred_lab = first_drop_rep$kms[[first_drop_ana$drop_to_index]]$cluster,
  true_lab = true_lab,
  dist = first_drop_rep$dist[first_drop_ana$drop_to_index],
  mis_rate = first_drop_rep$mis_rate[first_drop_ana$drop_to_index]
)

dev.off()


res <- list(
  rep_runs = rep_runs,
  analysis = analysis,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")