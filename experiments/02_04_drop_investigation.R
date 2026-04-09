rm(list = ls())
source("R/load_all.R")

params <- list(
  sigma = 1,
  n = 500,
  dist_min = 0,
  dist_max = 8,
  n_dist = 150,
  reps = 10,#200,
  kmeans_nstart = 20,
  drop_threshold = 0.05,
  conf = 0.95
)

seed <- 1

meta <- start_run(
  "exp_02_04_drop_investigation",
  params = params,
  seed = seed
)

runs <- run_kmeans_separation_sigmas(
  sigmas = params$sigma,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  n_dist = params$n_dist,
  kmeans_nstart = params$kmeans_nstart,
  reps = params$reps,
  keep_all = TRUE
)

#  misclassification summary statistics
mis_summary <- summarise_mis_ci_over_reps(
  runs,
  conf = params$conf
)

# Since sigma is fixed, extract the single run from each rep
runs <- lapply(runs, function(run) run[[1]])


# Info on drops df
drops_info <- drops_summary(runs,drop_threshold = params$drop_threshold)

drops_info_summary <- data.frame(
  n_total = nrow(drops_info),
  n_drop = sum(drops_info$has_drop, na.rm = TRUE),
  n_no_drop = sum(!drops_info$has_drop, na.rm = TRUE),
  prop_drop = mean(drops_info$has_drop, na.rm = TRUE)
)



# Histogram over drop locations
pdf(make_fig_path(meta, "drop_hist_d_over_sigma.pdf"), width = 8, height = 6)
hist(
  drops_info$dist_before[drops_info$has_drop],
  breaks = 20,
  main = "Drop locations",
  xlab = expression(d(mu[1], mu[2]) / sigma)
)
dev.off()

pdf(make_fig_path(meta, "drop_vs_no_drop_barplot.pdf"), width = 6, height = 5)
barplot(
  height = c(drops_info_summary$n_drop, drops_info_summary$n_no_drop),
  names.arg = c("Drop", "No drop"),
  ylab = "Number of runs",
  main = "Runs w/ and w/o drop"
)
dev.off()


# save rds
res <- list(
  params = params,
  runs = runs,
  mis_summary = mis_summary,
  drops_info = drops_info,
  drops_info_summary = drops_info_summary
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")

