rm(list = ls())
source("R/load_all.R")

params <- list(
  sigma = 1,
  R = 200,
  n = 300,
  side_length = 20,
  pull_max = 5,
  n_dist = 200,
  kmeans_nstart = 20
)

seed <- 1

meta <- start_run("exp_03_01_triangle_pull_bins", params = params, seed = seed)

out <- bin_clusterings_by_ari(
  sigma = params$sigma,
  R = params$R,
  efunc = function(sigma) e3_run_sigma(
    sigma = sigma,
    seed = sample.int(1e7, 1),
    n = params$n,
    side_length = params$side_length,
    pull_max = params$pull_max,
    n_dist = params$n_dist,
    kmeans_nstart = params$kmeans_nstart
  )
)

summary_tbl <- summarise_bins_over_pull(out)

pdf(draft_fig_path(meta, "n_bins_vs_pull.pdf"), width = 7, height = 5)
plot(summary_tbl$pull, summary_tbl$n_bins, type = "l", lwd = 2,
     xlab = "Pull", ylab = "Number of distinct clusterings",
     main = "Distinct clustering solutions vs pull")
dev.off()

pdf(draft_fig_path(meta, "dominant_fraction_vs_pull.pdf"), width = 7, height = 5)
plot(summary_tbl$pull, summary_tbl$dominant_frac, type = "l", lwd = 2,
     xlab = "Pull", ylab = "Dominant bin fraction",
     main = "Dominance of most frequent clustering vs pull")
dev.off()

save_result(meta, list(out = out, summary_tbl = summary_tbl, params = params))