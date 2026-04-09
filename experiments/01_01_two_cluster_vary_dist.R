# experiments/01_01_two_cluster_vary_dist.R
rm(list = ls())

source("R/load_all.R")

params <- list(
  sigmas = c(1, 2),
  n = 1000,
  dist_min = 0,
  dist_max = 3,
  n_dist = 1000,
  kmeans_nstart = 40,
  iter.max = 50,
  algorithm = "Hartigan-Wong"
)
seed <- 1

meta <- start_run("01_01_two_cluster_vary_dist", params = params, seed = seed)

runs <- run_kmeans_separation_sigmas(
  sigmas = params$sigmas,
  seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  n_dist = params$n_dist,
  kmeans_nstart = params$kmeans_nstart,
  iter.max = params$iter.max,
  algorithm = params$algorithm
)

# Save results
save_result(meta, list(params = params, runs = runs))

# Vary dist visualisation
pdf(make_fig_path(meta, "01_01_vary_dist_vis.pdf"),
    width = 10,
    height = 4.5
)

amount_to_plot <- 5
par(mfrow=c(1, 5))
run = runs[[1]][[1]]
for(i in floor(seq(1,length(run$Xs), length.out=amount_to_plot)))
{
  plot_kmeans(run$Xs[[i]],run$kms[[i]]$cluster, rep(1:nrow(run$kms[[i]]$centers), each = run$n), run$dist[[i]], run$mis_rate[[i]])
}

dev.off()

# Misclass. rate vs distance
pdf(make_fig_path(meta, "01_01_mis_vs_dist.pdf"),
  width = 10,
  height = 4.5
)
par(mfrow=c(1,1))
plot_mis_dist(runs[[1]])
dev.off()

# Collapse by dist/sigma
pdf(make_fig_path(meta, "01_01_mis_vs_dist_over_sigma.pdf"),
    width = 10,
    height = 4.5
)
par(mfrow=c(1,1))
plot_mis_dist_over_sigma(runs[[1]])
dev.off()

cat("Done:", meta$experiment, meta$run_id, "\n")