rm(list = ls())
source("R/load_all.R")

params <- list(
  n = 500,#3000,
  d = c(2,4,8),
  ks = 1:10,
  reps = 100,#50,
  kmeans_nstart = 20#50
)

seed <- 1

meta <- start_run(
  "exp_04_01_rim_fixed_d",
  params = params,
  seed = seed
)

runs <- run_rim_reps_over_d(
  n = params$n,
  ds = params$d,
  ks = params$ks,
  nstart = params$kmeans_nstart,
  reps = params$reps,
  seed = seed
)

summary_d1 <- rim_reps_sum(runs[[1]])
summary_d2 <- rim_reps_sum(runs[[2]])
summary_d3 <- rim_reps_sum(runs[[3]])
summary_over_d <- rim_reps_over_d_sum(runs)


pdf(make_fig_path(meta, "all_rim_curves.pdf"), width = 8, height = 6)
par(mfrow=c(1,3))
plot_all_rim_curves(runs[[1]])
plot_all_rim_curves(runs[[2]])
plot_all_rim_curves(runs[[3]])
dev.off()

pdf(make_fig_path(meta, "rim_ci.pdf"), width = 8, height = 6)
par(mfrow=c(1,1))
plot_rim_ci(summary_d1, ylim=c(0,0.5))
plot_rim_ci_add_line(summary_d2,  col=2)
plot_rim_ci_add_line(summary_d3, col = 3)
legend("topleft", legend=c(paste0("d = ", params$d[1]),paste0("d = ", params$d[2]),paste0("d = ", params$d[3])),
       col=c(1,2,3), lwd=2, bty="n")

dev.off()

pdf(make_fig_path(meta, "rim_sd.pdf"), width = 8, height = 6)

par(mfrow=c(1,1))
plot_rim_sd_vs_k(summary_d1, ylim=c(0,0.12))
plot_rim_sd_vs_k_add_line(summary_d2, col=2)
plot_rim_sd_vs_k_add_line(summary_d3, col = 3)
legend("topleft", legend=c(paste0("d = ", params$d[1]),paste0("d = ", params$d[2]),paste0("d = ", params$d[3])),
       col=c(1,2,3), lwd=2, bty="n")
dev.off()
res <- list(
  runs = runs,
  summary_d1 = summary_d1,
  summary_d2 = summary_d2,
  summary_d3 = summary_d3,
  summary_over_d = summary_over_d,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")