rm(list = ls())
source("R/load_all.R")

params <- list(
  n = 500,#1000,
  ds = c(2, 4,32),#, 16, 32, 64),
  ks = 1:10,
  reps = 100,
  kmeans_nstart = 10
  )

seed <- 1

meta <- start_run(
  "exp_04_02_rim_many_d",
  params = params,
  seed = seed
)


runs_by_d_s1 <- run_rim_reps_over_d_scaled(
    ds = params$ds,
    n = params$n,
    ks = params$ks,
    nstart = params$kmeans_nstart,
    reps = params$reps,
    seed = seed,
    scale_func=std_data_1
  )

runs_by_d_s2 <- run_rim_reps_over_d_scaled(
    ds = params$ds,
    n = params$n,
    ks = params$ks,
    nstart = params$kmeans_nstart,
    reps = params$reps,
    seed = seed,
    scale_func=std_data_2
  )

runs_by_d_s3 <- run_rim_reps_over_d_scaled(
  ds = params$ds,
  n = params$n,
  ks = params$ks,
  nstart = params$kmeans_nstart,
  reps = params$reps,
  seed = seed,
  scale_func=std_data_3
)

summary_by_d_s1 <- rim_reps_over_d_sum(runs_by_d_s1)
summary_across_d_s1 <- across_dimension_means_sum(summary_by_d_s1)
summary_by_d_s2 <- rim_reps_over_d_sum(runs_by_d_s2)
summary_across_d_s2 <- across_dimension_means_sum(summary_by_d_s2)
summary_by_d_s3 <- rim_reps_over_d_sum(runs_by_d_s3)
summary_across_d_s3 <- across_dimension_means_sum(summary_by_d_s3)

pdf(make_fig_path(meta, "mean_ci_scaled.pdf"), width = 8, height = 6)
par(mfrow=c(1,1))
plot_rim_ci(summary_by_d_s1[summary_by_d_s1$d==2,], ylim=c(0,0.5))
plot_rim_ci_add_line(summary_by_d_s2[summary_by_d_s2$d==2,],  col=2)
plot_rim_ci_add_line(summary_by_d_s3[summary_by_d_s3$d==2,], col = 3)

plot_rim_ci_add_line(summary_by_d_s1[summary_by_d_s1$d==4,],  col=1)
plot_rim_ci_add_line(summary_by_d_s2[summary_by_d_s2$d==4,],  col=2)
plot_rim_ci_add_line(summary_by_d_s3[summary_by_d_s3$d==4,], col = 3)

plot_rim_ci_add_line(summary_by_d_s1[summary_by_d_s1$d==32,],  col=1)
plot_rim_ci_add_line(summary_by_d_s2[summary_by_d_s2$d==32,],  col=2)
plot_rim_ci_add_line(summary_by_d_s3[summary_by_d_s3$d==32,], col = 3)

legend("topleft", legend=c("id","scale","normalize"),
       col=c(1,2,3), lwd=2, bty="n")
dev.off()

pdf(make_fig_path(meta, "rim_ci_across_d.pdf"), width = 8, height = 6)
par(mfrow=c(1,2))
plot_overall_dimension_mean_ci(summary_across_d_s2)
plot_overall_dimension_mean_ci(summary_across_d_s3)
dev.off()

res <- list(
  runs_by_d_s2 = runs_by_d_s2,
  runs_by_d_s3 = runs_by_d_s3,
  summary_by_d_s2 = summary_by_d_s2,
  summary_by_d_s3 = summary_by_d_s3,
  summary_across_d_s2 = summary_across_d_s2,
  summary_across_d_s3 = summary_across_d_s3,
  params = params
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
