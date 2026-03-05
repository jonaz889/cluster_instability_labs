# experiments/01_smoketest_logged.R
rm(list = ls())

source("R/hello_fun.R")
source("R/utils-runlog.R")

params <- list(n = 400, centers = 2, nstart = 10)
seed <- 42

meta <- start_run("smoketest", params = params, seed = seed)

res <- hello_kmeans(seed = seed, n = params$n)

# save result
rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")

# save draft figure
png(draft_fig_path(meta, "clusters.png"), width = 1200, height = 900)
plot(res$X, col = res$km$cluster, pch = 16, cex = 0.7,
     main = paste("Smoke test | seed =", seed, "| git =", meta$git_hash))
dev.off()