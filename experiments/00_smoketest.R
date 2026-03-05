# experiments/00_smoketest.R
rm(list = ls())

# Always run from project root
# (If you run via Rscript from root or RStudio project, this is fine)
source("R/hello_fun.R")

run_id <- format(Sys.time(), "%Y%m%d-%H%M%S")

res <- hello_kmeans(seed = 42, n = 400)

# --- save result object (ignored by git) ---
saveRDS(res, file = file.path("outputs", "rds", paste0("smoketest_", run_id, ".rds")))

# --- save a log (ignored by git) ---
writeLines(c(
  paste("run_id:", run_id),
  paste("seed:", res$seed),
  paste("n:", res$n),
  "",
  "sessionInfo():",
  paste(capture.output(sessionInfo()), collapse = "\n")
), con = file.path("outputs", "logs", paste0("smoketest_", run_id, "_session.txt")))

# --- save a draft figure (ignored by git) ---
png(file.path("outputs", "figs", "draft", paste0("smoketest_", run_id, ".png")),
    width = 1200, height = 900)
plot(res$X, col = res$km$cluster, pch = 16, cex = 0.7,
     main = paste("Smoke test kmeans | seed =", res$seed))
dev.off()

cat("OK: wrote outputs for run_id =", run_id, "\n")
