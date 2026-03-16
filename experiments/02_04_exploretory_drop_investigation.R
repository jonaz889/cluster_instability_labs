rm(list = ls())
source("R/load_all.R")

sigma <- 1
n <- 400
N <- 2*n
seed <- 2

set.seed(seed)
x_r <- rnorm(N)
y_r <- rnorm(N)

params <- list(
  sigma = sigma,
  n = n,
  dist_min = 0,
  dist_max = 5*sigma,
  exper_amount = 200,
  reps = 500,
  kmeans_nstart = 25,
  drop_threshold = 0.05
)

meta <- start_run("exp_02_04_investigating_drops", params=params, seed=seed)

drop_tbl <- screen_drop_runs(
  sigma = params$sigma,
  reps = params$reps,
  base_seed = seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  exper_amount = params$exper_amount,
  kmeans_nstart = params$kmeans_nstart,
  drop_threshold = 0.05
)

write.csv(
  drop_tbl,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_drop_table.csv")
  ),
  row.names = FALSE
)

# Freq drop vs no drop
n_total <- nrow(drop_tbl)
n_drop <- sum(drop_tbl$has_drop)
n_no_drop <- n_total - n_drop
prop_drop <- n_drop / n_total

cat("Total runs:", n_total, "\n")
cat("Runs with drop:", n_drop, "\n")
cat("Runs without drop:", n_no_drop, "\n")
cat("Proportion with drop:", round(prop_drop, 4), "\n")

freq_tbl <- data.frame(
  category = c("Drop", "No drop"),
  count = c(n_drop, n_no_drop),
  proportion = c(prop_drop, 1 - prop_drop)
)

write.csv(
  freq_tbl,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_drop_frequency.csv")
  ),
  row.names = FALSE
)

pdf(
  draft_fig_path(meta, "drop_vs_no_drop_barplot.pdf"),
  width = 6,
  height = 5
)

barplot(
  height = c(n_drop, n_no_drop),
  names.arg = c("Drop", "No drop"),
  ylab = "Number of runs",
  main = "Frequency of runs with and without drops"
)

dev.off()

# Where is the drops?
drop_locs <- drop_tbl$delta_after[drop_tbl$has_drop]

drop_summary <- data.frame(
  n_drop = length(drop_locs),
  mean_drop_location = mean(drop_locs, na.rm = TRUE),
  sd_drop_location = sd(drop_locs, na.rm = TRUE),
  median_drop_location = median(drop_locs, na.rm = TRUE),
  q25 = quantile(drop_locs, 0.25, na.rm = TRUE),
  q75 = quantile(drop_locs, 0.75, na.rm = TRUE)
)

print(drop_summary)

write.csv(
  drop_summary,
  file.path(
    "outputs", "tables",
    paste0(meta$experiment, "_", meta$run_id, "_drop_location_summary.csv")
  ),
  row.names = FALSE
)

pdf(
  draft_fig_path(meta, "drop_location_histogram.pdf"),
  width = 7,
  height = 5
)
  
mean_loc <- mean(drop_locs, na.rm = TRUE)

hist(
  drop_locs,
  breaks = 30,
  main = expression("Distribution of first drop locations"),
  xlab = expression(d(mu[1], mu[2]) / sigma),
  ylab = "Count",
  cex.main = 1.5, cex.lab = 1.3
)

abline(v = mean_loc, lty = 2, lwd = 2)

legend(
  "topright",
  legend = paste0("Mean = ", round(mean_loc, 3)),
  lty = 2,
  lwd = 2,
  bty = "n"
)
dev.off()


# Plot representative drop
drop_rows <- which(drop_tbl$has_drop)

if (length(drop_rows) == 0) {
  stop("No runs with drops were found.")
}

target <- median(drop_tbl$delta_after[drop_rows], na.rm = TRUE)
chosen_row <- drop_rows[which.min(abs(drop_tbl$delta_after[drop_rows] - target))]
chosen_seed <- drop_tbl$rep_seed[chosen_row]

cat("Chosen seed for representative drop-run:", chosen_seed, "\n")


drop_rows <- which(drop_tbl$has_drop)

if (length(drop_rows) <= 1) {
  stop("Not enough drop runs to skip one.")
}

target <- median(drop_tbl$delta_after[drop_rows], na.rm = TRUE)

# compute distance from median
dists <- abs(drop_tbl$delta_after[drop_rows] - target)

# remove the closest one (the bad one)
drop_rows <- drop_rows[order(dists)][-1]

# choose the next best candidate
chosen_row <- drop_rows[1]
chosen_seed <- drop_tbl$rep_seed[chosen_row]

cat("Chosen seed for representative drop-run:", chosen_seed, "\n")

run_full <- rerun_drop_rep(
  sigma = params$sigma,
  rep_seed = chosen_seed,
  n = params$n,
  dist_min = params$dist_min,
  dist_max = params$dist_max,
  exper_amount = params$exper_amount,
  kmeans_nstart = params$kmeans_nstart
)

drop <- find_first_drop(run_full, drop_threshold = params$drop_threshold)

if (is.null(drop)) {
  stop("Chosen run unexpectedly has no drop.")
}

i_from <- drop$drop_from_index
i_to   <- drop$drop_to_index

cat("Drop occurs from index", i_from, "to", i_to, "\n")
cat("Delta before:", drop$delta_before, "\n")
cat("Delta after :", drop$delta_after, "\n")
cat("Drop size   :", drop$drop_size, "\n")


X_before <- run_full$Xs[[i_from]]
X_after  <- run_full$Xs[[i_to]]

km_before <- run_full$kms[[i_from]]
km_after  <- run_full$kms[[i_to]]

pdf(
  draft_fig_path(meta, "before_after_drop.pdf"),
  width = 10,
  height = 5
)

    par(mfrow = c(1, 2))
    
    plot(
      X_before,
      col = km_before$cluster,
      pch = 16,
      main = paste("Before drop"),
      xlab = "",
      ylab = "",
      cex.main = 1.5, cex.lab=1.3,
      axes=FALSE
    )
    points(km_before$centers, pch = 4, cex = 1.5, lwd = 2,col=3)
    box()
    
    plot(
      X_after,
      col = km_after$cluster,
      pch = 16,
      main = "After drop",
      xlab = "",
      ylab = "",
      cex.main = 1.5, cex.lab=1.3,
      axes=FALSE
    )
    points(km_after$centers, pch = 4, cex = 1.5, lwd = 2,col=3)
    box()

dev.off()

pdf(
  draft_fig_path(meta, "representative_drop_curve.pdf"),
  width = 7,
  height = 5
)

par(mfrow=c(1,1))
plot(
  run_full$dist,
  run_full$mis_rate,
  type = "l",
  lwd = 2,
  xlab = expression(d(mu[1], mu[2]) / sigma),
  ylab = "Misclassification rate",
  main = "Misclassification vs distance/variance",
  cex.main = 1.5, cex.lab = 1.3,
)

abline(v = run_full$dist[i_from], lty = 2)
abline(v = run_full$dist[i_to], lty = 3)

dev.off()


res <- list(
  params = params,
  drop_tbl = drop_tbl,
  drop_summary = drop_summary,
  chosen_seed = chosen_seed,
  chosen_drop = drop
)

rds_path <- save_result(meta, res)
cat("Saved RDS:", rds_path, "\n")
