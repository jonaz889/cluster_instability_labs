
e1_plot_drop_mechanism <- function(res,
                                   sigma_index = 1,
                                   main_left = "Center orientation vs separation",
                                   main_right = "Misclassification with drop location") {
  run_obj <- res$runs[[sigma_index]]
  drop_idx <- res$drop_info$drop_indices[[sigma_index]]
  
  if (length(drop_idx) == 0) {
    stop("No drops found for this sigma.")
  }
  
  d <- drop_idx[1]
  
  par(mfrow = c(1, 2))
  
  plot(
    run_obj$dist,
    apply(res$drop_info$dist_to_center[[sigma_index]], 1, mean),
    type = "l", lwd = 2,
    xlab = expression(Delta),
    ylab = "Mean |y-coordinate| of centers",
    main = main_left,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  abline(v = run_obj$dist[d], lty = 2)
  
  plot(
    run_obj$dist,
    run_obj$mis_rate,
    type = "l", lwd = 2,
    xlab = expression(Delta),
    ylab = "Misclassification rate",
    main = main_right,
    cex.lab = 1.3,
    cex.main = 1.5
  )
  abline(v = run_obj$dist[d], lty = 2)
}


e1_plot_before_after_drop <- function(res,
                                      sigma_index = 1,
                                      main_before = "Before drop",
                                      main_after = "After drop") {
  run_obj <- res$runs[[sigma_index]]
  drop_idx <- res$drop_info$drop_indices[[sigma_index]]
  
  if (length(drop_idx) == 0) {
    stop("No drops found for this sigma.")
  }
  
  d <- drop_idx[1]
  
  if (d <= 1) {
    stop("Drop occurs at first grid point; cannot plot before/after.")
  }
  
  par(mfrow = c(1, 2))
  
  plot(
    run_obj$Xs[[d - 1]],
    col = run_obj$kms[[d - 1]]$cluster,
    pch = 16, cex = 0.5,
    xlab = "", ylab = "", axes = FALSE,
    main = main_before,
    cex.main = 1.5
  )
  points(run_obj$kms[[d - 1]]$centers, pch = 19, cex = 1.5, col = 5)
  box()
  
  plot(
    run_obj$Xs[[d]],
    col = run_obj$kms[[d]]$cluster,
    pch = 16, cex = 0.5,
    xlab = "", ylab = "", axes = FALSE,
    main = main_after,
    cex.main = 1.5
  )
  points(run_obj$kms[[d]]$centers, pch = 19, cex = 1.5, col = 5)
  box()
}


## Experiment 2C

e1_plot_drop_distribution <- function(drop_tbl,
                                      main = "Distribution of normalized drop location") {
  sigmas <- sort(unique(drop_tbl$sigma))
  
  boxplot(
    split(drop_tbl$r_star, drop_tbl$sigma),
    names = sigmas,
    xlab = expression(sigma),
    ylab = expression(r^"*" == Delta^"*" / sigma),
    main = main,
    cex.lab = 1.3,
    cex.main = 1.5
  )
}