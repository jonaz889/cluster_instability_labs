#R/e2_plots.R


plot_mis_dist <- function(run_sigma_list) {
  sigmas <- sapply(run_sigma_list, function(run) run$sigma)
  cols <- seq_along(run_sigma_list)
  
  xlim <- range(unlist(lapply(run_sigma_list, function(run) run$dist)))
  ylim <- c(0, max(unlist(lapply(run_sigma_list, function(run) run$mis_rate))))
  
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = expression(d(mu[1], mu[2])),
       ylab = "Misclassification rate",
       main = "Misclassification vs distance",
       cex.main = 1.5, cex.lab = 1.3,
       )
  
  for (i in seq_along(run_sigma_list)) {
    lines(run_sigma_list[[i]]$dist, run_sigma_list[[i]]$mis_rate, col = cols[i], lwd = 2)
  }
  
  legend("topright",
         legend = paste0("σ = ", sigmas),
         col = cols,
         lwd = 2,
         cex = 0.9,
         bty = "n")         
  }

plot_mis_dist_over_variance <- function(run_sigma_list, xlim=NULL) {
  sigmas <- sapply(run_sigma_list, function(run) run$sigma)
  cols <- seq_along(run_sigma_list)
  
  if(is.null(xlim)){xlim <- range(unlist(lapply(run_sigma_list, function(run) run$dist / run$sigma)))}
  ylim <- c(0, max(unlist(lapply(run_sigma_list, function(run) run$mis_rate))))
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = expression(d(mu[1], mu[2]) / sigma),
       ylab = "Misclassification rate",
       main = "Misclassification vs distance/variance",
       cex.main = 1.5, cex.lab = 1.3,
  )
  

  for (i in seq_along(run_sigma_list)) {
    lines(run_sigma_list[[i]]$dist / run_sigma_list[[i]]$sigma,
          run_sigma_list[[i]]$mis_rate, col = cols[i], lwd = 2)
  }
  
  
  legend("topright",
         legend = paste0("σ = ", sigmas),
         col = cols,
         lwd = 2,
         cex = 0.9,
         bty = "n")         
}



plot_single_sigma_ci <- function(sigma_run, xlim=c(0,7)) {
  
  plot(NULL,
       xlim = xlim,
       ylim = c(0,0.5),
       xlab = expression(d(mu[1], mu[2]) / sigma),
       ylab = "Misclassification rate",
       main = paste0("CI bands (σ = ", sigma_run$sigma[1], ")"),
       cex.lab = 1.3, cex.main = 1.5)
  
  # lower CI
  lines(sigma_run$dist_over_sigma,
        sigma_run$lower,
        lty = 1, col="red")
  
  # upper CI
  lines(sigma_run$dist_over_sigma,
        sigma_run$upper,
        lty = 1, col="blue")
}

