# Returns the misclassification rate of the alignment 
# with least mis classification. This is only for kmeans w/ k=2
# true_lab is the true labes, pred_lab is the kmean predicted labels.
mis_rate_aligned_k2 <- function(true_lab, pred_lab) {
  pred_swap_lab <- ifelse(pred_lab == 1, 2, 1)
  min(mean(pred_lab != true_lab), mean(pred_swap_lab != true_lab))
}


# Euclidian dist
dist_euclid <- function(a, b) sqrt(sum((a - b)^2))
