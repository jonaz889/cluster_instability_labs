
# Alignment with dist mat
# Aligns labels such that the distance between
# old and new centers are minimized
align_k_to_k1_leg <- function(km_k, km_k1) {
  # Get amount of centers
  k <- nrow(km_k$centers)
  
  # Initialize a k x k+1 distance matrix
  dist_mat <- matrix(NA_real_, nrow = k, ncol = k + 1)
  
  # Populate the distance matrix
  for (i in 1:k) {
    for (j in 1:(k + 1)) {
      dist_mat[i, j] <- dist_euclid(km_k$centers[i, ], km_k1$centers[j, ])
    }
  }
  
  # Get all possible permutations (size k) of 1:k
  perms <- gtools::permutations(k, k)
  
  # To store the minimal sum of distances, each time excluding cluster
  # i as the new cluster
  best_sum_distances <- numeric(k + 1)
  # Matrix where each row stores the best permutation excluding cluster i as
  # the new cluster
  best_perm_for_i <- matrix(NA_integer_, nrow = k + 1, ncol = k)
  
  # Loop, each excludes cluster i
  for (i in 1:(k + 1)) {
    # Get the distance matrix excluding cluster i
    dist_mat_wo_i <- as.matrix(dist_mat[, -i])
    
    # temp for computing best distance sum and storing best permutation
    best_sum <- Inf
    best_perm <- rep(NA_integer_, k)
    
    # Loop through each permutation, compute the sum of distances and 
    # check if it is the smallest observed this far
    for (perm in 1:nrow(perms)) {
      total <- 0
      for (j in 1:k) {
        total <- total + dist_mat_wo_i[j, perms[perm, j]]
      }
      
      if (total < best_sum) {
        best_sum <- total
        best_perm <- perms[perm, ]
      }
    }
    
    # Safe the permutation and minimal distance, when excluding new cluster i
    best_sum_distances[i] <- best_sum
    best_perm_for_i[i, ] <- best_perm
  }
  
  # Store the new cluster which when removed yields the minimal sum distance
  # as well as the permutation
  best_removed_i <- which.min(best_sum_distances)
  perm <- best_perm_for_i[best_removed_i, ]
  
  # Make a "map" which pairs the old clusters to the new clusters
  pair_old_to_new <- numeric(k)
  for (i in 1:k) {
    pair_old_to_new[i] <- (1:(k + 1))[-best_removed_i][perm[i]]
  }
  
  # Make a mapping which pairs new clusters to old clusters
  new_to_old_map <- numeric(k + 1) + (k + 1)
  for (i in 1:k) {
    new_to_old_map[pair_old_to_new[i]] <- i
  }
  
  # Return the new clustering labels aligned
  new_to_old_map[km_k1$cluster]
}


# Alignment with dist mat
# Aligns labels such that the distance between
# old and new centers are minimized
# This one uses the hungarian algorithm to align the labels as an
# assignment problem
align_k_to_k1 <- function(km_k, km_k1) {
  # Get amount of centers
  k <- nrow(km_k$centers)
  
  # Initialize a k x k+1 distance matrix
  dist_mat <- matrix(NA_real_, nrow = k, ncol = k + 1)
  
  # Populate the distance matrix
  for (i in 1:k) {
    for (j in 1:(k + 1)) {
      dist_mat[i, j] <- dist_euclid(km_k$centers[i, ], km_k1$centers[j, ])
    }
  }
  
  # assignment problem
  sol <- clue::solve_LSAP(dist_mat)
  
  # assignment[r] = column assigned to row r
  assigned_cols <- as.integer(sol)
  
  
  # Make a mapping which pairs new clusters to old clusters
  new_to_old_map <- numeric(k + 1) + (k + 1)

  # Real rows 1:k correspond to old clusters
  for (i in 1:k) {
    new_to_old_map[assigned_cols[i]] <- i
  }
  
  # Return the new clustering labels aligned
  new_to_old_map[km_k1$cluster]
}


# Function to compute the rim rate
rim_rate_between_k_and_k1 <- function(km_k, km_k1) {
  k <- nrow(km_k$centers)
  km_k1_aligned <- align_k_to_k1(km_k, km_k1)
  
  # Each rim point is set to 1 if rim and 0 if not rim
  # Points are counted as rimdata iff they do not stay in their own cluster
  # and do not move into the new cluster. I.e. if it moves between old clusters
  is_rim <- (km_k1_aligned != km_k$cluster) & (km_k1_aligned != k + 1)
  
  # The the proportion of rimdata
  mean(is_rim)
}
