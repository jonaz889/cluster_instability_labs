# Returns some information on drop if found else null
find_first_drop <- function(run, drop_threshold = 0.05) {
  drop_sizes <- -diff(run$mis_rate)
  drop_index <- which(drop_sizes >= drop_threshold)
  
  if (length(drop_index) == 0) {
    return(NULL)
  }
  
  i <- drop_index[1]
  
  list(
    drop_from_index = i,
    drop_to_index = i + 1,
    drop_size= drop_sizes[i],
    dist_before= run$dist[i],
    dist_after= run$dist[i + 1]
  )
}


# Summarizes drops across runs.¨
# Populate DF by adding each row to a list and lastly rowbinding
drops_summary <- function(rep_runs, drop_threshold = 0.05) {
  out <- lapply(seq_along(rep_runs), function(r) {
    
    drop_info <- find_first_drop(rep_runs[[r]], drop_threshold = drop_threshold)
    
    if (is.null(drop_info)) {
      data.frame(
        rep = r,
        has_drop = FALSE,
        drop_from_index = NA_integer_,
        drop_to_index   = NA_integer_,
        drop_size       = NA_real_,
        dist_before     = NA_real_,
        dist_after      = NA_real_
      )
    } else {
      data.frame(
        rep = r,
        has_drop = TRUE,
        drop_from_index = drop_info$drop_from_index,
        drop_to_index   = drop_info$drop_to_index,
        drop_size       = drop_info$drop_size,
        dist_before     = drop_info$dist_before,
        dist_after      = drop_info$dist_after
      )
    }
  })
  
  do.call(rbind, out)
}

