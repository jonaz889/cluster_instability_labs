# R/utils-runlog.R


# Creates a unique run ID based on the exact current time.
run_id_now <- function() { format(Sys.time(), "%Y%m%d-%H%M%S") }

# Return the short Git hash of the current repository.
git_hash <- function() {
  out <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA_character_)

  if (length(out) == 0) NA_character_ else out[1]
}

# Save session information for reproducibility
# I.e. just hard convert all output of sessionInfo() into a list and then
# write it all to path
write_session <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

# Save parameter values in a simple text format
write_params <- function(path, params) {
  # params should be a named list
  txt <- unlist(Map(function(nm, val) paste0(nm, ": ", paste(val, collapse = ",")),
                    names(params), params))
  writeLines(txt, con = path)
}

start_run <- function(experiment_name, params = list(), seed = NA_integer_) {
  rid <- run_id_now()
  gh  <- git_hash()
  
  meta <- list(
    run_id = rid,
    experiment = experiment_name,
    seed = seed,
    git_hash = gh,
    time = as.character(Sys.time()),
    params = params
  )
  
  # Log meta information
  writeLines(
    c(
      paste0("run_id: ", rid),
      paste0("experiment: ", experiment_name),
      paste0("time: ", meta$time),
      paste0("seed: ", ifelse(is.na(seed), "NA", seed)),
      paste0("git_hash: ", ifelse(is.na(gh), "NA", gh))
    ),
    con = file.path("outputs", "logs", paste0(experiment_name, "_", rid, "_meta.txt"))
  )
  
  # Log parameters 
  write_params(file.path("outputs", "logs", paste0(experiment_name, "_", rid, "_params.txt")),c(list(seed = seed), params))
  
  # Log session
  write_session(file.path("outputs", "logs", paste0(experiment_name, "_", rid, "_session.txt")))
  
  meta
}

# Return the path for logging to terminal!!!
# Saves R object so we can access later
# Remember to acess saved RDS objects, e.g: 
# res <- readRDS("outputs/rds/exp_02_02_kmeans_20260309_1425.rds")
save_result <- function(meta, obj) {
  path <- file.path("outputs", "rds", paste0(meta$experiment, "_", meta$run_id, ".rds"))
  saveRDS(obj, path)
  path
}

# Get the path to the output folder
# file.path is system independent!!! just use comma as "go into folder"
make_fig_path <- function(meta, filename) {
  file.path("outputs", "figs", paste0(meta$experiment, "_", meta$run_id, "_", filename))
}
