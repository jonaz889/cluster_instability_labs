# R/utils-runlog.R


# Creates a unique run ID based on the exact current time.
run_id_now <- function() {
  format(Sys.time(), "%Y%m%d-%H%M%S")
}

# Return the short Git hash of the current repository.
git_hash <- function() {
  out <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), 
                  error = function(e) NA_character_)
  if (length(out) == 0) NA_character_ else out[1]
}

# Create directories if they do not already exist
ensure_dirs <- function(paths) {
  for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

# Save session information for reproducibility
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

start_run <- function(experiment_name, params = list(), seed = NA_integer_,
                      outputs_dir = "outputs") {
  
  ensure_dirs(c(
    file.path(outputs_dir, "logs"),
    file.path(outputs_dir, "rds"),
    file.path(outputs_dir, "figs", "draft"),
    file.path(outputs_dir, "tables")
  ))
  
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
  
  # meta log
  writeLines(
    c(
      paste0("run_id: ", rid),
      paste0("experiment: ", experiment_name),
      paste0("time: ", meta$time),
      paste0("seed: ", ifelse(is.na(seed), "NA", seed)),
      paste0("git_hash: ", ifelse(is.na(gh), "NA", gh))
    ),
    con = file.path(outputs_dir, "logs", paste0(experiment_name, "_", rid, "_meta.txt"))
  )
  
  # params log
  write_params(
    file.path(outputs_dir, "logs", paste0(experiment_name, "_", rid, "_params.txt")),
    c(list(seed = seed), params)
  )
  
  # session log
  write_session(
    file.path(outputs_dir, "logs", paste0(experiment_name, "_", rid, "_session.txt"))
  )
  
  meta
}

save_result <- function(meta, obj, outputs_dir = "outputs") {
  path <- file.path(outputs_dir, "rds",
                    paste0(meta$experiment, "_", meta$run_id, ".rds"))
  saveRDS(obj, path)
  path
}

draft_fig_path <- function(meta, filename, outputs_dir = "outputs") {
  # filename like "mis.png"
  file.path(outputs_dir, "figs", "draft",
            paste0(meta$experiment, "_", meta$run_id, "_", filename))
}

paper_fig_path <- function(filename, outputs_dir = "outputs") {
  # stable filename, tracked by git if you want
  ensure_dirs(file.path(outputs_dir, "figs", "paper"))
  file.path(outputs_dir, "figs", "paper", filename)
}