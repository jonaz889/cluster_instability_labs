# R/load_all.R

# utilities
source("R/utils/log_run.R")

# metrics
source("R/analysis/metrics.R")

# gauss experiments
source("R/simulation/eq_var_gauss_sim.R")
source("R/analysis/eq_var_gauss_summary_statistics.R")
source("R/analysis/drop_analysis.R")

# Gaus plots 
source("R/plots/gauss_plots.R")
source("R/plots/collapse_plots.R")
source("R/plots/drop_plots.R")

# rim experiments 
source("R/analysis/rim_alignment.R")
source("R/simulation/rim_sim.R")
source("R/analysis/rim_summary_statistics.R")
source("R/simulation/rim_hd_same_data.R")
source("R/plots/rim_plots.R")


# equilateral triangles
source("R/simulation/e3_fun.R")