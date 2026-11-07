# Runs the full analysis from the top-level folder of this repository.
# Usage: Rscript run_all.R
dir.create("tables", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
sink("results/run_output.txt", split = TRUE)
source("r/fbd_model.R")
sink()
