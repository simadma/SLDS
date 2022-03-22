file_names <- c(
  "collapse",
  "filtering",
  "particle_filtering",
  "score",
  "simulate_slds"
)
purrr::walk(file_names, ~ source(paste0("R/", .x, ".R")))
rm(file_names)