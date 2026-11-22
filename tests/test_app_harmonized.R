# Drives the app's harmonized-matrix mode headless with shiny::testServer.
# Usage: Rscript tests/test_app_harmonized.R <test_data_dir> <output_dir>
args <- commandArgs(trailingOnly = TRUE); td <- normalizePath(args[1]); out <- normalizePath(args[2])
library(shiny)
app <- shinyAppDir(".")
testServer(app, {
  session$setInputs(input_mode = "harmonized", project_name = "AppHarm", platform = "NULISA", tissue = "Plasma",
                    expr_file = file.path(td, "cleanDat.NULISA.Plasma.combined.TAMPOR.csv"),
                    traits_file = file.path(td, "numericMeta.combined_preProc.csv"),
                    output_root = out, matrix_scale = "log2", fixes_file = "", run_outliers = TRUE, test_bimodal = TRUE,
                    max_missing = 0.5, outlier_sd = 3, outlier_distance_screen = FALSE, n_boot = 20, n_cores = 1, seed = 1)
  session$setInputs(load_traits = 1)
  session$setInputs(traits_id_col = "row_names", group_col = "Group", ctl_value = "CTL", age_col = "age_at_sample",
                    sex_col = "Sex", sex_male_value = "M", sex_female_value = "F",
                    batch_col = "Batch", pair_col = "PairID", duplicate_mode = "keep_one", keep_batch = "y2025")
  cat("---- check_text\n", output$check_text, "\n")
  s <- match_summary(); print(table(s$group[s$in_data]))
  stopifnot(length(s$not_in_traits) == 0)
  session$setInputs(run = 1); session$setInputs(blank_remove = 1)
  cat("---- run_text\n", rv$run, "\n")
  stopifnot(!is.null(rv$report), file.exists(rv$report), grepl("_harmonized$", dirname(rv$report)))
  y <- yaml::read_yaml(file.path(dirname(rv$report), "run_parameters.yaml"))
  stopifnot(y$input_mode == "harmonized", y$duplicate_mode == "keep_one", y$keep_batch == "y2025", y$traits_id_col == "row_names")
  cat("HARMONIZED APP TEST PASSED\n")
})
