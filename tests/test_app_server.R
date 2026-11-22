# Drives the app's server logic headless with shiny::testServer.
# Usage: Rscript tests/test_app_server.R <test_data_dir> <output_dir>
args <- commandArgs(trailingOnly = TRUE); td <- normalizePath(args[1]); out <- normalizePath(args[2])
library(shiny)
app <- shinyAppDir(".")
testServer(app, {
  session$setInputs(project_name = "AppTest", platform = "NULISA", tissue = "Plasma",
                    expr_file = file.path(td, "nulisa_MACPath.csv"), traits_file = file.path(td, "nulisa_MACPath_traits.csv"),
                    output_root = out, expr_id_col = "SampleName", fixes_file = "", fixed_lod_file = "",
                    olink_npx = "raw", use_fixed_lod = FALSE, lod_cutoff = 0.5, sn_crit = 0.3, drop_qc_fail = TRUE,
                    outlier_sd = 3, outlier_distance_screen = FALSE, plate_correct = FALSE, n_boot = 20, n_cores = 1, seed = 1)
  session$setInputs(load_traits = 1)
  session$setInputs(traits_id_col = "ID", group_col = "Dx", ctl_value = "CTL", age_col = "AGE", sex_col = "SEX",
                    sex_male_value = "M", sex_female_value = "F")
  cat("---- check_text\n", output$check_text, "\n")
  s <- match_summary(); print(table(s$group[s$in_data]))
  session$setInputs(run = 1)              # blank groups present -> modal
  session$setInputs(blank_remove = 1)     # choose 'remove' -> runs
  cat("---- run_text\n", rv$run, "\n")
  stopifnot(!is.null(rv$report), file.exists(rv$report))
  session$setInputs(run = 2); session$setInputs(blank_case = 1)   # second run: folder exists -> overwrite modal
  session$setInputs(go_overwrite = 1)
  cat("---- second run\n", rv$run, "\n")
  y <- yaml::read_yaml(file.path(dirname(rv$report), "run_parameters.yaml")); stopifnot(y$blank_group == "case")
  cat("APP SERVER TEST PASSED\n")
})
