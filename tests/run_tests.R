###############################################################################
# run_tests.R: render the report for each test configuration and check outputs.
# Usage: Rscript tests/run_tests.R <test_data_dir> <output_dir> [config names]
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
td <- normalizePath(args[1]); out <- args[2]; only <- args[-(1:2)]
dir.create(out, showWarnings = FALSE, recursive = TRUE); out <- normalizePath(out)
app <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
base <- list(group_col = "Dx", ctl_value = "CTL", age_col = "AGE", sex_col = "SEX", n_boot = 30, n_cores = 2,
             output_root = out, functions_file = file.path(app, "R", "preproc_functions.R"))
configs <- list(
  soma_plasma = list(project_name = "SomaTest", platform = "SomaScan", tissue = "Plasma",
                     expr_file = file.path(td, "soma_example.adat"), traits_file = file.path(td, "soma_traits.csv"),
                     traits_id_col = "SampleId"),
  soma_csf_options = list(project_name = "SomaTestOpts", platform = "SomaScan", tissue = "CSF",
                     expr_file = file.path(td, "soma_example.adat"), traits_file = file.path(td, "soma_traits.csv"),
                     traits_id_col = "SampleId", outlier_distance_screen = TRUE, plate_correct = TRUE, blank_group = "case"),
  olink_allftd = list(project_name = "OlinkALLFTD", platform = "OLink", tissue = "CSF",
                     expr_file = file.path(td, "olink_ALLFTD.parquet"), traits_file = file.path(td, "olink_ALLFTD_traits.csv")),
  olink_macwide = list(project_name = "OlinkMAC", platform = "OLink", tissue = "CSF", sex_male_value = "1", sex_female_value = "0",
                     expr_file = file.path(td, "olink_MACWide.parquet"), traits_file = file.path(td, "olink_MACWide_traits.csv"),
                     drop_qc_fail = FALSE),
  nulisa_allftd = list(project_name = "NulisaALLFTD", platform = "NULISA", tissue = "CSF", blank_group = "case",
                     expr_file = file.path(td, "nulisa_ALLFTD.csv"), traits_file = file.path(td, "nulisa_ALLFTD_traits.csv")),
  nulisa_macpath = list(project_name = "NulisaMACPath", platform = "NULISA", tissue = "Plasma",
                     expr_file = file.path(td, "nulisa_MACPath.csv"), traits_file = file.path(td, "nulisa_MACPath_traits.csv")),
  harmonized_keepall = list(template = "harmonized_report.Rmd", project_name = "ALLFTDHarm", platform = "NULISA", tissue = "Plasma",
                     expr_file = file.path(td, "cleanDat.NULISA.Plasma.combined.TAMPOR.csv"),
                     traits_file = file.path(td, "numericMeta.combined_preProc.csv"),
                     group_col = "Group", age_col = "age_at_sample", sex_col = "Sex", batch_col = "Batch", pair_col = "PairID"),
  harmonized_keepone = list(template = "harmonized_report.Rmd", project_name = "ALLFTDHarmOne", platform = "NULISA", tissue = "Plasma",
                     expr_file = file.path(td, "cleanDat.NULISA.Plasma.combined.TAMPOR.csv"),
                     traits_file = file.path(td, "numericMeta.combined_preProc.csv"), blank_group = "case",
                     group_col = "Group", age_col = "age_at_sample", sex_col = "Sex", sex_male_value = "M", sex_female_value = "F",
                     batch_col = "Batch", pair_col = "PairID", duplicate_mode = "keep_one", keep_batch = "y2025",
                     test_bimodal = TRUE, outlier_distance_screen = TRUE),
  harmonized_average = list(template = "harmonized_report.Rmd", project_name = "ALLFTDHarmAvg", platform = "NULISA", tissue = "Plasma",
                     expr_file = file.path(td, "cleanDat.NULISA.Plasma.combined.TAMPOR.csv"),
                     traits_file = file.path(td, "numericMeta.combined_preProc.csv"),
                     group_col = "Group", age_col = "age_at_sample", sex_col = "Sex", batch_col = "Batch", pair_col = "PairID",
                     duplicate_mode = "average", keep_batch = "y2024", run_outliers = FALSE)
)
`%||%` <- function(a, b) if (is.null(a)) b else a
if (length(only)) configs <- configs[only]
results <- list()
for (nm in names(configs)) {
  prm <- modifyList(base, configs[[nm]])
  tmpl <- prm$template %||% "preproc_report.Rmd"; prm$template <- NULL
  t0 <- Sys.time()
  res <- tryCatch({
    f <- rmarkdown::render(file.path(app, "templates", tmpl), params = prm, envir = new.env(),
                           output_dir = out, output_file = paste0(nm, ".html"), intermediates_dir = tempdir(),
                           knit_root_dir = tempdir(), quiet = TRUE)
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  results[[nm]] <- res
  cat(sprintf("%-18s %s (%.1f s)\n", nm, res, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
if (any(results != "OK")) quit(status = 1)
