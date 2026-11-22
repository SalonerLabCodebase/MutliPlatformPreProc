###############################################################################
# Proteomics pre-processing app (SomaScan, OLink Explore HT, NULISA)
#
# Run from RStudio: open this file and click "Run App", or
#   shiny::runApp("path/to/proteomicsPreProcApp")
#
# The app collects the run settings, checks the inputs, and renders
# templates/preproc_report.Rmd (single cohort, raw platform file) or
# templates/harmonized_report.Rmd (harmonized cleanDat + numericMeta). The
# report does all the processing and writes the HTML report, cleanDat files,
# numericMeta and QC figures to a per-run folder.
###############################################################################

# ---- Install missing packages at start-up (conflicts row 26) -----------------
if (!nzchar(Sys.getenv("PREPROC_SKIP_INSTALL"))) local({
  cran <- c("shiny", "rmarkdown", "knitr", "dplyr", "tidyr", "arrow", "readxl", "diptest", "mclust", "boot",
            "doParallel", "foreach", "ggplot2", "yaml", "SomaDataIO", "rstudioapi", "WGCNA")
  bioc <- c("limma", "impute", "preprocessCore", "GO.db", "AnnotationDbi")
  miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_bioc)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(miss_bioc, ask = FALSE, update = FALSE)
  }
  miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran)) install.packages(miss_cran)
})

library(shiny)
app_dir <- normalizePath(getwd())
fn_file <- file.path(app_dir, "R", "preproc_functions.R")
templates <- c(raw = file.path(app_dir, "templates", "preproc_report.Rmd"),
               harmonized = file.path(app_dir, "templates", "harmonized_report.Rmd"))
source(fn_file, local = TRUE)

has_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
default_cores <- max(1, min(7, parallel::detectCores() - 1, na.rm = TRUE))
id_defaults <- c(SomaScan = "SampleId", OLink = "SampleID", NULISA = "SampleName")

path_input <- function(id, label, placeholder = "", button = TRUE, dir = FALSE) {
  div(style = "display:flex; gap:6px; align-items:flex-end;",
      div(style = "flex:1;", textInput(id, label, placeholder = placeholder, width = "100%")),
      if (button && has_rstudio) actionButton(paste0(id, "_browse"), if (dir) "Folder..." else "File...",
                                              style = "margin-bottom:15px;"))
}

ui <- fluidPage(
  titlePanel("Proteomics pre-processing"),
  sidebarLayout(
    sidebarPanel(width = 5,
      h4("1. Data"),
      radioButtons("input_mode", "Input", c("Single cohort: raw platform file" = "raw",
                                             "Harmonized matrix: cleanDat + numericMeta" = "harmonized")),
      textInput("project_name", "Project name (used in the output folder name)", "MyProject"),
      fluidRow(column(6, radioButtons("platform", "Platform", c("SomaScan", "OLink", "NULISA"), inline = TRUE)),
               column(6, radioButtons("tissue", "Tissue", c("CSF", "Plasma"), inline = TRUE))),
      path_input("expr_file", "Expression file (.adat, .parquet or NULISA .csv)"),
      path_input("traits_file", "Traits file (.csv)"),
      path_input("output_root", "Output folder", dir = TRUE),
      conditionalPanel("input.input_mode == 'raw'",
        textInput("expr_id_col", "Sample ID column in the expression file", id_defaults[["SomaScan"]])),
      conditionalPanel("input.input_mode == 'harmonized'",
        radioButtons("matrix_scale", "Matrix scale", c("log2 (TAMPOR output, NPQ, log2 RFU/NPX)" = "log2", "Linear" = "linear"), inline = TRUE),
        helpText("Rows are features (first column = feature ID, e.g. 'Target|UniProt'); columns are samples named as in the numericMeta ID column.")),
      actionButton("load_traits", "Load traits and check sample matching", class = "btn-primary"),
      hr(),
      h4("2. Groups and covariates"),
      uiOutput("traits_ui"),
      hr(),
      h4("3. Options"),
      conditionalPanel("input.input_mode == 'harmonized'",
        helpText("LOD and signal-to-noise filtering were done before harmonization and are skipped."),
        checkboxInput("run_outliers", "Remove sample outliers (connectivity)", TRUE),
        checkboxInput("test_bimodal", "Re-test and remove bimodal features", FALSE),
        numericInput("max_missing", "Remove features missing in more than this fraction of samples", 0.5, min = 0, max = 1, step = 0.05)),
      conditionalPanel("input.input_mode == 'raw'",
      conditionalPanel("input.platform == 'OLink'",
        radioButtons("olink_npx", "OLink values", c("Raw NPX" = "raw", "PC-normalised NPX" = "PC"), inline = TRUE),
        checkboxInput("use_fixed_lod", "Use Olink fixed LOD (switches to PC-normalised NPX)", FALSE),
        conditionalPanel("input.use_fixed_lod", path_input("fixed_lod_file", "Olink fixed LOD file"))),
      numericInput("lod_cutoff", "Keep features with at most this fraction of samples below LOD", 0.75, min = 0.05, max = 1, step = 0.05),
      conditionalPanel("input.platform != 'NULISA'",
        numericInput("sn_crit", "2-pass signal-to-noise threshold", 0.3, min = 0, step = 0.05)),
      checkboxInput("drop_qc_fail", "Drop wells that failed platform sample QC", TRUE)),
      numericInput("outlier_sd", "Outlier cutoff: connectivity z below minus", 3, min = 1, step = 0.5),
      checkboxInput("outlier_distance_screen", "Also flag outliers by distance from centroid / top PCs (robust z > 6)", FALSE),
      conditionalPanel("input.input_mode == 'raw'",
        checkboxInput("plate_correct", "Regress out plate (control/case difference protected)", FALSE)),
      fluidRow(column(4, numericInput("n_boot", "Bootstrap resamples", 1000, min = 10, step = 100)),
               column(4, numericInput("n_cores", "Cores", default_cores, min = 1, max = parallel::detectCores())),
               column(4, numericInput("seed", "Seed", 8675309))),
      path_input("fixes_file", "Dataset-fixes script (optional; see dataset_fixes_TEMPLATE.R)"),
      hr(),
      actionButton("run", "Run pre-processing", class = "btn-success btn-lg", width = "100%")
    ),
    mainPanel(width = 7,
      h4("Input check"),
      verbatimTextOutput("check_text"),
      tableOutput("group_table"),
      h4("Run"),
      verbatimTextOutput("run_text"),
      uiOutput("open_report")
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(traits = NULL, pre = NULL, check = "Load the traits file to begin.", run = "", report = NULL)

  # File and folder pickers (RStudio only; otherwise paste paths)
  for (id in c("expr_file", "traits_file", "fixed_lod_file", "fixes_file")) local({
    i <- id
    observeEvent(input[[paste0(i, "_browse")]], {
      f <- rstudioapi::selectFile(caption = "Select file")
      if (!is.null(f)) updateTextInput(session, i, value = f)
    })
  })
  observeEvent(input$output_root_browse, {
    d <- rstudioapi::selectDirectory(caption = "Select output folder")
    if (!is.null(d)) updateTextInput(session, "output_root", value = d)
  })
  observeEvent(input$input_mode, {
    h <- input$input_mode == "harmonized"
    updateTextInput(session, "expr_file", label = if (h) "Harmonized cleanDat (.csv, features x samples)" else "Expression file (.adat, .parquet or NULISA .csv)")
    updateTextInput(session, "traits_file", label = if (h) "numericMeta (.csv)" else "Traits file (.csv)")
    rv$traits <- NULL; rv$pre <- NULL; rv$check <- "Load the traits file to begin."
  }, ignoreInit = TRUE)
  observeEvent(input$platform, {
    updateTextInput(session, "expr_id_col", value = id_defaults[[input$platform]])
    updateNumericInput(session, "lod_cutoff", value = default_lod_cutoff(input$platform))
  })

  observeEvent(input$load_traits, {
    rv$traits <- NULL; rv$pre <- NULL
    tryCatch({
      validate_file <- function(p, what) if (!nzchar(p) || !file.exists(path.expand(p))) stop(what, " not found: ", p)
      validate_file(input$traits_file, "Traits file")
      tr <- fix_blank_colnames(utils::read.csv(path.expand(input$traits_file), check.names = FALSE, colClasses = "character", na.strings = "NA"))
      rv$traits <- tr
      msg <- sprintf("Traits: %d rows, %d columns.", nrow(tr), ncol(tr))
      if (nzchar(input$expr_file)) {
        validate_file(input$expr_file, "Expression file")
        if (identical(input$input_mode, "harmonized")) {
          hdr <- colnames(utils::read.csv(path.expand(input$expr_file), nrows = 1, check.names = FALSE))
          rv$pre <- list(ids = trimws(hdr[-1]), panel = NULL)
          msg <- c(msg, sprintf("Harmonized matrix: %d sample columns.", length(rv$pre$ids)))
        } else {
          withProgress(message = "Reading sample IDs from the expression file", {
            rv$pre <- preflight_expression(path.expand(input$expr_file), input$platform, input$expr_id_col)
          })
          msg <- c(msg, sprintf("Expression file: %d sample wells.", length(rv$pre$ids)))
        }
      }
      rv$check <- paste(msg, collapse = "\n")
    }, error = function(e) rv$check <- paste("Error:", conditionMessage(e)))
  })

  output$traits_ui <- renderUI({
    tr <- rv$traits
    if (is.null(tr)) return(helpText("Load the traits file first."))
    cols <- colnames(tr)
    guess <- function(pat) { g <- grep(pat, cols, ignore.case = TRUE, value = TRUE); if (length(g)) g[1] else "" }
    tagList(
      selectInput("traits_id_col", "Traits column holding the sample ID", cols, selected = cols[1]),
      selectInput("group_col", "Diagnostic grouping column", cols, selected = guess("group|dx|diag|clin")),
      uiOutput("ctl_ui"),
      selectInput("age_col", "Age column (blank = no age regression)", c("", cols), selected = guess("^age")),
      selectInput("sex_col", "Sex column (blank = no sex regression)", c("", cols), selected = guess("^sex")),
      uiOutput("sex_ui"),
      if (identical(input$input_mode, "harmonized")) tagList(
        selectInput("batch_col", "Batch column (MDS colouring; blank = none)", c("", cols), selected = guess("^batch")),
        selectInput("pair_col", "Column shared by repeat measurements of one sample (blank = none)", c("", cols), selected = guess("^pair")),
        conditionalPanel("input.pair_col != ''",
          radioButtons("duplicate_mode", "Repeat measurements (e.g. GIS samples run in both batches)",
                       c("Keep all" = "keep_all", "Keep one per sample" = "keep_one", "Average them" = "average")),
          conditionalPanel("input.duplicate_mode != 'keep_all'", uiOutput("keep_batch_ui"))))
    )
  })

  output$keep_batch_ui <- renderUI({
    req(rv$traits, input$batch_col)
    vals <- sort(unique(trimws(rv$traits[[input$batch_col]])))
    selectInput("keep_batch", "Keep the measurement (and its metadata) from batch", vals, selected = rev(vals)[1])
  })

  output$ctl_ui <- renderUI({
    req(rv$traits, input$group_col)
    v <- rv$traits[[input$group_col]]
    tab <- table(ifelse(is.na(v) | trimws(v) == "", "(blank)", trimws(v)))
    choices <- setdiff(names(tab), "(blank)")
    sel <- if ("CTL" %in% choices) "CTL" else choices[1]
    selectInput("ctl_value", "Value meaning control (every other non-blank value becomes Case)",
                setNames(choices, paste0(choices, " (n = ", tab[choices], ")")), selected = sel)
  })

  output$sex_ui <- renderUI({
    req(rv$traits)
    if (is.null(input$sex_col) || !nzchar(input$sex_col)) return(NULL)
    v <- trimws(rv$traits[[input$sex_col]]); vals <- sort(unique(v[!is.na(v) & v != ""]))
    gm <- vals[toupper(vals) %in% c("M", "MALE", "1")][1]; gf <- vals[toupper(vals) %in% c("F", "FEMALE", "0", "2")][1]
    fluidRow(column(6, selectInput("sex_male_value", "Value meaning male (coded 1)", vals, selected = gm %||% vals[1])),
             column(6, selectInput("sex_female_value", "Value meaning female (coded 2)", vals, selected = gf %||% vals[min(2, length(vals))])))
  })

  match_summary <- reactive({
    req(rv$traits, input$traits_id_col, input$group_col, input$ctl_value)
    tr <- rv$traits
    ids <- trimws(tr[[input$traits_id_col]])
    g <- trimws(tr[[input$group_col]])
    out <- list(n_blank = sum(is.na(g) | g == ""),
                group = ifelse(is.na(g) | g == "", "(blank)", ifelse(g == input$ctl_value, "CTL", "Case")))
    if (!is.null(rv$pre)) {
      out$not_in_traits <- setdiff(rv$pre$ids, ids)
      out$not_in_data <- setdiff(ids, rv$pre$ids)
      out$in_data <- ids %in% rv$pre$ids
    }
    out
  })

  output$check_text <- renderText({
    s <- tryCatch(match_summary(), error = function(e) NULL)
    if (is.null(s)) return(rv$check)
    lines <- rv$check
    if (!is.null(s$not_in_traits)) lines <- c(lines,
      sprintf("Samples in the expression file with no traits row: %d %s", length(s$not_in_traits),
              if (length(s$not_in_traits)) paste0("(", paste(head(s$not_in_traits, 15), collapse = ", "), if (length(s$not_in_traits) > 15) ", ..." , ")") else ""),
      sprintf("Traits rows with no sample in the expression file: %d", length(s$not_in_data)))
    lines <- c(lines, sprintf("Blank group values: %d", s$n_blank))
    paste(lines, collapse = "\n")
  })

  output$group_table <- renderTable({
    s <- tryCatch(match_summary(), error = function(e) NULL)
    req(s)
    keep <- if (is.null(s$in_data)) TRUE else s$in_data
    as.data.frame(table(`Pre-processing group (matched samples)` = s$group[keep]))
  })

  # ---- Run -------------------------------------------------------------------
  collect_params <- function(blank_group) {
    if (identical(input$input_mode, "harmonized")) return(list(
         project_name = input$project_name, platform = input$platform, tissue = input$tissue,
         expr_file = path.expand(input$expr_file), traits_file = path.expand(input$traits_file),
         output_root = path.expand(input$output_root), traits_id_col = input$traits_id_col,
         matrix_scale = input$matrix_scale, group_col = input$group_col, ctl_value = input$ctl_value,
         blank_group = blank_group, age_col = input$age_col %||% "", sex_col = input$sex_col %||% "",
         sex_male_value = if (nzchar(input$sex_col %||% "")) input$sex_male_value else "",
         sex_female_value = if (nzchar(input$sex_col %||% "")) input$sex_female_value else "",
         batch_col = input$batch_col %||% "", pair_col = input$pair_col %||% "",
         duplicate_mode = if (nzchar(input$pair_col %||% "")) input$duplicate_mode %||% "keep_all" else "keep_all",
         keep_batch = input$keep_batch %||% "", run_outliers = input$run_outliers, outlier_sd = input$outlier_sd,
         outlier_distance_screen = input$outlier_distance_screen, test_bimodal = input$test_bimodal,
         max_missing = input$max_missing, n_boot = input$n_boot, n_cores = input$n_cores, seed = input$seed,
         fixes_file = path.expand(input$fixes_file), functions_file = fn_file))
    list(project_name = input$project_name, platform = input$platform, tissue = input$tissue,
         expr_file = path.expand(input$expr_file), traits_file = path.expand(input$traits_file),
         output_root = path.expand(input$output_root), traits_id_col = input$traits_id_col,
         expr_id_col = input$expr_id_col, group_col = input$group_col, ctl_value = input$ctl_value,
         blank_group = blank_group, age_col = input$age_col %||% "", sex_col = input$sex_col %||% "",
         sex_male_value = if (nzchar(input$sex_col %||% "")) input$sex_male_value else "",
         sex_female_value = if (nzchar(input$sex_col %||% "")) input$sex_female_value else "",
         drop_qc_fail = input$drop_qc_fail, olink_npx = input$olink_npx, use_fixed_lod = input$use_fixed_lod,
         fixed_lod_file = path.expand(input$fixed_lod_file), lod_cutoff = input$lod_cutoff,
         sn_crit = input$sn_crit, outlier_sd = input$outlier_sd,
         outlier_distance_screen = input$outlier_distance_screen, plate_correct = input$plate_correct,
         n_boot = input$n_boot, n_cores = input$n_cores, seed = input$seed,
         fixes_file = path.expand(input$fixes_file), functions_file = fn_file)
  }

  validate_inputs <- function() {
    errs <- character(0)
    if (!nzchar(input$project_name)) errs <- c(errs, "Project name is empty.")
    for (f in c("expr_file", "traits_file")) if (!file.exists(path.expand(input[[f]]))) errs <- c(errs, paste("File not found:", input[[f]]))
    if (!dir.exists(path.expand(input$output_root))) errs <- c(errs, paste("Output folder not found:", input$output_root))
    if (is.null(rv$traits) || is.null(input$ctl_value)) errs <- c(errs, "Load the traits file and choose the control value.")
    if (nzchar(input$fixes_file) && !file.exists(path.expand(input$fixes_file))) errs <- c(errs, "Dataset-fixes script not found.")
    if (identical(input$input_mode, "harmonized") && !identical(input$duplicate_mode %||% "keep_all", "keep_all") &&
        nzchar(input$pair_col %||% "") && !nzchar(input$batch_col %||% "")) errs <- c(errs, "Choose a batch column so the app knows which repeat to keep.")
    if (!identical(input$input_mode, "harmonized") && input$platform == "OLink" && input$use_fixed_lod && !file.exists(path.expand(input$fixed_lod_file))) errs <- c(errs, "Fixed LOD file not found.")
    if (nzchar(input$sex_col %||% "") && identical(input$sex_male_value, input$sex_female_value)) errs <- c(errs, "Male and female sex values are the same.")
    errs
  }

  pending_blank <- reactiveVal("remove")

  observeEvent(input$run, {
    errs <- validate_inputs()
    if (length(errs)) { rv$run <- paste(c("Cannot run:", errs), collapse = "\n"); return() }
    s <- match_summary()
    nb <- sum(s$group[if (is.null(s$in_data)) TRUE else s$in_data] == "(blank)")
    if (nb > 0) {
      showModal(modalDialog(title = "Samples with a blank group",
        sprintf("%d matched samples have a blank value in '%s'. By default they are removed. Remove them, or label them Case?", nb, input$group_col),
        footer = tagList(actionButton("blank_remove", "Remove them"), actionButton("blank_case", "Label them Case"), modalButton("Cancel"))))
    } else confirm_overwrite("remove")
  })
  observeEvent(input$blank_remove, { removeModal(); confirm_overwrite("remove") })
  observeEvent(input$blank_case, { removeModal(); confirm_overwrite("case") })

  confirm_overwrite <- function(blank_group) {
    pending_blank(blank_group)
    harm <- identical(input$input_mode, "harmonized")
    npx <- if (!harm && input$platform == "OLink") (if (input$use_fixed_lod) "PC" else input$olink_npx) else NULL
    rd <- run_dir_path(path.expand(input$output_root), input$project_name, input$platform, input$tissue,
                       panel = if (harm) "harmonized" else if (input$platform == "NULISA" && !is.null(rv$pre)) rv$pre$panel else NULL,
                       npx_type = npx)
    if (dir.exists(rd)) {
      showModal(modalDialog(title = "Output folder exists",
        paste("This run will overwrite files in", rd), footer = tagList(actionButton("go_overwrite", "Overwrite"), modalButton("Cancel"))))
    } else do_run(blank_group)
  }
  observeEvent(input$go_overwrite, { removeModal(); do_run(pending_blank()) })

  do_run <- function(blank_group) {
    prm <- collect_params(blank_group)
    rv$run <- "Running... the bootstrap regression can take a while. The app is busy until it finishes."
    env <- new.env()
    t0 <- Sys.time()
    res <- withProgress(message = "Pre-processing (see the R console for progress)", value = 0.5, {
      tryCatch(rmarkdown::render(templates[[input$input_mode %||% "raw"]], params = prm, envir = env, output_dir = tempdir(),
                                 output_file = paste0(input$project_name, "_report.html"),
                                 intermediates_dir = tempdir(), knit_root_dir = tempdir(), quiet = FALSE),
               error = function(e) e)
    })
    if (inherits(res, "error")) {
      rv$run <- paste("The run stopped with an error:\n", conditionMessage(res),
                      "\nNothing after that step was written. Fix the input (a dataset-fixes script handles dataset quirks) and run again.")
      rv$report <- NULL
      return()
    }
    final <- file.path(env$run_dir, basename(res))
    file.copy(res, final, overwrite = TRUE)
    rv$report <- final
    rv$run <- sprintf("Finished in %.1f minutes.\nOutputs: %s\nReport: %s",
                      as.numeric(difftime(Sys.time(), t0, units = "mins")), env$run_dir, final)
  }

  output$run_text <- renderText(rv$run)
  output$open_report <- renderUI({ req(rv$report); actionButton("open_html", "Open report") })
  observeEvent(input$open_html, utils::browseURL(rv$report))
}

shinyApp(ui, server)
