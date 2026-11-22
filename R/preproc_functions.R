###############################################################################
# preproc_functions.R
#
# Function library for the proteomics pre-processing app (SomaScan, OLink
# Explore HT, NULISA; CSF or plasma). Sourced by templates/preproc_report.Rmd.
#
# Consolidated from the 2024-2026 notebooks (files 1-16 in the audit). Every
# function takes its inputs as arguments and returns its outputs: nothing reads
# or writes the global environment, so the report can be rendered from Shiny
# in its own environment (audit B10, conflicts row 27).
#
# Design contract: whatever the platform, build_dataset() returns the same
# objects (see its header), so no step after it needs to ask which platform
# produced the data except where the science differs (LOD, S:N, scale).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a)) || identical(a, "")) b else a

log_msg <- function(...) {
  cat(paste0(...), "\n", sep = "")
  invisible(NULL)
}

row_median <- function(m) apply(m, 1, stats::median, na.rm = TRUE)
row_sd     <- function(m) apply(m, 1, stats::sd, na.rm = TRUE)
row_mad    <- function(m) apply(m, 1, stats::mad, na.rm = TRUE)

# Fill a features x wells matrix from long-format vectors without pivot_wider.
# Much faster than tidyr for the 10M-row Olink exports, and it cannot silently
# create list-columns: duplicated (feature, well) pairs are an error.
long_to_matrix <- function(feature, well, value, feature_levels = NULL, well_levels = NULL) {
  feature_levels <- feature_levels %||% unique(feature)
  well_levels    <- well_levels %||% unique(well)
  ri <- match(feature, feature_levels)
  ci <- match(well, well_levels)
  key <- ri + (ci - 1) * length(feature_levels)
  if (anyDuplicated(key)) {
    d <- which(duplicated(key))[1:min(5, sum(duplicated(key)))]
    stop("Duplicate measurements for the same feature and well, e.g. ",
         paste0(feature[d], " / ", well[d], collapse = "; "),
         ". Fix them in a dataset-fixes script before running.")
  }
  m <- matrix(NA_real_, nrow = length(feature_levels), ncol = length(well_levels),
              dimnames = list(feature_levels, well_levels))
  m[cbind(ri, ci)] <- value
  m
}

# Make labels unique by appending a suffix only to the repeated ones
dedupe_labels <- function(labels, suffix) {
  dup <- labels %in% labels[duplicated(labels)]
  labels[dup] <- paste0(labels[dup], "^", suffix[dup])
  labels
}

# -----------------------------------------------------------------------------
# Platform defaults
# -----------------------------------------------------------------------------

# Display labels for the four standard well types, per platform. The same
# vector feeds colouring and legends, so the two can never disagree
# (audit: 13/16 plate- vs sample-control legend swap).
well_type_labels <- function(platform) {
  switch(platform,
    SomaScan = c(Sample = "Sample", Buffer = "Buffer", Calibrator = "Calibrator", QC = "QC"),
    OLink    = c(Sample = "Sample", Buffer = "Negative control", Calibrator = "Plate control", QC = "Sample control"),
    NULISA   = c(Sample = "Sample", Buffer = "Negative control (NC)", Calibrator = "Inter-plate control (IPC)", QC = "Sample control (SC)"),
    stop("Unknown platform: ", platform))
}

default_lod_cutoff <- function(platform) if (platform == "NULISA") 0.50 else 0.75

# -----------------------------------------------------------------------------
# Reading inputs
# -----------------------------------------------------------------------------

#' A CSV written by write.csv() with row names has a blank first header
#' (numericMeta files do). Name it "row_names" so it can be chosen as the ID.
fix_blank_colnames <- function(df) {
  cn <- colnames(df)
  blank <- is.na(cn) | cn == ""
  if (any(blank)) {
    cn[blank] <- ifelse(which(blank) == 1, "row_names", paste0("V", which(blank)))
    colnames(df) <- cn
  }
  df
}

#' Read the traits (numericMeta) file.
#' id_col = "" uses the first column as the sample ID (the old row.names = 1
#' behaviour); otherwise the named column is used (conflicts row 17).
read_traits <- function(path, id_col = "") {
  tr <- utils::read.csv(path, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE, colClasses = "character",
                        na.strings = c("NA"))
  tr <- fix_blank_colnames(tr)
  if (identical(id_col, "") || is.null(id_col)) id_col <- colnames(tr)[1]
  if (!id_col %in% colnames(tr)) stop("Traits ID column '", id_col, "' not found in ", basename(path))
  ids <- trimws(tr[[id_col]])
  if (anyNA(ids) || any(ids == "")) stop("Traits ID column '", id_col, "' has blank IDs.")
  if (anyDuplicated(ids)) stop("Traits ID column '", id_col, "' has duplicated IDs: ",
                               paste(head(unique(ids[duplicated(ids)]), 10), collapse = ", "))
  # Convert columns that are fully numeric back to numeric (colClasses read all as text
  # so that IDs with leading zeros survive)
  for (cc in colnames(tr)) {
    if (cc == id_col) next
    v <- tr[[cc]]
    vv <- suppressWarnings(as.numeric(v))
    if (all(is.na(v) | v == "" | !is.na(vv))) tr[[cc]] <- vv
  }
  rownames(tr) <- ids
  attr(tr, "id_col") <- id_col
  tr
}

#' Read the expression file in its raw platform format.
read_expression_raw <- function(path, platform) {
  if (!file.exists(path)) stop("Expression file not found: ", path)
  switch(platform,
    SomaScan = SomaDataIO::read_adat(path),
    OLink = {
      cols <- names(arrow::open_dataset(path)$schema)
      want <- intersect(c("SampleID", "SampleType", "WellID", "PlateID", "OlinkID", "UniProt",
                          "Assay", "AssayType", "Panel", "Block", "Count", "NPX",
                          "PCNormalizedNPX", "SampleQC", "AssayQC", "LOD", "PCNormalizedLOD"),
                        cols)
      as.data.frame(arrow::read_parquet(path, col_select = dplyr::all_of(want)))
    },
    NULISA = utils::read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unknown platform: ", platform))
}

#' Source an optional dataset-fixes script (audit item 7). The script sees
#' `expr.raw`, `traits.raw`, `platform` and `tissue`, may modify `expr.raw`
#' and/or `traits.raw`, and nothing else is kept from it.
apply_dataset_fixes <- function(fixes_path, expr.raw, traits.raw, platform, tissue) {
  if (is.null(fixes_path) || identical(fixes_path, "") || !file.exists(fixes_path))
    return(list(expr.raw = expr.raw, traits.raw = traits.raw, applied = FALSE))
  env <- new.env(parent = globalenv())
  env$expr.raw <- expr.raw
  env$traits.raw <- traits.raw
  env$platform <- platform
  env$tissue <- tissue
  sys.source(fixes_path, envir = env)
  list(expr.raw = env$expr.raw, traits.raw = env$traits.raw, applied = TRUE)
}

#' Light pre-flight read for the app: sample IDs and panel only, without
#' building the full matrix (the Olink read selects 2-3 columns).
preflight_expression <- function(path, platform, sample_id_col = "") {
  if (platform == "SomaScan") {
    a <- SomaDataIO::read_adat(path)
    idc <- sample_id_col %||% "SampleId"
    ids <- unique(as.character(a[[idc]][a$SampleType == "Sample"]))
    return(list(ids = trimws(ids), panel = NULL, id_cols = setdiff(colnames(a), SomaDataIO::getAnalytes(a))))
  }
  if (platform == "OLink") {
    idc <- sample_id_col %||% "SampleID"
    cols <- names(arrow::open_dataset(path)$schema)
    d <- as.data.frame(arrow::read_parquet(path, col_select = dplyr::all_of(intersect(c(idc, "SampleType", "Panel"), cols))))
    return(list(ids = trimws(unique(d[[idc]][toupper(d$SampleType) == "SAMPLE"])), panel = NULL, id_cols = cols))
  }
  idc <- sample_id_col %||% "SampleName"
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  list(ids = trimws(unique(d[[idc]][toupper(d$SampleType) == "SAMPLE"])), panel = unique(d$Panel), id_cols = colnames(d))
}

# -----------------------------------------------------------------------------
# Building the common dataset object
# -----------------------------------------------------------------------------
# Every builder returns a list with:
#   exprMat0  features x wells. Linear scale for SomaScan (RFU) and OLink (2^NPX),
#             log2 NPQ for NULISA. Includes all control wells.
#   wells     data.frame, one row per column of exprMat0 (same order):
#             well_id, sample_id (NA for controls), well_type (Sample, Buffer,
#             Calibrator, QC, Other), plate, qc_fail (logical)
#   featData  data.frame, one row per row of exprMat0 (same order, same
#             rownames): label, symbol, uniprot, is_human_protein, plus
#             platform columns
#   lodMat    NULISA only: per-well target LOD (same dims as exprMat0); NULL otherwise
#   olinkLOD  OLink only, when present in the file: per-well Olink LOD (NPX scale)
#   scale     "linear" or "log2"
#   panel     panel name(s) for output naming (NULISA, OLink)
#   notes     character vector of messages for the report
# -----------------------------------------------------------------------------

# Symbols SomaLogic ships as "None" (repaired at the feature-ID step, conflicts row 22)
SOMA_NONE_SYMBOLS <- c(F5H284 = "PPIAL4D", K9M1U5 = "IFNL4", H3BQW9 = "FAM229A", C9JTQ0 = "ANKRD63",
                       I3L1E1 = "C19orf84", B7ZBB8 = "PPP1R3G", P30042 = "GATD3A", B0FP48 = "UPK3BL1",
                       Q8NFS9 = "GCNT2", C9JXX5 = "C11orf94", E0CX11 = "STMP1", B9ZVM9 = "TCP10L2",
                       B1AK53 = "ESPN", B2RUY7 = "VWC2L", B1AKI9 = "ISM1")

#' SomaScan feature IDs: "Symbol|UniProt", with "^SomaId@AptName" appended to
#' repeated targets. Symbols reported as "None" are repaired here so every
#' matrix inherits the corrected names.
build_soma_feature_ids <- function(featData) {
  fd <- featData
  nonprot <- fd$Type != "Protein"
  fd$UniProt[nonprot] <- fd$SomaId[nonprot]
  fd$EntrezGeneSymbol[nonprot] <- fd$AptName[nonprot]
  fd$EntrezGeneSymbol[fd$UniProt %in% "Q99LC4"] <- fd$AptName[fd$UniProt %in% "Q99LC4"]
  gene <- ifelse(is.na(fd$EntrezGeneSymbol) | fd$EntrezGeneSymbol == "", "NA", fd$EntrezGeneSymbol)
  uni  <- ifelse(is.na(fd$UniProt) | fd$UniProt == "", "NA", fd$UniProt)
  gene <- vapply(strsplit(gene, "|", fixed = TRUE), function(p) if (length(p)) p[1] else "NA", "")
  uni  <- vapply(strsplit(uni,  "|", fixed = TRUE), function(p) if (length(p)) p[1] else "NA", "")
  none <- gene == "None"
  fixable <- none & uni %in% names(SOMA_NONE_SYMBOLS)
  gene[fixable] <- SOMA_NONE_SYMBOLS[uni[fixable]]
  unrepaired <- unique(uni[none & !fixable])
  ids <- paste0(gene, "|", uni)
  dup <- ids %in% ids[duplicated(ids)]
  ids[dup] <- paste0(ids[dup], "^", fd$SomaId[dup], "@", fd$AptName[dup])
  if (anyDuplicated(ids)) stop("SomaScan feature IDs are not unique after disambiguation.")
  list(ids = ids, symbol = gene, uniprot = uni, unrepaired_none = unrepaired)
}

build_dataset_somascan <- function(adat, sample_id_col = "SampleId", drop_qc_fail = TRUE) {
  notes <- character(0)
  if (!SomaDataIO::is.soma_adat(adat)) stop("SomaScan input is not a soma_adat object.")
  analytes <- SomaDataIO::getAnalytes(adat)
  meta_cols <- setdiff(colnames(adat), analytes)
  notes <- c(notes, sprintf("Adat: %d wells, %d analytes, %d metadata columns (found by name, not position).",
                            nrow(adat), length(analytes), length(meta_cols)))
  if (!sample_id_col %in% meta_cols) stop("Sample ID column '", sample_id_col, "' not in the adat metadata. Available: ",
                                          paste(meta_cols, collapse = ", "))
  meta <- as.data.frame(adat[, meta_cols, drop = FALSE])
  info <- as.data.frame(SomaDataIO::getAnalyteInfo(adat))
  info <- info[match(analytes, info$AptName), , drop = FALSE]
  ids <- build_soma_feature_ids(info)
  if (length(ids$unrepaired_none))
    notes <- c(notes, paste0("SomaScan symbols still 'None' (UniProt not in the repair table): ",
                             paste(ids$unrepaired_none, collapse = ", ")))

  st <- as.character(meta$SampleType)
  well_type <- ifelse(st %in% c("Sample", "Buffer", "Calibrator", "QC"), st, "Other")
  sid <- trimws(as.character(meta[[sample_id_col]]))
  sid[well_type != "Sample"] <- NA
  plate <- as.character(meta$PlateId %||% NA)
  pos   <- as.character(meta$PlatePosition %||% seq_len(nrow(meta)))
  well_id <- ifelse(well_type == "Sample", sid, paste(well_type, plate, pos, sep = "."))
  if (anyNA(well_id) || any(well_id == "")) stop("Some Sample wells have a blank '", sample_id_col, "'.")
  if (anyDuplicated(well_id)) stop("Duplicated sample IDs in the adat: ",
                                   paste(head(unique(well_id[duplicated(well_id)]), 10), collapse = ", "),
                                   ". Resolve them in a dataset-fixes script.")
  # SomaLogic's RowCheck FLAG is a normalisation-scale note, not a failed sample,
  # so it is reported but never used to drop wells.
  if ("RowCheck" %in% meta_cols)
    notes <- c(notes, sprintf("RowCheck FLAG on %d of %d wells (reported only; not used to drop wells).",
                              sum(meta$RowCheck == "FLAG", na.rm = TRUE), nrow(meta)))

  m <- t(as.matrix(as.data.frame(adat[, analytes, drop = FALSE])))
  storage.mode(m) <- "double"
  dimnames(m) <- list(ids$ids, well_id)

  featData <- data.frame(label = ids$ids, symbol = ids$symbol, uniprot = ids$uniprot,
                         is_human_protein = info$Type == "Protein" & info$Organism == "Human",
                         AptName = info$AptName, SeqId = info$SeqId, SomaId = info$SomaId,
                         TargetFullName = info$TargetFullName, Target = info$Target,
                         EntrezGeneSymbol = info$EntrezGeneSymbol, Type = info$Type,
                         Organism = info$Organism, ColCheck = info$ColCheck %||% NA,
                         row.names = ids$ids, stringsAsFactors = FALSE, check.names = FALSE)

  wells <- data.frame(well_id = well_id, sample_id = sid, well_type = well_type, plate = plate,
                      position = pos, qc_fail = FALSE, stringsAsFactors = FALSE)
  list(exprMat0 = m, wells = wells, featData = featData, lodMat = NULL, olinkLOD = NULL,
       scale = "linear", panel = NULL, notes = notes)
}

build_dataset_olink <- function(df, npx_type = c("raw", "PC"), sample_id_col = "SampleID",
                                drop_qc_fail = TRUE) {
  npx_type <- match.arg(npx_type)
  notes <- character(0)
  need <- c(sample_id_col, "SampleType", "PlateID", "WellID", "OlinkID", "Assay", "UniProt")
  miss <- setdiff(need, colnames(df))
  if (length(miss)) stop("OLink file is missing columns: ", paste(miss, collapse = ", "))
  val_col <- if (npx_type == "raw") "NPX" else "PCNormalizedNPX"
  if (!val_col %in% colnames(df)) stop("OLink file has no '", val_col, "' column for npx type '", npx_type, "'.")

  # Control assays out (AssayType where present, assay name otherwise)
  if ("AssayType" %in% colnames(df)) {
    n0 <- length(unique(df$OlinkID))
    df <- df[df$AssayType == "assay", , drop = FALSE]
    notes <- c(notes, sprintf("Removed %d control assays (AssayType != 'assay').", n0 - length(unique(df$OlinkID))))
  } else {
    keep <- !grepl("control", df$Assay, ignore.case = TRUE)
    notes <- c(notes, sprintf("No AssayType column; removed %d rows whose Assay name contains 'control'.", sum(!keep)))
    df <- df[keep, , drop = FALSE]
  }

  st <- toupper(df$SampleType)
  well_type <- dplyr::case_when(st == "SAMPLE" ~ "Sample", st == "NEGATIVE_CONTROL" ~ "Buffer",
                                st == "PLATE_CONTROL" ~ "Calibrator", st == "SAMPLE_CONTROL" ~ "QC",
                                TRUE ~ "Other")
  sid <- trimws(as.character(df[[sample_id_col]]))
  well_id <- ifelse(well_type == "Sample", sid, paste(well_type, df$PlateID, df$WellID, sep = "."))

  # Wells that failed Olink sample QC in any block are dropped as whole samples
  # (default), or kept with drop_qc_fail = FALSE.
  fail_wells <- character(0)
  if ("SampleQC" %in% colnames(df)) {
    fw <- unique(well_id[toupper(df$SampleQC) == "FAIL"])
    fail_wells <- fw
    notes <- c(notes, sprintf("Olink SampleQC FAIL in at least one block: %d wells (%d of them samples).",
                              length(fw), length(intersect(fw, well_id[well_type == "Sample"]))))
  }

  feat_levels <- unique(df$OlinkID)
  well_levels <- unique(well_id)
  m_npx <- long_to_matrix(df$OlinkID, well_id, df[[val_col]], feat_levels, well_levels)
  olinkLOD <- NULL
  lod_col <- if (npx_type == "raw") "LOD" else "PCNormalizedLOD"
  if (lod_col %in% colnames(df) && any(!is.na(df[[lod_col]]))) {
    olinkLOD <- long_to_matrix(df$OlinkID, well_id, df[[lod_col]], feat_levels, well_levels)
  } else {
    notes <- c(notes, sprintf("No usable '%s' column in this export; Olink's own LOD is not available.", lod_col))
  }

  fi <- df[!duplicated(df$OlinkID), c("OlinkID", "Assay", "UniProt",
                                       intersect(c("Panel", "Block"), colnames(df))), drop = FALSE]
  fi <- fi[match(feat_levels, fi$OlinkID), , drop = FALSE]
  label <- dedupe_labels(paste0(fi$Assay, "|", fi$UniProt), fi$OlinkID)
  rownames(m_npx) <- label
  if (!is.null(olinkLOD)) rownames(olinkLOD) <- label
  featData <- data.frame(label = label, symbol = fi$Assay, uniprot = fi$UniProt,
                         is_human_protein = TRUE, OlinkID = fi$OlinkID, Assay = fi$Assay,
                         Panel = fi$Panel %||% NA, Block = fi$Block %||% NA,
                         row.names = label, stringsAsFactors = FALSE)

  wi <- !duplicated(well_id)
  wells <- data.frame(well_id = well_id[wi], sample_id = ifelse(well_type[wi] == "Sample", sid[wi], NA),
                      well_type = well_type[wi], plate = as.character(df$PlateID[wi]),
                      position = as.character(df$WellID[wi]), stringsAsFactors = FALSE)
  wells <- wells[match(well_levels, wells$well_id), , drop = FALSE]
  wells$qc_fail <- wells$well_id %in% fail_wells
  dup_s <- duplicated(wells$sample_id) & !is.na(wells$sample_id)
  if (any(dup_s)) stop("Sample IDs appear in more than one well: ", paste(head(wells$sample_id[dup_s], 10), collapse = ", "))

  # Linear scale, following Dammer et al. 2022 (template convention)
  list(exprMat0 = 2^m_npx, wells = wells, featData = featData, lodMat = NULL, olinkLOD = olinkLOD,
       scale = "linear", panel = unique(featData$Panel), npx_type = npx_type, notes = notes)
}

build_dataset_nulisa <- function(df, sample_id_col = "SampleName", drop_qc_fail = TRUE) {
  notes <- character(0)
  need <- c(sample_id_col, "SampleType", "PlateID", "Target", "UniProtID", "NPQ", "targetLOD_NPQ")
  miss <- setdiff(need, colnames(df))
  if (length(miss)) stop("NULISA file is missing columns: ", paste(miss, collapse = ", "))
  st <- toupper(df$SampleType)
  well_type <- dplyr::case_when(st == "SAMPLE" ~ "Sample", st == "NC" ~ "Buffer",
                                st == "IPC" ~ "Calibrator", st == "SC" ~ "QC", TRUE ~ "Other")
  sid <- trimws(as.character(df[[sample_id_col]]))
  feat <- paste0(df$Target, "|", df$UniProtID)
  # Control wells: number repeated control names within each plate and target
  # in file order (the same order for every target in NULISA exports)
  base <- ifelse(well_type == "Sample", sid, paste(well_type, df$PlateID, sid, sep = "."))
  ord <- stats::ave(seq_along(base), feat, base, FUN = seq_along)
  well_id <- ifelse(well_type == "Sample", sid, paste0(base, ".", ord))
  if (any(well_type == "Sample" & ord > 1))
    stop("Sample IDs measured more than once for the same target: ",
         paste(head(unique(sid[well_type == "Sample" & ord > 1]), 10), collapse = ", "),
         ". Resolve them in a dataset-fixes script.")

  fail_wells <- character(0)
  if ("SampleQC" %in% colnames(df)) {
    fw <- unique(well_id[tolower(trimws(df$SampleQC)) != "passed"])
    fail_wells <- fw
    notes <- c(notes, sprintf("NULISA SampleQC not 'passed': %d wells (%d of them samples).",
                              length(fw), length(intersect(fw, well_id[well_type == "Sample"]))))
  }
  feat_levels <- unique(feat)
  well_levels <- unique(well_id)
  m   <- long_to_matrix(feat, well_id, df$NPQ, feat_levels, well_levels)
  lod <- long_to_matrix(feat, well_id, df$targetLOD_NPQ, feat_levels, well_levels)
  n_zero <- sum(lod == 0, na.rm = TRUE)
  if (n_zero) notes <- c(notes, sprintf("targetLOD_NPQ = 0 in %d target-well cells; kept as 0, so those values count as above LOD.", n_zero))

  fi <- df[!duplicated(feat), , drop = FALSE]
  fi <- fi[match(feat_levels, paste0(fi$Target, "|", fi$UniProtID)), , drop = FALSE]
  featData <- data.frame(label = feat_levels, symbol = fi$Target, uniprot = fi$UniProtID,
                         is_human_protein = TRUE, Target = fi$Target,
                         ProteinName = fi$ProteinName %||% NA, row.names = feat_levels,
                         stringsAsFactors = FALSE)
  wi <- !duplicated(well_id)
  wells <- data.frame(well_id = well_id[wi], sample_id = ifelse(well_type[wi] == "Sample", sid[wi], NA),
                      well_type = well_type[wi], plate = as.character(df$PlateID[wi]),
                      position = NA_character_, stringsAsFactors = FALSE)
  wells <- wells[match(well_levels, wells$well_id), , drop = FALSE]
  wells$qc_fail <- wells$well_id %in% fail_wells
  matrix_type <- if ("SampleMatrixType" %in% colnames(df)) unique(df$SampleMatrixType[well_type == "Sample"]) else NA
  list(exprMat0 = m, wells = wells, featData = featData, lodMat = lod, olinkLOD = NULL,
       scale = "log2", panel = unique(df$Panel %||% NA), matrix_type = matrix_type, notes = notes)
}

build_dataset <- function(expr.raw, platform, sample_id_col, npx_type = "raw", drop_qc_fail = TRUE) {
  ds <- switch(platform,
    SomaScan = build_dataset_somascan(expr.raw, sample_id_col %||% "SampleId", drop_qc_fail),
    OLink    = build_dataset_olink(expr.raw, npx_type, sample_id_col %||% "SampleID", drop_qc_fail),
    NULISA   = build_dataset_nulisa(expr.raw, sample_id_col %||% "SampleName", drop_qc_fail),
    stop("Unknown platform: ", platform))
  ds$platform <- platform
  check_dataset(ds)
  ds
}

#' Assert the build contract (catches a broken platform branch immediately).
check_dataset <- function(ds) {
  stopifnot(
    "exprMat0 must be a numeric matrix" = is.matrix(ds$exprMat0) && is.numeric(ds$exprMat0),
    "wells must have one row per exprMat0 column, in order" =
      nrow(ds$wells) == ncol(ds$exprMat0) && identical(ds$wells$well_id, colnames(ds$exprMat0)),
    "featData must have one row per exprMat0 row, in order" =
      nrow(ds$featData) == nrow(ds$exprMat0) && identical(rownames(ds$featData), rownames(ds$exprMat0)),
    "feature IDs must be unique" = !anyDuplicated(rownames(ds$exprMat0)),
    "well IDs must be unique" = !anyDuplicated(colnames(ds$exprMat0)),
    "at least one Sample well" = any(ds$wells$well_type == "Sample")
  )
  if (!is.null(ds$lodMat)) stopifnot(identical(dimnames(ds$lodMat), dimnames(ds$exprMat0)))
  if (!is.null(ds$olinkLOD)) stopifnot(identical(dimnames(ds$olinkLOD), dimnames(ds$exprMat0)))
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Traits: matching, grouping, covariates
# -----------------------------------------------------------------------------

#' Match Sample wells to traits rows and define the pre-processing groups.
#' Every sample in the data that has no traits row is listed (audit B5/B6).
#' The pre-processing group is two-level: CTL (rows whose group_col equals
#' ctl_value) and Case (every other non-blank value), as in file 15
#' (conflicts row 14). Blank groups are removed or labelled Case (row 15).
match_samples_to_traits <- function(ds, traits.raw, group_col, ctl_value,
                                    blank_group = c("remove", "case"),
                                    drop_qc_fail = TRUE) {
  blank_group <- match.arg(blank_group)
  if (!group_col %in% colnames(traits.raw)) stop("Group column '", group_col, "' not in the traits file.")
  w <- ds$wells
  is_s <- w$well_type == "Sample"
  in_traits <- w$sample_id %in% rownames(traits.raw)
  report <- list(
    samples_not_in_traits = sort(w$sample_id[is_s & !in_traits]),
    traits_not_in_data = sort(setdiff(rownames(traits.raw), w$sample_id[is_s]))
  )
  grp_raw <- as.character(traits.raw[w$sample_id, group_col])
  grp_raw[!is_s | !in_traits] <- NA
  blank <- is_s & in_traits & (is.na(grp_raw) | trimws(grp_raw) == "")
  report$blank_group <- sort(w$sample_id[blank])
  pp_group <- ifelse(trimws(grp_raw) == as.character(ctl_value), "CTL", "Case")
  pp_group[blank] <- if (blank_group == "case") "Case" else NA
  pp_group[!is_s | !in_traits] <- NA
  failed <- is_s & w$qc_fail
  report$qc_fail_samples <- sort(w$sample_id[failed & !is.na(pp_group)])
  if (drop_qc_fail) pp_group[failed] <- NA
  if (sum(pp_group == "CTL", na.rm = TRUE) < 3)
    stop("Fewer than 3 control samples after matching: check that '", ctl_value,
         "' is the control value of column '", group_col, "'.")
  w$group_raw <- grp_raw
  w$pp_group <- pp_group
  ds$wells <- w
  ds$match_report <- report
  ds
}

#' Recode sex to numeric (male = 1, female = 2), as the user chooses in the
#' GUI. With no mapping given: numeric columns are used as they are; two-value
#' columns that look like M/F or Male/Female are mapped automatically
#' (file 15 behaviour); anything else stops with a message.
recode_sex <- function(x, male_value = "", female_value = "") {
  if (is.numeric(x) && (identical(male_value, "") || is.null(male_value))) return(list(sex = x, note = "Sex column is numeric; used as is."))
  xs <- trimws(as.character(x))
  if (identical(male_value, "") || is.null(male_value)) {
    vals <- unique(toupper(xs[!is.na(xs) & xs != ""]))
    if (all(vals %in% c("M", "F"))) { male_value <- "M"; female_value <- "F" }
    else if (all(vals %in% c("MALE", "FEMALE"))) { male_value <- "MALE"; female_value <- "FEMALE" }
    else stop("Cannot tell which sex values are male and female (", paste(vals, collapse = ", "),
              "). Choose them in the app.")
    xs <- toupper(xs)
  }
  out <- ifelse(xs == male_value, 1, ifelse(xs == female_value, 2, NA_real_))
  other <- sum(!is.na(xs) & xs != "" & is.na(out))
  list(sex = out, note = sprintf("Sex recoded: '%s' = 1, '%s' = 2; %d other non-blank values set to missing.",
                                 male_value, female_value, other))
}

# -----------------------------------------------------------------------------
# LOD, signal-to-noise, group medians
# -----------------------------------------------------------------------------

#' Per-feature reference statistics from the Buffer wells (SomaScan buffers,
#' Olink negative controls, NULISA NC wells), group medians and S:N, on the
#' scale each platform uses for its LOD (conflicts rows 5, 7):
#'   SomaScan: linear RFU; LOD = buffer median + 3 SD
#'   OLink:    LOD on log2 (NPX) scale = NC median + 3 SD; S:N on linear values
#'   NULISA:   log2 NPQ; vendor per-well LOD; S:N as a log2 difference (reported only)
compute_reference_stats <- function(ds) {
  m <- ds$exprMat0
  w <- ds$wells
  buf <- w$well_type == "Buffer"
  if (sum(buf) < 2 && ds$platform != "NULISA") stop("Fewer than 2 buffer / negative-control wells; cannot compute LOD.")
  lodScale <- if (ds$platform == "OLink") log2(m) else m
  out <- data.frame(row.names = rownames(m))
  out$Buffer    <- if (sum(buf)) row_median(m[, buf, drop = FALSE]) else NA_real_
  out$BufferLOD <- if (sum(buf)) row_median(lodScale[, buf, drop = FALSE]) else NA_real_
  out$BufferSD  <- if (sum(buf) > 1) row_sd(lodScale[, buf, drop = FALSE]) else NA_real_
  out$BufferMAD <- if (sum(buf) > 1) row_mad(lodScale[, buf, drop = FALSE]) else NA_real_
  out$LOD.3SD   <- out$BufferLOD + 3.0 * out$BufferSD
  out$LOD.MAD   <- out$BufferLOD + 4.9 * out$BufferMAD
  if (ds$platform == "NULISA") {
    s <- w$well_type == "Sample"
    out$LOD.vendor <- row_median(ds$lodMat[, s, drop = FALSE])
  }
  for (g in c("CTL", "Case")) {
    cols <- which(w$pp_group == g)
    out[[g]] <- if (length(cols)) row_median(m[, cols, drop = FALSE]) else NA_real_
    out[[paste0("StoN.", g)]] <- if (ds$scale == "log2") out[[g]] - out$Buffer else (out[[g]] - out$Buffer) / out$Buffer
  }
  out$StoN.max <- pmax(out$StoN.CTL, out$StoN.Case, na.rm = TRUE)
  if (ds$platform == "SomaScan") {
    # SomaLogic S:N groups, with the parenthesis fix from file 3 (audit files 1-2)
    samp <- which(!is.na(w$pp_group))
    n <- length(samp)
    out$S0 <- ds$featData$Type == "Hybridization Control Elution"
    out$S1 <- !out$S0 & rowSums(m[, samp, drop = FALSE] < 2 * out$LOD.3SD, na.rm = TRUE) >= 0.9 * n
    out$S3 <- !out$S0 & rowSums(m[, samp, drop = FALSE] > 80000, na.rm = TRUE) > 0.5 * n
    out$S2 <- !(out$S0 | out$S1 | out$S3)
  }
  out
}

#' Below-LOD mask for the study samples: TRUE = below LOD.
below_lod_mask <- function(ds, ref, samples, use_fixed_lod = FALSE, fixed_lod = NULL) {
  m <- ds$exprMat0[, samples, drop = FALSE]
  if (ds$platform == "SomaScan") {
    below <- m < ref$LOD.3SD
  } else if (ds$platform == "OLink") {
    if (use_fixed_lod) {
      if (is.null(fixed_lod)) stop("Fixed LOD requested but no fixed LOD values were loaded.")
      below <- log2(m) < fixed_lod[rownames(m)]
    } else {
      below <- log2(m) < ref$LOD.3SD
    }
  } else {
    lod <- ds$lodMat[, samples, drop = FALSE]
    below <- !is.na(lod) & m < lod      # LOD of 0 stays 0 (user decision)
  }
  below[is.na(below)] <- FALSE          # NA abundance is not counted as below LOD (conflicts row 4)
  below
}

#' Signal-to-noise for the study samples (SomaScan, OLink): (x - buffer)/buffer on linear values.
signal_to_noise <- function(ds, ref, samples) {
  m <- ds$exprMat0[, samples, drop = FALSE]
  (m - ref$Buffer) / ref$Buffer
}

#' Olink fixed LODs (optional, PC-normalised NPX only), via OlinkAnalyze.
load_olink_fixed_lod <- function(expr.raw, fixed_lod_path, featData) {
  if (!requireNamespace("OlinkAnalyze", quietly = TRUE)) stop("OlinkAnalyze is needed for the fixed-LOD option.")
  d <- OlinkAnalyze::olink_lod(expr.raw, lod_file_path = fixed_lod_path, lod_method = "FixedLOD")
  lodcol <- intersect(c("PCNormalizedLOD", "LOD"), colnames(d))[1]
  agg <- stats::aggregate(d[[lodcol]], by = list(OlinkID = d$OlinkID), FUN = stats::median, na.rm = TRUE)
  setNames(agg$x[match(featData$OlinkID, agg$OlinkID)], featData$label)
}

# -----------------------------------------------------------------------------
# Sample outliers
# -----------------------------------------------------------------------------

#' Network-connectivity outliers (z of connectivity below -sdout, with the
#' precedence fix, audit B1). Optionally adds the distance-from-centroid /
#' top-PC robust-z screen from files 13 and 16 (conflicts row 9).
detect_outliers <- function(mat, sdout = 3, distance_screen = FALSE, robust_crit = 6, top_pcs = 3) {
  normadj <- (0.5 + 0.5 * WGCNA::bicor(mat, use = "pairwise.complete.obs")^2)
  ku <- WGCNA::fundamentalNetworkConcepts(normadj)$Connectivity
  z.ku <- (ku - mean(ku)) / sqrt(stats::var(ku))
  out_conn <- z.ku < -sdout
  names(out_conn) <- colnames(mat)
  res <- data.frame(sample = colnames(mat), z.connectivity = z.ku, outlier.connectivity = out_conn,
                    stringsAsFactors = FALSE)
  out_all <- out_conn
  if (distance_screen) {
    X <- mat[stats::complete.cases(mat), , drop = FALSE]
    robust_z <- function(x) {
      s <- stats::mad(x, na.rm = TRUE); if (s == 0) s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) return(rep(0, length(x)))
      (x - stats::median(x, na.rm = TRUE)) / s
    }
    z.dist <- robust_z(sqrt(colSums((X - rowMeans(X))^2)))
    pca <- stats::prcomp(t(X), center = TRUE, scale. = FALSE)
    npc <- min(top_pcs, ncol(pca$x))
    z.pc.max <- apply(abs(apply(pca$x[, seq_len(npc), drop = FALSE], 2, robust_z)), 1, max)
    out_dist <- abs(z.dist) > robust_crit | z.pc.max > robust_crit
    res$z.distance <- z.dist
    res$z.pc.max <- z.pc.max
    res$outlier.distance <- out_dist
    out_all <- out_conn | out_dist
  }
  res$outlier <- out_all
  list(outliers = setNames(out_all, colnames(mat)), table = res)
}

# -----------------------------------------------------------------------------
# Bimodal features
# -----------------------------------------------------------------------------

check_bimodality <- function(x) {
  x <- stats::na.omit(x)
  if (length(x) < 10 || stats::sd(x) == 0) return(list(bimodal = FALSE, dip_p_value = NA_real_, gmm_components = NA_integer_))
  dip_p <- diptest::dip.test(x)$p.value
  # Mclust only works with mclust attached (it calls mclustBIC unqualified);
  # detect_bimodal() attaches it. A failed fit is returned as NA and counted in
  # the report, never silently treated as "not bimodal".
  gmm <- tryCatch(suppressWarnings(mclust::Mclust(x, G = 1:2, verbose = FALSE)), error = function(e) NULL)
  g <- if (is.null(gmm)) NA_integer_ else gmm$G
  list(bimodal = isTRUE(dip_p < 0.05 && g == 2), dip_p_value = dip_p, gmm_components = g)
}

detect_bimodal <- function(mat) {
  suppressPackageStartupMessages(library(mclust))
  res <- lapply(seq_len(nrow(mat)), function(i) check_bimodality(mat[i, ]))
  data.frame(Variable = rownames(mat),
             Bimodal = vapply(res, `[[`, logical(1), "bimodal"),
             DipPValue = vapply(res, function(r) as.numeric(r$dip_p_value), numeric(1)),
             GMMComponents = vapply(res, function(r) as.numeric(r$gmm_components), numeric(1)),
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Plate regression (optional, conflicts row 20)
# -----------------------------------------------------------------------------

#' Per-feature OLS of y ~ group + plate; subtracts the mean-centred plate term
#' only, so the CTL/Case difference is protected (file 16 logic).
regress_plate <- function(mat, plate, group) {
  plate <- factor(plate); group <- stats::relevel(factor(group), ref = "CTL")
  if (anyNA(plate)) stop("Plate regression: some samples have no plate.")
  if (nlevels(plate) < 2) stop("Plate regression: fewer than 2 plates.")
  X <- stats::model.matrix(~ group + plate)
  pc <- grep("^plate", colnames(X))
  Xp <- scale(X[, pc, drop = FALSE], center = TRUE, scale = FALSE)
  Y <- t(mat)
  ok <- stats::complete.cases(Y)
  if (!all(ok)) stop("Plate regression needs complete data; ", sum(!ok), " samples have missing values.")
  beta <- qr.coef(qr(X), Y)
  bp <- beta[pc, , drop = FALSE]; bp[is.na(bp)] <- 0
  out <- t(Y - Xp %*% bp)
  dimnames(out) <- dimnames(mat)
  out
}

# -----------------------------------------------------------------------------
# Covariates: PC1 and blood eigenprotein
# -----------------------------------------------------------------------------

#' PC1 from control samples, projected onto every sample with the control
#' centring and scaling (file 3 fix, conflicts row 10). Uses the matrix it is
#' given (fix for 13/16, which ignored the argument). Features with missing
#' values or zero variance in controls are left out of the PCA.
compute_pc1_scores <- function(mat, ctl_samples) {
  ctl <- mat[, ctl_samples, drop = FALSE]
  use <- stats::complete.cases(mat) & apply(ctl, 1, stats::sd) > 0
  if (sum(use) < 2) stop("PC1: fewer than 2 usable features.")
  pca <- stats::prcomp(t(ctl[use, , drop = FALSE]), center = TRUE, scale. = TRUE)
  load1 <- pca$rotation[, 1]
  all_s <- scale(t(mat[use, , drop = FALSE]), center = pca$center, scale = pca$scale)
  list(scores = setNames(as.vector(all_s %*% load1), colnames(mat)),
       n_features = sum(use),
       var_explained = summary(pca)$importance[2, 1])
}

BLOOD_MARKERS <- c("HBA1", "HBB", "HBG1", "FGA", "FGB", "FGG", "HP", "CAT", "CA1")

#' Blood-contamination eigenprotein (CSF only). Reports which of the 9 markers
#' were used and which were removed earlier as bimodal (audit B8).
blood_eigenprotein <- function(mat, featData, bimodal_features = character(0), markers = BLOOD_MARKERS) {
  sym_mat <- featData[rownames(mat), "symbol"]
  used_rows <- rownames(mat)[sym_mat %in% markers]
  sym_all_bimodal <- featData[intersect(bimodal_features, rownames(featData)), "symbol"]
  report <- list(markers_used = sort(unique(sym_mat[sym_mat %in% markers])),
                 rows_used = used_rows,
                 markers_removed_bimodal = sort(unique(intersect(markers, sym_all_bimodal))),
                 markers_absent = sort(setdiff(markers, c(sym_mat, sym_all_bimodal))))
  if (length(used_rows) < 2) return(list(ep = NULL, report = report))
  x <- t(mat[used_rows, , drop = FALSE])
  x <- x[, colSums(is.na(x)) == 0, drop = FALSE]
  if (ncol(x) < 2) return(list(ep = NULL, report = report))
  me <- WGCNA::moduleEigengenes(x, rep("red", ncol(x)), impute = FALSE)$eigengenes$MEred
  list(ep = setNames(me, colnames(mat)), report = report)
}

# -----------------------------------------------------------------------------
# Covariate regression (control-only bootstrap)
# -----------------------------------------------------------------------------

#' Fit y ~ covariates in control samples for every feature with an ordinary
#' nonparametric bootstrap (median coefficient over n_boot resamples), the
#' same estimator and seed as the notebooks. Returns a features x coefficients
#' matrix with NA coefficients set to 0.
fit_ctl_bootstrap <- function(mat, covars, ctl_samples, n_boot = 1000, seed = 8675309, n_cores = 1) {
  stopifnot(all(ctl_samples %in% rownames(covars)), all(ctl_samples %in% colnames(mat)))
  cv <- covars[ctl_samples, , drop = FALSE]
  Y <- mat[, ctl_samples, drop = FALSE]
  vars <- colnames(cv)
  form <- stats::as.formula(paste("thisexp ~", paste(vars, collapse = " + ")))
  bs <- function(data, indices, formula) stats::coef(stats::lm(formula, data = data[indices, , drop = FALSE]))
  one <- function(i) {
    set.seed(seed)
    d <- stats::na.omit(data.frame(thisexp = as.numeric(Y[i, ]), cv))
    if (nrow(d) < length(vars) + 2) return(rep(NA_real_, length(vars) + 1))
    b <- boot::boot(data = d, statistic = bs, R = n_boot, formula = form)
    apply(b$t, 2, stats::median, na.rm = TRUE)
  }
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)   # clusters are always stopped (conflicts row 27)
    parallel::clusterExport(cl, c("Y", "cv", "form", "bs", "seed", "n_boot", "vars"), envir = environment())
    res <- parallel::parLapply(cl, seq_len(nrow(Y)), one)
  } else {
    res <- lapply(seq_len(nrow(Y)), one)
  }
  coefmat <- do.call(rbind, res)
  colnames(coefmat) <- c("(Intercept)", vars)
  rownames(coefmat) <- rownames(mat)
  coefmat[is.na(coefmat)] <- 0
  coefmat
}

#' Centre every covariate over all samples; missing values become 0 (the
#' mean) with a note (conflicts rows 11, 12).
center_covariates <- function(covars_all) {
  notes <- character(0)
  out <- lapply(colnames(covars_all), function(v) {
    x <- covars_all[[v]] - mean(covars_all[[v]], na.rm = TRUE)
    if (anyNA(x)) { notes <<- c(notes, sprintf("'%s' missing for %d samples; set to the mean (0 after centring).", v, sum(is.na(x)))); x[is.na(x)] <- 0 }
    x
  })
  names(out) <- colnames(covars_all)
  list(centered = as.data.frame(out, row.names = rownames(covars_all)), notes = notes)
}

COVARIATE_LABELS <- c(bloodEP = "Blood", pc1score = "PC", age = "Age", sex = "Sex")

#' All regressed matrices, as a named list (audit B10). For CSF, Blood is
#' included in every combination (conflicts row 19); for plasma there is no
#' Blood covariate and every non-empty combination is produced.
generate_regressed <- function(mat, coefmat, centered, mandatory = character(0), base_name = "cleanDat") {
  optional <- setdiff(colnames(centered), mandatory)
  combos <- unlist(lapply(seq_along(optional), function(k) utils::combn(optional, k, simplify = FALSE)), recursive = FALSE)
  if (length(mandatory)) combos <- c(list(character(0)), combos)
  out <- list()
  for (cmb in combos) {
    vars <- c(mandatory, cmb)
    vars <- vars[order(match(vars, names(COVARIATE_LABELS)))]
    corr <- Reduce(`+`, lapply(vars, function(v) outer(coefmat[, v], centered[colnames(mat), v])))
    nm <- paste(c(base_name, COVARIATE_LABELS[vars], "Reg"), collapse = ".")
    out[[nm]] <- mat - corr
  }
  out
}

# -----------------------------------------------------------------------------
# Feature filters and winsorization
# -----------------------------------------------------------------------------

#' Keep features whose fraction of samples below LOD is at most `cutoff`
#' (75% SomaScan/OLink, 50% NULISA by default; conflicts row 1), counted over
#' the post-outlier samples (audit B2).
lod_filter <- function(below_mask, cutoff) {
  frac <- rowMeans(below_mask)
  list(keep = names(frac)[frac <= cutoff], frac_below = frac)
}

#' 2-pass: features passing the LOD filter AND with S:N below `sn_crit` in at
#' most `cutoff` of samples (SomaScan and OLink only; conflicts row 7).
sn_filter <- function(sn, cutoff, sn_crit = 0.3) {
  frac <- rowMeans(sn < sn_crit, na.rm = TRUE)
  list(keep = names(frac)[frac <= cutoff], frac_below = frac)
}

#' Scaled winsorization at mean +/- n_sd SD (5 SD by default; conflicts row 8).
#' Values beyond a bound are pulled toward it by factor k, preserving their
#' order (no clipping). NA-safe (audit B3).
scaled_winsorize <- function(x, k = 0.25, n_sd = 5) {
  mu <- mean(x, na.rm = TRUE); s <- stats::sd(x, na.rm = TRUE)
  lo <- mu - n_sd * s; hi <- mu + n_sd * s
  x <- ifelse(!is.na(x) & x < lo, lo + k * (x - lo), x)
  ifelse(!is.na(x) & x > hi, hi + k * (x - hi), x)
}

winsorize_matrix <- function(mat, k = 0.25, n_sd = 5) {
  out <- mat
  for (i in seq_len(nrow(mat))) out[i, ] <- scaled_winsorize(mat[i, ], k = k, n_sd = n_sd)
  out
}

# -----------------------------------------------------------------------------
# QC table
# -----------------------------------------------------------------------------

pass_fail <- function(fail) ifelse(is.na(fail), NA, ifelse(fail, "Fail", "Pass"))

#' One row per feature in the expression matrix, all flags as Pass/Fail
#' (conflicts row 25); group columns found by name, not typed in (audit, 16 line 2339).
build_master_qc <- function(ds, ref, bimodal_tbl, lod_res, sn_res, platform) {
  fd <- ds$featData
  q <- data.frame(Feature = rownames(fd), Symbol = fd$symbol, UniProt = fd$uniprot,
                  HumanProtein = fd$is_human_protein, stringsAsFactors = FALSE, check.names = FALSE)
  extra <- switch(platform, SomaScan = c("AptName", "SeqId", "SomaId", "TargetFullName", "Type", "Organism", "ColCheck"),
                  OLink = c("OlinkID", "Panel", "Block"), NULISA = c("ProteinName"))
  for (e in extra) q[[e]] <- fd[[e]]
  unit <- switch(platform, SomaScan = "RFU", OLink = "NPX", NULISA = "NPQ")
  tolog <- function(x) if (platform == "OLink") log2(x) else x
  q[[paste0("Median ", unit, " ", well_type_labels(platform)[["Buffer"]])]] <- tolog(ref$Buffer)
  q[[paste0("Median ", unit, " CTL")]]  <- tolog(ref$CTL)
  q[[paste0("Median ", unit, " Case")]] <- tolog(ref$Case)
  q$StoN.CTL <- ref$StoN.CTL; q$StoN.Case <- ref$StoN.Case; q$StoN.max <- ref$StoN.max
  q$LOD.3SD <- ref$LOD.3SD; q$LOD.MAD <- ref$LOD.MAD
  if (platform == "NULISA") q$`LOD vendor (median)` <- ref$LOD.vendor
  if (platform == "SomaScan") q$`SomaLogic S:N group` <- ifelse(ref$S0, "S0", ifelse(ref$S1, "S1", ifelse(ref$S3, "S3", "S2")))
  bi <- bimodal_tbl[match(q$Feature, bimodal_tbl$Variable), ]
  q$DipPValue <- bi$DipPValue; q$GMMComponents <- bi$GMMComponents
  q$`Bimodal Flag` <- pass_fail(bi$Bimodal)
  q$`%BelowLOD` <- lod_res$frac_below[q$Feature]
  q$`1-Pass Flag` <- pass_fail(!(q$Feature %in% lod_res$keep))
  if (!is.null(sn_res)) {
    q$`%SN<crit` <- sn_res$frac_below[q$Feature]
    q$`2-Pass Flag` <- pass_fail(!(q$Feature %in% intersect(lod_res$keep, sn_res$keep)))
  }
  q$`QC Status` <- pass_fail(q$`Bimodal Flag` %in% "Fail" | q$`1-Pass Flag` %in% "Fail")
  q
}

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

#' Per-run output folder: <output_root>/<project>_<platform>_<tissue>[_<panel>][_<npx>]
run_dir_path <- function(output_root, project, platform, tissue, panel = NULL, npx_type = NULL) {
  clean <- function(x) gsub("[^A-Za-z0-9.-]+", "", x)
  parts <- c(clean(project), platform, tissue)
  if (!is.null(panel) && length(panel) && !all(is.na(panel))) parts <- c(parts, clean(paste(panel[!is.na(panel)], collapse = "-")))
  if (platform == "OLink" && !is.null(npx_type)) parts <- c(parts, paste0(npx_type, "NPX"))
  file.path(path.expand(output_root), paste(parts, collapse = "_"))
}

make_run_dir <- function(output_root, project, platform, tissue, panel = NULL, npx_type = NULL) {
  d <- run_dir_path(output_root, project, platform, tissue, panel, npx_type)
  for (s in c("-cleanDat", "-traits", "-pre-proc", "QC")) dir.create(file.path(d, s), recursive = TRUE, showWarnings = FALSE)
  d
}

#' Write a named list of matrices to <run_dir>/-cleanDat, human proteins only
#' (conflicts row 21).
write_cleanDats <- function(mat_list, run_dir, featData) {
  for (nm in names(mat_list)) {
    m <- mat_list[[nm]]
    keep <- featData[rownames(m), "is_human_protein"] %in% TRUE
    utils::write.csv(m[keep, , drop = FALSE], file.path(run_dir, "-cleanDat", paste0(nm, ".csv")))
  }
  invisible(names(mat_list))
}

#' Draw a QC figure in the report and save the same figure to QC/<name>.pdf
#' (replaces dev.new(), conflicts row 27).
qc_figure <- function(name, draw, qc_dir, width = 8, height = 7) {
  grDevices::pdf(file.path(qc_dir, paste0(name, ".pdf")), width = width, height = height)
  tryCatch(draw(), finally = grDevices::dev.off())
  draw()
  invisible(NULL)
}

#' MDS plot coloured by well type / group. Controls are subsampled to at most
#' `max_per_type` wells per type with a fixed seed (16 lines 1274-1282).
plot_mds_wells <- function(ds, title, include_controls = TRUE, samples = NULL, max_per_type = 20, seed = 1) {
  w <- ds$wells
  lab <- well_type_labels(ds$platform)
  keep <- if (!is.null(samples)) w$well_id %in% samples else (!is.na(w$pp_group) | (include_controls & w$well_type %in% c("Buffer", "Calibrator", "QC")))
  set.seed(seed)
  for (t in c("Buffer", "Calibrator", "QC")) {
    idx <- which(keep & w$well_type == t)
    if (length(idx) > max_per_type) keep[setdiff(idx, sample(idx, max_per_type))] <- FALSE
  }
  m <- ds$exprMat0[, keep, drop = FALSE]
  if (ds$scale == "linear") m <- log2(m)
  m <- m[stats::complete.cases(m) & apply(m, 1, stats::sd) > 0, , drop = FALSE]
  ww <- w[keep, ]
  cat_ <- ifelse(!is.na(ww$pp_group), ww$pp_group, lab[ww$well_type])
  pal <- c(CTL = "chartreuse3", Case = "goldenrod2")
  pal[lab[["Buffer"]]] <- "grey50"; pal[lab[["Calibrator"]]] <- "turquoise3"; pal[lab[["QC"]]] <- "purple"
  cat_[is.na(cat_)] <- "Other"; pal["Other"] <- "black"
  limma::plotMDS(m, col = pal[cat_], pch = 16, main = title)
  present <- intersect(names(pal), unique(cat_))
  graphics::legend("topright", legend = present, col = pal[present], pch = 16, cex = 0.8, bg = "white")
}

# =============================================================================
# Harmonized / pre-processed input (templates/harmonized_report.Rmd)
# -----------------------------------------------------------------------------
# Input: a cleanDat CSV (features x samples, first column = feature IDs such as
# "Target|UniProt") that has already been LOD-filtered and harmonized (e.g. the
# TAMPOR output), plus a numericMeta CSV whose ID column matches the cleanDat
# column names. Platform-specific steps (wells, LOD, S:N, 1-/2-pass) were done
# upstream and are skipped; everything from sample outliers onward is shared
# with the single-cohort pipeline.
# =============================================================================

#' Read a harmonized cleanDat CSV into a numeric log2 matrix.
read_harmonized_matrix <- function(path, scale = c("log2", "linear")) {
  scale <- match.arg(scale)
  if (!file.exists(path)) stop("Harmonized matrix not found: ", path)
  d <- fix_blank_colnames(utils::read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE))
  feats <- trimws(as.character(d[[1]]))
  if (anyNA(feats) || any(feats == "")) stop("The first column (feature IDs) has blanks.")
  if (anyDuplicated(feats)) stop("Duplicated feature IDs: ", paste(head(unique(feats[duplicated(feats)]), 10), collapse = ", "))
  raw <- d[, -1, drop = FALSE]
  m <- suppressWarnings(vapply(raw, as.numeric, numeric(nrow(raw))))
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(raw))
  bad <- sum(is.na(m) & !(is.na(as.matrix(raw)) | as.matrix(raw) == ""))
  if (bad) stop(bad, " non-numeric values in the matrix (is the first column the feature ID and every other column a sample?).")
  dimnames(m) <- list(feats, trimws(colnames(raw)))
  if (anyDuplicated(colnames(m))) stop("Duplicated sample columns: ", paste(head(unique(colnames(m)[duplicated(colnames(m))]), 10), collapse = ", "))
  notes <- character(0)
  if (scale == "linear") {
    nonpos <- sum(m <= 0, na.rm = TRUE)
    if (nonpos) notes <- c(notes, sprintf("%d values <= 0 set to missing before log2.", nonpos))
    m[m <= 0] <- NA
    m <- log2(m)
  } else if (stats::median(m, na.rm = TRUE) > 100) {
    warning("The matrix median is ", round(stats::median(m, na.rm = TRUE)), "; it looks linear, not log2. Check the 'Matrix scale' setting.")
  }
  list(mat = m, notes = notes)
}

#' Feature annotation from IDs like "Target|UniProt" (NULISA harmonization)
#' or "Symbol|UniProt[^...]" (SomaScan). Inputs are already human-only.
harmonized_featdata <- function(labels) {
  parts <- strsplit(labels, "|", fixed = TRUE)
  symbol <- vapply(parts, function(p) p[1], "")
  uniprot <- vapply(parts, function(p) if (length(p) >= 2) sub("\\^.*$", "", p[2]) else NA_character_, "")
  uniprot[grepl("^noUniProt", uniprot)] <- NA
  data.frame(label = labels, symbol = symbol, uniprot = uniprot, is_human_protein = TRUE,
             row.names = labels, stringsAsFactors = FALSE)
}

#' Match matrix columns to numericMeta rows and define CTL/Case.
match_harmonized <- function(mat, traits.raw, group_col, ctl_value, blank_group = c("remove", "case")) {
  blank_group <- match.arg(blank_group)
  if (!group_col %in% colnames(traits.raw)) stop("Group column '", group_col, "' not in the numericMeta file.")
  ids <- colnames(mat)
  report <- list(samples_not_in_traits = sort(setdiff(ids, rownames(traits.raw))),
                 traits_not_in_data = sort(setdiff(rownames(traits.raw), ids)))
  ids <- ids[ids %in% rownames(traits.raw)]
  g <- trimws(as.character(traits.raw[ids, group_col]))
  blank <- is.na(g) | g == ""
  report$blank_group <- sort(ids[blank])
  pp <- ifelse(g == as.character(ctl_value), "CTL", "Case")
  pp[blank] <- if (blank_group == "case") "Case" else NA
  keep <- ids[!is.na(pp)]
  meta <- traits.raw[keep, , drop = FALSE]
  meta$PreProcGroup <- pp[!is.na(pp)]
  if (sum(meta$PreProcGroup == "CTL") < 3)
    stop("Fewer than 3 control samples after matching: check that '", ctl_value, "' is the control value of '", group_col, "'.")
  list(mat = mat[, keep, drop = FALSE], meta = meta, report = report)
}

#' Samples that share a pair ID (e.g. the GIS bridge samples, measured once
#' per batch). Always reports the pairs and whether their group and sex agree;
#' mode decides what is kept:
#'   keep_all  every measurement stays (the old post-harmonization notebook)
#'   keep_one  one measurement per pair, from batch `keep_value` when present
#'   average   the pair's values are averaged into the `keep_value` sample's column
#' Only pair IDs that occur more than once are touched; singletons are unchanged.
resolve_duplicates <- function(mat, meta, pair_col = "", mode = c("keep_all", "keep_one", "average"),
                               batch_col = "", keep_value = "", sex_col = "") {
  mode <- match.arg(mode)
  if (!nzchar(pair_col %||% "")) return(list(mat = mat, meta = meta, pairs = NULL, notes = character(0)))
  if (!pair_col %in% colnames(meta)) stop("Pair column '", pair_col, "' not in the numericMeta file.")
  pid <- trimws(as.character(meta[[pair_col]]))
  dups <- sort(unique(pid[!is.na(pid) & pid != "" & duplicated(pid)]))
  if (!length(dups)) return(list(mat = mat, meta = meta, pairs = NULL, notes = sprintf("No repeated values in '%s'.", pair_col)))
  batch <- if (nzchar(batch_col %||% "") && batch_col %in% colnames(meta)) as.character(meta[[batch_col]]) else rep(NA_character_, nrow(meta))
  sexv <- if (nzchar(sex_col %||% "") && sex_col %in% colnames(meta)) as.character(meta[[sex_col]]) else rep(NA_character_, nrow(meta))
  rows <- lapply(dups, function(p) {
    i <- which(pid == p)
    pref <- i[batch[i] %in% keep_value][1]
    if (is.na(pref)) pref <- i[1]
    data.frame(pair = p, n = length(i), samples = paste(rownames(meta)[i], collapse = "; "),
               batches = paste(batch[i], collapse = "; "), groups = paste(meta$PreProcGroup[i], collapse = "; "),
               group_agrees = length(unique(meta$PreProcGroup[i])) == 1,
               sex_agrees = length(unique(sexv[i][!is.na(sexv[i])])) <= 1,
               kept = rownames(meta)[pref], stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, rows)
  notes <- character(0)
  if (any(!pairs$group_agrees)) notes <- c(notes, sprintf("%d repeated samples have a different group in different rows: %s",
                                                          sum(!pairs$group_agrees), paste(head(pairs$pair[!pairs$group_agrees], 20), collapse = ", ")))
  if (any(!pairs$sex_agrees)) notes <- c(notes, sprintf("%d repeated samples have a different sex in different rows.", sum(!pairs$sex_agrees)))
  if (mode == "keep_all") return(list(mat = mat, meta = meta, pairs = pairs, notes = notes))
  drop <- character(0)
  meta$n_measurements_averaged <- 1
  for (k in seq_len(nrow(pairs))) {
    members <- rownames(meta)[pid == pairs$pair[k]]
    kept <- pairs$kept[k]
    if (mode == "average") {
      v <- rowMeans(mat[, members, drop = FALSE], na.rm = TRUE)
      v[is.nan(v)] <- NA
      mat[, kept] <- v
      meta[kept, "n_measurements_averaged"] <- length(members)
    }
    drop <- c(drop, setdiff(members, kept))
  }
  if (mode == "keep_one") meta$n_measurements_averaged <- NULL
  keep <- setdiff(colnames(mat), drop)
  list(mat = mat[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE], pairs = pairs, notes = notes)
}

#' Fraction of samples missing per feature; keep features at or under max_missing.
missing_filter <- function(mat, max_missing = 0.5) {
  frac <- rowMeans(is.na(mat))
  list(keep = names(frac)[frac <= max_missing], frac_missing = frac)
}

#' MDS of samples coloured by a metadata column.
plot_mds_samples <- function(mat, meta, color_col, title) {
  m <- mat[stats::complete.cases(mat), , drop = FALSE]
  m <- m[apply(m, 1, stats::sd) > 0, , drop = FALSE]
  f <- if (nzchar(color_col %||% "") && color_col %in% colnames(meta)) factor(as.character(meta[colnames(m), color_col])) else factor(rep("all", ncol(m)))
  pal <- if (identical(levels(f), c("Case", "CTL"))) c(Case = "goldenrod2", CTL = "chartreuse3") else
    setNames(grDevices::hcl.colors(max(3, nlevels(f)), "Dark 3")[seq_len(nlevels(f))], levels(f))
  cols <- pal[as.character(f)]; cols[is.na(cols)] <- "black"
  limma::plotMDS(m, col = cols, pch = 16, main = title)
  graphics::legend("topright", legend = names(pal), col = pal, pch = 16, cex = 0.8, bg = "white", title = color_col)
}

#' One row per input feature: annotation, medians, missingness, bimodal test.
build_master_qc_harmonized <- function(featData, mat, meta, miss_res, bimodal_tbl = NULL) {
  q <- data.frame(Feature = rownames(featData), Symbol = featData$symbol, UniProt = featData$uniprot,
                  stringsAsFactors = FALSE, check.names = FALSE)
  for (g in c("CTL", "Case")) {
    cols <- rownames(meta)[meta$PreProcGroup == g]
    q[[paste("Median", g)]] <- if (length(cols)) row_median(mat[q$Feature, cols, drop = FALSE]) else NA_real_
  }
  q$FractionMissing <- miss_res$frac_missing[q$Feature]
  q$`Missingness Flag` <- pass_fail(!(q$Feature %in% miss_res$keep))
  if (!is.null(bimodal_tbl)) {
    bi <- bimodal_tbl[match(q$Feature, bimodal_tbl$Variable), ]
    q$DipPValue <- bi$DipPValue; q$GMMComponents <- bi$GMMComponents
    q$`Bimodal Flag` <- pass_fail(bi$Bimodal)
  }
  q$`QC Status` <- pass_fail(q$`Missingness Flag` %in% "Fail" | (q$`Bimodal Flag` %||% FALSE) %in% "Fail")
  q
}
