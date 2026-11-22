###############################################################################
# make_test_data.R: synthetic inputs for testing the app on every platform.
# No real data: SomaScan comes from SomaDataIO's public example adat; OLink and
# NULISA are simulated in the layouts of the real exports. Traits are random.
# Usage: Rscript tests/make_test_data.R <output_dir>
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "test_data"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(42)
suppressPackageStartupMessages({ library(SomaDataIO); library(arrow) })

make_traits <- function(ids, n_extra = 2, id_name = "ID", sex_style = c("MF", "MaleFemale", "num")) {
  sex_style <- match.arg(sex_style)
  n <- length(ids)
  grp <- sample(c("CTL", "FTD", "AD", "PSP"), n, replace = TRUE, prob = c(.35, .3, .2, .15))
  grp[sample(n, 2)] <- ""                                  # blank groups
  sex <- sample(c("M", "F"), n, replace = TRUE)
  sex <- switch(sex_style, MF = sex, MaleFemale = ifelse(sex == "M", "Male", "Female"), num = ifelse(sex == "M", 1, 0))
  age <- round(runif(n, 45, 85), 1); age[sample(n, 1)] <- NA
  tr <- data.frame(ids, Dx = grp, AGE = age, SEX = sex, Site = sample(c("A", "B"), n, TRUE), stringsAsFactors = FALSE)
  names(tr)[1] <- id_name
  extra <- data.frame(paste0("NOTINDATA", seq_len(n_extra)), "CTL", 60, sex[1], "A")
  names(extra) <- names(tr)
  rbind(tr, extra)
}

## SomaScan: public example adat (5k plasma, 34 metadata columns) -------------
d <- SomaDataIO::example_data
write_adat(d, file.path(out, "soma_example.adat"))
sids <- d$SampleId[d$SampleType == "Sample"]
tr <- make_traits(sids[-1], id_name = "SampleId")         # first sample has no traits row
write.csv(tr, file.path(out, "soma_traits.csv"), row.names = FALSE)

## OLink Explore HT, two export layouts ----------------------------------------
sim_olink <- function(layout = c("ALLFTD", "MACWide"), n_assay = 300, n_plate = 3, n_samp = 60) {
  layout <- match.arg(layout)
  syms <- c("HBA1", "HBB", "HBG1", "FGA", "FGB", "FGG", "HP", "CAT", "CA1", paste0("PROT", seq_len(n_assay - 9)))
  assays <- data.frame(OlinkID = sprintf("OID%05d", seq_len(n_assay)), Assay = syms,
                       UniProt = sprintf("P%05d", seq_len(n_assay)), AssayType = "assay",
                       Block = as.character(rep(1:8, length.out = n_assay)), stringsAsFactors = FALSE)
  ctrl <- data.frame(OlinkID = sprintf("OIDC%03d", 1:24), Assay = paste(rep(c("Extension control", "Amplification control", "Incubation control"), each = 8), rep(1:8, 3)),
                     UniProt = rep(c("EXT1", "AMP1", "INC1"), each = 8), AssayType = rep(c("ext_ctrl", "amp_ctrl", "inc_ctrl"), each = 8),
                     Block = as.character(rep(1:8, 3)), stringsAsFactors = FALSE)
  feats <- rbind(assays, ctrl)
  base <- rnorm(nrow(feats), 4, 1.5); bimodal <- sample(n_assay, 5)
  wells <- list(); k <- 0
  for (pl in seq_len(n_plate)) {
    pid <- if (layout == "ALLFTD") sprintf("Run842_PL%d", pl) else sprintf("r34p%d.plate.layout", pl)
    sids <- if (layout == "ALLFTD") sprintf("A%07d", (pl - 1) * n_samp + seq_len(n_samp)) else as.character(1000000 + (pl - 1) * n_samp + seq_len(n_samp))
    if (layout == "ALLFTD") sids[1:2] <- paste0("Ref-", pl, "-", 1:2)
    ty <- c(rep("SAMPLE", n_samp), rep("NEGATIVE_CONTROL", 4), rep("PLATE_CONTROL", 5), rep("SAMPLE_CONTROL", 3))
    nm <- c(sids, if (layout == "ALLFTD") c(paste0("NEG", 1:4, "_", pid), paste0("PC", 1:5, "_", pid), paste0("SC", 1:3, "_", pid)) else c(paste0("NC", 1:4), paste0("PC", 1:5), paste0("SC", 1:3)))
    wells[[pl]] <- data.frame(SampleID = nm, SampleType = ty, PlateID = pid, WellID = sprintf("W%02d", seq_along(nm)), stringsAsFactors = FALSE)
  }
  wells <- do.call(rbind, wells)
  long <- merge(wells, feats, by = NULL)
  grp_eff <- ifelse(long$SampleType == "SAMPLE", rnorm(nrow(long), 0, 0.6), 0)
  npx <- base[match(long$OlinkID, feats$OlinkID)] + grp_eff
  npx[long$SampleType == "NEGATIVE_CONTROL"] <- base[match(long$OlinkID, feats$OlinkID)][long$SampleType == "NEGATIVE_CONTROL"] - 3 + rnorm(sum(long$SampleType == "NEGATIVE_CONTROL"), 0, 0.3)
  bm <- long$OlinkID %in% feats$OlinkID[bimodal] & long$SampleType == "SAMPLE" & runif(nrow(long)) < 0.5
  npx[bm] <- npx[bm] + 3
  npx[long$AssayType != "assay"] <- 0
  long$Panel <- "Explore_HT"; long$Count <- as.integer(round(2^npx * 100))
  long$NPX <- npx; long$PCNormalizedNPX <- npx + 0.1
  long$SampleQC <- "PASS"; long$SampleQC[long$SampleID == long$SampleID[long$SampleType == "SAMPLE"][5] & long$Block == "3"] <- "FAIL"
  long$AssayQC <- "NA"; long$DataAnalysisRefID <- "D10002"; long$Normalization <- "Intensity"; long$ExtNPX <- 0
  if (layout == "ALLFTD") {
    long$LOD <- base[match(long$OlinkID, feats$OlinkID)] - 2; long$PCNormalizedLOD <- long$LOD + 0.1
    long$ExploreVersion <- "7.4.2"; long$SubjectID <- ifelse(long$SampleType == "SAMPLE", long$SampleID, NA)
    long$InterCV <- NA_real_; long$IntraCV <- NA_real_
  } else {
    long$SoftwareVersion <- "1.3.0"; long$SoftwareName <- "NPX Map"; long$IntraCV <- NaN; long$InterCV <- NaN
    long$SampleBlockQCWarn <- 1L; long$SampleBlockQCFail <- 1L; long$BlockQCFail <- 1L; long$AssayQCWarn <- 0L
  }
  long[sample(nrow(long)), ]
}
for (lay in c("ALLFTD", "MACWide")) {
  ol <- sim_olink(lay)
  write_parquet(ol, file.path(out, paste0("olink_", lay, ".parquet")))
  sids <- unique(ol$SampleID[ol$SampleType == "SAMPLE" & !grepl("^Ref", ol$SampleID)])
  write.csv(make_traits(sids[-1], sex_style = if (lay == "ALLFTD") "MF" else "num"),
            file.path(out, paste0("olink_", lay, "_traits.csv")), row.names = FALSE)
}

## NULISA, two CSV variants ----------------------------------------------------
sim_nulisa <- function(variant = c("ALLFTD", "MACPath"), n_target = 120, n_plate = 2, n_samp = 70) {
  variant <- match.arg(variant)
  tg <- data.frame(Target = c("HBB", "FGA", paste0("TGT", 1:(n_target - 2))), UniProtID = sprintf("Q%05d", 1:n_target),
                   ProteinName = paste("Protein", 1:n_target), stringsAsFactors = FALSE)
  base <- rnorm(n_target, 12, 2)
  rows <- list()
  for (pl in seq_len(n_plate)) {
    pid <- sprintf("P%s_Plate%d_CNS.xml", if (variant == "ALLFTD") "029" else "078", pl)
    sn <- if (variant == "ALLFTD") paste0("CSF", (pl - 1) * n_samp + seq_len(n_samp)) else as.character(1900000 + (pl - 1) * n_samp + seq_len(n_samp))
    wells <- data.frame(SampleName = c(sn, rep("IPC", 3), rep("NC", 4), rep("SC", 3)),
                        SampleType = c(rep("Sample", n_samp), rep("IPC", 3), rep("NC", 4), rep("SC", 3)), stringsAsFactors = FALSE)
    lod <- base - 4 + rnorm(n_target, 0, .5); lod[1] <- if (pl == 2) 0 else lod[1]   # a LOD of 0, as in the ALLFTD sheet
    for (t in seq_len(n_target)) {
      npq <- base[t] + rnorm(nrow(wells), 0, 1); npq[wells$SampleType == "NC"] <- lod[t] - 1 + rnorm(4, 0, .2)
      if (t > n_target - 15) npq[wells$SampleType == "Sample"] <- lod[t] - 0.5 + rnorm(n_samp, 0, 0.6)   # mostly below LOD
      rows[[length(rows) + 1]] <- data.frame(Panel = "CNS Disease Panel", PlateID = pid, wells,
        SampleMatrixType = if (variant == "ALLFTD") "CSF" else "PLASMA", Target = tg$Target[t], UniProtID = tg$UniProtID[t],
        ProteinName = tg$ProteinName[t], SampleQC = if (variant == "ALLFTD") "passed" else "Passed",
        targetLOD_NPQ = lod[t], NPQ = npq, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  df$SampleQC[df$SampleName == df$SampleName[1]] <- "failed"
  df[, c("Panel", "PlateID", "SampleName", "SampleType", "SampleMatrixType", "Target", "UniProtID", "ProteinName", "SampleQC", "targetLOD_NPQ", "NPQ")]
}
for (v in c("ALLFTD", "MACPath")) {
  nu <- sim_nulisa(v)
  write.csv(nu, file.path(out, paste0("nulisa_", v, ".csv")), row.names = FALSE)
  sids <- unique(nu$SampleName[nu$SampleType == "Sample"])
  write.csv(make_traits(sids[-1], sex_style = if (v == "ALLFTD") "MaleFemale" else "MF"),
            file.path(out, paste0("nulisa_", v, "_traits.csv")), row.names = FALSE)
}
## Harmonized NULISA plasma: the layout harmonize.ALLFTD.Plasma.NULISA.Rmd writes --
# cleanDat.NULISA.Plasma.combined.TAMPOR.csv: "Target|UniProt" x "yBATCH.DID_cycle", log2 ratio near 0
# numericMeta.combined_preProc.csv: row names = NewName, blank first header
sim_harmonized <- function() {
  n24 <- 70; n25_new <- 50; n_gis <- 30
  did24 <- sprintf("DID%04d", seq_len(n24))
  gis <- did24[seq_len(n_gis)]
  did25 <- c(gis, sprintf("DID%04d", 5000 + seq_len(n25_new)))
  k24 <- data.frame(DID = did24, cycle = "1", Batch = "y2024")
  k25 <- data.frame(DID = did25, cycle = "1", Batch = "y2025")
  k25 <- rbind(k25, data.frame(DID = sample(did25, 15), cycle = "2", Batch = "y2025"),
               data.frame(DID = sample(did25, 4), cycle = "Collection", Batch = "y2025"))
  key <- rbind(k24, k25)
  key$NewName <- paste0(key$Batch, ".", key$DID, "_", key$cycle)
  key$PairID <- paste0(key$DID, "_", key$cycle)
  key$GIS <- ifelse(key$DID %in% gis & key$cycle == "1", "GIS", NA)
  person <- unique(key$DID)
  pgrp <- setNames(sample(c("CTL", "Case"), length(person), TRUE, prob = c(.4, .6)), person)
  psex <- setNames(sample(c("M", "F"), length(person), TRUE), person)
  page <- setNames(round(runif(length(person), 40, 85), 1), person)
  key$Group <- pgrp[key$DID]; key$Sex <- psex[key$DID]
  key$age_at_sample <- page[key$DID] + ifelse(key$Batch == "y2025", 1.2, 0) + suppressWarnings(as.numeric(key$cycle) - 1)
  key$age_at_sample[is.na(key$age_at_sample)] <- page[key$DID[is.na(key$age_at_sample)]]
  flip <- which(key$Batch == "y2025" & key$GIS %in% "GIS")[1:4]            # group mismatches, like the real data
  key$Group[flip] <- ifelse(key$Group[flip] == "CTL", "Case", "CTL")
  key$Group[sample(nrow(key), 2)] <- ""                                     # blank groups
  key$SampleName <- paste0("NULISA_", seq_len(nrow(key)))
  key$Repeats <- ifelse(key$GIS %in% "GIS", 2, NA)
  targets <- c("NEFL", "GFAP", "APOE4", "pTau217", "pTau181", "Aβ40", "Aβ42", "UCHL1", "CHI3L1", "IL6",
               sprintf("TGT%02d", 1:70))
  uni <- c("P07196", "P14136", "P02649", "P10636", "P10636", "P05067", "P05067", "P09936", "P36222", "P05231",
           sprintf("Q%05d", 1:70))
  feats <- paste(targets, uni, sep = "|")
  ns <- nrow(key)
  person_eff <- setNames(rnorm(length(person)), person)
  m <- sapply(seq_len(ns), function(j) {
    d <- key$DID[j]
    0.4 * person_eff[d] * runif(length(feats), 0.3, 1) + 0.5 * (key$Group[j] == "Case") * (seq_along(feats) %% 5 == 0) +
      0.02 * (key$age_at_sample[j] - 60) + 0.3 * (key$Sex[j] == "M") * (seq_along(feats) %% 7 == 0) + rnorm(length(feats), 0, 0.35)
  })
  carriers <- setNames(runif(length(person)) < 0.3, person)
  m[3, ] <- ifelse(carriers[key$DID], 3, 0) + rnorm(ns, 0, 0.2)            # APOE4: genotype-driven, bimodal
  m[, ns] <- m[, ns] * 0.1 + rnorm(length(feats), 0, 3)                     # one outlier sample
  m[sample(length(m), 40)] <- NA
  m[80, sample(ns, round(0.6 * ns))] <- NA                                  # one mostly-missing feature
  dimnames(m) <- list(feats, key$NewName)
  meta <- key[, c("SampleName", "DID", "cycle", "Batch", "age_at_sample", "Sex", "Group", "NewName", "GIS", "PairID", "Repeats")]
  rownames(meta) <- meta$NewName
  list(mat = m, meta = meta)
}
h <- sim_harmonized()
write.csv(h$mat, file.path(out, "cleanDat.NULISA.Plasma.combined.TAMPOR.csv"))
write.csv(h$meta, file.path(out, "numericMeta.combined_preProc.csv"))
cat("Test data written to", normalizePath(out), "\n")
