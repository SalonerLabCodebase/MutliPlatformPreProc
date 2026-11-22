# Proteomics Pre-processing App

A Shiny app for pre-processing fluid proteomics data from **SomaScan**, **OLink Explore HT** and **NULISA**, for **CSF** or **plasma**. You choose the run settings in the app. The app checks your inputs and then renders an R Markdown report that does all of the processing. Each run writes an HTML report, cleaned expression matrices (`cleanDat`), sample metadata (`numericMeta`), QC tables and QC figures to its own output folder.

The app has two input modes:

| Mode | Input | Use it for |
|---|---|---|
| **Single cohort** | A raw platform file (`.adat`, OLink `.parquet` or NULISA `.csv`) plus a traits `.csv` | Processing one cohort from the raw platform export |
| **Harmonized matrix** | A harmonized `cleanDat` `.csv` (features x samples) plus a `numericMeta` `.csv` | Processing a matrix that was already LOD-filtered and harmonized across batches or cohorts (e.g. TAMPOR output) |

> **Does your dataset have quirks?** For example, sample IDs in a non-standard column, known bad wells, duplicated measurements, or traits IDs that don't match the expression file. Don't edit the app, the functions or the report templates. Write a small **dataset-fixes script** instead. See [Dataset-fixes scripts](#dataset-fixes-scripts) below.

---

## Contents

- [Installation](#installation)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Dataset-fixes scripts](#dataset-fixes-scripts)
- [Input files](#input-files)
- [What the pipeline does](#what-the-pipeline-does)
- [Settings reference](#settings-reference)
- [Outputs](#outputs)
- [Troubleshooting](#troubleshooting)
- [Data privacy](#data-privacy)

---

## Installation

1. Clone or download this repository:
   ```bash
   git clone <REPO_URL>
   ```
2. Open R (version 4.0 or later) or RStudio.
3. Start the app (see [Quick start](#quick-start)). **The first time it starts, the app installs any missing packages**, from CRAN and Bioconductor:

   - CRAN: `shiny`, `rmarkdown`, `knitr`, `dplyr`, `tidyr`, `arrow`, `readxl`, `diptest`, `mclust`, `boot`, `doParallel`, `foreach`, `ggplot2`, `yaml`, `SomaDataIO`, `rstudioapi`, `WGCNA`
   - Bioconductor: `limma`, `impute`, `preprocessCore`, `GO.db`, `AnnotationDbi` (and `BiocManager`, if you don't have it)

   This can take several minutes. To turn off the automatic install, for example on a managed server, set an environment variable before you launch the app:
   ```r
   Sys.setenv(PREPROC_SKIP_INSTALL = "1")
   ```
4. **Optional:** the OLink **fixed LOD** option needs `OlinkAnalyze`, which is **not** installed automatically:
   ```r
   install.packages("OlinkAnalyze")
   ```

---

## Repository layout

```
proteomicsPreProcApp/
├── app.R                          # The Shiny app (UI + input checks + launches the report)
├── R/
│   └── preproc_functions.R        # All processing functions, used by both reports
├── templates/
│   ├── preproc_report.Rmd         # Report for single-cohort (raw platform file) runs
│   └── harmonized_report.Rmd      # Report for harmonized cleanDat + numericMeta runs
├── dataset_fixes_TEMPLATE.R       # Starting point for your own dataset-fixes script
└── README.md
```

The app finds `R/` and `templates/` relative to its own folder, so keep these files together.

---

## Quick start

1. **Launch the app.** In RStudio, open `app.R` and click **Run App**. Or, from any R session, run:
   ```r
   shiny::runApp("path/to/proteomicsPreProcApp")
   ```
2. **Section 1: Data**
   - Choose the input mode, platform and tissue, and enter a project name. The project name is used in the output folder name.
   - Enter the expression file, the traits file and an output folder. In RStudio, the **File...** and **Folder...** buttons open a file browser. Outside RStudio, paste the full paths.
   - *Single cohort only:* check the **Sample ID column** in the expression file. It defaults to `SampleId` for SomaScan, `SampleID` for OLink and `SampleName` for NULISA.
   - Click **Load traits and check sample matching**.
3. **Read the Input check panel** on the right. It shows:
   - how many sample wells (or, in harmonized mode, sample columns) were found
   - expression samples that have no traits row (the first 15 IDs are listed)
   - traits rows that have no sample in the expression file
   - how many samples have a blank group value
   - a table of matched samples by group (CTL / Case / blank)

   This check reads the files **before** any dataset-fixes script runs (see the note under [Dataset-fixes scripts](#dataset-fixes-scripts)).
4. **Section 2: Groups and covariates.** Pick the traits ID column, the diagnostic grouping column and the value that means *control*. Every other non-blank value becomes **Case**. You can also choose age and sex columns for covariate regression.
5. **Section 3: Options.** Adjust the filters, outlier settings and bootstrap settings (see [Settings reference](#settings-reference)). If you have a dataset-fixes script, select it here.
6. Click **Run pre-processing**.
   - If any matched samples have a blank group, the app asks whether to **remove** them (the default) or **label them Case**.
   - If the output folder for this run already exists, the app asks before it overwrites the folder.
7. When the run finishes, the **Run** panel shows the output folder and an **Open report** button.

> The app is busy until the run finishes. The control-only bootstrap regression is the slow step: it runs `Bootstrap resamples` × number of features linear models. Progress messages appear in the **R console**, not in the app.

---

## Dataset-fixes scripts

**This is how you handle dataset-specific problems.** A dataset-fixes script is an optional R script that runs **once, on the raw inputs, before any other step**. It replaces the old practice of hand-editing a `manual fix` chunk in the notebook for each dataset. **The full text of the script is printed in the report**, under *Dataset fixes applied*, so every run records exactly what was patched.

### How to use one

1. Copy `dataset_fixes_TEMPLATE.R` into the folder with your data, and give it a descriptive name (e.g. `MyCohort_CSF_fixes.R`).
2. Edit the copy: uncomment or adapt the examples, or write your own fixes.
3. In the app, under **Section 3: Options**, set **Dataset-fixes script** to your file.

Leave that field blank if your data need no fixes.

### Objects available in the script

| Object | Single cohort mode | Harmonized mode |
|---|---|---|
| `expr.raw` | The expression file as it was read. **SomaScan:** a `soma_adat` (wells x columns). **OLink:** a long `data.frame`, one row per well x assay. **NULISA:** a long `data.frame`, one row per well x target. | The harmonized matrix: a numeric **matrix**, features x samples, on the **log2** scale (a linear input has already been log2-transformed). Row names are feature IDs and column names are sample IDs. |
| `traits.raw` | The traits `data.frame`. | The numericMeta `data.frame`. |
| `platform` | `"SomaScan"`, `"OLink"` or `"NULISA"` | Same (in this mode, used only for labels and the folder name) |
| `tissue` | `"CSF"` or `"Plasma"` | Same |

In both modes, **the row names of `traits.raw` are the sample IDs**, taken from the traits ID column you chose. Samples are matched on these row names, so changing them changes the matching. Columns whose values are all numbers have already been converted to numeric. Every other column is text.

### Rules

- **Change `expr.raw` and/or `traits.raw` in place. Nothing else is kept.** Any other objects you create in the script are discarded afterwards.
- **Put each fix behind a check**, so the fix does nothing when the problem isn't there. An unguarded fix can do real damage. In an earlier version of the pipeline, an unconditional "drop samples with missing values" line emptied the whole matrix when there was nothing to drop. Use `if (...)` and `if (length(bad))` guards, as in the examples below.
- **Guard by data shape as well as by platform when a script could be used in both modes.** `platform` is set in harmonized mode too, so a fix guarded only by `platform == "OLink"` will also run on a harmonized OLink matrix, where `expr.raw$SampleID` fails. Add `is.data.frame(expr.raw)` (long OLink/NULISA tables) or `is.matrix(expr.raw)` (harmonized input) to the check.
- **Single cohort mode:** the platform data object is built from `expr.raw` *after* your script runs. If your script moves the sample IDs into a different column, set **Sample ID column** in the app to the column that holds them after the fix.
- The script runs in a clean environment. It can use base R and installed packages (call them as `pkg::fun()`; `dplyr` is attached), but **not** the helper functions in `preproc_functions.R`.
- A script can be reused across runs and shared with collaborators. Because the report prints the script, each report stays self-documenting.

### Examples (from `dataset_fixes_TEMPLATE.R`)

**SomaScan sample IDs stored in a `DID` column instead of `SampleId`:**
```r
if (platform == "SomaScan" && "DID" %in% names(expr.raw)) {
  expr.raw$SampleId <- ifelse(expr.raw$SampleType == "Sample", expr.raw$DID, expr.raw$SampleId)
}
```

**Drop samples that are known to be bad:**
```r
bad <- c("A0001234", "A0005678")
if (platform == "OLink" && is.data.frame(expr.raw)) expr.raw <- expr.raw[!expr.raw$SampleID %in% bad, ]
```

**Drop samples with missing values in many assays (written safely):**
```r
if (platform == "OLink" && is.data.frame(expr.raw)) {
  na_by_sample <- tapply(is.na(expr.raw$NPX), expr.raw$SampleID, mean)
  bad <- names(na_by_sample)[na_by_sample > 0.10]
  if (length(bad)) expr.raw <- expr.raw[!expr.raw$SampleID %in% bad, ]
}
```

**Make traits IDs match the expression file (e.g. hyphens vs dots):**
```r
rownames(traits.raw) <- gsub("-", ".", rownames(traits.raw))
```

**Harmonized mode: drop sample columns:**
```r
bad <- c("S0123", "S0456")
if (is.matrix(expr.raw) && any(colnames(expr.raw) %in% bad))
  expr.raw <- expr.raw[, !colnames(expr.raw) %in% bad, drop = FALSE]
```

> **Tip:** The app's sample-matching check (step 3 of the Quick start) reads the files directly, so it does **not** apply your fixes script. If you are fixing ID mismatches with a script, you'll still see those mismatches in the Input check panel. The report's *Match samples* section shows the matching after your fixes have been applied, so check it there.

---

## Input files

### Single cohort mode

| File | Format and required columns |
|---|---|
| SomaScan `.adat` | Read with `SomaDataIO::read_adat()`. Metadata needs `SampleType` and the sample ID column. `PlateId` and `PlatePosition` are used when present. |
| OLink Explore HT `.parquet` | Needs the sample ID column plus `SampleType`, `PlateID`, `WellID`, `OlinkID`, `Assay`, `UniProt`, and `NPX` (raw) or `PCNormalizedNPX` (PC-normalised). Uses `AssayType`, `SampleQC`, `LOD`/`PCNormalizedLOD`, `Panel` and `Block` when present. |
| NULISA `.csv` | Needs the sample ID column plus `SampleType`, `PlateID`, `Target`, `UniProtID`, `NPQ` and `targetLOD_NPQ`. Uses `SampleQC`, `Panel`, `ProteinName` and `SampleMatrixType` when present. |
| Traits `.csv` | One row per sample. The ID column must have no blank and no duplicated values. All columns are read as text first, so IDs with leading zeros are kept. The string `NA` means missing. |

Each sample may appear in only one well. Duplicated sample IDs stop the run with a message pointing you to a dataset-fixes script.

### Harmonized matrix mode

| File | Format |
|---|---|
| Harmonized `cleanDat` | `.csv`, **features x samples**. The first column holds the feature IDs, which must be unique (e.g. `Target|UniProt`). The remaining column names are sample IDs. Every other cell must be numeric or blank. |
| `numericMeta` | `.csv`, one row per sample. One ID column must match the `cleanDat` column names. A file written by `write.csv()` with row names has a blank first header. The app names that column `row_names`, so you can select it as the ID column. |

Set **Matrix scale** to say whether the matrix is **log2** (TAMPOR output, NPQ, log2 RFU/NPX) or **linear**. If a matrix declared as log2 has a median above 100, the report warns that it looks linear.

**Feature IDs in harmonized mode:** the text before the first `|` is used as the gene symbol. For **CSF**, the blood-contamination covariate is built from blood markers found by symbol (`HBA1`, `HBB`, `HBG1`, `FGA`, `FGB`, `FGG`, `HP`, `CAT`, `CA1`). Feature IDs therefore need to start with gene symbols for that step to work. All features are assumed to be human proteins.

---

## What the pipeline does

### Single cohort (`templates/preproc_report.Rmd`)

1. **Load** the traits and the raw expression file, then **apply the dataset-fixes script**.
2. **Build a common dataset**: a features x wells matrix with well types (sample, buffer/negative control, calibrator/plate control, QC/sample control) and feature annotation. Values are linear for SomaScan (RFU) and OLink (2^NPX), and log2 NPQ for NULISA. OLink control assays are removed. SomaScan symbols that are listed as `None` are repaired, where possible.
3. **Match samples to traits** and assign the pre-processing groups (CTL / Case). The run stops if there are fewer than 3 controls.
4. **Reference statistics and LOD**:
   - SomaScan: buffer median + 3 SD (RFU), plus the SomaLogic S:N groups
   - OLink: negative-control median + 3 SD on the NPX scale, or Olink's fixed LOD
   - NULISA: the vendor's per-well `targetLOD_NPQ`
5. **Log2 transform** the study samples.
6. **Remove sample outliers** by network connectivity, optionally adding the distance/PC screen.
7. **Remove bimodal features**: dip test p < 0.05 *and* a 2-component Gaussian mixture.
8. *Optional:* **regress out plate**, with the CTL/Case difference protected.
9. **Covariates**: PC1 of the controls, projected onto all samples; for CSF, a blood eigenprotein; age and sex if chosen. A covariate that doesn't vary among the controls is dropped.
10. **Control-only bootstrap regression**: per-feature `lm` in the controls, using the median coefficient over the bootstrap resamples. The fitted effects are then removed from all samples.
11. **Feature filters**:
    - 1-pass: at most the LOD-cutoff fraction of samples below LOD
    - 2-pass (SomaScan, OLink): 1-pass *and* S:N below the threshold in at most the same fraction of samples
12. **Winsorize**: scaled winsorization at ±5 SD. Values beyond a bound are pulled toward it, and their order is kept.
13. **Write the outputs.**

### Harmonized (`templates/harmonized_report.Rmd`)

The well, LOD, signal-to-noise and 1-/2-pass steps are **skipped**, because they were done before harmonization. The steps are:

1. Load the matrix and apply the fixes script.
2. Match samples to numericMeta.
3. Resolve repeat measurements.
4. Remove sample outliers (optional).
5. Remove features with too much missing data.
6. Re-test bimodality (optional).
7. Build the covariates.
8. Run the control-only bootstrap regression.
9. Winsorize.
10. Write the outputs.

---

## Settings reference

### Groups and covariates (both modes)

| Setting | Notes |
|---|---|
| Traits sample ID column | The column that matches the expression sample IDs |
| Diagnostic grouping column | The app tries to guess it from names like `group`, `dx`, `diag` or `clin` |
| Control value | This value becomes **CTL**. Every other non-blank value becomes **Case**. |
| Age column | Leave blank to skip age regression. Non-numeric values become missing. |
| Sex column + male/female values | Leave blank to skip sex regression. The chosen values are recoded to male = 1 and female = 2. Any other value becomes missing. |

Covariate values that are missing are set to the mean (0 after centring), and the report notes this.

**Harmonized mode only:**

| Setting | Notes |
|---|---|
| Batch column | Used to colour the MDS plots. Also needed to choose which repeat to keep (see below). |
| Repeat-measurement column | A column whose value is shared by repeat measurements of the same sample (e.g. bridging samples run in both batches). Only values that occur more than once are affected. |
| Repeat measurements | **Keep all**: every measurement stays. The report notes when control repeats count more than once in the control fit. **Keep one per sample**: the measurement from the chosen batch is kept, along with its metadata. **Average them**: the values are averaged into the chosen batch's column, and `n_measurements_averaged` is added to numericMeta. In every mode, the report lists repeats whose group or sex disagree. |

### Options

| Setting | Mode | Default | Notes |
|---|---|---|---|
| OLink values | Single (OLink) | Raw NPX | Raw NPX or PC-normalised NPX |
| Use Olink fixed LOD | Single (OLink) | Off | Needs an Olink fixed-LOD file and the `OlinkAnalyze` package. Switches the values to PC-normalised NPX. |
| LOD cutoff | Single | 0.75 (NULISA: 0.50) | Keep features with at most this fraction of samples below LOD |
| 2-pass signal-to-noise threshold | Single (SomaScan, OLink) | 0.3 | |
| Drop wells that failed platform sample QC | Single | On | Applies to OLink (`SampleQC` = FAIL in any block) and NULISA (`SampleQC` not "passed"). SomaScan `RowCheck` flags are reported only and are never used to drop wells. |
| Remove sample outliers (connectivity) | Harmonized | On | Always run in single cohort mode |
| Re-test and remove bimodal features | Harmonized | Off | Always run in single cohort mode |
| Maximum fraction missing per feature | Harmonized | 0.5 | Features missing in more samples than this are removed |
| Outlier cutoff | Both | 3 | A sample is an outlier when its connectivity z-score is below −cutoff. Connectivity is computed from bicor adjacency. |
| Distance-based outlier screen | Both | Off | Also flags samples whose distance from the centroid, or whose score on any of the top 3 PCs, has a robust z > 6 |
| Regress out plate | Single | Off | Needs at least 2 plates and no missing values in the matrix |
| Bootstrap resamples | Both | 1000 | |
| Cores | Both | detected cores − 1 (max 7) | Used for the bootstrap regression |
| Seed | Both | 8675309 | Set the same seed to reproduce a run |
| Dataset-fixes script | Both | none | See [Dataset-fixes scripts](#dataset-fixes-scripts) |

**Fixed defaults that are not in the app:** the winsorization pull factor (`winsor_k = 0.25`), the winsorization bound (`winsor_sd = 5`) and the robust-z cutoff for the distance screen (`outlier_robust_crit = 6`). To change them, edit the `params:` block at the top of the report template.

---

## Outputs

Each run writes to its own folder:

```
<output folder>/<Project>_<Platform>_<Tissue>[_<Panel>][_<raw|PC>NPX]/
```

- `<Project>` is the project name with any character other than letters, digits, `.` and `-` removed.
- `_<Panel>` is added for NULISA (from the file's `Panel` column) and is `harmonized` for harmonized runs.
- `_rawNPX` or `_PCNPX` is added for single-cohort OLink runs.

If you rerun with the same settings, the app writes to the same folder. It asks before it overwrites anything.

```
<run folder>/
├── <Project>_report.html             # Full report, including the dataset-fixes script
├── run_parameters.yaml               # Every setting used (plus the LOD cutoff / NPX type actually applied)
├── -cleanDat/                        # Output matrices (see below)
├── -traits/
│   ├── numericMeta.PostSampOut.csv   # Sample metadata after outlier removal, + PreProcGroup, PC1score,
│   │                                 #   bloodEigenprotein (CSF), sex.numeric, PlateID.expr (single cohort)
│   └── input_<traits file>.csv       # Copy of the input traits file (the original is never moved)
├── -pre-proc/
│   ├── masterQC.csv                  # One row per feature: annotation, medians, LOD/S:N, bimodal, Pass/Fail flags
│   ├── sample_outliers.csv
│   ├── bimodal_test.csv              # when bimodality was tested
│   ├── repeated_samples.csv          # harmonized mode, when repeats were found
│   ├── regression_coefficients.csv
│   ├── sessionInfo.txt
│   └── final_workspace.RData         # Every object from the report, for follow-up work
└── QC/                               # PDF copies of every figure in the report
```

### cleanDat file names

Every matrix is features x samples, on the log2 scale, with the post-outlier samples only. In single cohort mode, the files contain **human proteins only**; SomaScan non-human and control analytes are left out.

The file names are built from these parts:

- **Base:** `cleanDat` is the matrix before covariate regression.
- **Regressed versions:** `cleanDat.<covariates>.Reg`, where the covariate labels are `Blood`, `PC`, `Age` and `Sex`.
  - For **CSF**, Blood is part of every regressed version, e.g. `cleanDat.Blood.Reg`, `cleanDat.Blood.PC.Reg` and `cleanDat.Blood.PC.Age.Sex.Reg`.
  - For **plasma**, every non-empty combination of the remaining covariates is produced.
- **Feature filters** (single cohort only): each matrix above also comes as `.1pass.filtered` and, for SomaScan and OLink, `.2pass.filtered`.
- **Winsorized copies:** every one of those matrices also has a copy with `.5SDWinsor` added to the name.

For example, a SomaScan CSF run with age and sex writes 9 regression variants × 3 filter levels × 2 (raw and winsorized) = **54 files**.

---

## Troubleshooting

| Problem | What to do |
|---|---|
| "Cannot run: ..." | The message lists what's missing: files or folders that don't exist, traits not loaded, no control value chosen, identical male and female values, a missing fixed-LOD file, or a missing batch column when you've asked to keep or average repeats. |
| "The run stopped with an error" | Nothing after the failed step was written. Read the error; the R console has the full log. Dataset-specific problems are usually best fixed with a [dataset-fixes script](#dataset-fixes-scripts). |
| "Duplicate measurements for the same feature and well" / "Duplicated sample IDs" / "Sample IDs appear in more than one well" | The same sample was measured in more than one well (or row). Decide which measurement to keep, and drop or rename the others in a dataset-fixes script. |
| "Sample ID column '...' not in the adat metadata" | The error lists the available columns. Set **Sample ID column** to the right one, or create it in a fixes script. |
| "Fewer than 3 control samples after matching" | The control value or grouping column is wrong, or the IDs don't match. Check the *Match samples* section of the report. |
| "Cannot tell which sex values are male and female" | Choose the male and female values in the app. |
| "Covariate '...' does not vary among controls" (warning) | That covariate was dropped from the regression, e.g. because all controls are the same sex. |
| "Fewer than 2 blood markers available" (CSF) | No blood regression was done. In harmonized mode, check that the feature IDs start with gene symbols. |
| "It looks linear, not log2" (harmonized) | Set **Matrix scale** to **Linear**. |
| "OlinkAnalyze is needed for the fixed-LOD option" | Run `install.packages("OlinkAnalyze")`. |
| No File... / Folder... buttons | They only appear when you run the app in RStudio. Elsewhere, paste the full paths. |
| Package install fails at start-up | Install the listed packages by hand, then relaunch the app with `PREPROC_SKIP_INSTALL=1`. |

---

## Data privacy

**Don't commit data to this repository.** Expression files, traits/metadata, dataset-fixes scripts that contain sample IDs, and output folders should stay out of the repository. This `.gitignore` covers the usual cases. It does not catch `dataset_fixes_TEMPLATE.R`, so the template stays in the repo.

```
*.adat
*.parquet
*.csv
*_fixes.R
*.html
*.RData
.Rhistory
.Rproj.user/
```

---

## Citation / contact

<!-- TODO: license and how to cite -->
Maintained by Joshua D. Downer, UCSF Memory and Aging Center.
