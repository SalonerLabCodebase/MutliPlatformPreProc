###############################################################################
# Dataset-fixes script (optional). Copy this file next to your data, edit it,
# and choose it in the app under "Dataset-fixes script".
#
# It replaces the old `adat.manual.fix` chunk. It runs once, on the raw
# inputs, before anything else, and its full text is printed in the report so
# every run records what was patched.
#
# Available objects:
#   expr.raw    the expression file as read:
#                 SomaScan: a soma_adat (wells x columns)
#                 OLink:    a long data.frame (one row per well x assay)
#                 NULISA:   a long data.frame (one row per well x target)
#   traits.raw  the traits data.frame; row names are the sample IDs
#   platform, tissue
#
# Modify `expr.raw` and/or `traits.raw` in place. Nothing else is kept.
# Keep each fix behind a check so the script does nothing when the problem is
# absent (the old file-16 line 517 emptied the matrix when there was nothing to fix).
###############################################################################

## Example: SomaScan sample IDs live in a DID column rather than SampleId (file 11)
# if (platform == "SomaScan" && "DID" %in% names(expr.raw)) {
#   expr.raw$SampleId <- ifelse(expr.raw$SampleType == "Sample", expr.raw$DID, expr.raw$SampleId)
# }

## Example: drop samples known to be bad
# bad <- c("A0001234", "A0005678")
# if (platform == "OLink") expr.raw <- expr.raw[!expr.raw$SampleID %in% bad, ]

## Example: drop samples with missing values in many assays (the file-16 hack, made safe)
# if (platform == "OLink") {
#   na_by_sample <- tapply(is.na(expr.raw$NPX), expr.raw$SampleID, mean)
#   bad <- names(na_by_sample)[na_by_sample > 0.10]
#   if (length(bad)) expr.raw <- expr.raw[!expr.raw$SampleID %in% bad, ]
# }

## Example: harmonise traits IDs to the expression file (e.g. hyphens vs dots)
# rownames(traits.raw) <- gsub("-", ".", rownames(traits.raw))
