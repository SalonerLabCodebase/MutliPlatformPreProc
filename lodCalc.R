## Filter npq.raw to only include samples found in the traits file ----
npq.raw <- npq.raw[npq.raw$SampleName %in% rownames(traits.raw), ]

## Create the Target | UniProtID identifier ----
npq.raw$TargetID <- paste(npq.raw$Target, 
                           npq.raw$UniProtID, 
                           sep = "|")

cleanDat <- pivot_wider(
  npq.raw,
  id_cols = TargetID,
  names_from = SampleName,
  values_from = NPQ
)

cleanDat<-data.frame(cleanDat)
rownames(cleanDat)<-cleanDat$TargetID
cleanDat$TargetID<-NULL
cleanDat.targetLOD_NPQ <- pivot_wider(
  npq.raw,
  id_cols = TargetID,
  names_from = SampleName,
  values_from = targetLOD_NPQ
)

above.lod.mat <- cleanDat >= cleanDat.targetLOD_NPQ[rownames(cleanDat), colnames(cleanDat)]
cleanDat.neat.softNA <- cleanDat
cleanDat.neat.softNA[!above.lod.mat]<-NA

missingnessCrit = floor(0.5*ncol(cleanDat.neat.softNA))
missing_mask <- is.na(cleanDat.neat.softNA)
num_missing <- rowSums(missing_mask)
keepRownames.1pass <- rownames(cleanDat.neat.softNA)[num_missing <= missingnessCrit]

keepRowIndices.1pass <- which(rownames(cleanDat.neat.PostSampOut.noBimodalAptamers) %in% keepRownames.1pass)

# Assays with <= 50% of samples below LOD are retained