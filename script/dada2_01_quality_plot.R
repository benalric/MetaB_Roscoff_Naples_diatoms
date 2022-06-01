# ---------------------------------------- #
# Inspect read quality profiles (step 2.1) #
# ---------------------------------------- #

## Some packages must be loaded
library(dada2)
library(magrittr)
library(ggplot2)
## Set up pathway to forward and reverse FASTQ files
pathF <- "/shared/projects/marechiara/Astan_18SV4/FWD" 
pathR <- "/shared/projects/marechiara/Astan_18SV4/REV"
## Create folders for quality profiles
dir.create(paste0(pathF, "/qualityplot"))
dir.create(paste0(pathR, "/qualityplot"))
## Set file paths where quality plots will be stored
filtpathF <- file.path(pathF, "qualityplot") 
filtpathR <- file.path(pathR, "qualityplot")
## Get a list of all FASTQ files in the work directory and separate FWD and REV 
## Forward and reverse FASTQ files have format: 
## SAMPLENAME_Cut1_trimmed.fastq.gz and SAMPLENAME_Cut2_trimmed.fastq.gz
fastqFs <- sort(list.files(pathF, pattern = "fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern = "fastq.gz"))
## Select file with a size above to 1000 bytes
fastqFs <- fastqFs[file.size(file.path(pathF, fastqFs)) > 1000]
fastqRs <- fastqRs[file.size(file.path(pathR, fastqRs)) > 1000]
if (length(fastqFs) != length(fastqRs)) {
  stop("Forward and reverse files do not match.")
}
## Identify the name of the run (here one run: GYJ)
runs <- sub("^RA\\d+_\\d+_([^_]+).+$", "\\1", fastqFs) %>% unique
## Quality plots to the sample level
# FWD
pdf(file.path(filtpathF, "indiv_F_Qplots.pdf"))
for (i in file.path(pathF, fastqFs)[1]) {
  print(plotQualityProfile(i))
}
dev.off()
# REV
pdf(file.path(filtpathR, "indiv_R_Qplots.pdf"))
for (i in file.path(pathR, fastqRs)) {
  print(plotQualityProfile(i))
}
dev.off()
## Version of packages used to build this document
sessionInfo()