# ----------------------------------------- #
# Denoising, and merging (steps 2.4 to 2.7) #
# ----------------------------------------- #

## Some packages must be installed and loaded
library(dada2)
library(data.table)
library(doParallel)
library(foreach)
library(magrittr)
registerDoParallel(cores = 6)
## Load the table of dada2 parameters
path <- "/shared/projects/marechiara/Astan_18SV4"
dada2_param <- read.csv2(paste0(path, "/dada2_pipeline_parameters.csv"), header = TRUE, stringsAsFactors = FALSE)
## Define the following path variable so that it points to the extracted directory
pathF <- "/shared/projects/marechiara/Astan_18SV4/FWD" 
pathR <- "/shared/projects/marechiara/Astan_18SV4/REV"
pathErr <- "/shared/projects/marechiara/Astan_18SV4/dada2/log"
## Arguments for denoising and merging
args <- dada2_param[dada2_param$step == "step02", ]
min_read_nb <- args[args$vairable == "MIN_READ_NUM", 3]
## Track reads through the pipeline
## As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
getA <- function(x) length(getUniques(x))
## selection based on the number of reads
tmp <- fread(paste0(path, "/dada2/log/filter.csv"))
tmpR <- tmpF <- tmp[reads.out >= min_read_nb, sample]
## Identify the name of the run (here one run: PHYTOPORT)
runs <- "GYJ"
runs <- c(paste(runs, "Cut1", sep = "_"), paste(runs, "Cut2", sep = "_"))
## File parsing
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered")
filtFsall <- paste(filtpathF, tmpF, sep = "/")
filtRsall <- paste(filtpathR, tmpR, sep = "/")
sample.names.all <- sub("_trimmed.+$", "", basename(filtFsall))
sample.namesR.all <- sub("_trimmed.+$" ,"", basename(filtRsall))
if(!identical(sample.names.all, sample.namesR.all)) {
  stop("Forward and reverse files do not match.")
}
names(filtFsall) <- sample.names.all
names(filtRsall) <- sample.names.all
filtErr <- sort(list.files(pathErr, pattern = ".rda"))
## set seed to ensure that randomized steps are replicatable
set.seed(100)
## Loop allowing to determinate learn error rates and infer sample composition
foreach(i = runs, .packages = c("data.table","dada2")) %dopar% {
  filtFs <- filtFsall[grep(i, names(filtFsall))]
  filtRs <- filtRsall[grep(i, names(filtRsall))]
  sample.names <- sub("_trimmed.+$","", basename(filtFs))
  # Select the FWD and REV error files
  err <- filtErr[grep(i, filtErr)]
  load(paste(pathErr, err[grep("errF", err)], sep = "/"))
  load(paste(pathErr, err[grep("errR", err)], sep = "/"))
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  track <- vector("list", length(sample.names))
  names(track) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    # Dereplication of forward reads+
    derepF <- derepFastq(filtFs[[sam]])
    # Infer sample composition of forward reads
    ddF <- dada(derepF, err = errF, multithread = FALSE, pool = TRUE)
    # Dereplication of revers reads
    derepR <- derepFastq(filtRs[[sam]])
    # Infer sample composition of reverse reads
    ddR <- dada(derepR, err = errR, multithread = FALSE, pool = TRUE)
    # Merge forward/reverse reads
    merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE, verbose = TRUE)
    mergers[[sam]] <- merger
    track[[sam]] <- data.table(sample = sam, 
                               denoisedF.read = getN(ddF), denoisedR.read = getN(ddR), merged.read = getN(merger),
                               denoisedF.seq = getA(ddF), denoisedR.seq = getA(ddR), merged.seq = getA(merger))
  }
  rm(derepF); rm(derepR)
  # Construct sequence table
  track <- rbindlist(track)
  seqtab <- makeSequenceTable(mergers)
  saveRDS(track, paste0(path, "/dada2/log/track_",i,".rds")) 
  saveRDS(seqtab, paste0(path, "/dada2/seqtab_",i,".rds")) 
  rm(track)
  rm(seqtab)
  for(i in 1:length(mergers)) {
    mergers[[i]] <- data.frame(samples = names(mergers)[i], mergers[[i]])
  }
  mergers <- rbindlist(mergers)
  saveRDS(mergers, paste0(path, "/dada2/log/mergers_",i,".rds"))
}
## Version of packages used to build this document
sessionInfo()