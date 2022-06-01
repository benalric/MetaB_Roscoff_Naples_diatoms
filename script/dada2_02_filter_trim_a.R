# --------------------------------- #
# Filtering and trimming (step 2.2) #
# --------------------------------- #

## Some packages must be installed and loaded
library(dada2)
library(data.table)
library(magrittr)
## Define the following path variable so that it points to the extracted directory
path <- "/shared/projects/marechiara/Astan_18SV4"
pathF <- "/shared/projects/marechiara/Astan_18SV4/FWD" 
pathR <- "/shared/projects/marechiara/Astan_18SV4/REV"
## Create folders for dada2 results
dir.create((paste0(path, "/dada2/log")))
## Load the table of dada2 parameters
dada2_param <- read.csv2(paste0(path, "/dada2_pipeline_parameters.csv"), header = TRUE, stringsAsFactors = FALSE)
## Read in the names of the fastq files 
## Forward and reverse fastq files have format: SAMPLENAME_Cut1_trimmed.fastq.gz and SAMPLENAME_Cut1_trimmed.fastq.gz
fastqFs <- sort(list.files(pathF, pattern = "fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern = "fastq.gz"))
## Select file with a size above to 1000 bytes
fastqFs <- fastqFs[file.size(file.path(pathF, fastqFs)) > 1000]
fastqRs <- fastqRs[file.size(file.path(pathR, fastqRs)) > 1000]
if (length(fastqFs) != length(fastqRs)) {
  stop("Forward and reverse files do not match.")
}
## Arguments for filter and trim
args <- dada2_param[dada2_param$step == "step01", ]
trunc_fwd <- args[args$variable == "TRUNC_FWD", 3] %>% as.numeric
trunc_rev <- args[args$variable == "TRUNC_REV", 3] %>% as.numeric
maxee <- args[args$variable == "MAXEE", 3] %>% as.numeric
## Create path towards directory where filtered results will be stored
pathF <- "/shared/projects/marechiara/Astan_18SV4/FWD" 
pathR <- "/shared/projects/marechiara/Astan_18SV4/REV"
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered")
## Filtering and trimming
for (i in fastqFs) {
  dada2::filterAndTrim(fwd = file.path(pathF, i), filt = file.path(filtpathF, i),
                                     rev = file.path(pathR, i), filt.rev = file.path(filtpathR, i), 
                       truncLen = c(trunc_fwd, trunc_rev),
                       maxEE = maxee, maxN = 0, compress = TRUE, verbose = TRUE, multithread = TRUE) -> out
  
  # Save filtered and trimmed table in the folder of dada2 results
  data.table(out, keep.rownames = TRUE) %>% 
    setnames(names(data.table(out, keep.rownames = TRUE))[1], c("sample")) %>% 
    fwrite(paste0(path, "/dada2/log/filter_", sub("_trimmed.+$", "", basename(i)), ".csv"), sep = ";", col.names = TRUE)
rm(out)
}
## Version of packages used to build this document
sessionInfo()