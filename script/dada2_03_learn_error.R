# ---------------------------- #
# Learn error rates (step 2.3) #
# ---------------------------- #

## Some packages must be installed and loaded
library(dada2)
library(data.table)
library(foreach)
library(magrittr)
## Load the table of dada2 parameters
path <- "/shared/projects/marechiara/Astan_18SV4"
dada2_param <- read.csv2(paste0(path, "/dada2_pipeline_parameters.csv"), header = TRUE, stringsAsFactors = FALSE)
## Define the following path variable so that it points to the extracted directory
pathF <- "/shared/projects/marechiara/Astan_18SV4/FWD" 
pathR <- "/shared/projects/marechiara/Astan_18SV4/REV"
## Arguments for denoising and merging
args <- dada2_param[dada2_param$step == "step02", ]
min_read_nb <- args[args$variable == "MIN_READ_NUM", 3]
## Track reads through the pipeline
## As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
getA <- function(x) length(getUniques(x))
## selection based on the number of reads
tmp <- fread(paste0(path, "/dada2/log/filter.csv"))
tmpR <- tmpF <- tmp[reads.out >= min_read_nb, sample]
## Identify the name of the run (here one run: MC)
runs <- "MC"
runs <- c(paste(runs, "Cut1", sep = "_"))
## File parsing
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered")
filtFsall <- paste(filtpathF, tmpF, sep = "/")
filtRsall <- paste(filtpathR, tmpR, sep = "/")
sample.names.all <- sub("_trimmed.+$", "", basename(filtFsall))
sample.namesR.all <- sub("_trimmed.+$" , "", basename(filtRsall))
if(!identical(sample.names.all, sample.namesR.all)) {
  stop("Forward and reverse files do not match.")
}
names(filtFsall) <- sample.names.all
names(filtRsall) <- sample.names.all
## set seed to ensure that randomized steps are replicable
set.seed(100)
## Loop allowing to determinate learn error rates and infer sample composition
filtFs <- filtFsall
filtRs <- filtRsall
sample.names <- sub("_trimmed.+$", "", basename(filtFs))
## Learn forward error rates
errF <- learnErrors(filtFs, nbases = 1e8, multithread = FALSE)
pdf(paste0("dada2/log/errF_", runs, ".pdf"))
print(plotErrors(errF, nominalQ = TRUE))
dev.off()
## Learn reverse error rates
errR <- learnErrors(filtRs, nbases = 1e8, multithread = FALSE)
pdf(paste0("dada2/log/errR_", runs, ".pdf"))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()
save(errF, file = paste0(path, "/dada2/log/errF_", runs,".rda")) 
save(errR, file = paste0(path, "/dada2/log/errR_", runs,".rda")) 
## Version of packages used to build this document
sessionInfo()