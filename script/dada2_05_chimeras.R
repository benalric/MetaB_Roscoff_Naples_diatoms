# -------------------------- #
# Remove Chimeras (step 2.8) #
# -------------------------- #

## Some packages must be installed and loaded
library(dada2); packageVersion("dada2")
library(data.table)
library(magrittr)
library(digest)
## Define the following path variable so that it points to the extracted directory
pathFilt <- "/shared/projects/marechiara/Astan_18SV4/log"
pathSave <- "/shared/projects/marechiara/Astan_18SV4"
## Merge multiple runs (if necessary)
x <- grep("seqtab_.+\\.rds$", dir("/shared/projects/marechiara/Astan_18SV4/dada2/"), value = TRUE)
x <- paste0("/shared/projects/marechiara/Astan_18SV4/dada2/", x)
st.all <- mergeSequenceTables(tables = x)
## Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method = "pooled", multithread = TRUE, minFoldParentOverAbundance = 2)
seqtab.nochim2 <- seqtab.nochim %>% t %>% data.table
## load statistics
x <- grep("track_.+\\.rds$",dir("/shared/projects/marechiara/Astan_18SV4/log/"), value = TRUE)
x <- paste0("/shared/projects/marechiara/Astan_18SV4/log/", x)
stattab <- lapply(x, readRDS) %>%
  rbindlist
stattab[, nochim.read := sapply(sample, function(X){
  sum(seqtab.nochim2[, get(X)])
})]
stattab[, nochim.seq := sapply(sample, function(X){
  sum(seqtab.nochim2[, get(X)] != 0)
})]
stattab <- stattab[, list(sample, denoisedF.read, denoisedR.read, merged.read, 
                          nochim.read, denoisedF.seq, denoisedR.seq, merged.seq, nochim.seq)]
filter <- read.csv2(paste0(pathFilt, "/filter.csv"), header = TRUE, stringsAsFactors = FALSE)
filter$sample <- sub("_trimmed.+$","", filter$sample)
stattab <- merge(stattab, filter, by.x = "sample", by.y = "sample")
stattab <- stattab[, c("sample", "reads.in", "reads.out",
                       "denoisedF.read", "denoisedR.read", "merged.read", "nochim.read",
                       "denoisedF.seq", "denoisedR.seq", "merged.seq", "nochim.seq")]
fwrite(stattab, "/shared/projects/marechiara/Astan_18SV4/log/statistics.tsv", sep = "\t")
rm(stattab, seqtab.nochim2)
## Construct sequence table
tmp <- seqtab.nochim %>% t %>% data.table(keep.rownames = TRUE)
setnames(tmp, "rn", "sequence")
saveRDS(seqtab.nochim, paste0(pathSave, "/dada2/seqtab.nochim",".rds"))
fwrite(seqtab.nochim, "/shared/projects/marechiara/Astan_18SV4/dada2/seqtab.nochim.tsv", sep = "\t")
rm(seqtab.nochim)
tmp <- data.table(amplicon = sapply(tmp[, sequence], digest, algo = "sha1"), tmp)
fwrite(tmp, "/shared/projects/marechiara/Astan_18SV4/dada2/seqtab_all.tsv", sep = "\t")
## Version of packages used to build this document
sessionInfo()