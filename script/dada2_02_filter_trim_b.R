# --------------------------------- #
# Filtering and trimming (step 2.2) #
# --------------------------------- #

## Some packages must be installed and loaded
library(data.table)
library(magrittr)
## Merge all files in one table
path <- "/shared/projects/marechiara/Astan_18SV4/dada2/log"
myfiles <- list.files(path = path, pattern = "*.csv", full.names = TRUE)
filter <- lapply(myfiles, read.csv2)
filter <- do.call('rbind', filter)
filter %>% 
  data.table() %>%
  fwrite(paste0(path, "/filter.csv"), sep = ";", col.names = TRUE)
## Version of packages used to build this document
sessionInfo()