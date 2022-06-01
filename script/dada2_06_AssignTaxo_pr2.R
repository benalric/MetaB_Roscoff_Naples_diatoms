# ----------------------------- #
# Taxonomic assignment (step 3) #
# ----------------------------- #

## Some packages must be installed and loaded
library(dada2); packageVersion("dada2")
## Define the following path variable so that it points to the extracted directory
pathSeq <- "/shared/projects/marechiara/Astan_18SV4"
pathRef <- "/shared/projects/marechiara/Astan_18SV4/refdb"
## Load sequence table
seqtab.nochim1 <- readRDS(paste0(pathSeq, "/dada2/seqtab.nochim.rds"))
## Pathway to load PR2 reference database
pr2_file <- paste0(pathRef, "/pr2_version_4.13.0_18S_dada2.fasta.gz")
## Fix the taxanomic levels
## PR2 database present specific taxonomic levels
PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", 
                    "Genus", "Species")
## Taxonomic assginment
taxa <- assignTaxonomy(seqtab.nochim1, refFasta = pr2_file, taxLevels = PR2_tax_levels,
                       minBoot = 0, outputBootstraps = TRUE, verbose = TRUE)
save(taxa, file = paste0(pathSeq, "/dada2/seqtab_all_assignTaxo_taxopr2.rda"))
## Version of packages used to build this document
sessionInfo()