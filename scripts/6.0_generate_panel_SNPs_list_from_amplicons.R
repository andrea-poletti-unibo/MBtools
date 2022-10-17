
library(biomaRt)
library(data.table)
library(tidyverse)

mart <- useEnsembl("ENSEMBL_MART_SNP","hsapiens_snp", "https://feb2014.archive.ensembl.org")
biomaRt::listAttributes(mart)
biomaRt::listFilters(mart)

my_attributes <- c("refsnp_id","chr_name","chrom_start","allele","chrom_strand","validated","minor_allele","minor_allele_freq","minor_allele_count","refsnp_source")

# example <- getBM(my_attributes, "snp_filter", "rs232423", mart)

# getBM(attributes = c("refsnp_id","chr_name","chrom_start"),
#       filters = "snp_filter",
#       values = "rs232423",
#       mart =  mart)

# ex2 <- getBM(attributes = my_attributes,
#              filters = "chromosomal_region",
#              values = "1:20000:21000,3:500000:510000",
#              mart =  mart)



panel <- fread("C:/Users/andre/Dropbox (Personal)/Boston_DFCI_Broad_work/Missionbio_scDNA/Data_from_JB_Ankit/Cell_lines_panel_3074/3074.amplicons")
names(panel) <- c("chr", "start", "end", "ampl_id")
panel$chr <- panel$chr %>% str_remove("chr")

filtervalue <- paste0(panel$chr,":",panel$start,":",panel$end) %>% paste0(collapse = ",")

query <- getBM(attributes = my_attributes,
             filters = "chromosomal_region",
             values = filtervalue,
             mart =  mart)

query$refsnp_id %>% duplicated() %>% table

write_tsv(query, "data/SNPs_list_in_panel.txt")

