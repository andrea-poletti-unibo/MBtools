
library(biomaRt)
library(data.table)
library(tidyverse)

mart <- useEnsembl("ENSEMBL_MART_SNP","hsapiens_snp", "https://feb2014.archive.ensembl.org")
biomaRt::listAttributes(mart)
biomaRt::listFilters(mart)

my_attributes <- c("refsnp_id", "chr_name", "chrom_start",
                   "allele", "chrom_strand", "validated",
                   "minor_allele", "minor_allele_freq",
                   "minor_allele_count", "refsnp_source")


panel <- fread("C:/Users/andre/Dropbox (Personal)/Boston_DFCI_Broad_work/Missionbio_scDNA/Data_from_JB_Ankit/Cell_lines_panel_3074/3074.amplicons")
names(panel) <- c("chr", "start", "end", "ampl_id")
panel$chr <- panel$chr %>% str_remove("chr")

filtervalue <- paste0(panel$chr,":",panel$start,":",panel$end) %>% paste0(collapse = ",")

query <- getBM(attributes = my_attributes,
             filters = "chromosomal_region",
             values = filtervalue,
             mart =  mart)

query$refsnp_id %>% duplicated() %>% table

query$snpID <- paste0("chr",query$chr_name,":",query$chrom_start ,":", query$allele)
query$snpID %>% duplicated() %>% table

query$var_pos <- query$snpID %>% str_extract("chr[0-9]+:[0-9]+")
query$var_pos %>% duplicated() %>% table


write_tsv(query, "data/SNPs_list_in_panel.txt")


###############################################################################

# BENCHMARK analysis on ~1/3 of amplicons (150 amplicons, randomly chosen)

panel_noTx <- panel %>% filter(chr %in% 1:22)

set.seed(123)
panelBench <- panel_noTx[sample(1:519, 150, replace = F),]

filtervalueBench <- paste0(panelBench$chr,":",panelBench$start,":",panelBench$end) %>% paste0(collapse = ",")

queryB <- getBM(attributes = my_attributes,
               filters = "chromosomal_region",
               values = filtervalueBench,
               mart =  mart)

queryB$refsnp_id %>% duplicated() %>% table

queryB$snpID <- paste0("chr",queryB$chr_name,":",queryB$chrom_start ,":", queryB$allele)
queryB$snpID %>% duplicated() %>% table

queryB$var_pos <- queryB$snpID %>% str_extract("chr[0-9]+:[0-9]+")
queryB$var_pos %>% duplicated() %>% table

queryB %>% group_by(var_pos) %>% filter(n()>1) %>% View

queryB_uniqueSNPpos <- queryB[!duplicated(queryB$var_pos),]
queryB_uniqueSNPpos$var_pos %>% duplicated() %>% table
queryB_uniqueSNPpos$snpID %>% duplicated() %>% table

write_tsv(queryB_uniqueSNPpos, "data/SNPs_list_in_panelBenchmark.txt")


