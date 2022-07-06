
#=========== Get the GC percentage content of genomic regions =============

# #___________________ first approach ________________
# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome)
# gr <- GenomicRanges::GRanges(seqnames = c("chr1", "chr2"),
#                              ranges = IRanges(start = c(123456,123490),
#                                               end = c(500020, 600020)))
# GetGC <- function(bsgenome, gr){
#   seqs <- BSgenome::getSeq(bsgenome, gr)
#   return(as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)))
# }
#
# GetGC(bsgenome = BSgenome.Hsapiens.UCSC.hg19, gr = gr)

#__________________ second approach _______________

# BiocManager::install("Repitools")
require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome)
library(Repitools)

# Granges<- gr
# gc<- gcContentCalc(Granges , organism=Hsapiens, verbose=TRUE)


#================= LOCALIZE DATA ==================

library(rhdf5)
library(tidyverse)

H5_files_dir <- "<H5_FILES_PATH>"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )

regions <- data.frame(chr=paste0("chr",names$CHROM),
                      start=names$start_pos,
                      end=names$end_pos,
                      ID_amplicon=names$id)

# filter out tx amplicons
regions_f <- regions %>% filter(!str_detect(chr, "tx"))


Gregions <- makeGRangesFromDataFrame(regions_f, keep.extra.columns = T)


############################## GC CONTENT ###############################


gc_Gregions <- gcContentCalc(Gregions , organism=Hsapiens, verbose=TRUE)


regions_f$GC_content <- gc_Gregions

regions_f %>%
  mutate(sort_ID= factor(ID_amplicon, levels = ID_amplicon[order(regions_f$GC_content)], ordered = T )) %>%
  ggplot(aes(sort_ID %>% as.numeric(), GC_content)) + geom_point()


############################## MAPPABILITY ###############################

# BiocManager::install("WGSmapp")
# BiocManager::install("SCOPE")

library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SCOPE)

mapp <- get_mapp(Gregions, hgref = "hg19")

mapp


regions_f$mappability_score <- mapp


#================== export list of amplicons with GC and Mapp info ====================

write_tsv(regions_f, "data/amplicons_list_GC_mapp.txt")
