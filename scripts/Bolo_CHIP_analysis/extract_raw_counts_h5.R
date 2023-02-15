# BiocManager::install("rhdf5",update = F)

library(rhdf5)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(RColorBrewer)

#================= LOCALIZE DATA ==================

H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/RUN5/"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# sample

for(i in 1:length(list_h5_files)){

  h5filePath <- list_h5_files[i]

  samp_num <- str_extract(h5filePath, "[1-9]+(?=\\.dna)")

  rhdf5::h5ls(h5filePath)


  #============ LOAD AND PREPARE RAW READ COUNT ===========

  rc <- h5read(h5filePath, "/assays/dna_read_counts/layers")

  df <- rc$read_counts %>% t %>% as.data.frame()

  # extract and set colnames(ampl) and rownames(cells)
  names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )
  cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")

  rownames(df) <- cells$barcode
  colnames(df) <- names$id

  exp <- rownames_to_column(df, var = "read_counts")


  write_csv(exp, paste0("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/python_notebook_CNA/Spyder_analysis/raw_reads/sample", samp_num, "_raw_reads.csv"))



  #===== reproduce mosaic normalization approach in R ==========

  mat <- df %>% as.matrix()

  #_____ find good cells ____
  cell_sum <- rowSums(mat)
  good <- cell_sum > sort(cell_sum, decreasing = T)[11]/10
  good %>% table

  #_____ Row/cell normalizarion______
  # divide each row/cell for it's mean across cols/ampl , then add 1

  row_means <- apply(mat, 1, mean ) + 1

  m2 <- mat %>% sweep(1, row_means, `/` )

  # pheatmap(m2, cluster_rows = T, cluster_cols = F,
  #          show_rownames = F, show_colnames = F,
  #          annotation_col = chr_ann )



  #_____ Column/amplicon normalizarion______
  # divide each column/amplicon for it's median across rows/cells (but removing cells not GOOD!), then add 0.05

  m2[good,] %>% dim
  col_medians <- apply(m2[good,], 2, median ) + 0.05

  m3 <- m2 %>% sweep(2, col_medians, `/` )
  dim(m3)

  #_______ multiply per 2 - for diploidy__________
  m4 <- (m3 * 2) %>% round(6)

  dim(m4)

  exp2 <- m4 %>% as.data.frame() %>% rownames_to_column( var = "read_counts")

  write_csv(exp2, paste0("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/python_notebook_CNA/Spyder_analysis/normalized_reads/sample", samp_num, "_norm_reads.csv"))


}
