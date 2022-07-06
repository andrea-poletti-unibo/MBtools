library(rhdf5)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

#================= LOCALIZE DATA ==================
list_h5_files <- list.files("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/", pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

all_df <- list()
for(h5filePath in list_h5_files){
  
  message(h5filePath)
  
  sample <- h5filePath %>% str_extract("(?<=scDNA_).*(?=_reference)")
  
  rc <- h5read(h5filePath, "/assays/dna_read_counts/layers")
  
  df <- rc$read_counts %>% t %>% as.data.frame() 
  
  # extract and set colnames(ampl) and rownames(cells)
  amplicons <- h5read(h5filePath, "/assays/dna_read_counts/ca" ) %>% as.data.frame()
  cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")
  
  amplicons_f <- amplicons %>% select(chr=CHROM, 
                                      start=start_pos, 
                                      end= end_pos)
  
  
  rownames(df) <- paste0(sample,"/", cells$barcode)
  colnames(df) <- paste0(amplicons$id,
                         "/chr", amplicons$CHROM,":", 
                         amplicons$start_pos, "-", 
                         amplicons$end_pos)
  
  all_df[[sample]] <- df
}


df <- Reduce(rbind, all_df)


#_____ load amplicons table with stat ______
ampl <- data.table::fread("amplicons_list_allStats.txt")
ampl$cov_quality2 <- ifelse(ampl$ampl_median <300 ,"ok", "high_cov_M300")

ampl$FILTER <- ifelse(ampl$mappability_score<0.8 | ampl$cov_quality != "ok", 1,0)
ampl$FILTER %>% table

okampl <- ampl %>% filter(FILTER==0)

ampl_to_keep <- okampl %>% .$ID


# exclude bad amplicons
df3 <- df %>% select(matches(ampl_to_keep))


#================= FILTERING CELLS =================

row_sd <- apply(df3, 1, sd )
good_cells <- row_sd < 300
df3 <- df3[good_cells,]
df3[df3==0] <- 1


#================= SELECTING TX AMPLICONS =================

# only tx amplicons
dft <- df %>% select(matches("tx"))
dft <- dft[good_cells,]



################################ CELL CLUSTERING ####################################

# create annotate chr df for heatmap
sample_ann <- data.frame(sample=umap_df$sample %>% as.factor(), 
                        row.names = rownames(df3))


dft %>% as.matrix() %>% as.vector() %>% summary

dft %>% as.matrix() %>% density() %>% plot
dft %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.001))

breaksList = seq(0,50, 1)

sampleGap <- sample_ann$sample %>% table %>% cumsum()

ann_colors = list(sample = c("KMM1"= '#7fc97f', "KMS11"= '#beaed4', "Mix"='#fdc086', "MM1S"= '#ffff99', "OPM2"= '#386cb0', "RPMI8226"="grey"))

pheatmap(dft, 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = T,
         gaps_row = sampleGap,
         breaks = breaksList,
         filename = "C:/Users/andre/Desktop/all_myClusters_TX.png", width = 15, height = 15,
         annotation_row = sample_ann,
         annotation_colors = ann_colors)

dev.off()

pheatmap(dft_c, 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = T,
         scale = "row",
         gaps_row = sampleGap,
         filename = "C:/Users/andre/Desktop/all_myClusters_TX_rowScaled.png", width = 15, height = 15,
         annotation_row = sample_ann,
         annotation_colors = ann_colors)

dev.off()
