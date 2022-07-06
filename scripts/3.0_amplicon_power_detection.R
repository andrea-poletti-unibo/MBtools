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



#_____ load amplicons table with stat ______
ampl <- data.table::fread("amplicons_list_allStats.txt")
# ampl$cov_quality2 <- ifelse(ampl$ampl_median <300 ,"ok", "high_cov_M300")

ampl$ampl_median %>% summary

ampl$FILTER <- ifelse(ampl$mappability_score<0.8 | ampl$cov_quality != "ok", 1,0)
ampl$FILTER %>% table

okampl <- ampl %>% filter(FILTER==0)

ampl_to_keep <- okampl %>% .$ID

okampl$ampl_median %>% plot


# add chromosome arm info to amplicons: Popeye function
ampl$chr <- parse_number(ampl$chr)
ampl <- ampl %>% rename(ID_amplicon="ID")

source("C:/Users/andre/Desktop/MM/R/PROJECTS/MM_workspace/[workshop]/Popeye2.R")
ampl_anno <- ampl %>% Popeye2(refGenome = "hg19", removeXY = T)

# ================= merge all samples in one single data.frame ==================

df <- Reduce(rbind, all_df)

#================= FILTERING AMPLICONS =================

# exclude tx amplicons
df2 <- df %>% select(!matches("tx"))

# rename amplicons with chrarm info

names(df2) <- paste0(df2 %>% names,"/", ampl_anno$chrarm)



# exclude bad amplicons
df3 <- df2 %>% select(matches(ampl_to_keep))


#================= FILTERING CELLS =================
df3 %>% rowSums() %>% plot
df3 %>% rowSums() %>% summary
df3 %>% rowSums() %>% quantile(probs = seq(0, 1, 0.1))

df3 %>% rowSums() %>% density %>% plot

df3 %>% rowSums() %>% `<`(30000) %>% table

row_means <- apply(df3, 1, mean )
row_means %>% density() %>%  plot 

row_sd <- apply(df3, 1, sd )

row_sd %>% density() %>%  plot + abline(v = 300)

(row_sd > 300) %>% table

good_cells <- row_sd < 300

df4 <- df3[good_cells,]

df4[df4==0] <- 1


#========== select normal cells based on cell lines CN annotation excel file =============

# load cell lines anno file

cl_CN <- xlsx::read.xlsx("MM_CELL_LINES_CN_annotation.xlsx", sheetIndex = 1)

# amplicon power = median of amplicon X / median of all amplicons medians

# compute medians of every single amplicon BUT based only on normal diploid cells (info from cell lines annot)


df_single_amps <- data.frame()

i= names(df4)[16]
for ( i in names(df4)){
  message(i)
  
  chrarm <- i %>% str_split("/") %>% unlist %>% .[3]

  diploid_cl <-  cl_CN %>% filter(CHR_ARM==chrarm) %>% select(where(~.==0)) %>% names %>% paste0(collapse = "|")
  
  selection_diploid_amplicon <- df4 %>% select(i) %>% filter(grepl(diploid_cl, row.names(df4) ))
  
  selection_diploid_amplicon %>% row.names() %>% str_extract("^.*/") %>% unique # sanity check
  
  n_cl <- selection_diploid_amplicon %>% row.names() %>% str_extract("^.*/") %>% unique %>% length()
  
  med <- selection_diploid_amplicon %>% unlist() %>% median
  
  res <- data.frame(i, med, n_cl, diploid_cl)
  
  df_single_amps <- rbind(df_single_amps, res)
}


df_single_amps %>% ggplot(aes(n_cl %>% as.factor(), med))  + geom_boxplot() + geom_jitter(height = 0, width = 0.2)

# exclude ampl with 6 cell lines (in every cell line there is a CN event - no diploid reference cells)
df_single_amps2 <- df_single_amps %>% filter(n_cl != 6) # 433 ampl remaining

overall_Median <- df_single_amps2$med %>% median

df_single_amps2$power <- df_single_amps2$med/ overall_Median

df_single_amps2$arm <- df_single_amps2$i %>%  str_extract("[0-9]+[pq]")

df_single_amps2 %>% ggplot(aes(arm, power)) + geom_boxplot() + geom_jitter(height = 0, width = 0.2)
df_single_amps2 %>% ggplot(aes(n_cl %>% as.factor(), power)) + geom_boxplot() + geom_jitter(height = 0, width = 0.2)

aov(n_cl ~ power, data = df_single_amps2) %>% summary

write_tsv(df_single_amps2, "amplicon_empirical_powers.txt")














