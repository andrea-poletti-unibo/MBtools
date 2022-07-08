
install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update = F)

library(tidyverse)

df <- data.table::fread("data/h5_files_stats.tsv", skip = 1)
df$proportion_mappedToCell_total_read_pairs <- df$n_read_pairs_mapped_to_cells %>% as.numeric() / df$n_read_pairs %>% as.numeric()
df$proportion_validCell_total_read_pairs <- df$n_read_pairs_valid_cell_barcodes %>% as.numeric() / df$n_read_pairs %>% as.numeric()
df$proportion_mappedInsert_total_reads <- df$n_reads_mapped_insert %>% as.numeric() / (df$n_read_pairs %>% as.numeric() *2)
df$proportion_bases_r1_q30 <- df$n_bases_r1_q30 %>% as.numeric() / df$n_bases_r1 %>% as.numeric()
df$proportion_bases_r2_q30 <- df$n_bases_r2_q30 %>% as.numeric() / df$n_bases_r2 %>% as.numeric()
df$proportion_bases_cellBarcode_q30 <- df$n_cell_barcode_bases_q30 %>% as.numeric() / df$n_cell_barcode_bases %>% as.numeric()

df$n_cell_barcode_bases_q30/1000000000
df$n_cell_barcode_bases/1000000000


melt <- reshape2::melt(df, id.vars="sample")



dfB1 <- data.table::fread("data/h5_files_stats_RUN1.tsv", skip = 1)
dfB1$proportion_mappedToCell_total_read_pairs <- dfB1$n_read_pairs_mapped_to_cells %>% as.numeric() / dfB1$n_read_pairs %>% as.numeric()
dfB1$proportion_validCell_total_read_pairs <- dfB1$n_read_pairs_valid_cell_barcodes %>% as.numeric() / dfB1$n_read_pairs %>% as.numeric()
dfB1$proportion_mappedInsert_total_reads <- dfB1$n_reads_mapped_insert %>% as.numeric() / (dfB1$n_read_pairs %>% as.numeric() * 2)
dfB1$proportion_bases_r1_q30 <- dfB1$n_bases_r1_q30 %>% as.numeric() / dfB1$n_bases_r1 %>% as.numeric()
dfB1$proportion_bases_r2_q30 <- dfB1$n_bases_r2_q30 %>% as.numeric() / dfB1$n_bases_r2 %>% as.numeric()
dfB1$proportion_bases_cellBarcode_q30 <- dfB1$n_cell_barcode_bases_q30 %>% as.numeric() / dfB1$n_cell_barcode_bases %>% as.numeric()

meltB1 <- reshape2::melt(dfB1, id.vars="sample")


dfB2 <- xlsx::read.xlsx("data/h5_files_stats.xlsx", sheetIndex = 1)
dfB2$proportion_mappedToCell_total_read_pairs <- dfB2$n_read_pairs_mapped_to_cells %>% as.numeric() / dfB2$n_read_pairs %>% as.numeric()
dfB2$proportion_validCell_total_read_pairs <- dfB2$n_read_pairs_valid_cell_barcodes %>% as.numeric() / dfB2$n_read_pairs %>% as.numeric()
dfB2$proportion_mappedInsert_total_reads <- dfB2$n_reads_mapped_insert %>% as.numeric() / (dfB2$n_read_pairs %>% as.numeric() * 2)
dfB2$proportion_bases_r1_q30 <- dfB2$n_bases_r1_q30 %>% as.numeric() / dfB2$n_bases_r1 %>% as.numeric()
dfB2$proportion_bases_r2_q30 <- dfB2$n_bases_r2_q30 %>% as.numeric() / dfB2$n_bases_r2 %>% as.numeric()
dfB2$proportion_bases_cellBarcode_q30 <- dfB2$n_cell_barcode_bases_q30 %>% as.numeric() / dfB2$n_cell_barcode_bases %>% as.numeric()

meltB2 <- reshape2::melt(dfB2, id.vars="sample")



allmelt <- rbind(melt, meltB1, meltB2)

interest_vars <- df %>% select_if(is.numeric) %>% names
interest_vars

allmelt$ID <- allmelt$sample %>% str_remove_all("scDNA_|_reference")

allmelt$ID <- allmelt$ID %>% factor(levels = str_sort(allmelt$ID %>% unique, numeric = T), ordered = T)



allmelt %>% filter(variable %in% interest_vars) %>%
  mutate(val=as.numeric(value),
         dataset=ifelse(str_detect(sample,"Sample"), "Bologna_run2",
                        ifelse(str_detect(sample,"sample"), "Bologna_run1","Ghobrial"))) %>%
  ggplot(aes(ID, val, colour=dataset)) +
  geom_point() +
  facet_wrap(~ variable, ncol = 5, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")


allmelt %>% filter(variable %in% interest_vars) %>%
  filter(sample %>% str_detect("DNA")) %>%
  mutate(val=as.numeric(value)) %>%
  ggplot(aes(ID, val)) +
  geom_point() +
  facet_wrap(~ variable, ncol = 5, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")


# group variable exploration

# CUTOFFS
allmelt %>% filter( variable %>% str_detect("cutoff")) %>%
  mutate(val=as.numeric(value),
         dataset=ifelse(str_detect(sample,"Sample"), "Bologna_run2",
                        ifelse(str_detect(sample,"sample"), "Bologna_run1","Ghobrial"))) %>%
  ggplot(aes(ID, val, colour=dataset)) +
  geom_point() +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")



# important
allmelt %>% filter( variable %>% str_detect("amplicons|avg_panel_uniformity|ado_rate|n_cells|n_passing_variants|read_pairs$|reads$")) %>%
  mutate(val=as.numeric(value),
         dataset=ifelse(str_detect(sample,"Sample"), "Bologna_run2",
                        ifelse(str_detect(sample,"sample"), "Bologna_run1","Ghobrial"))) %>%
  ggplot(aes(ID, val, colour=dataset)) +
  geom_point() +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")



# per base stats
allmelt %>% filter( variable %>% str_detect("proportion_bases")) %>%
  mutate(val=as.numeric(value),
         dataset=ifelse(str_detect(sample,"Sample"), "Bologna_run2",
                        ifelse(str_detect(sample,"sample"), "Bologna_run1","Ghobrial"))) %>%
  ggplot(aes(ID, val, colour=dataset)) +
  geom_point() +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
