#################### SNP score matching cell-lines samples ###################

library(tidyverse)
library(data.table)
library(rhdf5)
library(gridExtra)


H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"
list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

snpdata <- data.table::fread("data/SNPs_list_in_panel.txt")

profiles_files <- list.files("data/SNP_profiles/", full.names = T)
profiles_list <- lapply(profiles_files,fread)
names(profiles_list) <- profiles_files %>% str_remove("data/SNP_profiles/") %>% str_remove("_SNP_profile.txt")


dir.create("plots/SNP_analysis/SNP_demultplexing/cell-lines_samples/")

#=============== loop over samples ============

list_h5_files_noMix <- list_h5_files[-3]

SAMP=list_h5_files_noMix[1]
for(SAMP in list_h5_files_noMix){

  h5filePath <- SAMP

  filename <- h5filePath %>% str_split("/") %>% .[[1]] %>% tail(1)
  samplename <- filename %>% str_extract("(?<=scDNA_).*(?=_reference)")
  filename %>% message



  # extract barcodes (cells)
  barcodes <- h5read(h5filePath, "/assays/dna_variants/ra")
  bc_df <- Reduce(cbind, barcodes) %>% as.data.frame()
  colnames(bc_df) <- names(barcodes)

  bc_df$cell_idx <- paste0(bc_df$sample_name[1] %>% str_extract("(?<=scDNA_).*(?=_reference)"),"_", 1:nrow(bc_df))
  bc_df$cell_sample_id <- paste0(bc_df$cell_idx, "_", bc_df$barcode)

  # extract variants
  variants <- h5read(h5filePath, "/assays/dna_variants/ca")
  var_df <- Reduce(cbind, variants) %>% as.data.frame()
  colnames(var_df) <- names(variants)

  var_df <- var_df %>% dplyr::select(id, CHROM, POS, REF, ALT, QUAL, ado_gt_cells, ado_rate, amplicon, filtered)

  IDX <- var_df$filtered == "00"
  var_pass <- var_df[IDX,]

  h5f = H5Fopen(h5filePath)

  #_____ AF _______
  m_AF <- h5f$"assays/dna_variants/layers/AF" %>% as.data.frame()
  m_AF_f <- m_AF[IDX,]
  names(m_AF_f) <- bc_df$cell_sample_id
  rownames(m_AF_f) <- var_pass$id
  rm(m_AF)

  #_____ NGT _______
  m_NGT <- h5f$"assays/dna_variants/layers/NGT" %>% as.data.frame()
  m_NGT_f <- m_NGT[IDX,]
  names(m_NGT_f) <- bc_df$cell_sample_id
  rownames(m_NGT_f) <- var_pass$id
  rm(m_NGT)

  #===== long versions =====
  m_AF_fl <- m_AF_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "VAF")
  m_NGT_fl <- m_NGT_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "NGT")


  #====== merge all =======
  m_ALL_fl <- left_join(m_AF_fl, m_NGT_fl)


  #======= merge with SNP data =========
  m_ALL_fl$var_pos <- m_ALL_fl$varID %>% str_extract("chr[0-9]+:[0-9]+")
  snpdata$var_pos <- snpdata$snpID %>% str_extract("chr[0-9]+:[0-9]+")


  # m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by=c("varID"= "snpID"))
  m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by="var_pos")

  m_ALL_onlySNP_fl <- m_ALL_SNP_fl %>% filter(!is.na(refsnp_id))
  m_ALL_onlySNP_fl$cell_n <- m_ALL_onlySNP_fl$cellID %>% str_extract("(?<=_)[0-9]+(?=_)") %>% as.numeric()


  ########################### matching ################################

  total_cells <- m_ALL_onlySNP_fl$cellID %>% unique() %>% length

  i=1

  cells_res_df <- data.frame()
  for (i in 1:total_cells){

    message(i)
    cell_name <- m_ALL_onlySNP_fl %>% filter(m_ALL_onlySNP_fl$cell_n == i) %>% .$cellID %>% unique %>% as.character()

    celldf <- m_ALL_onlySNP_fl %>% filter(cellID== cell_name)

    j=1

    for( j in 1:length(profiles_list)){

      prof <- profiles_list[[j]]
      sample_name <- names(profiles_list)[j]

      df <- left_join(celldf, prof, by = "snpID")

      df$concordance <- df$NGT == df$callNGT

      res_tab <- data.frame(id= cell_name,
                            False= sum(df$concordance==F, na.rm = T),
                            True= sum(df$concordance==T, na.rm = T),
                            Na= sum(is.na(df$concordance)),
                            sample_match= sample_name)

      cells_res_df <- rbind(cells_res_df, res_tab)

    }

  }

  cells_res_df$tot <- cells_res_df$False + cells_res_df$True + cells_res_df$Na

  #________ create score1: T/F ratio ________

  cells_res_df$ratio_T_F <- cells_res_df$True/cells_res_df$False

  table(cells_res_df$ratio_T_F > 1) # set threshold at 1
  cutoff <- 1

  mytab <- table(cells_res_df$ratio_T_F > cutoff) %>% as.data.frame() %>% setNames(c("matched","n_cells"))

  cells_res_df %>% arrange(ratio_T_F) %>%
    ggplot(aes(1:nrow(cells_res_df), ratio_T_F, colour=sample_match)) +
    geom_point() + ggtitle("Knee plot for score1 (True/False)") +
    geom_hline(yintercept = cutoff, linetype=2) +
    annotation_custom(tableGrob(mytab), xmin=100, xmax=3000, ymin=2, ymax=3)


  cells_res_df$matched_score1 <- ifelse(cells_res_df$ratio_T_F > cutoff, "match", "no_match")


  # SAVE data matching
  dir.create("data/cells_cell-lines_SNP_matching/", showWarnings = F)
  write_tsv(cells_res_df, paste0("data/cells_cell-lines_SNP_matching/",samplename,"_cells_matching_SNPprofiles.txt"))


  #_____________ faceted-per-cell-line Knee plot of score1 _____________

  table(cells_res_df$sample_match, cells_res_df$matched_score1)

  mytab <- table(cells_res_df$sample_match, cells_res_df$matched_score1) %>% as.data.frame.matrix()

  a <- cells_res_df %>% arrange(ratio_T_F) %>%
    ggplot(aes(1:nrow(cells_res_df), ratio_T_F)) + geom_hline(yintercept = 1, linetype=2, size = 0.5) +
    geom_point(aes(colour=matched_score1), size=0.8) +
    ggtitle(label = samplename, subtitle = "Knee plot for score1 (True/False) faceted per cell-line") +
    facet_wrap(~sample_match) + theme(legend.position = "bottom")

  b <- ggplot() + annotation_custom(tableGrob(mytab), xmin=0, xmax=1, ymin=0, ymax=1)

  grid.arrange(a,b,ncol=2,widths=c(3/4,1/4))

  g <- arrangeGrob(a,b,ncol=2,widths=c(3/4,1/4))
  ggsave(plot = g, filename = paste0("plots/SNP_analysis/SNP_demultplexing/cell-lines_samples/",samplename,"_knee_plot_score1_facet_PER_CELL-LINE.png"),
         width = 12, height = 8)

  #__________ single cell analysis __________


  nsc <- cells_res_df %>% group_by(id) %>% summarise()

  sc <- cells_res_df %>% group_by(id) %>% group_map(~ sort(.x$ratio_T_F, decreasing = T))

  sc_df <- Reduce(rbind, sc) %>% as.data.frame(row.names = nsc$id) %>% rownames_to_column(var = "id")

  names(sc_df) <- c("id", paste0("score_best_",1:5))

  sc_df$diff2best <- sc_df$score_best_1 - sc_df$score_best_2

  sc_df$diff2best %>% sort %>% plot

  sct <- cells_res_df %>% select(id, sample_match, ratio_T_F) %>% reshape2::dcast(id~sample_match)

  ranks <- apply(sct[,-1], 1, rank) %>% t %>% as.data.frame()

  ranks$KMM1 %>% table
  ranks$KMS11 %>% table
  ranks$MM1S %>% table
  table(ranks)

  ranks_m <- reshape2::melt(ranks)

  ranks_df <- table(ranks_m$variable, ranks_m$value) %>% as.data.frame.matrix() %>% setNames(paste0("rank_",5:1))

  pheatmap::pheatmap(ranks_df, cluster_rows = F, cluster_cols = F, display_numbers = ranks_df,
                     main = paste0(samplename," rank score cell assignments" ),
                     color=colorRampPalette(c("white","orange", "red"))(50), angle_col = 45,
                     filename = paste0("plots/SNP_analysis/SNP_demultplexing/cell-lines_samples/",samplename,"_rankTable.png"),
                     width = 5, height = 4)

}
