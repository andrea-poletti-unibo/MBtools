
################ Please Note!: BENCHMARK version of the demult analysis on 1/3 of amplicons ####################

library(tidyverse)
library(data.table)
library(rhdf5)
library(gridExtra)


H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"
list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# USE THE BENCHMARK PANEL SNP list!
snpdataB <- data.table::fread("data/SNPs_list_in_panelBenchmark.txt")



#=============== select MIX sample ============

h5filePath <- list_h5_files[3]
h5filePath

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

m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdataB, by="var_pos")

m_ALL_onlySNP_fl <- m_ALL_SNP_fl %>% filter(!is.na(refsnp_id))
m_ALL_onlySNP_fl$cell_n <- m_ALL_onlySNP_fl$cellID %>% str_extract("(?<=_)[0-9]+(?=_)") %>% as.numeric()



########################### matching ################################

profiles_files <- list.files("data/SNP_profiles/", full.names = T)
profiles_list <- lapply(profiles_files,fread)
names(profiles_list) <- profiles_files %>% str_remove("data/SNP_profiles/") %>% str_remove("_SNP_profile.txt")

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


#________ create score 1: T/F ratio ________

cells_res_df$ratio_T_F <- cells_res_df$True/cells_res_df$False

table(cells_res_df$ratio_T_F > 1.5) # set threshold at 1.5
cutoff <- 1.5


mytab <- table(cells_res_df$ratio_T_F > cutoff) %>% as.data.frame() %>% setNames(c("matched","n_cells"))

cells_res_df %>% arrange(ratio_T_F) %>%
  ggplot(aes(1:nrow(cells_res_df), ratio_T_F)) +
  geom_point() + ggtitle("Knee plot for score1 (True/False) - BENCHMARK only 150 amplicons") +
  geom_hline(yintercept = cutoff, linetype=2) +
  annotation_custom(tableGrob(mytab), xmin=100, xmax=3000, ymin=2, ymax=3)

ggsave("plots/SNP_analysis/SNP_demultplexing/BMK_knee_plot_score1_T-F.png", width = 6, height = 6)

cells_res_df$matched_score1 <- ifelse(cells_res_df$ratio_T_F > cutoff, "match", "no_match")


# SAVE
write_tsv(cells_res_df, "data/BMK_cells_matching_SNPprofiles.txt")



#============== knee plots faceted per cell-line =========

cells_res_df %>% arrange(ratio_T_F) %>%
  ggplot(aes(1:nrow(cells_res_df), ratio_T_F)) + geom_hline(yintercept = cutoff, linetype=2, size = 0.5) +
  geom_point(aes(colour=matched_score1), size=0.8) + ggtitle("Knee plot for score1 (True/False) faceted per cell-line - BMK") +
  facet_wrap(~sample_match)

ggsave("plots/SNP_analysis/SNP_demultplexing/BMK_knee_plot_score1_facet_PER_CELL-LINE.png", width = 6, height = 6)


#================ validate assignments with UMAP on CNA data ===================

#____ score 1 _____

cells_matched <- cells_res_df %>% filter(matched_score1 == "match")
cells_matched$cell_id <- cells_matched$id %>% str_remove("Mix_[0-9]+_")

solution <- fread("data/umap_deconvolution_clusters.txt")
solution_mix <- solution %>% filter(sample == "Mix")

solution_mix$cell_id <- solution_mix$ID %>% str_remove("Mix/")


merge <- left_join(solution_mix, cells_matched)
merge$samples_in_clust <- merge$sample_clust %>% str_remove(" - clust.")

table(merge$sample_match == merge$samples_in_clust)
1317/(1317+202) # 86.7 % concordance

table(merge$sample_match, merge$samples_in_clust)

tb <- table(merge$sample_match, merge$samples_in_clust) %>% as.data.frame.matrix()
merge %>% ggplot(aes(UMAP1, UMAP2, colour=sample_match)) + geom_point(size=0.8, alpha=0.5) +
  ggtitle("clusters of cells colored by SNP profile matching (scoreT/F > 1.5) - BMK analysis",
          subtitle = "(86.7% concordance)") +
  annotation_custom(tableGrob(tb,theme=ttheme_default(base_size = 10)), xmin=3, xmax=11, ymin=-8, ymax=-2)


ggsave("plots/SNP_analysis/SNP_demultplexing/SNP_assignment_in_CNA_clusters.png", width = 11, height = 11)


#================ single cell level analysis ==================

nsc <- cells_res_df %>% group_by(id) %>% summarise()

sc <- cells_res_df %>% group_by(id) %>% group_map(~ sort(.x$ratio_T_F, decreasing = T))
sc_df <- Reduce(rbind, sc) %>% as.data.frame(row.names = nsc$id) %>% rownames_to_column(var = "id")


names(sc_df) <- c("id", paste0("score_best_",1:5))

sc_df$diff2best <- sc_df$score_best_1 - sc_df$score_best_2

cut_doub <- 0.2
sc_df %>% arrange(diff2best) %>%  ggplot(aes(1:nrow(sc_df), diff2best)) +
  geom_point() + geom_hline(yintercept = cut_doub)


sc_df$doublet_suspect <- ifelse(sc_df$diff2best< cut_doub, "doublet_suspect","ok")

sc_df$doublet_suspect %>% table()

sc_df$cell_id <- sc_df$id %>% str_remove("Mix_[0-9]+_")


merge2 <- left_join(merge, sc_df, by="cell_id")

merge2 %>% ggplot(aes(UMAP1, UMAP2, colour=doublet_suspect)) + geom_point(size=0.8, alpha=0.5) +
  ggtitle("clusters of cells colored by doublet suspect - BENCHMARK ANALYSIS")

merge2$cell_concordance <- ifelse(merge2$sample_match == merge2$samples_in_clust, "concordant", "discordant")


table(merge2$cell_concordance, merge2$doublet_suspect)
merge2 %>% ggplot(aes(UMAP1, UMAP2, colour=doublet_suspect, shape=cell_concordance)) + geom_point(size=0.8, alpha=0.5) +
  ggtitle("clusters of cells colored by SNP profile matching - score1 > 1.5")


# cluster plot clean from suspect doublets
merge2_c <- merge2 %>% filter(doublet_suspect=="ok")
table(merge2_c$sample_match == merge2_c$samples_in_clust )
1288/(1288+176)

merge2_c %>% ggplot(aes(UMAP1, UMAP2, colour=sample_match)) + geom_point(size=0.8, alpha=0.5) +
  ggtitle("clusters of cells colored by SNP profile matching (scoreT/F > 1.5) - BMK analysis",
          subtitle = "CLENANED from doublet suspect cells: 88.0% concordance")

ggsave("plots/SNP_analysis/SNP_demultplexing/BMK_SNP_assignment_in_CNA_clusters_CLEAN_DOUB.png", width = 7, height = 7)


mytab <- table(merge2$cell_concordance, merge2$doublet_suspect) %>% as.data.frame.matrix() %>% setNames(c("doublet suspect", "no doublet"))
merge2 %>% arrange(diff2best) %>%  ggplot(aes(1:nrow(merge2), diff2best, colour= cell_concordance)) +
  geom_point() + geom_hline(yintercept = cut_doub, linetype=2) +
  annotation_custom(tableGrob(mytab), xmin=100, xmax=1400, ymin=1.5, ymax=2.5)

ggsave("plots/SNP_analysis/SNP_demultplexing/BMK_cell_concordance_and_doublet_score.png", width = 8, height = 8)


