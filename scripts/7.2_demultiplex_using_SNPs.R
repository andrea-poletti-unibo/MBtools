
library(tidyverse)
library(data.table)
library(rhdf5)


H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"
list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

snpdata <- data.table::fread("data/SNPs_list_in_panel.txt")
snpdata$snpID <- paste0("chr",snpdata$chr_name,":",snpdata$chrom_start ,":", snpdata$allele)


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
snpdata$var_pos <- snpdata$snpID %>% str_extract("chr[0-9]+:[0-9]+")


# m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by=c("varID"= "snpID"))
m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by="var_pos")

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


cells_res_df$ratio_T_F <- cells_res_df$True/cells_res_df$False

cells_res_df$ratio_T_F %>% plot

table(cells_res_df$ratio_T_F > 1) # set threshold at 1

cells_res_df$matched <- ifelse(cells_res_df$ratio_T_F > 1, "match", "no_match")


write_tsv(cells_res_df, "cells_matching_SNPprofiles.txt")


#================ validate assignments with UMAP on CNA data ===================

cells_matched <- cells_res_df %>% filter(matched == "match")

cells_matched$cell_id <- cells_matched$id %>% str_remove("Mix_[0-9]+_")

solution <- fread("data/umap_deconvolution_clusters.txt")
solution_mix <- solution %>% filter(sample == "Mix")

solution_mix$cell_id <- solution_mix$ID %>% str_remove("Mix/")


merge <- left_join(solution_mix, cells_matched)

merge$samples_in_clust <- merge$sample_clust %>% str_remove(" - clust.")

table(merge$sample_match == merge$samples_in_clust)
1588/(1588+57) # 96.5 % concordance!!!!

table(merge$sample_match, merge$samples_in_clust)

merge %>% ggplot(aes(UMAP1, UMAP2, colour=sample_match)) + geom_point()











