
library(rhdf5)
library(tidyverse)


H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# one sample
h5filePath <- list_h5_files[1]
h5filePath

rhdf5::h5ls(h5filePath)

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

var_df <- var_df %>% select(id, CHROM, POS, REF, ALT, QUAL, ado_gt_cells, ado_rate, amplicon, filtered)

var_df %>% group_by(filtered) %>% summarise(mean= mean(QUAL %>% as.numeric()),
                                            median=median(QUAL %>% as.numeric()))

var_df$filtered %>% table

IDX <- var_df$filtered == "00"

var_pass <- var_df[IDX,]

h5f = H5Fopen(h5filePath)

i <- h5f&"assays/dna_variants/layers"
i

m_AF <- h5f$"assays/dna_variants/layers/AF" %>% as.data.frame()
m_AF_f <- m_AF[IDX,]
names(m_AF_f) <- bc_df$cell_sample_id
rownames(m_AF_f) <- var_pass$id
rm(m_AF)

m_DP <- h5f$"assays/dna_variants/layers/DP" %>% as.data.frame()
m_DP_f <- m_DP[IDX,]
names(m_DP_f) <- bc_df$cell_sample_id
rownames(m_DP_f) <- var_pass$id
rm(m_DP)


m_GQ <- h5f$"assays/dna_variants/layers/GQ" %>% as.data.frame()
m_GQ_f <- m_GQ[IDX,]
names(m_GQ_f) <- bc_df$cell_sample_id
rownames(m_GQ_f) <- var_pass$id
rm(m_GQ)


m_NGT <- h5f$"assays/dna_variants/layers/NGT" %>% as.data.frame()
m_NGT_f <- m_NGT[IDX,]
names(m_NGT_f) <- bc_df$cell_sample_id
rownames(m_NGT_f) <- var_pass$id
rm(m_NGT)


m_RGQ <- h5f$"assays/dna_variants/layers/RGQ" %>% as.data.frame()
m_RGQ_f <- m_RGQ[IDX,]
names(m_RGQ_f) <- bc_df$cell_sample_id
rownames(m_RGQ_f) <- var_pass$id
rm(m_RGQ)


#===== merge all =====

m_AF_fl <- m_AF_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "VAF")
m_DP_fl <- m_DP_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "DP")
m_GQ_fl <- m_GQ_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "GQ")
m_NGT_fl <- m_NGT_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "NGT")
m_RGQ_fl <- m_RGQ_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "RGQ")


m_ALL_fl <- left_join(m_AF_fl, m_DP_fl) %>% left_join(m_GQ_fl) %>% left_join(m_NGT_fl) %>% left_join(m_RGQ_fl)


snpdata <- data.table::fread("C:/Users/andre/Dropbox (Personal)/Boston_DFCI_Broad_work/Missionbio_scDNA/SNP_selection/SNP6probeset_dbSNP151_annot.txt")

snpdata$snpID <- paste0(snpdata$chrom,":",snpdata$start + 1,":", snpdata$observed)

SNP_idx <- rownames(m_AF_f) %in% snpdata$snpID

SNP_idx %>% table


m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by=c("varID"= "snpID"))

m_ALL_onlySNP_fl <- m_ALL_SNP_fl %>% filter(!is.na(rsID))
m_ALL_onlySNP_fl$cell_n <- m_ALL_onlySNP_fl$cellID %>% str_extract("(?<=_)[0-9]+(?=_)") %>% as.numeric()


SNPs_VAF_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(cell_n, VAF)) + geom_point(alpha=0.3) + facet_wrap(~varID)
SNPs_VAF_plot

SNPs_VAFdens_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(VAF)) + geom_density() + facet_wrap(~varID, scales = "free")
SNPs_VAFdens_plot







############################ loop all samples ##################################

library(rhdf5)
library(tidyverse)

H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"
list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# # FROM SNP5 LIST
# snpdata2 <- data.table::fread("C:/Users/andre/Dropbox (Personal)/Boston_DFCI_Broad_work/Missionbio_scDNA/SNP_selection/SNP6probeset_dbSNP151_annot.txt")
# snpdata2$snpID <- paste0(snpdata2$chrom,":",snpdata2$start + 1,":", snpdata2$observed)

# QUERY BIOMART
snpdata <- data.table::fread("data/SNPs_list_in_panel.txt")
snpdata$snpID <- paste0("chr",snpdata$chr_name,":",snpdata$chrom_start ,":", snpdata$allele)

# # sanity check
# j <- inner_join(snpdata, snpdata2, by= c("refsnp_id"= "rsID" ))
# (j$snpID.x == j$snpID.y) %>% table


# one sample test
h5filePath <- list_h5_files[1]
h5filePath

for( SAMP in list_h5_files){

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

  #______approach no filter___
  # var_df_snp <- inner_join(var_df, snpdata, by = c("id" = "snpID"))
  # var_df_snp$filtered %>% table
  # var_df_snp$QUAL <- as.numeric(var_df_snp$QUAL)
  # all_SNP_ids <- var_df_snp$id
  #____________________________

  var_df %>% group_by(filtered) %>% summarise(mean= mean(QUAL %>% as.numeric()),
                                              median=median(QUAL %>% as.numeric()))

  var_df$filtered %>% table

  IDX <- var_df$filtered == "00"

  var_pass <- var_df[IDX,]

  h5f = H5Fopen(h5filePath)

  i <- h5f&"assays/dna_variants/layers"
  i

  #_____ AF _______
  m_AF <- h5f$"assays/dna_variants/layers/AF" %>% as.data.frame()
  m_AF_f <- m_AF[IDX,]
  names(m_AF_f) <- bc_df$cell_sample_id
  rownames(m_AF_f) <- var_pass$id

  # #approach no filter
  # names(m_AF) <- bc_df$cell_sample_id
  # m_AF$id <- var_df$id
  # var_df_snp2 <- inner_join(var_df_snp, m_AF)
  # var_df_snp2_l <- var_df_snp2 %>% reshape2::melt(id.vars=names(var_df_snp), variable.name="cellID", value.name = "VAF" )

  rm(m_AF)

  #_____ DP _______
  m_DP <- h5f$"assays/dna_variants/layers/DP" %>% as.data.frame()
  m_DP_f <- m_DP[IDX,]
  names(m_DP_f) <- bc_df$cell_sample_id
  rownames(m_DP_f) <- var_pass$id

  # #approach no filter
  # names(m_DP) <- bc_df$cell_sample_id
  # m_DP$id <- var_df$id
  # var_df_snp4 <- inner_join(var_df_snp, m_DP)
  # var_df_snp4_l <- var_df_snp4 %>% reshape2::melt(id.vars=names(var_df_snp), variable.name="cellID", value.name = "DP" )

  rm(m_DP)

  #_____ GQ _______
  m_GQ <- h5f$"assays/dna_variants/layers/GQ" %>% as.data.frame()
  m_GQ_f <- m_GQ[IDX,]
  names(m_GQ_f) <- bc_df$cell_sample_id
  rownames(m_GQ_f) <- var_pass$id
  rm(m_GQ)

  #_____ NGT _______
  m_NGT <- h5f$"assays/dna_variants/layers/NGT" %>% as.data.frame()
  m_NGT_f <- m_NGT[IDX,]
  names(m_NGT_f) <- bc_df$cell_sample_id
  rownames(m_NGT_f) <- var_pass$id

  # #approach no filter
  # names(m_NGT) <- bc_df$cell_sample_id
  # m_NGT$id <- var_df$id
  # var_df_snp3 <- inner_join(var_df_snp, m_NGT)
  # var_df_snp3_l <- var_df_snp3 %>% reshape2::melt(id.vars=names(var_df_snp), variable.name="cellID", value.name = "NGT" )

  rm(m_NGT)

  #_____ RGQ _______
  m_RGQ <- h5f$"assays/dna_variants/layers/RGQ" %>% as.data.frame()
  m_RGQ_f <- m_RGQ[IDX,]
  names(m_RGQ_f) <- bc_df$cell_sample_id
  rownames(m_RGQ_f) <- var_pass$id
  rm(m_RGQ)


  #===== long versions =====

  m_AF_fl <- m_AF_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "VAF")
  m_DP_fl <- m_DP_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "DP")
  m_GQ_fl <- m_GQ_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "GQ")
  m_NGT_fl <- m_NGT_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "NGT")
  m_RGQ_fl <- m_RGQ_f %>% rownames_to_column(var = "varID") %>% reshape2::melt(id.vars="varID", variable.name="cellID", value.name = "RGQ")


  #====== merge all =======

  m_ALL_fl <- left_join(m_AF_fl, m_DP_fl) %>% left_join(m_GQ_fl) %>% left_join(m_NGT_fl) %>% left_join(m_RGQ_fl)


  # # approach no filter
  # var_df_snp_all_l <- left_join(var_df_snp2_l, var_df_snp3_l) %>% left_join(var_df_snp4_l)
  # var_df_snp_all_l$cell_n <- var_df_snp_all_l$cellID %>% str_extract("(?<=_)[0-9]+(?=_)") %>% as.numeric()
  # var_df_snp_all_l$filtered <- var_df_snp_all_l$filtered %>% dplyr::recode("01" = "filtered", "00" = "pass")

  #======= merge with SNP data =========

  m_ALL_fl$var_pos <- m_ALL_fl$varID %>% str_extract("chr[0-9]+:[0-9]+")
  snpdata$var_pos <- snpdata$snpID %>% str_extract("chr[0-9]+:[0-9]+")


  # m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by=c("varID"= "snpID"))
  m_ALL_SNP_fl <- left_join(m_ALL_fl, snpdata, by="var_pos")

  m_ALL_onlySNP_fl <- m_ALL_SNP_fl %>% filter(!is.na(refsnp_id))
  m_ALL_onlySNP_fl$cell_n <- m_ALL_onlySNP_fl$cellID %>% str_extract("(?<=_)[0-9]+(?=_)") %>% as.numeric()

  m_ALL_onlySNP_fl$varID %>% unique()


  # #======== PLOTS ============
  #
  # dir.create("plots/SNP_analysis",showWarnings = F)
  #
  #
  # SNPs_VAF_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(cell_n, VAF, colour= NGT %>% as.factor())) +
  #   geom_point(alpha=0.3) + facet_wrap(~varID)+ theme(legend.position = "bottom") + ggtitle(samplename, subtitle = "SNPs VAF plot")
  # SNPs_VAF_plot
  # ggsave(paste0("plots/SNP_analysis/",samplename, "SNP_VAF_plot.png" ), width = 30, height = 30)
  #
  #
  # SNPs_VAFdens_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(VAF)) +
  #   geom_density(alpha=0.5) + facet_wrap(~varID, scales = "free_y") + theme(legend.position = "bottom") + ggtitle(samplename,  subtitle = "SNPs VAF density plot")
  # SNPs_VAFdens_plot
  # ggsave(paste0("plots/SNP_analysis/",samplename, "SNP_VAFdens_plot.png" ), width = 30, height = 30)
  #
  #
  # SNPs_VAFhist_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(VAF, fill=NGT %>% as.factor())) +
  #   geom_histogram() + facet_wrap(~varID, scales = "free_y") + theme(legend.position = "bottom") + ggtitle(samplename,  subtitle = "SNPs VAF histogram plot")
  # SNPs_VAFhist_plot
  # ggsave(paste0("plots/SNP_analysis/",samplename, "SNP_VAFhist_plot.png" ), width = 30, height = 30)


  #approach no filter

  # SNPs_VAF_plot <- var_df_snp_all_l %>% ggplot(aes(cell_n, VAF, colour= NGT %>% as.factor())) +
  #   geom_point(alpha=0.3) + facet_wrap(filtered~id)+ theme(legend.position = "bottom")
  # SNPs_VAF_plot
  # ggsave(paste0("plots/SNP_analysis/",samplename, "SNP_VAF_plot.png" ), width = 12, height = 12)
  #
  # SNPs_VAFdens_plot <- m_ALL_onlySNP_fl %>% ggplot(aes(VAF, fill=NGT %>% as.factor())) +
  #   geom_density(alpha=0.5) + facet_wrap(filtered~id, scales = "free") + theme(legend.position = "bottom")
  # SNPs_VAFdens_plot
  # ggsave(paste0("plots/SNP_analysis/",samplename, "SNP_VAFdens_plot.png" ), width = 12, height = 12)




  #__________ generate SAMPLE profile to use as a ground truth for SNP demultiplexing __________

  summ_profile <- m_ALL_onlySNP_fl %>% group_by(snpID) %>% summarise(
    Het= sum(NGT==1),
    Hom= sum(NGT==2),
    Ref= sum(NGT==0),
    Na = sum(NGT==3),
    tot_cell=n())

  summ_profile$Prop_het <- (summ_profile$Het / summ_profile$tot_cell) %>% round(4)
  summ_profile$Prop_het %>% sort %>%  plot

  summ_profile$call_HET <- ifelse(summ_profile$Prop_het>0.8,1,0)
  summ_profile$call_HET %>% table


  summ_profile$Prop_Hom <- (summ_profile$Hom / summ_profile$tot_cell) %>% round(4)
  summ_profile$Prop_ref <- (summ_profile$Ref / summ_profile$tot_cell) %>% round(4)
  summ_profile$callNGT <- case_when(summ_profile$Prop_het > 0.8 ~ 1,
                                    summ_profile$Prop_Hom > 0.8 ~ 2,
                                    summ_profile$Prop_ref > 0.6 ~ 0,
                                    T ~ 3 )



  dir.create("data/SNP_profiles/", showWarnings = F)
  write_tsv(summ_profile, paste0("data/SNP_profiles/", samplename, "_SNP_profile.txt" ))

}


##############################################################################################



example <- m_ALL_onlySNP_fl %>% filter(cell_n==1)

example %>% ggplot(aes(snpID, VAF, colour=NGT %>% as.factor)) + geom_point()


prof <- fread("data/SNP_profiles/RPMI8226_SNP_profile.txt")
prof$Prop_Hom <- (prof$Hom / prof$tot_cell) %>% round(4)
prof$Prop_ref <- (prof$Ref / prof$tot_cell) %>% round(4)
prof$callNGT <- case_when(prof$Prop_het > 0.8 ~ 1, prof$Prop_Hom > 0.8 ~ 2, prof$Prop_ref > 0.6 ~ 0, T ~ 3 )
df <- left_join(example, prof, by = "snpID")

table(df$NGT == df$callNGT)

# dfHet <- df %>% filter(NGT==1)
# table(dfHet$NGT == dfHet$call_HET)



prof <- fread("data/SNP_profiles/KMM1_SNP_profile.txt")
prof$Prop_Hom <- (prof$Hom / prof$tot_cell) %>% round(4)
prof$Prop_ref <- (prof$Ref / prof$tot_cell) %>% round(4)
prof$callNGT <- case_when(prof$Prop_het > 0.8 ~ 1, prof$Prop_Hom > 0.8 ~ 2, prof$Prop_ref > 0.6 ~ 0, T ~ 3 )
df <- left_join(example, prof, by = "snpID")

table(df$NGT == df$callNGT)

# dfHet <- df %>% filter(NGT==1)
# table(dfHet$NGT == dfHet$call_HET)




prof <- fread("data/SNP_profiles/KMS11_SNP_profile.txt")
prof$Prop_Hom <- (prof$Hom / prof$tot_cell) %>% round(4)
prof$Prop_ref <- (prof$Ref / prof$tot_cell) %>% round(4)
prof$callNGT <- case_when(prof$Prop_het > 0.8 ~ 1, prof$Prop_Hom > 0.8 ~ 2, prof$Prop_ref > 0.6 ~ 0, T ~ 3 )
df <- left_join(example, prof, by = "snpID")

table(df$NGT == df$callNGT)

# dfHet <- df %>% filter(NGT==1)
# table(dfHet$NGT == dfHet$call_HET)




prof <- fread("data/SNP_profiles/MM1S_SNP_profile.txt")
prof$Prop_Hom <- (prof$Hom / prof$tot_cell) %>% round(4)
prof$Prop_ref <- (prof$Ref / prof$tot_cell) %>% round(4)
prof$callNGT <- case_when(prof$Prop_het > 0.8 ~ 1, prof$Prop_Hom > 0.8 ~ 2, prof$Prop_ref > 0.6 ~ 0, T ~ 3 )
df <- left_join(example, prof, by = "snpID")

table(df$NGT == df$callNGT)

# dfHet <- df %>% filter(NGT==1)
# table(dfHet$NGT == dfHet$call_HET)





prof <- fread("data/SNP_profiles/OPM2_SNP_profile.txt")
prof$Prop_Hom <- (prof$Hom / prof$tot_cell) %>% round(4)
prof$Prop_ref <- (prof$Ref / prof$tot_cell) %>% round(4)
prof$callNGT <- case_when(prof$Prop_het > 0.8 ~ 1, prof$Prop_Hom > 0.8 ~ 2, prof$Prop_ref > 0.6 ~ 0, T ~ 3 )
df <- left_join(example, prof, by = "snpID")

table(df$NGT == df$callNGT)

# dfHet <- df %>% filter(NGT==1)
# table(dfHet$NGT == dfHet$call_HET)

