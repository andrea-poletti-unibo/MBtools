library(rhdf5)
library(tidyverse)

H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna(1)/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

# extract barcodes (cells)
barcodes <- h5read(h5filePath, "/assays/dna_variants/ra")
bc_df <- Reduce(cbind, barcodes) %>% as.data.frame()
colnames(bc_df) <- names(barcodes)

bc_df$cell_idx <- 1:nrow(bc_df)
bc_df$sample <- paste0(bc_df$sample_name[1] %>% str_extract("(?<=scDNA_).*(?=_reference)"))
bc_df$cell_sample_id <- paste0(bc_df$sample, "/", bc_df$barcode)



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


# LOAD 4 mutations features sparse matrices

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




######################### DEMULTIPLEXING ###############################

# import correct samples assignments
cell_annot <- data.table::fread("data/umap_deconvolution_clusters.txt")
cell_annot$sample_clust <- recode(cell_annot$cluster,
                                  "1"="KMM1 - clust1",
                                  "2"="OPM2 - clust2",
                                  "3"="KMS11 - clust3",
                                  "4"="RPMI8226 - clust4",
                                  "5"="MM1S - clust5")

write_tsv(cell_annot, "data/umap_deconvolution_clusters.txt")

cell_annot_mix <- cell_annot %>% filter(sample=="Mix") %>% select(ID, sample_clust)


#========= explore AF: use only AF info for demult =========== very good structure

m_input_deconv <- rbind(m_AF_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 825 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using all PASS Mutations AF")



#========= explore DP: use only DP info for demult =========== some structure

m_input_deconv <- rbind(m_DP_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 825 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using all PASS Mutations DP")



#========= explore NGT: use only NGT info for demult =========== very good structure

m_input_deconv <- rbind(m_NGT_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 825 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using all PASS Mutations NGT")



#========= explore GQ: use only GQ info for demult =========== good structure

m_input_deconv <- rbind(m_GQ_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 825 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using PASS Mutations GQ")




#========= approach 1: use both AF and DP info for demult ===========

m_input_deconv <- rbind(m_AF_f, m_DP_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 1650 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using PASS Mutations AF & DP features")




#========= approach 2: use all features AF + DP + GQ + NGT info for demult ===========

m_input_deconv <- rbind(m_AF_f, m_DP_f, m_GQ_f, m_NGT_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 3300 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using all PASS Mutations features (AF + DP + GQ + NGT)")




#========= approach 3: use three features AF + DP + GQ info for demult ===========

m_input_deconv <- rbind(m_AF_f, m_DP_f, m_GQ_f) %>% as.matrix() %>% t
m_input_deconv %>% dim # 2575 features

umap_res <- umap::umap(m_input_deconv)

umap_df <- data.frame(ID= rownames(m_input_deconv), UMAP1=umap_res$layout[,1], UMAP2=umap_res$layout[,2])

umap_df_annot <- left_join(umap_df, cell_annot_mix)

umap_df_annot %>% ggplot(aes(x = UMAP1, y = UMAP2, colour=sample_clust)) +
  geom_point(size=0.8, alpha=0.5) + ggtitle("Demultiplexing using three PASS Mutations features: AF + DP + GQ")


