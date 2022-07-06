
library(rhdf5)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(RColorBrewer)

#================= LOCALIZE DATA ==================

H5_files_dir <- "<H5_FILES_PATH>"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

# missionBio default 4-cell-lines test data
# h5filePath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/test_samples_missionbio/4-cell-lines-AML-CNV.dna.h5"


#============ LOAD AND PREPARE RAW READ COUNT ===========

rc <- h5read(h5filePath, "/assays/dna_read_counts/layers")

df <- rc$read_counts %>% t %>% as.data.frame()

# extract and set colnames(ampl) and rownames(cells)
names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )
cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")

rownames(df) <- cells$barcode
colnames(df) <- paste0("chr",names$CHROM,":", names$start_pos, "-", names$end_pos)

# create annotate chr df for heatmap
chr_ann <- data.frame(chr=names$CHROM %>% factor(levels = str_sort(names$CHROM %>% unique,numeric = T), ordered = T), row.names = colnames(df))


# downsample
set.seed(1)
number_downsampled_cells <- 500
df <- df[sample(1:nrow(df), number_downsampled_cells),]


# pheatmap(df, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
#          # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
#          annotation_col = chr_ann )


#_____ load normalized reads by mosaic (custom missionbio approach)____

norm_reads <- data.table::fread("data/scDNA_Mix_reference_norm_reads.csv")


#========== MISSION BIO normalization / workflow copied from mosaic "normalize_reads" function ============

# ____________ Python code from Mosaic: _____________

# def normalize_reads(self):
#   """
#         Normalize the read counts.
#
#         This normalization is used to determine
#         the ploidy of the cells.
#         """
#
# dat = self.layers["read_counts"]
#
# # Find good barcodes: those with at least 1/10th the reads of the 10th highest barcode
# dat_sum = dat.sum(axis=1)
# s = np.sort(dat_sum)[::-1]
# good = dat_sum > s[10] / 10
#
# normal_counts = deepcopy(dat)
# normal_counts = normal_counts / (normal_counts.mean(axis=1)[:, None] + 1) # dividere per la media
# normal_counts = normal_counts / (np.median(normal_counts[good, :], axis=0)[None, :] + 0.05)
# normal_counts = normal_counts * 2  # For diploid
#
# self.add_layer(NORMALIZED_READS, normal_counts)



#===== reproduce this approach in R ==========

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


#_______ multiply per 2 - for diploidy__________
m4 <- (m3 * 2) %>% round(6)

# comparison with official Mosaic normalized data
nrm <- norm_reads[,-1] %>% as.matrix() %>% round(6)
table(m4 == nrm) # they are exactly the same!


# plot with heatmap
m4 %>% density() %>% plot
m4 %>% as.vector() %>% summary()
m4 %>% as.vector() %>% quantile(probs = seq(0, 1, 0.05))

# breaksList = seq(0, 4, by = 0.2)
# pheatmap(m4, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          filename = "C:/Users/andre/Desktop/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann )
# dev.off()

ndf <- m4 %>% as.data.frame()



# plot with heatmap
m3 %>% density() %>% plot
m3 %>% as.vector() %>% summary()
m3 %>% as.vector() %>% quantile(probs = seq(0, 1, 0.05))



breaksList = seq(0, 2, by = 0.2)
# pheatmap(m3, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          # filename = "C:/Users/andre/Desktop/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann )
# dev.off()


m3l <- apply(m3+0.01, c(1,2), log2)

m3l %>% density() %>% plot
m3l %>% as.vector() %>% summary()
m3l %>% as.vector() %>% quantile(probs = seq(0, 1, 0.05))

breaksList = seq(-3, 3, by = 0.1)
# pheatmap(m3l, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          # filename = "plots/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann )



################################ CELL CLUSTERING ####################################

input_mat <- m3l

PCA <- prcomp(input_mat)

fviz_eig(PCA, ncp = 30)
max_dim <- 15

post_pca <- PCA$x[,1:max_dim]

library(umap)

umap_data <- umap(post_pca)

umap_df <- umap_data$layout %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)


dbscan_clusters <- dbscan::dbscan(umap_df, eps = 0.3)

umap_df$cluster <- dbscan_clusters$cluster %>% as.factor()

umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(colour=cluster))



# create annotate chr df for heatmap
clust_ann <- data.frame(cluster=umap_df$cluster %>% as.factor(), row.names = rownames(df))

# pheatmap(m3l, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          # filename = "plots/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann,
#          annotation_row = clust_ann)

m3l_1 <- m3l[rownames(m3l) %in% rownames(subset(clust_ann, subset = cluster==1)),]

m3l_1 %>% dim

# pheatmap(m3l_1, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          # filename = "plots/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann,
#          annotation_row = clust_ann)




melt_1 <- m3l_1 %>% as.data.frame() %>% rownames_to_column("cell") %>%  reshape2::melt(id.vars="cell")
melt_1s <- melt_1 %>% group_by(variable) %>% summarise(mean= mean(value))
melt_1s$chr <- melt_1s$variable %>% str_extract("chr[1-9]+")
melt_1s$start <- melt_1s$variable %>% str_extract("(?<=:)[0-9]+(?=-)")



melt_1s %>% ggplot(aes(variable, mean, colour=chr)) + geom_point() + facet_grid(~chr)
