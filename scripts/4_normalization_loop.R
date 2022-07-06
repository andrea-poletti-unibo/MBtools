library(rhdf5)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

#================= LOCALIZE DATA ==================
H5_files_dir <- "<H5_FILES_PATH>"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)


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
ampl <- data.table::fread("data/amplicons_list_allStats.txt")
ampl$cov_quality2 <- ifelse(ampl$ampl_median <300 ,"ok", "high_cov_M300")

ampl$FILTER <- ifelse(ampl$mappability_score<0.8 | ampl$cov_quality != "ok", 1,0)
ampl$FILTER %>% table

okampl <- ampl %>% filter(FILTER==0)

ampl_to_keep <- okampl %>% .$ID

okampl$ampl_median %>% plot



# add chromosome arm info to amplicons: Popeye function
ampl$chr <- parse_number(ampl$chr)
ampl <- ampl %>% rename(ID_amplicon="ID")

# use this function to add ARM-LEVEL annotation (Popeye2)
source("scripts/Popeye2.R")
ampl_anno <- ampl %>% Popeye2(refGenome = "hg19", removeXY = T)

ampl_anno_ok <- ampl_anno %>% filter(ID_amplicon %in% ampl_to_keep)

###################### analyze ############################

# (optional) select only single sample
# df <- all_df$Mix

df <- Reduce(rbind, all_df)


#================= FILTERING AMPLICONS =================

# exclude tx amplicons
df2 <- df %>% select(!matches("tx"))

# exclude bad quality amplicons
df3 <- df2 %>% select(matches(ampl_to_keep))


#================= FILTERING CELLS =================
df3 %>% rowSums() %>% plot
df3 %>% rowSums() %>% summary
df3 %>% rowSums() %>% quantile(probs = seq(0, 1, 0.1))

df3 %>% rowSums() %>% density %>% plot


df3 %>% rowSums() %>% `<`(30000) %>% table

row_sd <- apply(df3, 1, sd )

row_sd %>% density() %>%  plot + abline(v = 300) # threshold of 300 by visual inspection of density profile

(row_sd > 300) %>% table

# filter cells based on standard deviation of read counts
good_cells <- row_sd < 300
df3 <- df3[good_cells,]


df3[df3==0] <- 1


#=================== normalization ====================

# #_____ amplicon GC correction ______ NOT RUN - included in amplicon power correction
# col_adjustGC <- okampl$corr_fact
# df4 <- df3 %>% sweep(2, col_adjustGC, `*` )


#_____ amplicon power computation______
ampl_power_df <- data.table::fread("data/amplicon_empirical_powers.txt")

ampl_power <- ampl_power_df$power

df3 %>% names
names(df3) <- paste0(names(df3),"/",ampl_anno_ok$chrarm)

# exclude amplicons on 1q (no reference cells)
df4 <- df3 %>% select(!matches("/1q"))

df5 <- df4 %>% sweep(2, ampl_power, `/` )



#_____ log2 scale conversion ______
df5[df5==0] <- 1
df6 <- log2(df5)


#_____ Row/cell normalizarion______

row_means <- apply(df6, 1, mean )
df7 <- df6 %>% sweep(1, row_means, `-` )

row_sd <- apply(df7, 1, sd )
df8 <- df7 %>% sweep(1, row_sd, `/` )


#===================== Visualization =====================

# df8 <- df

chrGap <- colnames(df8) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()
chr_ann <- data.frame(chr=ampl_anno$chr %>% factor(levels = str_sort(amplicons$CHROM %>% unique,numeric = T), ordered = T),
                      row.names = paste0(colnames(df2),"/",ampl_anno$chrarm))


df8 %>% as.matrix() %>% density() %>% plot + abline(v=c(-2,0,2)) + abline(v=0.2, col="red")
df8 %>% as.matrix() %>% as.vector() %>% summary
df8 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.05))

corr_factor <- 0.2

breaksList = seq(-2, 2, by =0.1) + corr_factor

set.seed(1)
test <- df8[sample(1:nrow(df8),300),]

pheatmap(test,
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         annotation_col = chr_ann )




################################ CELL CLUSTERING ####################################

input_mat <- df8 %>% as.matrix()

# RUN PCA pre UMAP
PCA <- prcomp(input_mat)

# look at scree plot and select the number of components on the knee
factoextra::fviz_eig(PCA, ncp = 30)
max_dim <- 15 # <----------------------select here!
post_pca <- PCA$x[,1:max_dim]

library(umap)

# run UMAP for dimensionality reduction
umap_data <- umap(post_pca)
umap_df <- umap_data$layout %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(alpha=0.33, size=0.5)


# run DBscan for clustering
dbscan_clusters <- dbscan::dbscan(umap_df, eps = 0.5)
umap_df$cluster <- dbscan_clusters$cluster %>% as.factor()
umap_df$index <- 1:nrow(umap_df)


# PLOT UMAP + CLUSTERING RESULT
umap_df %>% ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(colour=cluster), alpha=0.33, size=0.5) +
  ggtitle("UMAP + DBscan clustering results")
ggsave("plots/UMAP_demultiplexing_clusters.png", height = 8, width = 10)

# check concordance for cluster and samples
umap_df$sample <-  rownames(umap_df) %>% str_extract(".*(?=/)")
table(umap_df$sample, umap_df$cluster) # perfect concordance

umap_df %>% ggplot(aes(UMAP1, UMAP2)) +
  geom_point(aes(colour=sample), alpha=0.33, size=0.5) +
  ggtitle("UMAP with samples annotated")
ggsave("plots/UMAP_demultiplexing_samples.png", height = 8, width = 10)


#============ HEATMAPS CN visualization ===========

IDX_clusterd <- umap_df %>% arrange(cluster) %>% .$index

# create annotate chr df for heatmap
clust_ann <- data.frame(cluster=umap_df$cluster %>% as.factor(),
                        sample=umap_df$sample %>% as.factor(),
                        row.names = rownames(df8))

chrGap <- colnames(df8) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()
chr_ann <- data.frame(chr=ampl_anno$chr %>% factor(levels = str_sort(amplicons$CHROM %>% unique,numeric = T), ordered = T),
                      row.names = paste0(colnames(df2),"/",ampl_anno$chrarm))

df8 %>% as.matrix() %>% as.vector() %>% summary
adj <- 0.2

df8 %>% as.matrix() %>% density() %>% plot + abline(v = c(-1.5,0,1.5) + adj, col="blue")
df8 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.1))

breaksList = seq(-1.5,1.5, 0.1) + adj

# set.seed(1)
# pheatmap(df8[sample(1:nrow(df8),500),],
#          cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          cutree_rows = 5,
#          gaps_col = chrGap,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          annotation_row = clust_ann,
#          annotation_col = chr_ann )

df8_c <- df8[IDX_clusterd,]

clusterGap <- clust_ann$cluster %>% table %>% cumsum()
breaksList = seq(-1.5,1.5, 0.1) + adj

ann_colors = list(
  cluster = c("1"="#e41a1c","2"="#377eb8","3"="#4daf4a","4"="#984ea3","5"="#ff7f00"),
  sample = c("KMM1"= '#7fc97f',"KMS11"= '#beaed4', "Mix"='#fdc086',"MM1S"= '#ffff99',"OPM2"= '#386cb0', "RPMI8226"="grey" ))
pheatmap(df8_c,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap,
         gaps_row = clusterGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         filename = "plots/all_myClusters.png",  width = 15, height = 15,
         annotation_row = clust_ann,
         annotation_col = chr_ann,
         annotation_colors = ann_colors)

dev.off()

clust_ann$sample %>% table


# Visualize with gaps for chr arms

cag <- colnames(df8) %>% str_extract("(?<=/)[0-9]+[pq]")
cag2 <- factor(cag, levels = unique(cag) %>% str_sort(numeric = T))

chrArmGap <- cag2 %>% table %>% as.vector() %>% cumsum()

pheatmap(df8_c,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrArmGap,
         gaps_row = clusterGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         filename = "plots/all_myClusters_arms.png",  width = 18, height = 15,
         annotation_row = clust_ann,
         annotation_col = chr_ann,
         annotation_colors = ann_colors)

dev.off()



######################## single samples CN analysis ########################

all_samples <- clust_ann$sample %>% unique %>% as.character()

i <- all_samples[1]

for (i in all_samples){

  message(i)

  df8_1 <- df8 %>% filter(rownames(df8) %in% rownames(clust_ann %>% filter(sample==i)))

  probes <- df8_1 %>% apply(2, median)

  probes %>% plot

  c <- df8_1 %>% colnames() %>% str_extract("chr[0-9]+")
  a <- df8_1 %>% colnames() %>% str_extract("[0-9]+[pq]")
  probes_1_df <- data.frame(probe= df8_1 %>% colnames(),
                            signal=probes,
                            chr= c %>% factor(levels = c %>% unique() %>% str_sort(numeric = T)),
                            arm= a %>% factor(levels = a %>% unique() %>% str_sort(numeric = T)))

  arm_1_df <- probes_1_df %>% group_by(arm) %>%
    summarise(arm_median=median(signal),
              baseline=median(probes_1_df$signal)) %>%
    mutate(adj_signal = arm_median - baseline) %>%
    mutate(status=ifelse(adj_signal >0.1, "amp",
                         ifelse(adj_signal < -0.1, "del", "normal")))


  probes_1_df %>% left_join(arm_1_df, by = "arm") %>%
    ggplot(aes(arm, signal)) +
    geom_boxplot(outlier.shape = NA, alpha=0.5, aes(fill=status)) +
    geom_jitter(height = 0 , width = 0.1, alpha=0.7, size=1) +
    geom_hline(yintercept = median(probes_1_df$signal)) +
    geom_hline(yintercept = c(median(probes_1_df$signal)- 0.1,
                              median(probes_1_df$signal)+ 0.1), linetype=2) +
    ylim(-2,2) +
    ggtitle(i)

  ggsave(paste0("plots/",i,"_CN_profile_normalized.png"), width = 12 , height = 5)

}









####################################################################
#################### MISSION BIO approach ##########################
####################################################################

# MissionBio normalization

#=============== MISSION BIO normalization ===============

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


#_______ multiply per 2 - for diploidy __________
m4 <- (m3 * 2) %>% round(6)


m4 %>% as.matrix() %>% as.vector() %>% summary
adj <- 0.14

m4 %>% as.matrix() %>% density() %>% plot #+ abline(v = c(-1.5,0,1.5) + adj, col="blue")
m4 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.1))

breaksList = seq(0,4,0.1) + adj



chrGap <- colnames(m4) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()
chr_ann <- data.frame(chr=amplicons$CHROM %>% factor(levels = str_sort(amplicons$CHROM %>% unique,numeric = T), ordered = T), row.names = colnames(df))

pheatmap(m4,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap,
         # gaps_row = clusterGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         filename = "C:/Users/andre/Desktop/all_MissionBio_norm.png",  width = 15, height = 15,
         # annotation_row = clust_ann,
         annotation_col = chr_ann)

dev.off()


#================= CELL CLUSTERING ====================

input_mat <- m4 %>% as.matrix()

# RUN PCA pre UMAP
PCA <- prcomp(input_mat)

# look at scree plot and select the number of components on the knee
factoextra::fviz_eig(PCA, ncp = 30)
max_dim <- 15 # <----------------------select here!
post_pca <- PCA$x[,1:max_dim]

library(umap)

# run UMAP for dimensionality reduction
umap_data <- umap(post_pca)
umap_df <- umap_data$layout %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

# run DBscan for clustering
dbscan_clusters <- dbscan::dbscan(umap_df, eps = 0.5)
umap_df$cluster <- dbscan_clusters$cluster %>% as.factor()
umap_df$index <- 1:nrow(umap_df)


# PLOT UMAP + CLUSTERING RESULT
umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(colour=cluster), alpha=0.33, size=0.5)


# check concordance for cluster and samples
umap_df$sample <-  rownames(umap_df) %>% str_extract(".*(?=/)")
table(umap_df$sample, umap_df$cluster) # perfect concordance

umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(colour=sample), alpha=0.33, size=0.5)


IDX_clusterd <- umap_df %>% arrange(cluster) %>% .$index

# create annotate chr df for heatmap
clust_ann <- data.frame(cluster=umap_df$cluster %>% as.factor(),
                        sample=umap_df$sample %>% as.factor(),
                        row.names = rownames(df8))

chrGap <- colnames(df8) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()


df8 %>% as.matrix() %>% as.vector() %>% summary
adj <- 0.14

df8 %>% as.matrix() %>% density() %>% plot + abline(v = c(-1.5,0,1.5) + adj, col="blue")
df8 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.1))

breaksList = seq(-1.5,1.5, 0.1) + adj

set.seed(1)
pheatmap(df8[sample(1:nrow(df8),500),],
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         cutree_rows = 5,
         gaps_col = chrGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         annotation_row = clust_ann,
         annotation_col = chr_ann )



df8_c <- df8[IDX_clusterd,]

clusterGap <- clust_ann$cluster %>% table %>% cumsum()
breaksList = seq(-1.5,1.5, 0.1) + adj

ann_colors = list(
  cluster = c("1"="#e41a1c","2"="#377eb8","3"="#4daf4a","4"="#984ea3","5"="#ff7f00"),
  sample = c("KMM1"= '#7fc97f',"KMS11"= '#beaed4', "Mix"='#fdc086',"MM1S"= '#ffff99',"OPM2"= '#386cb0', "RPMI8226"="grey" ))
pheatmap(df8_c,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap,
         gaps_row = clusterGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         filename = "C:/Users/andre/Desktop/all_myClusters.png",  width = 15, height = 15,
         annotation_row = clust_ann,
         annotation_col = chr_ann,
         annotation_colors = ann_colors)

dev.off()

clust_ann$sample %>% table


##### singleclust ######

df8_1 <- df8 %>% filter(rownames(df8) %in% rownames(clust_ann %>% filter(cluster==2)))

probes <- df8_1 %>% apply(2, median)

probes %>% plot







