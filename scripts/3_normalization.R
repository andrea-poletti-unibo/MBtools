library(rhdf5)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)



#================= LOCALIZE DATA ==================
list_h5_files <- list.files("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/", pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath


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


#_____ load amplicons table with stat ______
ampl <- data.table::fread("amplicons_list_allStats.txt")
ampl$cov_quality2 <- ifelse(ampl$ampl_median <300 ,"ok", "high_cov_M300")

ampl$FILTER <- ifelse(ampl$mappability_score<0.8 | ampl$cov_quality != "ok" | ampl$cov_quality2 != "ok", 1,0)
ampl$FILTER %>% table

okampl <- ampl %>% filter(FILTER==0)

ampl_to_keep <- okampl %>% .$ID

okampl$ampl_median %>% plot


###################### analyze sample 3 (mix) ############################

#================= FILTERING AMPLICONS =================

# exclude tx amplicons
df2 <- df %>% select(!matches("tx"))

# exclude bad amplicons
df3 <- df2 %>% select(matches(ampl_to_keep))

#================= FILTERING CELLS =================
df3 %>% rowSums() %>% plot
df3 %>% rowSums() %>% summary
df3 %>% rowSums() %>% quantile(probs = seq(0, 1, 0.1))

df3 %>% rowSums() %>% `<`(30000) %>% table

row_sd <- apply(df3, 1, sd )

row_sd %>% density() %>%  plot

(row_sd > 200) %>% table

good_cells <- row_sd < 200

df3 <- df3[good_cells,]


df3[df3==0] <- 1



#========= visualize raw unnormalized data ============
# create chromosome annotation
chr_ann <- data.frame(chr=names$CHROM %>% factor(levels = str_sort(names$CHROM %>% unique,numeric = T), ordered = T), row.names = colnames(df))

chrGap <- colnames(df3) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()

df3 %>% as.matrix() %>% density() %>% plot
df3 %>% as.matrix() %>% as.vector() %>% summary
df3 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.05))

df3[df3==0] <- 1


df3l <- apply(df3, c(1,2), log)

df3l %>% as.matrix() %>% density() %>% plot
df3l %>% as.matrix() %>% as.vector() %>% summary
df3l %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.05))


breaksList = seq(-2, 2, by =0.1)
pheatmap(df3l[sample(1:nrow(df3l),300),], 
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         # filename = "C:/Users/andre/Desktop/hm.png", width = 15, height = 10,
         annotation_col = chr_ann )


#=================== normalization ====================

#_____ amplicon GC correction ______
col_adjustGC <- okampl$corr_fact

df4 <- df3 %>% sweep(2, col_adjustGC, `*` )


#_____ amplicon ability computation______
ampl_medians <- apply(df4, 2, median)
ampl_ability <- ampl_medians/ median(ampl_medians)

df5 <- df4 %>% sweep(2, ampl_ability, `/` )


#_____ log2 scale conversion ______
df5[df5==0] <- 1
df6 <- log2(df5)


#_____ Row/cell normalizarion______

row_means <- apply(df6, 1, mean )
df7 <- df6 %>% sweep(1, row_means, `-` )

row_sd <- apply(df7, 1, sd )
df8 <- df7 %>% sweep(1, row_sd, `/` )


################################ CELL CLUSTERING ####################################

input_mat <- df8 %>% as.matrix()

PCA <- prcomp(input_mat) 

factoextra::fviz_eig(PCA, ncp = 30)
max_dim <- 15

post_pca <- PCA$x[,1:max_dim]

library(umap)

umap_data <- umap(post_pca)

umap_df <- umap_data$layout %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)


dbscan_clusters <- dbscan::dbscan(umap_df, eps = 0.5)

umap_df$cluster <- dbscan_clusters$cluster %>% as.factor()

umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(colour=cluster))

umap_df$index <- 1:nrow(umap_df)

IDX_clusterd <- umap_df %>% arrange(cluster) %>% .$index

# create annotate chr df for heatmap
clust_ann <- data.frame(cluster=umap_df$cluster %>% as.factor(), row.names = rownames(df8))

# pheatmap(m3l, cluster_rows = T, cluster_cols = F,
#          show_rownames = F, show_colnames = F,
#          color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
#          # filename = "C:/Users/andre/Desktop/hm.png", width = 15, height = 10,
#          annotation_col = chr_ann,
#          annotation_row = clust_ann)


chrGap <- colnames(df8) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()

df8 %>% as.matrix() %>% density() %>% plot + abline(v = c(-1.5,0.25,1.5), col="blue")



df8 %>% as.matrix() %>% as.vector() %>% summary
df8 %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.1))

breaksList = seq(-1.5,1.5, 0.1)

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
breaksList = seq(-1.5,1.2, 0.1)
pheatmap(df8_c, 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap, 
         gaps_row = clusterGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         filename = "C:/Users/andre/Desktop/myClusters2.png",  width = 15, height = 10,
         annotation_row = clust_ann,
         annotation_col = chr_ann)
dev.off()









#__________ single cluster visualization _____________
df8_1 <- df8[rownames(df8) %in% rownames(subset(clust_ann, subset = cluster==1)),]

pheatmap(df8_1, 
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         annotation_row = clust_ann,
         annotation_col = chr_ann )


ampl_track_1 <- df8_1 %>% apply(2, median)

library(DNAcopy)

CNAdat <- CNA(ampl_track_1, 
              chrom = okampl$chr,
              maploc = okampl$start)

seg <- segment(CNAdat, alpha = 0.001)

plotSample(seg, xmaploc = T, col=c("red","blue"))

