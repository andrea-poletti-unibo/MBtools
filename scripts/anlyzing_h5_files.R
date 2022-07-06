# BiocManager::install("rhdf5")


library(rhdf5)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(RColorBrewer)



list_h5_files <- list.files("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/", pattern = ".*h5$", full.names = T)

outpath <- "C:/Users/andre/Dropbox (Personale)/Boston_DFCI_Broad_work/Missionbio_scDNA/plots_custom_Rscript/"

dir.create(outpath,showWarnings = F)


i <- list_h5_files[3]

for(i in list_h5_files){
  
  message(which(i==list_h5_files)," - ", i)

  h5filePath <- i
  
  sampleName <- str_extract(i, "scDNA_.*_reference")
  
  rc <- h5read(h5filePath, "/assays/dna_read_counts/layers" )
  
  df <- rc$read_counts %>% as.data.frame()
  
  names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )
  
  cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")
  
  colnames(df) <- cells$barcode
  row.names(df) <- paste0(names$CHROM,"_", names$start_pos)
  
  row_ann <- data.frame(chr=names$CHROM %>% factor(levels = str_sort(names$CHROM %>% unique,numeric = T), ordered = T), row.names = rownames(df))
  scaled_df <-  scale(df)
  pheatmap::pheatmap(scaled_df, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F,
                     annotation_row = row_ann, filename = paste0(outpath, sampleName,".pdf") )
  
}


i <- list_h5_files[3]

h5filePath <- i


h5filePath <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/test_samples_missionbio/4-cell-lines-AML-CNV.dna.h5"

h5ls(h5filePath)

rc <- h5read(h5filePath, "/assays/dna_read_counts/layers")


df <- rc$read_counts %>% t %>% as.data.frame() 

set.seed(1)
df <- df[sample(1:nrow(df), 400),]


names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )
cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")

rownames(df) <- cells$barcode
colnames(df) <- paste0("chr",names$CHROM,"_", names$start_pos)




pheatmap::pheatmap(df, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                   # breaks = breaksList,
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )




# normlize rows (cells) dividing by rowsum (total reads)
m <- sweep(df, 1 ,rowSums(df),`/`)



# normalize cols (amplicons) dividing by their median
medAmp <- apply(m, 2, mean)
m2 <- sweep(m, 2, medAmp, `/`)


m2 <- scale(m)


m2 %>% density(na.rm = T) %>% plot


m3 <- apply(as.matrix(m2) + 1, c(1,2), log2)

breaksList = seq(0, 4, by = 0.2)
pheatmap::pheatmap(m2, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                   breaks = breaksList,
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )



#================ Z score filtering ==================
z_cells <- (rowSums(df) - mean(rowSums(df)))/sd(rowSums(df))

z_cells %>% plot

outlier_cells <- abs(z_cells) > 3 
outlier_cells %>% table


z_amp <- (colSums(df) - mean(colSums(df)))/sd(colSums(df))

z_amp %>% plot

outlier_amp <- abs(z_amp) > 3 
outlier_amp %>% table

cdf <- df[!outlier_cells, !outlier_amp]

cdf %>% as.matrix() %>% density() %>% plot

chr_ann <- data.frame(chr=names$CHROM[!outlier_amp] %>% factor(levels = str_sort(names$CHROM[!outlier_amp] %>% unique,numeric = T), ordered = T), row.names = colnames(cdf))


pheatmap::pheatmap(cdf, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                   # breaks = breaksList,
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )

dev.off()



breaksList = seq(0, 1500, by = 10)
pheatmap::pheatmap(cdf, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   # color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                   # breaks = breaksList,
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )




# normalize cols (amplicons) dividing by their median
medAmp <- apply(cdf, 2, mean, trim=0.0025)
m2 <- sweep(cdf, 2, medAmp, `/`)


# normlize rows (cells) dividing by rowsum (total reads)
m <- sweep(cdf, 1 ,rowSums(cdf),`/`)





breaksList = seq(0, 20, by = 1)
pheatmap::pheatmap(m2, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                   breaks = breaksList,
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )




# log scaling
df1 <- df
df1[df1==0] <- 1
log2df1 <- df1 %>% apply(c(1,2), log2)


pheatmap::pheatmap(log2df1, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   # filename = "C:/Users/andre/Desktop/heatmap_raw.png", width = 15, height = 10,
                   annotation_col = chr_ann )


log2df1 %>% apply(2, function(x){  })


log2df1[,1]



# normalize amplicons (columns)
scaled_log2df1 <-  scale(log2df1)
scaled_log2df1 %>% dim


# normalize cells (rows)
scaled2_log2df1 <- scale(t(scaled_log2df1)) %>% t
scaled2_log2df1 %>% dim

scaled2_log2df1 %>% as.matrix() %>% density(na.rm=T) %>% plot

  
pheatmap::pheatmap(scaled2_log2df1, cluster_rows = T, cluster_cols = F, 
                   show_rownames = F, show_colnames = F,
                   # filename = "C:/Users/andre/Desktop/heatmap_normalized.png",  width = 15, height = 10,
                   annotation_col = chr_ann  )




log2_sc_df <- apply(scaled_df, c(1,2), log2)

pheatmap::pheatmap(log2_sc_df, cluster_rows = F, cluster_cols = F, 
                   show_rownames = F, show_colnames = F,
                   # filename = "C:/Users/andre/Desktop/heatmap_normalized.png",  width = 15, height = 10,
                   annotation_col = chr_ann  )



dimRedMat %>% colSums() %>% unname()

dev.off()

dimRedMat <- m2[,-63]


dimRedMat %>% is.na() %>% table

PCA <- prcomp(dimRedMat %>% t) 

screeplot(PCA, npcs = 40, type = "l")

fviz_eig(PCA, ncp = 30)


post_pca <- PCA$x[,1:15]

library(umap)

umap_data <- umap(post_pca)

umap_df <- umap_data$layout %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)


dbscan_clusters <- dbscan::dbscan(umap_df, eps = 0.3)

umap_df$cluster <- dbscan_clusters$cluster %>% as.factor()

umap_df %>% ggplot(aes(UMAP1, UMAP2)) + geom_point(aes(colour=cluster))






# ======================== VARIANTS =======================

var <- h5read(h5filePath, "/assays/dna_variants/layers")





