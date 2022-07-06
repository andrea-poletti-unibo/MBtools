
library(rhdf5)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

#================= LOCALIZE DATA ==================
list_h5_files <- list.files("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/", pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

#============ LOAD AND PREPARE RAW READ COUNT ===========

all_samp_rc <- data.frame()

s <- list_h5_files[1]
for(s in list_h5_files ){

  h5filePath <- s
  sample <- s %>% str_extract("(?<=scDNA_).*(?=_reference)")

  message(sample)
  rc <- h5read(h5filePath, "/assays/dna_read_counts/layers")

  df <- rc$read_counts %>% t %>% as.data.frame()

  # extract and set colnames(ampl) and rownames(cells)
  names <- h5read(h5filePath, "/assays/dna_read_counts/ca" )
  cells <- h5read(h5filePath, "/assays/dna_read_counts/ra")

  rownames(df) <- paste0(sample,"/", cells$barcode)
  colnames(df) <- paste0(names$id,"/chr",names$CHROM,":", names$start_pos, "-", names$end_pos)


  all_samp_rc <- rbind(all_samp_rc, df)

}

table( str_extract( rownames(all_samp_rc), ".*(?=/)") )


#_____ load missionbio mosaic normalized reads for mix sample (for comparison) ____
norm_reads_mix <- data.table::fread("data/scDNA_Mix_reference_norm_reads.csv")


#_____ load amplicon stats ____
ampl <- data.table::fread("data/amplicons_list_GC_mapp.txt")


AMPL_tot_reads <- colSums(all_samp_rc %>% select(!matches("tx")))

AMPL_tot_reads %>% plot


ampl$tot_reads <- AMPL_tot_reads


ampl %>% ggplot(aes(tot_reads, GC_content)) + geom_point()

ampl$GC_content_cat <- ampl$GC_content %>% cut(breaks=c(seq(0,1,0.05))) %>% as.numeric() %>% `/`(20)
ampl$GC_content_cat %>% table

ampl$mappability_cat <- ampl$mappability_score %>% cut(breaks=c(seq(0,1,0.1))) %>% as.numeric() %>% `/`(10)
ampl$mappability_cat %>% table

ampl %>% ggplot(aes(GC_content_cat %>% as.factor, tot_reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.05)



ampl %>% ggplot(aes(mappability_cat %>% as.factor(), tot_reads)) +
  geom_boxplot(outlier.shape = NA) +
  geom_violin() +
  geom_jitter(height = 0, width = 0.05)



gg <- ampl %>% ggplot(aes(GC_content_cat, tot_reads)) +
  geom_jitter(alpha=0.5, width = 0.005, height = 0) +
  geom_smooth() +
  ylim(0,1.5*10^7)

gg



MED_OVERALL <- median(ampl$tot_reads)

GCcatTab <- ampl %>% group_by(GC_content_cat) %>% summarise(median_tot_reads_GCcat=median(tot_reads))

GCcatTab$median_overall <- MED_OVERALL
GCcatTab$corr_fact <- GCcatTab$median_overall / GCcatTab$median_tot_reads_GCcat


ampl2 <- left_join(ampl, GCcatTab, by = "GC_content_cat")

ampl2$tot_reads_GCcorr <- ampl2$tot_reads * ampl2$corr_fact

ampl2%>% ggplot(aes(GC_content_cat, tot_reads_GCcorr)) +
  geom_point() +
  geom_smooth() +
  ylim(0,1.5*10^7)


check <- data.frame(
  amplicon= rep(ampl$ID_amplicon, 2),
  GC_content_cat= rep(ampl$GC_content_cat, 2),
  totReads= c(ampl$tot_reads, ampl2$tot_reads_GCcorr),
  phase=c(rep("1-pre", nrow(ampl)), rep("2-post", nrow(ampl))) )


check %>% ggplot(aes(GC_content_cat %>% as.factor, totReads)) +
  geom_point() +
  geom_boxplot() + facet_grid(~phase) +
  ylim(0,1.5*10^7)



#####################

samp <- all_samp_rc %>% filter(str_detect(rownames(all_samp_rc), "Mix"))




col_adjustGC <- c( ampl2$corr_fact, rep(1, 13))



#_____ Row/cell normalizarion______
# divide each row/cell for it's mean across cols/ampl

row_means <- apply(samp, 1, mean )

samp_norm <- samp %>% sweep(1, row_means, `/` )

samp[1:3, 1:2]
samp_norm[1:3, 1:2]


#_____ Column/amplicon normalizarion______

samp_GC <- samp_norm %>% sweep(2, col_adjustGC, `*` )

samp_GC_cell_log <- apply(samp_GC_cell, c(1,2), log2)

samp_GC_cell_log %>% as.matrix() %>% density() %>% plot

samp_GC_cell_log[samp_GC_cell_log==-Inf] <- -10


chr_ann <- data.frame(chr=names$CHROM %>% factor(levels = str_sort(names$CHROM %>% unique,numeric = T), ordered = T), row.names = colnames(df))
chrGap <- chr_ann$chr %>% table %>% as.vector() %>% cumsum()


breaksList = seq(-5, 5, by = 0.2)
pheatmap(samp_GC_cell_log, cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         # filename = "plots/hm.png", width = 15, height = 10,
         annotation_col = chr_ann )


##################################################
# remove bad amplicons

all_samp_rc[,1:519] %>% apply( 2, median) %>% `>`(5) %>% table
idx_low <- (all_samp_rc[,1:519] %>% apply( 2, median) %>% `<`(5))


all_samp_rc[,1:519] %>% apply( 2, median) %>% `>`(500) %>% table
idx_high <- (all_samp_rc[,1:519] %>% apply( 2, median) %>% `>`(500))

ampl2$ampl_median <- all_samp_rc[,1:519] %>% apply( 2, median)

ampl2$cov_quality <- ifelse(idx_low==T, "low_cov_m5", ifelse(idx_high==T, "high_cov_M500", "ok"))

ampl2$quality %>% table


ampl2 %>% ggplot(aes(ID_amplicon, ampl_median %>% log, colour=quality)) + geom_point()

keep_ampl <- ampl2$quality == "ok"

ampl2 %>% filter(cov_quality=="ok") %>% .$ampl_median %>% sort %>% plot

write_tsv(ampl2, "data/amplicons_list_allStats.txt")

ampl3 <- ampl2 %>% filter(cov_quality=="ok")


#====================

s <- all_samp_rc %>% filter(str_detect(rownames(all_samp_rc), "Mix")) %>% select(1:519)
s <- s[,keep_ampl]


#_____ Column/amplicon ability computation______

ampl_medians <- apply(s, 2, median)
ampl_medians %>% sort %>% plot

ampl_ability <- ampl_medians/ median(ampl_medians)

#_____ Row/cell normalizarion______

row_means <- apply(s, 1, mean )

s2 <- s %>% sweep(1, row_means, `/` )


#_____ Column/amplicon normalizarion______

s3 <- s2 %>% sweep(2, ampl_ability, `/` )


#_____ Log2 transformation ________
s3l <- apply(s3, c(1,2), log2)

s3l %>% as.matrix() %>%
  density() %>% plot + abline(v=c(-1.8,-0.8,0))

s3l %>% as.matrix() %>% median

s3l[s3l==-Inf] <- -10

seq(-2, 2, by =1) -0.8

chr_ann <- data.frame(chr=names$CHROM %>% factor(levels = str_sort(names$CHROM %>% unique,numeric = T), ordered = T), row.names = colnames(df))

breaksList = seq(-2.8, 1.2, by =1)
pheatmap(s3l, cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         # filename = "plots/hm.png", width = 15, height = 10,
         annotation_col = chr_ann )


s3lt <- s3l %>% as.matrix()
s3lt[s3l < -1.8] <- 0
s3lt[s3l>0] <- 2
s3lt[s3l<0 & s3l> -1.8] <- 1

s3lt %>% as.matrix() %>% density() %>% plot

breaksList = seq(-0.5, 2.5, by =1)
pheatmap(s3lt, cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         color = c("navy", "white", "firebrick3"), breaks = breaksList,
         # filename = "plots/hm.png", width = 15, height = 10,
         annotation_col = chr_ann )



#============= experiments ===============

sub <- s3[sample(1:nrow(s3),400),]

sub[sub>1] <- 1
sub %>% as.matrix() %>% density() %>% plot
sub %>% as.matrix() %>% quantile()


chrGap <- colnames(sub) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()

breaksList = seq(-2, 2, by =0.05)
pheatmap(sub, cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         cutree_rows = 5, scale = "row",
         gaps_col = chrGap,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         # filename = "plots/hm.png", width = 15, height = 10,
         annotation_col = chr_ann )


subS <- apply(sub, 1, scale)


