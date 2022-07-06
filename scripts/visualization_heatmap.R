
TEST <- df8

chrGap <- colnames(TEST) %>% str_extract("(?<=chr)[0-9]+") %>% as.numeric() %>% table %>% as.vector() %>% cumsum()


TEST %>% as.matrix() %>% density() %>% plot
TEST %>% as.matrix() %>% as.vector() %>% summary
TEST %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.1))

breaksList = seq(-2, 2, by =0.1)

breaksList = c(-2,-0.85,0.85,2)

set.seed(1)
pheatmap(TEST[sample(1:nrow(TEST),300),], 
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         cutree_rows = 5,
         gaps_col = chrGap, 
         color = c("navy", "white", "firebrick3"), breaks = breaksList,
         # color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         annotation_col = chr_ann )



#____________________________ LOG SCALE _____________________________


TEST[TEST==0] <- 1

TESTl <- log2(TEST)

TESTl %>% as.matrix() %>% density() %>% plot
TESTl %>% as.matrix() %>% as.vector() %>% summary
TESTl %>% as.matrix() %>% quantile(probs = seq(0, 1, 0.05))


breaksList = seq(-2, 2 , by =0.1)
set.seed(1)
pheatmap(TESTl[sample(1:nrow(TESTl),300),], 
         cluster_rows = T, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         gaps_col = chrGap, 
         cutree_rows = 5,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)), breaks = breaksList,
         annotation_col = chr_ann )
