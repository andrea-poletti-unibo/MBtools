library(tidyverse)
library(data.table)

files <- list.files("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/python_notebook_CNA/Spyder_analysis/raw_reads/", full.names = T)

i=9

res <- data.frame()
for( i in 1:length(files)){

  message(i)
  df <- fread(files[i])

  samp_name <- files[i] %>%  str_extract("[Ss]ample.*(?=_raw)")
  message(samp_name)

  df$sample <- samp_name

  res <- rbind(res, df)

}

res$sample %>% table



res2 <- res %>% group_by(sample) %>% group_map(~ apply(.x[,-1],2,mean))

# res2 <- res %>% group_by(sample) %>% group_map(~ colSums(.x[,-c(1)]))

grouped <- group_by(res, sample)
names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))


res3 <- res2 %>% unlist

dat <- data.frame( ampl= names(res3), val= res3, sample= rep(names,each=127))


dat %>% ggplot(aes(ampl, val)) + geom_point() + geom_boxplot() + ylim(NA, 500) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



dat$amplicon <- dat$ampl %>% str_remove_all("AML_v2_")


odf <- dat %>% group_by(amplicon) %>% summarise(median=median(val)) %>% arrange(desc(median))

order <- odf$amplicon

dat$ampl_order <- dat$amplicon %>% factor(levels = order)

dat$gene <- dat$amplicon %>% str_extract(".*(?=_[0-9]+$)")


dat %>% ggplot(aes(ampl_order, val)) +
  geom_boxplot(outlier.shape = NA, alpha=0.6, aes(fill=gene)) +
  geom_jitter(height = 0, width = 0.03, alpha=0.3, size=0.8) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dat %>% ggplot(aes(ampl_order, val)) +
  stat_summary(fun = median, geom = "bar", fill="dodgerblue4") + # plot mean (column)
  coord_cartesian(ylim = c(0, 150)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
# Create a ggplot with 18 colors
# Use scale_fill_manual

dat %>% ggplot(aes(ampl_order, val)) +
  geom_boxplot(outlier.shape = NA, alpha=0.8, aes(fill=gene)) +
  # geom_jitter(height = 0, width = 0.03, alpha=0.5, size=0.8) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
  scale_fill_manual(values = mycolors) +
  ggtitle("Amplicon median coverage per cell ") + xlab("Amplicon") + ylab("Median coverage")

palette2 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

palette3 <- sample(palette2,20)

dat %>% ggplot(aes(ampl_order, val)) +
  geom_boxplot(outlier.shape = NA, alpha=0.8, aes(fill=gene)) +
  # geom_jitter(height = 0, width = 0.03, alpha=0.5, size=0.8) +
  coord_cartesian(ylim = c(0, 400)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)) +
  scale_fill_manual(values = palette3) +
  ggtitle("Amplicon median coverage per cell ") + xlab("Amplicon") + ylab("Median coverage")


ggsave("C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/OUTPUT_PLOTS/Amplicon_coverage_plot_EnricA.png", width = 12, height = 8)

