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

samples

res2 <- res %>% group_by(sample) %>% group_map(~ apply(.x[,-1],2,mean))

# res2 <- res %>% group_by(sample) %>% group_map(~ colSums(.x[,-c(1)]))

grouped <- group_by(res, sample)
names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))


res3 <- res2 %>% unlist

dat <- data.frame( ampl= names(res3), val= res3, sample= rep(names,each=127))


dat %>% ggplot(aes(ampl, val)) + geom_point() + geom_boxplot() + ylim(NA, 500)


