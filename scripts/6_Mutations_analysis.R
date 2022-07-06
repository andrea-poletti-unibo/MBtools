
library(rhdf5)
library(tidyverse)

H5_files_dir <- "C:/Users/andre/Alma Mater Studiorum Università di Bologna/PROJECT Single-Cell - Documenti/DFCI_Ghobrial_RUN/"

list_h5_files <- list.files(H5_files_dir, pattern = ".*h5$", full.names = T)

# mix sample
h5filePath <- list_h5_files[3]
h5filePath

rhdf5::h5ls(h5filePath)

variants <- h5read(h5filePath, "/assays/dna_variants/ca")

summary(variants)

names(variants)

var_df <- Reduce(cbind, variants) %>% as.data.frame()
colnames(var_df) <- names(variants)


var_df$filtered %>% table

var_pass <- var_df %>% filter(filtered=="00")

