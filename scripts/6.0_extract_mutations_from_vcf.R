

# BiocManager::install("vcfR", update = F)

vcf.gz_files <- list.files("~/data/VCF_files/",pattern = "vcf.gz", full.names = T)
vcf.gz_files

library(vcfR)

for(i in vcf.gz_files){
  i=1

  message(vcf.gz_files[i])

  vcf <- read.vcfR(vcf.gz_files[i])

  tidyVCF<- vcfR2tidy(vcf, single_frame = T, dot_is_NA=T)

  df <- tidyVCF$dat

  readr::write_tsv(df, paste0("~/data/VCF_files/extracted_mutations_tables/",i,"mutations.tsv"))
}
