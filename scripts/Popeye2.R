
Popeye2 <- function(segments, refGenome = "hg19", removeXY=T) {

  options(scipen=999)

  require(GenomicRanges)
  require(plyranges)
  require(data.table)
  require(tidyverse)

  # load chr arm table with arms start/end info - hg19 or hg38
  if(refGenome == "hg19") {
    bands <- data.table::fread("data/cytoBandIdeo_hg19.txt")
    } else if (refGenome =="hg38") {
    bands <- data.table::fread("data/cytoBandIdeo_hg38.txt")
    } else { stop("Invalid reference genome. Please choose hg19 or hg38")}

  names(bands) <- c("chr", "start", "end", "band", "giemsa")
  bands$arm <- str_extract(bands$band, "[pq]")
  bands$chrarm <- paste0(str_remove(bands$chr,"chr"), bands$arm)
  bands2 <- bands %>% filter(!is.na(arm))
  chrtab <- bands2 %>% group_by(chrarm) %>%
    summarise(chr= unique(chr) %>% str_remove("chr"),
              start=min(start),
              end=max(end),
              arm=str_extract(chrarm,"[pq]") %>% unique,
              ref=refGenome) %>%
    select(chr, start, end, chrarm, arm, ref)

  chrtabGR <- makeGRangesFromDataFrame(chrtab, keep.extra.columns = T)

  # optional: remove X Y chr
  if(removeXY == T){
    segments <- segments %>% filter(!chr %in% c("X", "Y"))
  }

  segmentsGR <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)

  # plyranges "join" the two granges
  segmentsGR_w_arms <- join_overlap_intersect(chrtabGR, segmentsGR)

  # arrange per decreasing chr and per decreasing start pos - rename the seqnames column
  segments_anno <- segmentsGR_w_arms %>% as.data.frame() %>%
    dplyr::rename("chr"="seqnames") %>%
    arrange(chr, start)

  segments_anno
}
