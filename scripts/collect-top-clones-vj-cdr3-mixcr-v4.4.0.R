library(tidyverse)
library(RColorBrewer)

projdir <- "C:/Users/Lily.Zheng/OneDrive - FDA/Documents/projects/BCR_SEQC/50-cell-lines"
mixcrdir <- file.path(projdir, "results-mixcr-v4.4.0")
plotdir <- file.path(mixcrdir, "plots")
dir.create(plotdir)

################################################################################################
# Functions
################################################################################################
readMixcrResults <- function(file, samp_field = 5) {
  bn <- basename(file)
  assay <- basename(dirname(file))
  sample <- gsub("-IGH", "", gsub("-IGL", "", strsplit(bn, "_|\\.")[[1]][samp_field]))
  chain <- ifelse(grepl("heavy", file), "Heavy", "Light")

  res <- read.delim(file, quote = NULL) %>%
    as_tibble() %>%
    mutate(Chain = chain,
           Sample = sample,
           Assay = assay)
  
  # check if table has columns bestVGene, bestJGene
  hasGeneCols <- any(any(grepl("bestVGene", colnames(res))),
                     any(grepl("bestJGene", colnames(res))))
  
  if (!hasGeneCols) {
    getName <- function(x) {
      gsub("\\*.*", "", x)
    }
    res$bestVGene <- getName(res$allVHitsWithScore)
    res$bestJGene <- getName(res$allJHitsWithScore)
  }
  
  return(res)
}

# collapse clones by V, J, CDR3 nucleotide seq 
collapseClones <- function(mixcr_tb) {
  mixcr_tb$nSeqCDR3[is.na(mixcr_tb$nSeqCDR3)] <- ""
  s <- unique(mixcr_tb$Sample) 
  ch <- unique(mixcr_tb$Chain) 
  assay <- unique(mixcr_tb$Assay)
  
  if (any(grepl("uniqueMolecule", colnames(mixcr_tb)))) {
    # use unique molecule 
    collapsed_x <- aggregate(cbind(uniqueMoleculeCount, uniqueMoleculeFraction) ~ bestVGene + bestJGene + nSeqCDR3 + aaSeqCDR3, data = mixcr_tb, FUN=sum)
    new_sum <- sum(collapsed_x$uniqueMoleculeCount)
    collapsed_x$cloneFraction <- collapsed_x$uniqueMoleculeCount / new_sum
    collapsed_x$uniqueMoleculeFraction <- NULL
    collapsed_x <- collapsed_x %>% 
      as_tibble() %>% 
      rename(tagCount = uniqueMoleculeCount) %>% 
      mutate(tagType = "UMI")
  } else {
    # use read count
    collapsed_x <- aggregate(cbind(readCount, readFraction) ~ bestVGene + bestJGene + nSeqCDR3 + aaSeqCDR3, data = mixcr_tb, FUN=sum)
    new_sum <- sum(collapsed_x$readCount)
    collapsed_x$cloneFraction <- collapsed_x$readCount / new_sum
    collapsed_x$readFraction <- NULL
    collapsed_x <- collapsed_x %>% 
      as_tibble() %>% 
      rename(tagCount = readCount) %>% 
      mutate(tagType = "read")
  }
  
  collapsed_x <- collapsed_x %>% 
    as_tibble() %>% 
    rename(V = bestVGene, J = bestJGene, CDR3 = nSeqCDR3, CDR3_AA = aaSeqCDR3) %>% 
    mutate(Sample = s, Chain = ch, Assay = assay) %>% 
    arrange(-cloneFraction) %>% 
    select(Sample, Chain, Assay, V, J, CDR3, CDR3_AA, tagCount, cloneFraction, tagType) %>%
    mutate(Percentage = round(tagCount / sum(tagCount) * 100, 1)) %>% 
    mutate(cloneNum = 1:n())
  return(collapsed_x)
}

readInAssayRes <- function(mixcrdir, assay, samp_field = 5) {
  assay_files <- list.files(file.path(mixcrdir, assay))
  assay_res <- lapply(assay_files, 
                      function(x) readMixcrResults(file.path(mixcrdir, assay, x),
                                                   samp_field = samp_field))
  assay_res <- assay_res[sapply(assay_res, function(x) nrow(x) != 0)] 
  assay_res_collapse <- lapply(assay_res, collapseClones) %>% 
    bind_rows 
  return(assay_res_collapse)
}

################################################################################################

all_assays <- c("NEB", "TakaraBio", "RNA-seq", "WGS-Illumina", "WGS-MGI")
sf <- c(5, 5, 4, 5, 4)

all_res <- mapply(function(assay, field) readInAssayRes(mixcrdir, assay, field),
                  assay = all_assays,
                  field = sf,
                  SIMPLIFY = F) %>%
  bind_rows()

top10 <- all_res %>% 
  filter(cloneNum <= 10) %>% 
  mutate(cloneNum = paste0("Clone ", cloneNum)) %>%
  mutate(V = gsub("D", "", V))

write.table(top10, file.path(mixcrdir, "50-cell-lines-mixcr-v4.4.0-top10-clones.csv"),
            sep = ",", quote = F, row.names = F, col.names = T)


