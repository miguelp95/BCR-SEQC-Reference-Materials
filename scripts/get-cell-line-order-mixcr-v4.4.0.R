library(tidyverse)

projdir <- "C:/Users/Lily.Zheng/OneDrive - FDA/Documents/projects/BCR_SEQC/50-cell-lines"
mixcrdir <- file.path(projdir, "results-mixcr-v4.4.0")
plotdir <- file.path(mixcrdir, "plots")

################################################################################################
# read in data
all_res <- read.csv(file.path(mixcrdir,"50-cell-lines-mixcr-v4.4.0-top10-clones.csv")) %>%
  as_tibble() %>%
  mutate(clone = factor(cloneNum, levels = paste0("Clone ", 10:1))) %>% 
  mutate(Percentage = ifelse(Chain == "Light", -Percentage, Percentage)) %>% 
  mutate(Assay = gsub("TakaraBio", "Takara Bio", Assay))


################################################################################################
# order cell lines by consistency 

all_comb <- all_res %>%
  filter(Assay %in% c("NEB", "Takara Bio", "RNA-seq")) %>% 
  mutate(Assay = factor(Assay, levels = c("NEB", "Takara Bio", "RNA-seq", "WGS-Illumina", "WGS-MGI")))

###   Collect top 2 clone info   ###############################################
clone2 <- all_comb %>% 
  filter(clone == "Clone 1" | clone == "Clone 2") %>% 
  select(Sample, Chain, Assay, Percentage, cloneNum) %>% 
  group_by(Sample, Chain, Assay) %>% 
  pivot_wider(names_from = cloneNum, values_from = Percentage) %>% 
  ungroup()

###   Classify clonality   #####################################################
top <- clone2 %>%
  mutate(Clonality = ifelse(`Clone 1` >= 90 & `Clone 2` < 5, "Mono",
                            ifelse(`Clone 1` + `Clone 2` >= 90, "Bi", "Multi")))

top_wider <- top %>%
  pivot_wider(id_cols = c(Sample, Chain),
              names_from = Assay,
              values_from = Clonality)

###   Sort by NAs   ############################################################
top_wider <- top_wider %>%
  mutate(num_NA = apply(top_wider, 1,
                        function(x) sum(sapply(x[3:5], is.na)))) 

###   sort by agreement    #####################################################
top_wider <- top_wider %>%
  mutate(num_classifications =  apply(top_wider, 1,
                                      function(x) length(unique(x[3:5][!is.na(x[3:5])])))) %>%
  arrange(num_NA, num_classifications)
top_wider

###   determine top classification across assays   #############################
getTopClass <- function(x) {
  x <- x[!is.na(x)]
  top_class <- names(which.max(table(x)))
  return(factor(top_class, levels = c("Mono", "Bi", "Multi")))
}

light_class <- top_wider %>%
  filter(Chain == "Light")
colnames(light_class)[3:7] <- paste0("light_", colnames(light_class)[3:7])
light_class <- light_class %>%
  mutate(top_light_class = apply(light_class, 1,
                                 function(x) getTopClass(x[3:5])))

heavy_class <- top_wider %>%
  filter(Chain == "Heavy")
colnames(heavy_class)[3:7] <- paste0("heavy_", colnames(heavy_class)[3:7])
heavy_class <- heavy_class %>%
  mutate(top_heavy_class = apply(heavy_class, 1,
                                 function(x) getTopClass(x[3:5])))

cell_line_class <- full_join(heavy_class, light_class, by = "Sample") 

###   get order of top clone percentage    #####################################
top_clone1_tb <- top %>%
  pivot_wider(id_cols = c(Sample, Chain),
              values_from = `Clone 1`, 
              names_from = Assay)

top_clone1_h <- top_clone1_tb %>%
  filter(Chain == "Heavy")
colnames(top_clone1_h)[3:5] <- paste0("c1_heavy_", colnames(top_clone1_h)[3:5])
top_clone1_h$Chain <- NULL
top_clone1_l <- top_clone1_tb %>%
  filter(Chain == "Light")
colnames(top_clone1_l)[3:5] <- paste0("c1_light_", colnames(top_clone1_l)[3:5])
top_clone1_l$Chain <- NULL

top_clone1 <- full_join(top_clone1_h, top_clone1_l)

cell_line_class <- left_join(cell_line_class, top_clone1)

###   final ordering of cell lines    ##########################################
cell_line_class_ordered <- cell_line_class %>%
  arrange(heavy_num_NA, light_num_NA, 
          heavy_num_classifications, light_num_classifications,
          top_heavy_class, top_light_class,
          desc(c1_heavy_NEB), desc(c1_light_NEB), 
          desc(`c1_heavy_Takara Bio`), desc(`c1_light_Takara Bio`), 
          desc(`c1_heavy_RNA-seq`), desc(`c1_light_RNA-seq`))  %>% 
  select(Sample, 
         heavy_NEB, `heavy_RNA-seq`, `heavy_Takara Bio`, 
         light_NEB, `light_RNA-seq`, `light_Takara Bio`)

cell_line_order <- cell_line_class_ordered$Sample

colnames(cell_line_class_ordered)[2:ncol(cell_line_class_ordered)] <- 
  gsub(" ", "", paste0(colnames(cell_line_class_ordered)[2:ncol(cell_line_class_ordered)], "_MiXCR_clonality"))
order_tb <- bind_cols(tibble(index = 1:length(cell_line_order)),
                      cell_line_class_ordered)
write.table(order_tb, file.path(mixcrdir, "top-clones", "tables", "mixcr-cell-line-order-clonality.csv"),
            sep = ",", quote = F, row.names = F, col.names = T)