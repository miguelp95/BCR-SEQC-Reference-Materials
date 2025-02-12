library(tidyverse)
projdir <- "C:/Users/Miguel.Pacheco/OneDrive - FDA/Documents/BCR-SEQC/resources"

order_tb <- read.csv(file.path(projdir, "mixcr-cell-line-order-clonality.csv"),
                     stringsAsFactors = F) %>% 
  as_tibble()

# Read all TPM data
all_tpm <- read.csv(file.path(projdir, "BCR-SEQC_50_cell_lines_IG_TPM_reads.csv")) %>%
  as_tibble()

# Filter for heavy chain only and calculate total heavy chain TPM for each cell line
total_tpm <- all_tpm %>%
  group_by(cell_line) %>%
  summarise(total_tpm = sum(tpm, na.rm = TRUE)) %>%
  arrange(desc(total_tpm))

# Reorder the cell_line factor based on total TPM in descending order
all_tpm <- all_tpm %>%
  mutate(cell_line = factor(cell_line, levels = total_tpm$cell_line))



pdf(file.path(projdir, "50-cell-line-tpm-bar_sorted_desc_rotlabels.pdf"),
    height = 6, width = 10)
ggplot(all_tpm, aes(cell_line, tpm, fill = chain, group = cell_line)) +
  geom_bar(stat = "identity") +
  theme_light() + 
  ylab("IG Expression Abundance (TPM)") + xlab("Cell Line") +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
dev.off()
