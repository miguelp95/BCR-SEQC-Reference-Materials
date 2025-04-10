# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Input file path (modify this as needed)
projdir <- "C:/Users/Miguel.Pacheco/OneDrive - FDA/Documents/BCR-SEQC/resources"
plotdir <- "C:/Users/Miguel.Pacheco/OneDrive - FDA/Documents/BCR-SEQC/resources/plots"
input_file <- "cell_line_ploidy.csv"

# Read the data
data <- read.csv(file.path(projdir,input_file))

# Read all TPM data
all_tpm <- read.csv(file.path(projdir, "BCR-SEQC_50_cell_lines_IG_TPM_reads.csv")) %>%
  as_tibble()

# Filter for heavy chain only and calculate total heavy chain TPM for each cell line
total_tpm <- all_tpm %>%
  group_by(cell_line) %>%
  summarise(total_tpm = sum(tpm, na.rm = TRUE)) %>%
  arrange(desc(total_tpm))

# Check the structure of the data
if (!all(c("cell_line", "ploidy") %in% colnames(data))) {
  stop("The input file must contain 'cell_line' and 'ploidy' columns.")
}

# Rename the cell_line entries to keep only the portion between the first and second underscores
data <- data %>%
  mutate(cell_line = gsub("^[^_]+_([^_]+)_.*$", "\\1", cell_line))

# read in order 
order_tb <- read.csv(file.path(projdir, "mixcr-cell-line-order-clonality.csv"),
                     stringsAsFactors = F) %>% 
  as_tibble()
cell_line_order <- order_tb$Sample


# Arrange the cell lines horizontally
#data <- data %>%
#  mutate(index = match(cell_line, order_tb$Sample)) %>%
#  arrange(index) %>%
#  mutate(cell_line = factor(cell_line, levels = unique(cell_line)))

# Reorder the cell_line factor based on total TPM in descending order
data <- data %>%
  mutate(cell_line = factor(cell_line, levels = total_tpm$cell_line))

# Generate the heatmap
heatmap_plot <- ggplot(data, aes(x = cell_line, y = 1, fill = ploidy)) +
  geom_tile(color = "white", height = 2) + # Increase cell height
  scale_fill_gradientn(colors = c("royalblue2", "white", "deeppink"),#c("royalblue2", "white", "deeppink"),
                       limits = c(1, 3.5),
                       values = scales::rescale(c(1, 2, 3.5)),
                       name = "Ploidy") +
  labs(title = "Ploidy Heatmap",
       x = "Cell Line",
       y = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),               # Remove Y-axis text
    axis.ticks.y = element_blank(),              # Remove Y-axis ticks
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), # Cell line names at 45 degrees
    axis.ticks.x = element_line(color = "lightgrey"), # Black tick marks below each cell
    panel.border = element_rect(color = "lightgrey", fill = NA), # Light grey border
    panel.grid = element_blank(),                # Remove grid lines
    legend.position = "right",                    # Legend to the right
    legend.direction = "horizontal",
    legend.text = element_text(angle=45, vjust=1, hjust=1)
  ) +
  coord_fixed(ratio = 0.7) # Adjust aspect ratio for taller cells

# Save the heatmap as a PDF
pdf(file.path(plotdir, "ploidy_heatmap_for_fig2.pdf"), width = 10, height = 5)
print(heatmap_plot)
dev.off()