library(tidyverse)
library(RColorBrewer)

projdir <- "C:/Users/Miguel.Pacheco/OneDrive - FDA/Documents/BCR-SEQC/resources"
plotdir <- file.path(projdir, "plots")
dir.create(plotdir)

################################################################################################
# read in data 
all_res <- read.csv(file.path(projdir,"50-cell-lines-mixcr-v4.4.0-top10-clones.csv")) %>%
  as_tibble() %>%
  mutate(clone = factor(cloneNum, levels = paste0("Clone ", 10:1))) %>% 
  mutate(Percentage = ifelse(Chain == "Light", -Percentage, Percentage)) %>% 
  mutate(Assay = gsub("TakaraBio", "Takara Bio", Assay))

# read in order 
order_tb <- read.csv(file.path(projdir, "mixcr-cell-line-order-clonality.csv"),
                     stringsAsFactors = F) %>% 
  as_tibble()
cell_line_order <- order_tb$Sample

################################################################################################
# plot 

all_comb_ordered <- all_res %>%
  filter(clone != "Clone 10") %>%
  mutate(Assay = factor(Assay, levels = c("NEB", "Takara Bio", "RNA-seq", "WGS-Illumina", "WGS-MGI"))) %>%
  mutate(Sample = factor(Sample, levels = cell_line_order)) %>%
  group_by(Sample, Assay, Chain) %>%
  mutate(label = ifelse(max(tagCount) <= 10 | 
                          ((Assay == "RNA-seq") & max(tagCount) <= 1000), 
                        "*", ""),
         nudge = ifelse(Chain == "Heavy", 0.2, -0.2)) %>%
  ungroup() 

top_2 <- all_comb_ordered %>% 
  filter(clone %in% c("Clone 1", "Clone 2")) %>%
  mutate(clone = factor(clone, 
                        levels = c("Other", "Clone 2", "Clone 1")))

other <- top_2 %>%
  mutate(Percentage = abs(Percentage)) %>%
  group_by(Sample, Chain, Assay) %>%
  summarize(Percentage = 100 - sum(Percentage)) %>%
  ungroup() %>%
  mutate(Percentage = ifelse(Chain == "Heavy", Percentage, -Percentage)) %>%
  distinct() %>%
  mutate(clone = factor("Other", 
                        levels = c("Other", "Clone 2", "Clone 1")))

top_2_other <- bind_rows(top_2, other) %>%
  arrange(clone)

# Modify ggplot to use different color palettes based on Chain
all_distbn_plot <- ggplot(top_2_other, aes(Sample, Percentage, fill = interaction(Chain, clone))) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Clone percent") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_abline(slope = 0, intercept = 0) +
  geom_abline(slope = 0, intercept = 90, lty = "dashed") +
  geom_abline(slope = 0, intercept = -90, lty = "dashed") +
  facet_wrap(~Assay, ncol = 1, strip.position = "right") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  geom_text(data = filter(top_2_other, Chain == "Heavy"), 
            aes(y = 5, label = label), size = 8) +
  geom_text(data = filter(top_2_other, Chain == "Light"), 
            aes(y = -20, label = label), size = 8) +
  scale_fill_manual(
    values = c(
      "Light.Other" = "#FF9FD4",#brewer.pal(5, "Set2")[5],
      "Light.Clone 2" = "lightgoldenrod1",#brewer.pal(5, "Set2")[3],
      "Light.Clone 1" = "lightskyblue",#brewer.pal(5, "Set2")[4],
      "Heavy.Other" = "deeppink",#brewer.pal(5, "Dark2")[5],
      "Heavy.Clone 2" = "darkgoldenrod1",#brewer.pal(5, "Dark2")[3],
      "Heavy.Clone 1" = "royalblue2"#brewer.pal(5, "Dark2")[4]
    ),
    limits = c(
      "Heavy.Clone 1", "Heavy.Clone 2", "Heavy.Other", 
      "Light.Clone 1", "Light.Clone 2", "Light.Other"
    )  # Explicit legend order
  ) +
  guides(fill = guide_legend(reverse = FALSE))


pdf(file.path(plotdir, "50-cell-lines-mixcr-clone-distrbn-20250211_cbf.pdf"),
    height = 12, width = 15)
all_distbn_plot
dev.off()

out_tb <- all_comb_ordered %>% 
  filter(clone == "Clone 1") %>%
  mutate(Pipeline = "MiXCR") %>%
  select(Assay, `Cell line` = Sample, Chain, V, J, CDR3, label)

write.table(out_tb, 
           file.path(projdir, "50-cell-lines-mixcr-top-clone-20250211.csv"),
            sep = ",", quote = F, col.names = T, row.names = F)
