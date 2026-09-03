# ==============================================================================
# MICROBIOME TAXONOMIC ABUNDANCE VISUALIZATION (Kraken2 / Bracken)
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(scales)

if (!dir.exists("plots")) dir.create("plots")

# 1. Ingest Bracken aggregated output
data_file <- ifelse(file.exists("bracken_output/ALL_SAMPLES_MICROBIOME.tsv"),
                    "bracken_output/ALL_SAMPLES_MICROBIOME.tsv",
                    "ALL_SAMPLES_MICROBIOME.tsv")

df <- read.delim(data_file, header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. STACKED BAR PLOT (Top 10 Taxa + "Other")
# ------------------------------------------------------------------------------
top_taxa <- df %>%
  group_by(Taxon) %>%
  summarise(Mean_Fraction = mean(Fraction_Total, na.rm = TRUE)) %>%
  arrange(desc(Mean_Fraction)) %>%
  slice(1:10) %>%
  pull(Taxon)

df_bar <- df %>%
  mutate(Taxon_Grouped = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
  group_by(Sample, Taxon_Grouped) %>%
  summarise(Fraction_Total = sum(Fraction_Total), .groups = "drop")

# Reorder factor so 'Other' sits at the top of the stack
df_bar$Taxon_Grouped <- factor(df_bar$Taxon_Grouped, levels = c(top_taxa, "Other"))

colors <- c(brewer.pal(min(length(top_taxa), 12), "Set3"), "#D3D3D3")

p_bar <- ggplot(df_bar, aes(x = Sample, y = Fraction_Total, fill = Taxon_Grouped)) +
  geom_bar(stat = "identity", position = "fill", color = "black", size = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = colors) +
  labs(
    title = "Microbiome Composition: Relative Species Abundance",
    x = "Sample Identifier",
    y = "Relative Abundance (%)",
    fill = "Taxon"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p_bar)
ggsave("plots/01_microbiome_relative_abundance.png", plot = p_bar, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3. TOP 20 TAXA HEATMAP
# ------------------------------------------------------------------------------
mat_df <- df %>%
  select(Sample, Taxon, Fraction_Total) %>%
  pivot_wider(names_from = Sample, values_from = Fraction_Total, values_fill = 0)

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$Taxon

# Retain top 20 taxa by overall average relative abundance
top20_taxa <- names(sort(rowMeans(mat), decreasing = TRUE))[1:min(20, nrow(mat))]
mat_top <- mat[top20_taxa, , drop = FALSE]

# Export Heatmap
png("plots/02_microbiome_taxa_heatmap.png", width = 2400, height = 2000, res = 300)
pheatmap(
  mat_top,
  color = colorRampPalette(c("#F7FCF5", "#41AB5D", "#00441B"))(100),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  display_numbers = FALSE,
  main = "Top 20 Species Abundance Heatmap (Fraction Total)"
)
dev.off()

cat("Visualization complete! High-resolution plots saved to 'plots/'\n")