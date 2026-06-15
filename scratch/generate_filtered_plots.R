library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(data.table)

cat("Loading final filtered biomarkers...\n")
degs_filtered <- readRDS("data_clean/sig_degs_hspc_filtered.rds")
cat("Loaded", nrow(degs_filtered), "filtered biomarkers.\n")

final_up <- degs_filtered %>% filter(direction == "Up") %>% pull(gene_symbol)
final_down <- degs_filtered %>% filter(direction == "Down") %>% pull(gene_symbol)

cat("Final up-regulated:", length(final_up), "\n")
cat("Final down-regulated:", length(final_down), "\n")

# Load full differential expression lists
cat("Loading full DEG lists...\n")
degs_d_all <- fread("data_clean/discovery_degs_all.csv")
degs_v_all <- fread("data_clean/gse13159_degs_all.csv")

# Set thresholds
adjp_thresh <- 0.01
logfc_d_thresh <- 1.5
logfc_v_thresh <- 0.15

# ------------------------------------------------------------------------------
# 1. Discovery Volcano Plot (Filtered)
# ------------------------------------------------------------------------------
cat("Generating Discovery Volcano Plot (Filtered)...\n")
temp_d <- degs_d_all %>%
  mutate(
    change = case_when(
      adj.P.Val < adjp_thresh & logFC > logfc_d_thresh & gene_symbol %in% final_up ~ "Final Robust Up",
      adj.P.Val < adjp_thresh & logFC < -logfc_d_thresh & gene_symbol %in% final_down ~ "Final Robust Down",
      adj.P.Val < adjp_thresh & logFC > logfc_d_thresh ~ "Discovery Up Only",
      adj.P.Val < adjp_thresh & logFC < -logfc_d_thresh ~ "Discovery Down Only",
      TRUE ~ "Not Significant"
    )
  )

# Labels for top 10 up and top 10 down from final filtered list
top_up_labels <- degs_filtered %>% filter(direction == "Up") %>% arrange(desc(logFC_discovery)) %>% head(10) %>% pull(gene_symbol)
top_down_labels <- degs_filtered %>% filter(direction == "Down") %>% arrange(logFC_discovery) %>% head(10) %>% pull(gene_symbol)
label_genes_d <- temp_d %>% filter(gene_symbol %in% c(top_up_labels, top_down_labels))

p_vol_d <- ggplot(temp_d, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c(
    "Final Robust Up" = "#C00000",
    "Final Robust Down" = "#2E75B6",
    "Discovery Up Only" = "#FFC000",
    "Discovery Down Only" = "#9BC2E6",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-logfc_d_thresh, logfc_d_thresh), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(adjp_thresh), linetype = "dashed", color = "darkgrey") +
  geom_text_repel(
    data = label_genes_d,
    aes(label = gene_symbol),
    size = 3,
    color = "black",
    max.overlaps = 20
  ) +
  theme_bw() +
  labs(
    title = paste0("Discovery Cohort Volcano (Safety Filtered, adj.P < ", adjp_thresh, ")"),
    x = "log2FC",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/discovery_volcano_hspc_filtered.png", plot = p_vol_d, width = 8, height = 6.5, dpi = 300)
cat("Discovery volcano saved.\n")

# ------------------------------------------------------------------------------
# 2. Validation Volcano Plot (Filtered)
# ------------------------------------------------------------------------------
cat("Generating Validation Volcano Plot (Filtered)...\n")
temp_v <- degs_v_all %>%
  mutate(
    change = case_when(
      adj.P.Val < adjp_thresh & logFC > logfc_v_thresh & gene_symbol %in% final_up ~ "Final Robust Up",
      adj.P.Val < adjp_thresh & logFC < -logfc_v_thresh & gene_symbol %in% final_down ~ "Final Robust Down",
      adj.P.Val < adjp_thresh & logFC > logfc_v_thresh ~ "Validation Up Only",
      adj.P.Val < adjp_thresh & logFC < -logfc_v_thresh ~ "Validation Down Only",
      TRUE ~ "Not Significant"
    )
  )

label_genes_v <- temp_v %>% filter(gene_symbol %in% c(top_up_labels, top_down_labels))

p_vol_v <- ggplot(temp_v, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c(
    "Final Robust Up" = "#C00000",
    "Final Robust Down" = "#2E75B6",
    "Validation Up Only" = "#FFC000",
    "Validation Down Only" = "#9BC2E6",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-logfc_v_thresh, logfc_v_thresh), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(adjp_thresh), linetype = "dashed", color = "darkgrey") +
  geom_text_repel(
    data = label_genes_v,
    aes(label = gene_symbol),
    size = 3,
    color = "black",
    max.overlaps = 20
  ) +
  theme_bw() +
  labs(
    title = paste0("Validation Cohort Volcano (Safety Filtered, adj.P < ", adjp_thresh, ")"),
    x = "logFC",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/validation_volcano_hspc_filtered.png", plot = p_vol_v, width = 8, height = 6.5, dpi = 300)
cat("Validation volcano saved.\n")

# ------------------------------------------------------------------------------
# 3. Discovery Heatmap (Filtered)
# ------------------------------------------------------------------------------
cat("Generating Discovery Heatmap (Filtered)...\n")
# Load expression data
discovery_expr <- readRDS("data_clean/discovery_expr_log2tpm.rds")
discovery_sample_info <- readRDS("data_clean/discovery_sample_info.rds")

# Sort by absolute logFC in discovery and take top 50
top50_filtered <- degs_filtered %>%
  arrange(desc(abs(logFC_discovery))) %>%
  head(50)

top50_expr_d <- discovery_expr[top50_filtered$ensembl_id, ]
rownames(top50_expr_d) <- top50_filtered$gene_symbol

sample_annot <- data.frame(
  Group = discovery_sample_info$group,
  Source = discovery_sample_info$source,
  row.names = discovery_sample_info$sample_id
)

annot_colors <- list(
  Group = c("Normal" = "#2E75B6", "AML" = "#C00000"),
  Source = c("TCGA" = "#ED7D31", "GTEx" = "#70AD47")
)

pheatmap(
  top50_expr_d,
  scale = "row",
  annotation_col = sample_annot,
  annotation_colors = annot_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  filename = "results/heatmap_hspc_filtered.png",
  width = 10,
  height = 8
)
cat("Heatmap saved.\n")
cat("Finished creating final filtered plots!\n")
