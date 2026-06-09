# scripts/05_diff_expression_analysis.R

# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ------------------------------------------------------------------------------
# [EN] Check for and install required BioConductor and CRAN packages, then load libraries.
# [ZH] 检查并自动安装所需的 BioConductor 和 CRAN 包，然后加载所需的 R 库。
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("limma")
required_cran <- c("data.table", "dplyr", "ggplot2", "pheatmap", "ggrepel")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(limma)
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Create output folders
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 1: Load Datasets | 第 1 部分：读取输入数据
# ------------------------------------------------------------------------------
# [EN] Load Discovery cohort (TCGA+GTEx) and Validation cohort (GSE13159) data.
# [ZH] 读取发现集（TCGA+GTEx）与独立验证集（GSE13159）的表达矩阵及临床元数据。
# ==============================================================================

cat("Loading datasets...\n")
# Discovery Cohort
discovery_expr <- readRDS("data_clean/discovery_expr_log2tpm.rds")
discovery_sample_info <- readRDS("data_clean/discovery_sample_info.rds")

# Validation Cohort
gse_expr <- readRDS("data_clean/gse13159_expr_clean.rds")
gse_sample_info <- readRDS("data_clean/gse13159_sample_info.rds")

cat("Discovery matrix dimensions:", dim(discovery_expr)[1], "genes x", dim(discovery_expr)[2], "samples\n")
cat("Validation matrix dimensions:", dim(gse_expr)[1], "genes x", dim(gse_expr)[2], "samples\n")


# ==============================================================================
# SECTION 2: PCA Diagnostic of Discovery Cohort | 第 2 部分：发现集 PCA 降维诊断
# ------------------------------------------------------------------------------
# [EN] Run PCA on the raw merged matrix to visualize platform batch and biological separation.
# [ZH] 对合并后的发现集运行 PCA，观察样本在测序来源与疾病分组上的分布情况。
# ==============================================================================

cat("Running PCA on Discovery Cohort...\n")
# Filter out zero variance genes for PCA
var_genes <- apply(discovery_expr, 1, var)
expr_pca_input <- discovery_expr[var_genes > 0, ]

pca_res <- prcomp(t(expr_pca_input), scale. = TRUE)
pca_df <- data.frame(
  sample_id = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  source = discovery_sample_info$source,
  group = discovery_sample_info$group
)

percent_var <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 2)

# Save PCA plot colored by biological group
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = source)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  scale_color_manual(values = c("Normal" = "#2E75B6", "AML" = "#C00000")) +
  labs(
    title = "PCA of Merged Discovery Cohort",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/discovery_pca.png", plot = p_pca, width = 7, height = 5.5, dpi = 300)
cat("Discovery PCA plot saved to results/discovery_pca.png\n")


# ==============================================================================
# SECTION 3: Differential Expression in Discovery Cohort (limma) | 第 3 部分：发现集差异表达分析
# ------------------------------------------------------------------------------
# [EN] Run limma directly on the uncorrected cohort. (Note: ComBat is bypassed
#      because batch (TCGA vs GTEx) is 100% confounded with biology (AML vs Normal).
#      Instead, robust biological markers are identified via cross-platform intersection).
# [ZH] 直接在未校正的合并发现集上运行 limma 差异分析。（注：由于测序批次与疾病分组
#      100% 完全混杂，使用 ComBat 强行去批次会抹杀所有生物学差异。我们将通过后续跨平台验证集
#      求交集来过滤掉技术批次假阳性）。
# ==============================================================================

cat("Running limma on Discovery Cohort...\n")
discovery_sample_info$group <- factor(discovery_sample_info$group, levels = c("Normal", "AML"))

# Construct design matrix (Normal is baseline)
mod_d <- model.matrix(~ group, data = discovery_sample_info)

# Fit linear model and run Empirical Bayes shrinkage
fit_d <- lmFit(discovery_expr, mod_d)
fit2_d <- eBayes(fit_d)

# Extract complete results
degs_d_all <- topTable(fit2_d, coef = 2, number = Inf, adjust.method = "BH")
degs_d_all$ensembl_id <- rownames(degs_d_all)

# Read gene mapping file to add Gene Symbols
gene_map <- read.csv("data_raw/gtex/gtex_gene_info.csv", stringsAsFactors = FALSE)
gene_map_clean <- gene_map %>%
  select(ensembl_id_clean, gene_symbol) %>%
  distinct()

degs_d_all <- degs_d_all %>%
  inner_join(gene_map_clean, by = c("ensembl_id" = "ensembl_id_clean")) %>%
  select(ensembl_id, gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Save complete table
write.csv(degs_d_all, file = "data_clean/discovery_degs_all.csv", row.names = FALSE)

# Filter discovery DEGs (adjp < 0.05, |log2FC| > 1.5)
degs_d_sig <- degs_d_all %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1.5)

cat("Discovery significant DEGs (|log2FC| > 1.5, adj.P < 0.05):", nrow(degs_d_sig), "\n")


# ==============================================================================
# SECTION 4: Differential Expression in Validation Cohort (limma) | 第 4 部分：独立验证集差异表达分析
# ------------------------------------------------------------------------------
# [EN] Filter GSE13159 for AML and Normal samples, and identify validation DEGs.
# [ZH] 筛选 GSE13159 中的 AML 与 Normal 正常对照样本，并使用 limma 计算验证集差异表达基因。
# ==============================================================================

cat("Running limma on GSE13159 Validation Cohort...\n")

# Filter for AML vs Normal samples in GSE13159
valid_samples <- gse_sample_info %>%
  filter(group %in% c("AML", "Normal")) %>%
  pull(sample_id)

gse_expr_subset <- gse_expr[, valid_samples]
gse_sample_subset <- gse_sample_info %>%
  filter(sample_id %in% valid_samples)

gse_sample_subset$group <- factor(gse_sample_subset$group, levels = c("Normal", "AML"))

cat("Validation subset size:", nrow(gse_expr_subset), "genes x", ncol(gse_expr_subset), "samples\n")

# Fit linear model
mod_v <- model.matrix(~ group, data = gse_sample_subset)
fit_v <- lmFit(gse_expr_subset, mod_v)
fit2_v <- eBayes(fit_v)

# Extract complete results
degs_v_all <- topTable(fit2_v, coef = 2, number = Inf, adjust.method = "BH")
degs_v_all$gene_symbol <- rownames(degs_v_all)
degs_v_all <- degs_v_all %>%
  select(gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# Save complete table
write.csv(degs_v_all, file = "data_clean/gse13159_degs_all.csv", row.names = FALSE)

# Filter validation DEGs
# Note: since the GSE13159 matrix is normalized to range [0, 1],
# the absolute differences are small. We use |logFC| > 0.15 and adj.P < 0.05 as thresholds.
degs_v_sig <- degs_v_all %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 0.15)

cat("Validation significant DEGs (|logFC| > 0.15, adj.P < 0.05):", nrow(degs_v_sig), "\n")


# ==============================================================================
# SECTION 5: Intersect Robust Biomarkers (Cross-Platform Overlap) | 第 5 部分：取交集提取稳定生物标志物
# ------------------------------------------------------------------------------
# [EN] Find genes showing consistent significant differential expression across both cohorts.
# [ZH] 求取两个数据集显著差异基因的交集，过滤掉测序批次效应假阳性，锁定最稳定的核心靶点。
# ==============================================================================

cat("Intersecting DEGs to identify robust biomarkers...\n")

# Make sure direction of fold-change is consistent (both up-regulated or both down-regulated)
discovery_up <- degs_d_sig %>% filter(logFC > 0) %>% pull(gene_symbol)
discovery_down <- degs_d_sig %>% filter(logFC < 0) %>% pull(gene_symbol)

validation_up <- degs_v_sig %>% filter(logFC > 0) %>% pull(gene_symbol)
validation_down <- degs_v_sig %>% filter(logFC < 0) %>% pull(gene_symbol)

overlap_up <- intersect(discovery_up, validation_up)
overlap_down <- intersect(discovery_down, validation_down)
overlap_all <- c(overlap_up, overlap_down)

cat("Robust Up-regulated genes (overlap):", length(overlap_up), "\n")
cat("Robust Down-regulated genes (overlap):", length(overlap_down), "\n")
cat("Total robust overlapping biomarkers:", length(overlap_all), "\n")

# Combine results for overlapping genes
degs_d_overlap <- degs_d_sig %>%
  filter(gene_symbol %in% overlap_all) %>%
  rename(logFC_discovery = logFC, adj.P.Val_discovery = adj.P.Val) %>%
  select(ensembl_id, gene_symbol, logFC_discovery, adj.P.Val_discovery)

degs_v_overlap <- degs_v_sig %>%
  filter(gene_symbol %in% overlap_all) %>%
  rename(logFC_validation = logFC, adj.P.Val_validation = adj.P.Val) %>%
  select(gene_symbol, logFC_validation, adj.P.Val_validation)

robust_biomarkers <- degs_d_overlap %>%
  inner_join(degs_v_overlap, by = "gene_symbol") %>%
  mutate(
    direction = ifelse(gene_symbol %in% overlap_up, "Up", "Down")
  ) %>%
  arrange(desc(abs(logFC_discovery)))

# Save robust overlapping biomarkers
write.csv(robust_biomarkers, file = "data_clean/sig_degs_overlap.csv", row.names = FALSE)
saveRDS(robust_biomarkers, file = "data_clean/sig_degs_overlap.rds")
cat("Robust overlapping biomarkers saved to data_clean/sig_degs_overlap.csv / rds\n")


# ==============================================================================
# SECTION 6: Visualization: Volcano Plots | 第 6 部分：可视化之火山图
# ------------------------------------------------------------------------------
# [EN] Generate volcano plots for both discovery and validation cohorts, highlighting overlapping genes.
# [ZH] 分别绘制发现集与验证集的火山图，并高亮标注交集共有的核心差异基因。
# ==============================================================================

cat("Generating Volcano Plots...\n")

# 1. Discovery Volcano Plot
degs_d_all <- degs_d_all %>%
  mutate(
    change = case_when(
      gene_symbol %in% overlap_up ~ "Robust Up",
      gene_symbol %in% overlap_down ~ "Robust Down",
      adj.P.Val < 0.05 & logFC > 1.5 ~ "Discovery Up Only",
      adj.P.Val < 0.05 & logFC < -1.5 ~ "Discovery Down Only",
      TRUE ~ "Not Significant"
    )
  )

# Select top 10 robust up and down genes to label
top_up_labels <- robust_biomarkers %>% filter(direction == "Up") %>% head(10) %>% pull(gene_symbol)
top_down_labels <- robust_biomarkers %>% filter(direction == "Down") %>% head(10) %>% pull(gene_symbol)
label_genes_d <- degs_d_all %>% filter(gene_symbol %in% c(top_up_labels, top_down_labels))

p_vol_d <- ggplot(degs_d_all, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c(
    "Robust Up" = "#C00000",
    "Robust Down" = "#2E75B6",
    "Discovery Up Only" = "#FFC000",
    "Discovery Down Only" = "#9BC2E6",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_text_repel(
    data = label_genes_d,
    aes(label = gene_symbol),
    size = 3,
    color = "black",
    box_padding = 0.4,
    point_padding = 0.3,
    max.overlaps = 20
  ) +
  theme_bw() +
  labs(
    title = "Discovery Cohort (TCGA + GTEx) Volcano Plot",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/discovery_volcano.png", plot = p_vol_d, width = 8, height = 6.5, dpi = 300)

# 2. Validation Volcano Plot
degs_v_all <- degs_v_all %>%
  mutate(
    change = case_when(
      gene_symbol %in% overlap_up ~ "Robust Up",
      gene_symbol %in% overlap_down ~ "Robust Down",
      adj.P.Val < 0.05 & logFC > 0.15 ~ "Validation Up Only",
      adj.P.Val < 0.05 & logFC < -0.15 ~ "Validation Down Only",
      TRUE ~ "Not Significant"
    )
  )

label_genes_v <- degs_v_all %>% filter(gene_symbol %in% c(top_up_labels, top_down_labels))

p_vol_v <- ggplot(degs_v_all, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c(
    "Robust Up" = "#C00000",
    "Robust Down" = "#2E75B6",
    "Validation Up Only" = "#FFC000",
    "Validation Down Only" = "#9BC2E6",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_text_repel(
    data = label_genes_v,
    aes(label = gene_symbol),
    size = 3,
    color = "black",
    box_padding = 0.4,
    point_padding = 0.3,
    max.overlaps = 20
  ) +
  theme_bw() +
  labs(
    title = "Validation Cohort (GSE13159) Volcano Plot",
    x = "logFC",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/validation_volcano.png", plot = p_vol_v, width = 8, height = 6.5, dpi = 300)
cat("Volcano plots saved successfully!\n")


# ==============================================================================
# SECTION 8: Visualization: Heatmap of Robust Overlapping Genes | 第 8 部分：可视化之共有差异基因热图
# ------------------------------------------------------------------------------
# [EN] Generate heatmap for the top 50 robust overlapping genes across discovery samples.
# [ZH] 在发现集样本中绘制 Top 50 个共有稳定差异基因的表达分化热图。
# ==============================================================================

cat("Generating Heatmap...\n")
# Select top 50 robust biomarkers ordered by significance in discovery
top50_overlap <- robust_biomarkers %>% head(50)
top50_expr_d <- discovery_expr[top50_overlap$ensembl_id, ]
rownames(top50_expr_d) <- top50_overlap$gene_symbol

# Create annotation
sample_annot <- data.frame(
  Group = discovery_sample_info$group,
  Source = discovery_sample_info$source,
  row.names = discovery_sample_info$sample_id
)

annot_colors <- list(
  Group = c("Normal" = "#2E75B6", "AML" = "#C00000"),
  Source = c("TCGA" = "#ED7D31", "GTEx" = "#70AD47")
)

# Plot Heatmap
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
  filename = "results/heatmap_robust_degs.png",
  width = 10,
  height = 8
)
cat("Heatmap saved to results/heatmap_robust_degs.png\n")

cat("Differential expression analysis completed successfully!\n")
