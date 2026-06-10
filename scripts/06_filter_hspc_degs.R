# scripts/06_filter_hspc_degs.R

# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ------------------------------------------------------------------------------
# [EN] Check for and install required BioConductor and CRAN packages, then load libraries.
# [ZH] 检查并自动安装所需的 BioConductor 和 CRAN 包，然后加载所需的 R 库。
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("GEOquery", "Biobase")
required_cran <- c("data.table", "dplyr", "stringr", "ggplot2")

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

library(GEOquery)
library(Biobase)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# Create output folders
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)
dir.create("data_raw/gse24759", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 1: Load Robust Biomarkers (DEGs) | 第 1 部分：读取强健的重叠差异表达基因
# ------------------------------------------------------------------------------
# [EN] Load the 131 robust biomarkers identified in Week 3.
# [ZH] 读取在第三周分析中筛选出的 131 个强健重叠生物标志物（差异表达基因）。
# ==============================================================================

cat("Loading 131 robust overlapping biomarkers...\n")
degs_file <- "data_clean/sig_degs_overlap.rds"
if (!file.exists(degs_file)) {
  stop("Input file 'data_clean/sig_degs_overlap.rds' not found. Please run scripts/05_diff_expression_analysis.R first.")
}
degs_overlap <- readRDS(degs_file)
cat("Loaded", nrow(degs_overlap), "robust biomarkers.\n")


# ==============================================================================
# SECTION 2: Load or Download GSE24759 Dataset | 第 2 部分：读取或下载 GSE24759 数据集
# ------------------------------------------------------------------------------
# [EN] Load GSE24759 from the local directory or download it from GEO if missing.
# [ZH] 读取本地的 GSE24759 数据集，如果不存在则自动从 GEO 下载。
# ==============================================================================

options(timeout = 100000)
gse_matrix_file <- "data_raw/gse24759/GSE24759_series_matrix.txt.gz"

if (!file.exists(gse_matrix_file)) {
  cat("GSE24759 local file not found. Downloading from GEO...\n")
  gse_list <- getGEO(
    GEO = "GSE24759",
    GSEMatrix = TRUE,
    getGPL = FALSE,
    destdir = "data_raw/gse24759"
  )
  gse <- gse_list[[1]]
} else {
  cat("Loading GSE24759 local Series Matrix file...\n")
  gse <- getGEO(filename = gse_matrix_file, getGPL = FALSE)
}

expr_matrix <- exprs(gse)
pheno <- pData(gse)

cat("GSE24759 expression matrix dimensions:", dim(expr_matrix)[1], "probes x", dim(expr_matrix)[2], "samples\n")


# ==============================================================================
# SECTION 3: Identify Normal HSPC Samples | 第 3 部分：识别正常造血干/祖细胞样本
# ------------------------------------------------------------------------------
# [EN] Select primitive cell samples (HSC, CMP, GMP, MEP) from normal controls in GSE24759.
# [ZH] 从 GSE24759 中筛选出代表原始造血干细胞/祖细胞（HSC, CMP, GMP, MEP）的样本。
# ==============================================================================

hspc_patterns <- "Hematopoietic stem cell|Common myeloid progenitor|Granulocyte/monocyte progenitor|Megakaryocyte/ erythroid progenitor"
hspc_idx <- which(grepl(hspc_patterns, pheno$title, ignore.case = TRUE))

if (length(hspc_idx) == 0) {
  stop("Error: No normal HSPC samples identified in GSE24759. Please check phenotype titles.")
}

hspc_samples <- pheno$geo_accession[hspc_idx]
hspc_titles <- pheno$title[hspc_idx]

cat("Identified", length(hspc_samples), "normal HSPC samples in GSE24759.\n")


# ==============================================================================
# SECTION 4: Map Probes to Gene Symbols | 第 4 部分：建立探针与基因名映射表并聚合
# ------------------------------------------------------------------------------
# [EN] Map GSE24759 probe IDs to Gene Symbols using local GPL570 annotation,
#      then aggregate expression matrix to gene-level values using mean.
# [ZH] 使用本地 GPL570 芯片注释文件将 GSE24759 的探针映射到基因名，然后求均值聚合成基因级别的表达矩阵。
# ==============================================================================

gpl_file <- "data_raw/gse13159/GPL570.annot.gz"
if (!file.exists(gpl_file)) {
  cat("Local GPL570 annotation not found at data_raw/gse13159/GPL570.annot.gz. Downloading...\n")
  dir.create("data_raw/gse13159", recursive = TRUE, showWarnings = FALSE)
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    destfile = gpl_file,
    mode = "wb"
  )
}

cat("Parsing GPL570 annotation file...\n")
gpl_data <- parseGEO(gpl_file)
gpl_table <- Table(gpl_data)

probe_gene_mapping <- data.frame(
  probe_id = as.character(gpl_table$ID),
  gene_symbol_raw = as.character(gpl_table$`Gene symbol`),
  stringsAsFactors = FALSE
) %>%
  mutate(
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1],
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---", gene_symbol != "NA")

cat("Mapping probes to Gene Symbols in GSE24759...\n")
expr_dt <- as.data.table(expr_matrix, keep.rownames = "probe_id")
probe_gene_dt <- as.data.table(probe_gene_mapping)

# Merge
expr_gene <- merge(expr_dt, probe_gene_dt[, .(probe_id, gene_symbol)], by = "probe_id")

# Aggregate duplicated genes by mean
cat("Aggregating multi-probe expressions by mean...\n")
expr_gene_mean <- expr_gene[, lapply(.SD, mean), by = gene_symbol, .SDcols = colnames(expr_matrix)]
cat("Aggregated matrix dimensions:", nrow(expr_gene_mean), "genes x", ncol(expr_gene_mean) - 1, "samples\n")

# Convert to data.frame
expr_gene_df <- as.data.frame(expr_gene_mean)
rownames(expr_gene_df) <- expr_gene_df$gene_symbol
expr_gene_df$gene_symbol <- NULL


# ==============================================================================
# SECTION 5: Filter Robust Biomarkers against HSPC baseline | 第 5 部分：根据干细胞基准过滤标志物
# ------------------------------------------------------------------------------
# [EN] Calculate mean HSPC expression for biomarkers, apply threshold filter (7.0),
#      and handle unmeasured genes.
# [ZH] 计算标志物在正常干细胞中的均值表达量，应用阈值（7.0）过滤，并妥善处理未检测基因。
# ==============================================================================

# Extract mean expression in normal HSPC samples
overlap_genes_in_gse <- intersect(degs_overlap$gene_symbol, rownames(expr_gene_df))
cat("Found", length(overlap_genes_in_gse), "out of", nrow(degs_overlap), "robust biomarkers in GSE24759.\n")

hspc_expr <- expr_gene_df[overlap_genes_in_gse, hspc_samples, drop = FALSE]
gene_hspc_mean <- rowMeans(hspc_expr, na.rm = TRUE)

# Set the expression threshold
threshold <- 7.0
cat("Applying normal HSPC expression threshold:", threshold, "(log2 scale)\n")

# Map values back to the full 131 DEGs list
degs_hspc_evaluated <- degs_overlap %>%
  mutate(
    hspc_mean_expr = ifelse(gene_symbol %in% overlap_genes_in_gse, gene_hspc_mean[gene_symbol], NA),
    hspc_status = case_when(
      is.na(hspc_mean_expr) ~ "Not-evaluated",
      hspc_mean_expr >= threshold ~ "HSPC-high",
      hspc_mean_expr < threshold ~ "HSPC-low"
    )
  )

# Separate up and down genes for detailed prints
up_degs <- degs_hspc_evaluated %>% filter(direction == "Up")
down_degs <- degs_hspc_evaluated %>% filter(direction == "Down")

cat("\n--- Up-regulated Candidate Biomarkers Evaluation ---\n")
print(up_degs %>% select(gene_symbol, logFC_discovery, hspc_mean_expr, hspc_status) %>% arrange(desc(hspc_mean_expr)))

cat("\nEvaluation Summary for Up-regulated Genes:\n")
print(table(up_degs$hspc_status))

cat("\nEvaluation Summary for Down-regulated Genes:\n")
print(table(down_degs$hspc_status))

# Filter biomarkers list
# Keep: "HSPC-low" (safe targets) AND "Not-evaluated" (unmeasured, retained by default)
# Exclude: "HSPC-high" (highly expressed in normal stem cells, high toxicity risk)
degs_filtered <- degs_hspc_evaluated %>%
  filter(hspc_status != "HSPC-high")

cat("\nFiltering summary:\n")
cat("Original robust biomarkers:", nrow(degs_overlap), "\n")
cat("Filtered out (HSPC-high):", sum(degs_hspc_evaluated$hspc_status == "HSPC-high"), "\n")
cat("Retained biomarkers:", nrow(degs_filtered), "\n")
cat("  - Safe (HSPC-low):", sum(degs_filtered$hspc_status == "HSPC-low"), "\n")
cat("  - Unmeasured (Not-evaluated):", sum(degs_filtered$hspc_status == "Not-evaluated"), "\n")


# ==============================================================================
# SECTION 6: Save Filtered Outputs | 第 6 部分：保存过滤后的输出结果
# ------------------------------------------------------------------------------
# [EN] Save full and filtered datasets to CSV and RDS formats.
# [ZH] 将完整和过滤后的数据集保存为 CSV 和 RDS 格式。
# ==============================================================================

# Save complete evaluated list (131 genes)
write.csv(degs_hspc_evaluated, file = "data_clean/sig_degs_overlap_hspc_expr.csv", row.names = FALSE)
saveRDS(degs_hspc_evaluated, file = "data_clean/sig_degs_overlap_hspc_expr.rds")

# Save final filtered list
write.csv(degs_filtered, file = "data_clean/sig_degs_hspc_filtered.csv", row.names = FALSE)
saveRDS(degs_filtered, file = "data_clean/sig_degs_hspc_filtered.rds")

cat("\nSaved output files in data_clean/:\n")
cat("  - sig_degs_overlap_hspc_expr.csv / .rds\n")
cat("  - sig_degs_hspc_filtered.csv / .rds\n")


# ==============================================================================
# SECTION 7: Generate Diagnostic Visualization | 第 7 部分：生成诊断可视化图表
# ==============================================================================

cat("\nGenerating diagnostic visualization for Up-regulated genes...\n")

# Prepare data for plotting
plot_df <- up_degs %>%
  mutate(
    # Set a dummy height for Not-evaluated genes so they appear on x-axis
    plot_expr = ifelse(hspc_status == "Not-evaluated", 0, hspc_mean_expr),
    # Sort gene symbols by expression level (unmeasured at the end)
    gene_symbol = factor(gene_symbol, levels = gene_symbol[order(is.na(hspc_mean_expr), hspc_mean_expr, decreasing = FALSE)])
  )

p_hspc <- ggplot(plot_df, aes(x = gene_symbol, y = plot_expr, fill = hspc_status)) +
  geom_col(color = "black", width = 0.7) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "#C00000", linewidth = 1) +
  # Add "N/A" text label on top of the x-axis for not-evaluated genes
  geom_text(
    data = filter(plot_df, hspc_status == "Not-evaluated"),
    aes(y = 0.3, label = "N/A"),
    color = "grey30",
    fontface = "bold",
    size = 3.5
  ) +
  coord_flip() +
  scale_fill_manual(
    name = "HSPC Status",
    values = c(
      "HSPC-low" = "#2E75B6",        # Safe (Blue)
      "HSPC-high" = "#C00000",       # Toxic / High expression (Red)
      "Not-evaluated" = "#D9D9D9"    # Unmeasured / Retained (Grey)
    )
  ) +
  labs(
    title = "Expression of Up-regulated AML Biomarkers in Normal HSPCs",
    subtitle = "Excluding high-risk targets (Normal HSPC Expression >= 7.0)",
    x = "Candidate Gene Symbol",
    y = "Mean Expression in Normal HSPCs (GSE24759, log2 scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic", color = "grey30"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave("results/hspc_expression_up_degs.png", plot = p_hspc, width = 7.5, height = 6, dpi = 300)
cat("Diagnostic plot saved to results/hspc_expression_up_degs.png\n")

cat("\nPipeline script 06 completed successfully!\n")
