# scripts/07_survival_analysis.R

# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ------------------------------------------------------------------------------
# [EN] Load required packages and create output directories.
# [ZH] 加载所需的 R 包并创建输出目录。
# ==============================================================================

library(data.table)
library(dplyr)

dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 1: Load Datasets | 第 1 部分：读取基因表达矩阵与安全特征基因列表
# ------------------------------------------------------------------------------
# [EN] Load TCGA-LAML expression matrix, sample info, and the 106 safety-filtered DEGs.
# [ZH] 读取 TCGA-LAML 表达矩阵、样本元数据以及过滤后的 106 个安全差异表达基因列表。
# ==============================================================================

cat("Loading TCGA-LAML expression matrix and sample metadata...\n")
tcga_expr <- readRDS("data_clean/tcga_laml_expr_log2tpm.rds")
tcga_sample_info <- read.csv("data_clean/tcga_laml_sample_info.csv", stringsAsFactors = FALSE)

cat("Loading 106 safety-filtered biomarkers...\n")
degs_filtered <- read.csv("data_clean/sig_degs_hspc_filtered.csv", stringsAsFactors = FALSE)

cat("TCGA Expression matrix dims:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")
cat("Filtered biomarkers count:", nrow(degs_filtered), "\n")


# ==============================================================================
# SECTION 2: Extract Expression and Perform Median-Based Grouping | 第 2 部分：提取表达值并进行中位数分组
# ------------------------------------------------------------------------------
# [EN] For each of the 106 genes, calculate the median expression across all samples,
#      and classify each sample into High (>= median) or Low (< median) group.
# [ZH] 针对 106 个基因中的每一个，计算其在所有样本中的表达量中位数，并分类为高表达组（>=中位数）或低表达组（<中位数）。
# ==============================================================================

# Align and subset expression matrix to the 106 safety-filtered genes
common_genes <- intersect(degs_filtered$ensembl_id, rownames(tcga_expr))
cat("Found", length(common_genes), "out of", nrow(degs_filtered), "filtered genes in the TCGA expression matrix.\n")

# Subset expression matrix and map row names to Gene Symbols for readability
expr_subset <- tcga_expr[common_genes, , drop = FALSE]
gene_map <- setNames(degs_filtered$gene_symbol, degs_filtered$ensembl_id)
rownames(expr_subset) <- gene_map[rownames(expr_subset)]

# Transpose so that samples are rows and genes are columns
expr_df <- as.data.frame(t(expr_subset))
expr_df$sample_id <- rownames(expr_df)

# Create an empty dataframe to store the High/Low grouping results
# Rows: 151 TCGA samples; Columns: sample_id + 106 genes
group_df <- data.frame(sample_id = expr_df$sample_id, stringsAsFactors = FALSE)

# Loop through each gene to calculate median and assign High/Low groups
cat("Performing median-based grouping for each gene...\n")
for (gene in colnames(expr_df)) {
  if (gene == "sample_id") next
  
  vals <- expr_df[[gene]]
  med_val <- median(vals, na.rm = TRUE)
  
  # Assign High (>= median) and Low (< median)
  group_df[[gene]] <- ifelse(vals >= med_val, "High", "Low")
}

# Preview the grouping results
cat("\nPreview of the grouping matrix (First 5 samples and 5 genes):\n")
print(group_df[1:5, 1:6])


# ==============================================================================
# SECTION 3: Save Grouping Results | 第 3 部分：保存分组结果
# ------------------------------------------------------------------------------
# [EN] Save the grouping results to both CSV and RDS formats in data_clean/.
# [ZH] 将分组结果保存为 CSV 和 RDS 格式于 data_clean/ 目录。
# ==============================================================================

write.csv(group_df, file = "data_clean/tcga_laml_gene_groups.csv", row.names = FALSE)
saveRDS(group_df, file = "data_clean/tcga_laml_gene_groups.rds")

cat("\nGrouping completed successfully! Saved output files:\n")
cat("  - data_clean/tcga_laml_gene_groups.csv\n")
cat("  - data_clean/tcga_laml_gene_groups.rds\n")
