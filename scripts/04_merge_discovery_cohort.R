# ==============================================================================
# SECTION 0: Load Required Packages | 第 0 部分：加载所需的 R 包
# ------------------------------------------------------------------------------
# [EN] Load libraries for data merging and processing.
# [ZH] 加载进行数据合并和处理所需的 R 包。
# ==============================================================================

library(data.table)
library(dplyr)


# ==============================================================================
# SECTION 1: Load TCGA and GTEx Preprocessed Matrices | 第 1 部分：加载已清洗的 TCGA 和 GTEx 表达矩阵
# ------------------------------------------------------------------------------
# [EN] Read the RDS files for TCGA and GTEx expression and sample metadata.
# [ZH] 读取清洗好的 TCGA 和 GTEx 基因表达矩阵及对应的样本注释文件。
# ==============================================================================

tcga_expr <- readRDS("data_clean/tcga_laml_expr_log2tpm.rds")
tcga_sample_info <- read.csv("data_clean/tcga_laml_sample_info.csv", stringsAsFactors = FALSE)

gtex_expr <- readRDS("data_clean/gtex_normal_expr_log2tpm.rds")
gtex_sample_info <- readRDS("data_clean/gtex_normal_sample_info.rds")


# ==============================================================================
# SECTION 2: Intersect Shared Genes | 第 2 部分：求取共有基因交集
# ------------------------------------------------------------------------------
# [EN] Find the overlapping gene identifiers and subset both matrices.
# [ZH] 获取两个表达矩阵共同拥有的基因 ID，并将两个矩阵均筛选至仅包含共有基因。
# ==============================================================================

common_genes <- intersect(rownames(tcga_expr), rownames(gtex_expr))
cat("Number of common genes found:", length(common_genes), "\n")

tcga_expr_subset <- tcga_expr[common_genes, , drop = FALSE]
gtex_expr_subset <- gtex_expr[common_genes, , drop = FALSE]


# ==============================================================================
# SECTION 3: Merge Expression Matrices Column-wise | 第 3 部分：按列合并表达矩阵
# ------------------------------------------------------------------------------
# [EN] Bind columns of the subsetted matrices to create the combined discovery matrix.
# [ZH] 将筛选后的两个表达矩阵按列合并，构建统一的发现集表达矩阵。
# ==============================================================================

discovery_expr <- cbind(tcga_expr_subset, gtex_expr_subset)
cat("Combined discovery expression matrix dimensions:\n")
print(dim(discovery_expr))


# ==============================================================================
# SECTION 4: Harmonize and Merge Sample Metadata | 第 4 部分：对齐并合并样本元数据
# ------------------------------------------------------------------------------
# [EN] Standardize variables to match columns and bind them row-wise.
# [ZH] 对齐两个元数据表的变量列，并将其按行合并以获得完整的样本注释表格。
# ==============================================================================

tcga_sample_clean <- tcga_sample_info %>%
  transmute(
    sample_id = sample_id,
    source = source,
    group = group,
    detail = "Acute Myeloid Leukemia"
  )

gtex_sample_clean <- gtex_sample_info %>%
  transmute(
    sample_id = sample_id,
    source = source,
    group = group,
    detail = tissue_detail
  )

discovery_sample_info <- rbind(tcga_sample_clean, gtex_sample_clean)
cat("Combined discovery sample metadata dimensions:\n")
print(dim(discovery_sample_info))

# Add a strict sample order check and alignment
cat("Aligning sample metadata rows to match expression matrix columns...\n")
discovery_sample_info <- discovery_sample_info[
  match(colnames(discovery_expr), discovery_sample_info$sample_id),
]

if (any(is.na(discovery_sample_info$sample_id))) {
  stop("Error: Some expression samples were not found in sample metadata.")
}

if (!identical(colnames(discovery_expr), discovery_sample_info$sample_id)) {
  stop("Error: Sample order mismatch between expression matrix columns and sample metadata.")
}
cat("Sample alignment check passed successfully!\n")


# ==============================================================================
# SECTION 5: Save Merged Discovery Cohort Outputs | 第 5 部分：保存合并的发现集数据文件
# ------------------------------------------------------------------------------
# [EN] Export the final merged discovery matrix and sample metadata.
# [ZH] 将合并好的最终发现集基因表达矩阵与样本注释表导出为 CSV 与 RDS 格式。
# ==============================================================================

write.csv(
  discovery_expr,
  file = "data_clean/discovery_expr_log2tpm.csv"
)

saveRDS(
  discovery_expr,
  file = "data_clean/discovery_expr_log2tpm.rds"
)

write.csv(
  discovery_sample_info,
  file = "data_clean/discovery_sample_info.csv",
  row.names = FALSE
)

saveRDS(
  discovery_sample_info,
  file = "data_clean/discovery_sample_info.rds"
)

cat("Combined discovery cohort saved successfully!\n")


