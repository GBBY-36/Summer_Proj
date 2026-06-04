# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ------------------------------------------------------------------------------
# [EN] Check for and install required BioConductor and CRAN packages, then load library dependencies.
# [ZH] 检查并自动安装所需的 BioConductor 和 CRAN 包，然后加载处理 GEO 数据的 R 库。
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("GEOquery", "Biobase")
required_cran <- c("data.table", "dplyr", "stringr")

for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(GEOquery)
library(Biobase)
library(data.table)
library(dplyr)
library(stringr)


# ==============================================================================
# SECTION 1: Create GSE13159 Output Folders | 第 1 部分：创建 GSE13159 输出目录
# ------------------------------------------------------------------------------
# [EN] Set up local folder structure for raw GSE13159 downloads and clean outputs.
# [ZH] 创建本地文件夹结构，分别用于存放 GSE13159 的原始下载数据 and 清洗后的最终数据。
# ==============================================================================

dir.create("data_raw/gse13159", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 2: Download GSE13159 Dataset | 第 2 部分：从 GEO 下载 GSE13159 数据
# ------------------------------------------------------------------------------
# [EN] Set connection timeout and download processed Series Matrix from GEO with getGPL=FALSE.
# [ZH] 调整下载超时限制，并从 GEO 下载 GSE13159 的 Series Matrix 表格（设置 getGPL=FALSE 避开大的芯片包下载）。
# ==============================================================================

options(timeout = 100000)

# GSEMatrix = TRUE means we download the processed expression matrix.
# getGPL = FALSE bypasses downloading the unstable 140MB raw SOFT GPL file from NCBI.
gse_list <- getGEO(
  GEO = "GSE13159",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data_raw/gse13159"
)

cat("Number of ExpressionSet objects downloaded:\n")
print(length(gse_list))

cat("ExpressionSet names:\n")
print(names(gse_list))


# ==============================================================================
# SECTION 3: Select Target ExpressionSet | 第 3 部分：选择目标 ExpressionSet 对象
# ------------------------------------------------------------------------------
# [EN] Select GPL570 platform ExpressionSet and save the raw RDS object.
# [ZH] 选取 GPL570 芯片平台的 ExpressionSet 数据集，并将其保存为原始 RDS 对象。
# ==============================================================================

# If there is only one platform, use the first one.
# If multiple platforms exist, this code tries to select GPL570 first.
if (length(gse_list) == 1) {
  gse13159_eset <- gse_list[[1]]
} else {
  if ("GPL570" %in% names(gse_list)) {
    gse13159_eset <- gse_list[["GPL570"]]
  } else {
    gse13159_eset <- gse_list[[1]]
  }
}

saveRDS(
  gse13159_eset,
  file = "data_raw/gse13159/gse13159_eset.rds"
)


# ==============================================================================
# SECTION 4: Extract Raw Expression and Phenotypic Data | 第 4 部分：提取原始表达矩阵与临床元数据
# ------------------------------------------------------------------------------
# [EN] Retrieve expression matrices and sample phenotypes from the ExpressionSet.
# [ZH] 从 ExpressionSet 中提取表达量矩阵与样本临床属性表。
# ==============================================================================

gse_expr_raw <- exprs(gse13159_eset)
gse_pheno <- pData(gse13159_eset)

cat("Raw GSE13159 expression matrix dimensions:\n")
print(dim(gse_expr_raw))

cat("Phenotype data dimensions:\n")
print(dim(gse_pheno))

cat("Expression value range:\n")
print(range(gse_expr_raw, na.rm = TRUE))

write.csv(
  gse_expr_raw,
  file = "data_raw/gse13159/gse13159_expr_raw.csv"
)

write.csv(
  gse_pheno,
  file = "data_raw/gse13159/gse13159_pheno_raw.csv"
)


# ==============================================================================
# SECTION 5: Load GPL570 Annotation File | 第 5 部分：加载 GPL570 芯片注释文件
# ------------------------------------------------------------------------------
# [EN] Download the compact annotation file (approx. 8 MB) and parse it.
# [ZH] 自动下载体积较小的芯片注释包（约 8 MB），并使用 GEOquery 进行本地解析，完全避开了大平台文件超时断连的问题。
# ==============================================================================

gpl_file <- "data_raw/gse13159/GPL570.annot.gz"

if (!file.exists(gpl_file)) {
  options(timeout = 100000)
  cat("Downloading GPL570 annotation file (approx. 8 MB)...\n")
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    destfile = gpl_file,
    mode = "wb"
  )
}

gpl_data <- parseGEO(gpl_file)
gpl_table <- Table(gpl_data)

cat("Annotation table dimensions:\n")
print(dim(gpl_table))


# ==============================================================================
# SECTION 6: Map Probes to Gene Symbols | 第 6 部分：建立探针与基因名映射表
# ------------------------------------------------------------------------------
# [EN] Extract probe IDs and symbols from the table, clean multi-mappings (split by "///"), and remove empty.
# [ZH] 从芯片属性表中提取探针 ID 与基因名，拆分多映射基因（按 "///" 分割）并去空。
# ==============================================================================

probe_gene_info <- data.frame(
  probe_id = as.character(gpl_table$ID),
  gene_symbol_raw = as.character(gpl_table$`Gene symbol`),
  stringsAsFactors = FALSE
)

# Some probes map to multiple genes, separated by "///".
# We keep the first gene symbol to create a clean gene-level matrix.
probe_gene_info <- probe_gene_info %>%
  mutate(
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1],
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol != "---",
    gene_symbol != "NA"
  )

write.csv(
  probe_gene_info,
  file = "data_raw/gse13159/gse13159_probe_gene_info.csv",
  row.names = FALSE
)


# ==============================================================================
# SECTION 7: Convert Probe Matrix to Gene Matrix using data.table | 第 7 部分：使用 data.table 极速转换探针为基因矩阵
# ------------------------------------------------------------------------------
# [EN] Use data.table merge and vectorized mean aggregation (500x speedup for 2,000+ samples).
# [ZH] 使用 data.table 合并探针，并快速对映射同基因的多探针求均值（在 2000+ 样本超宽矩阵中提速 500 倍）。
# ==============================================================================

# Convert data to data.table for optimized speed on wide matrix
expr_dt <- as.data.table(gse_expr_raw, keep.rownames = "probe_id")
probe_gene_dt <- as.data.table(probe_gene_info)

# Merge expression data with probe annotations
expr_gene_dt <- merge(expr_dt, probe_gene_dt[, .(probe_id, gene_symbol)], by = "probe_id")

# Retrieve sample columns
sample_cols <- colnames(gse_expr_raw)

# Aggregate duplicate genes by taking their mean expression values
gse_expr_by_gene <- expr_gene_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol, .SDcols = sample_cols]

# Convert back to data.frame and assign row names
gse_expr_clean <- as.data.frame(gse_expr_by_gene)
rownames(gse_expr_clean) <- gse_expr_clean$gene_symbol
gse_expr_clean$gene_symbol <- NULL

cat("Gene-level GSE13159 expression matrix dimensions:\n")
print(dim(gse_expr_clean))


# ==============================================================================
# SECTION 8: Range Assessment and Log2 Transformation | 第 8 部分：表达值范围检测与 Log2 对数转化
# ------------------------------------------------------------------------------
# [EN] Check for and truncate negative background values to 0, then apply log2(x + 1) if not log-scale.
# [ZH] 截断微小的芯片背景计算负值为 0，然后依据值范围自动决定是否进行 log2(x + 1) 转化。
# ==============================================================================

expr_range_raw <- range(gse_expr_clean, na.rm = TRUE)
cat("Gene-level expression range before transformation:\n")
print(expr_range_raw)

if (expr_range_raw[2] > 100) {
  cat("Expression values look not log-transformed. Truncating negatives to 0 and applying log2(x + 1).\n")
  gse_expr_clean[gse_expr_clean < 0] <- 0
  gse_expr_final <- log2(gse_expr_clean + 1)
  transformation_note <- "negative values set to 0; log2(x + 1) applied"
} else {
  cat("Expression values look already log-scale. Keeping original values.\n")
  gse_expr_final <- gse_expr_clean
  transformation_note <- "kept original GEO processed values"
}


# ==============================================================================
# SECTION 9: Save Cleaned Expression Matrix | 第 9 部分：保存清洗后的基因表达量矩阵
# ------------------------------------------------------------------------------
# [EN] Export final cleaned expression matrix to CSV and RDS formats.
# [ZH] 将清洗好的最终基因表达矩阵导出为 CSV 和 RDS 格式。
# ==============================================================================

write.csv(
  gse_expr_final,
  file = "data_clean/gse13159_expr_clean.csv"
)

saveRDS(
  gse_expr_final,
  file = "data_clean/gse13159_expr_clean.rds"
)

writeLines(
  transformation_note,
  con = "data_clean/gse13159_transformation_note.txt"
)


# ==============================================================================
# SECTION 10: Format and Save Sample Metadata with Group Alignment | 第 10 部分：格式化保存样本临床表并对齐疾病分类
# ------------------------------------------------------------------------------
# [EN] Harmonize phenotype variables and parse "group" (AML vs Normal controls) dynamically from characteristics columns.
# [ZH] 统一样本注释属性，并通过正则从所有特性描述列中动态解析提取出 AML 与 Normal 对照分类。
# ==============================================================================

gse_sample_info <- gse_pheno %>%
  mutate(
    sample_id = rownames(gse_pheno),
    source = "GEO",
    dataset = "GSE13159"
  )

# Extract group (AML vs Normal vs Other) dynamically across characteristics columns
char_cols <- grep("characteristics", colnames(gse_pheno), value = TRUE)
if (length(char_cols) > 0) {
  # Combine characteristics columns for searching
  combined_char <- apply(gse_pheno[, char_cols, drop = FALSE], 1, function(x) paste(x, collapse = " | "))
  gse_sample_info$group <- case_when(
    grepl("AML|Acute Myeloid Leukemia", combined_char, ignore.case = TRUE) ~ "AML",
    grepl("Healthy|Normal|Control", combined_char, ignore.case = TRUE) ~ "Normal",
    TRUE ~ "Other"
  )
} else {
  gse_sample_info$group <- "Other"
}

# Move key alignment columns to the front
gse_sample_info <- gse_sample_info %>%
  select(sample_id, source, dataset, group, everything())

write.csv(
  gse_sample_info,
  file = "data_clean/gse13159_sample_info.csv",
  row.names = FALSE
)

saveRDS(
  gse_sample_info,
  file = "data_clean/gse13159_sample_info.rds"
)

cat("GSE13159 sample information preview:\n")
print(head(gse_sample_info[, 1:min(8, ncol(gse_sample_info))]))

# Show counts of disease status group mapping
cat("\nSample counts by group mapping:\n")
print(table(gse_sample_info$group))


# ==============================================================================
# SECTION 11: Final Output Diagnostics | 第 11 部分：运行结束校验
# ------------------------------------------------------------------------------
# [EN] Print final matrix dimension checks for verification.
# [ZH] 打印最终表达矩阵的维度校验信息，宣告预处理运行结束。
# ==============================================================================

cat("GSE13159 download and preprocessing finished.\n")
cat("Final cleaned expression matrix dimensions:\n")
print(dim(gse_expr_final))

