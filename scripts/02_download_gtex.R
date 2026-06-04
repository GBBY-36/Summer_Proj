# ==============================================================================
# SECTION 0: Load Required Packages | 第 0 部分：加载所需的 R 包
# ------------------------------------------------------------------------------
# [EN] Check for and install missing CRAN packages, then load library dependencies for GTEx.
# [ZH] 检查并自动安装缺失的 CRAN 包，然后加载处理 GTEx 数据所需的相应 R 库。
# ==============================================================================

required_cran <- c("data.table", "dplyr", "stringr")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(data.table)
library(dplyr)
library(stringr)


# ==============================================================================
# SECTION 1: Create GTEx Output Folders | 第 1 部分：创建 GTEx 输出目录
# ------------------------------------------------------------------------------
# [EN] Initialize local directory structures for raw GTEx downloads and clean output matrices.
# [ZH] 创建本地文件夹结构，分别用于存放 GTEx 的原始下载数据和清洗后的最终数据。
# ==============================================================================

dir.create("data_raw/gtex", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 2: Define GTEx V8 URLs and File Paths | 第 2 部分：定义 GTEx V8 下载链接与本地路径
# ------------------------------------------------------------------------------
# [EN] Specify URLs and local destinations for GTEx V8 expression data and sample attributes.
# [ZH] 指定 GTEx V8 表达量数据和样本属性数据的 Google Cloud 存储链接及对应的本地保存路径。
# ==============================================================================

gtex_tpm_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"

gtex_sample_url <- "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

gtex_tpm_file <- "data_raw/gtex/GTEx_gene_tpm.gct.gz"
gtex_sample_file <- "data_raw/gtex/GTEx_sample_attributes.txt"


# ==============================================================================
# SECTION 3: Download GTEx Dataset | 第 3 部分：执行 GTEx 数据下载
# ------------------------------------------------------------------------------
# [EN] Set connection timeout and download GTEx expression and annotation files if missing.
# [ZH] 调整下载超时限制，并自动下载缺损的 GTEx 表达矩阵和样本属性文件。
# ==============================================================================

options(timeout = 100000)

download_if_missing <- function(url, destfile) {
  if (!file.exists(destfile)) {
    cat("Downloading:", destfile, "\n")
    download.file(
      url = url,
      destfile = destfile,
      mode = "wb"
    )
  } else {
    cat("File already exists:", destfile, "\n")
  }
}

download_if_missing(gtex_sample_url, gtex_sample_file)
download_if_missing(gtex_tpm_url, gtex_tpm_file)


# ==============================================================================
# SECTION 4: Load and Inspect Sample Attributes | 第 4 部分：加载并查看样本属性注释
# ------------------------------------------------------------------------------
# [EN] Read sample attributes table with data.table and export all unique tissue names.
# [ZH] 使用 data.table 快速读取样本注释表格，并保存数据集中所有可用的组织类型名称。
# ==============================================================================

gtex_sample_attr <- fread(gtex_sample_file)

cat("GTEx sample annotation dimensions:\n")
print(dim(gtex_sample_attr))

cat("GTEx sample annotation columns:\n")
print(colnames(gtex_sample_attr))

# Save available tissue names for checking
available_tissues <- sort(unique(gtex_sample_attr$SMTSD))

write.csv(
  data.frame(tissue = available_tissues),
  file = "data_raw/gtex/gtex_available_tissues.csv",
  row.names = FALSE
)

cat("Available GTEx tissues were saved to data_raw/gtex/gtex_available_tissues.csv\n")


# ==============================================================================
# SECTION 5: Filter Normal Control Tissues | 第 5 部分：筛选对照组正常组织
# ------------------------------------------------------------------------------
# [EN] Match target tissues (Whole Blood and Bone Marrow) and retrieve corresponding sample IDs.
# [ZH] 匹配目标对照组织（全血与骨髓），并提取其在 GTEx 数据库中对应的所有样本 ID。
# ==============================================================================

target_tissues <- c("Whole Blood", "Bone Marrow")

available_target_tissues <- intersect(target_tissues, available_tissues)
missing_target_tissues <- setdiff(target_tissues, available_tissues)

cat("Target tissues found in GTEx:\n")
print(available_target_tissues)

cat("Target tissues not found in GTEx:\n")
print(missing_target_tissues)

if (length(available_target_tissues) == 0) {
  stop("None of the target tissues were found in GTEx sample annotation.")
}

gtex_target_sample_info <- gtex_sample_attr %>%
  filter(SMTSD %in% available_target_tissues)

target_sample_ids <- gtex_target_sample_info$SAMPID

cat("Number of selected GTEx samples:\n")
print(length(target_sample_ids))

cat("Selected tissue counts:\n")
print(table(gtex_target_sample_info$SMTSD))


# ==============================================================================
# SECTION 6: Read Selected Samples from TPM Matrix | 第 6 部分：筛选并加载目标样本表达矩阵
# ------------------------------------------------------------------------------
# [EN] Scan TPM matrix header to match sample IDs and load only relevant expression columns.
# [ZH] 预读取大表达矩阵头部以匹配样本 ID，仅加载指定正常样本的表达量数据以节省内存。
# ==============================================================================

# GCT files have two metadata rows, so we use skip = 2.
# First, read only the header to avoid loading the full large matrix.
gtex_header <- fread(gtex_tpm_file, skip = 2, nrows = 0)

gct_columns <- colnames(gtex_header)

sample_ids_in_gct <- intersect(target_sample_ids, gct_columns)

if (length(sample_ids_in_gct) == 0) {
  stop("No selected GTEx samples were found in the TPM matrix.")
}

keep_cols <- c("Name", "Description", sample_ids_in_gct)

cat("Number of GTEx samples found in TPM matrix:\n")
print(length(sample_ids_in_gct))

# Read only gene ID, gene symbol, and selected tissue samples
gtex_tpm_raw <- fread(
  gtex_tpm_file,
  skip = 2,
  select = keep_cols
)

cat("GTEx selected TPM matrix dimensions before cleaning:\n")
print(dim(gtex_tpm_raw))

# Save raw selected TPM data
fwrite(
  gtex_tpm_raw,
  file = "data_raw/gtex/gtex_selected_normal_tpm_raw.csv"
)


# ==============================================================================
# SECTION 7: Clean Gene Annotation and ID Formats | 第 7 部分：整理与清洗基因注释标识符
# ------------------------------------------------------------------------------
# [EN] Map raw gene names, clean Ensembl version suffixes, and save gene metadata table.
# [ZH] 重命名基因标识列，提取去版本号的 Ensembl ID，并保存清洗后的基因元数据表格。
# ==============================================================================

setnames(
  gtex_tpm_raw,
  old = c("Name", "Description"),
  new = c("ensembl_id", "gene_symbol")
)

gtex_gene_info <- gtex_tpm_raw %>%
  select(ensembl_id, gene_symbol) %>%
  mutate(
    ensembl_id_clean = str_remove(ensembl_id, "\\..*$")
  ) %>%
  distinct()

write.csv(
  gtex_gene_info,
  file = "data_raw/gtex/gtex_gene_info.csv",
  row.names = FALSE
)


# ==============================================================================
# SECTION 8: Aggregate Duplicate Genes using data.table | 第 8 部分：使用 data.table 合并重复基因
# ------------------------------------------------------------------------------
# [EN] Strip version suffixes and merge duplicate rows (e.g., PAR genes) by taking their mean TPM.
# [ZH] 去除 Ensembl ID 版本号后缀，并通过 data.table 对同名重复行（如拟常染色体基因区域）求均值合并。
# ==============================================================================

expr_cols <- setdiff(colnames(gtex_tpm_raw), c("ensembl_id", "gene_symbol"))

# Add unversioned Ensembl ID column to the matrix
gtex_tpm_raw[, ensembl_id_clean := str_remove(ensembl_id, "\\..*$")]

# Filter out empty/NA Ensembl IDs (if any)
gtex_tpm_clean <- gtex_tpm_raw[!is.na(ensembl_id_clean) & ensembl_id_clean != ""]

# Aggregate duplicate Ensembl IDs by taking the mean of TPM values
gtex_tpm_by_gene <- gtex_tpm_clean[, lapply(.SD, mean, na.rm = TRUE), by = ensembl_id_clean, .SDcols = expr_cols]

cat("GTEx TPM matrix dimensions after Ensembl ID aggregation:\n")
print(dim(gtex_tpm_by_gene))


# ==============================================================================
# SECTION 9: Transform Matrix to log2(TPM + 1) | 第 9 部分：构建 log2(TPM + 1) 表达量矩阵
# ------------------------------------------------------------------------------
# [EN] Apply log2(TPM + 1) conversion and save the final clean GTEx expression matrix.
# [ZH] 对 TPM 表达值进行 log2(TPM + 1) 对数转化，并保存为清洗后的格式（CSV/RDS）。
# ==============================================================================

gtex_expr_df <- as.data.frame(gtex_tpm_by_gene)

rownames(gtex_expr_df) <- gtex_expr_df$ensembl_id_clean
gtex_expr_df$ensembl_id_clean <- NULL

gtex_log2_tpm <- log2(gtex_expr_df + 1)

write.csv(
  gtex_log2_tpm,
  file = "data_clean/gtex_normal_expr_log2tpm.csv"
)

saveRDS(
  gtex_log2_tpm,
  file = "data_clean/gtex_normal_expr_log2tpm.rds"
)

cat("GTEx log2(TPM + 1) matrix dimensions:\n")
print(dim(gtex_log2_tpm))



# ==============================================================================
# SECTION 10: Format and Save Sample Metadata | 第 10 部分：格式化并保存样本临床元数据
# ------------------------------------------------------------------------------
# [EN] Extract subject IDs and harmonize sample variables for downstream integration.
# [ZH] 提取患者 ID 并统一样本变量，为下游的数据合并做好临床注释信息的准备。
# ==============================================================================

gtex_sample_info_clean <- gtex_target_sample_info %>%
  filter(SAMPID %in% sample_ids_in_gct) %>%
  transmute(
    sample_id = SAMPID,
    subject_id = str_extract(SAMPID, "^GTEX-[^-]+"),
    source = "GTEx",
    group = "Normal",
    tissue_major = SMTS,
    tissue_detail = SMTSD
  )

write.csv(
  gtex_sample_info_clean,
  file = "data_clean/gtex_normal_sample_info.csv",
  row.names = FALSE
)

saveRDS(
  gtex_sample_info_clean,
  file = "data_clean/gtex_normal_sample_info.rds"
)

# Save used tissues mapping to document which tissues were available and processed
write.csv(
  data.frame(available_target_tissues = available_target_tissues),
  file = "data_clean/gtex_used_tissues.csv",
  row.names = FALSE
)

cat("GTEx sample information preview:\n")
print(head(gtex_sample_info_clean))

cat("GTEx normal tissue download and preprocessing finished.\n")

