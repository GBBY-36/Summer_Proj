# scripts/01_download_tcga_laml.R

# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ------------------------------------------------------------------------------
# [EN] Check for required BioConductor and CRAN packages, install if missing, and load libraries.
# [ZH] 检查所需的 BioConductor 和 CRAN 依赖包，若本地未安装则进行自动安装，并加载相应的 R 库。
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("TCGAbiolinks", "SummarizedExperiment")
required_cran <- c("data.table", "dplyr", "stringr", "jsonlite", "httr", "plyr")

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

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(stringr)


# ==============================================================================
# SECTION 1: In-Memory Namespace Patches for TCGAbiolinks | 第 1 部分：TCGAbiolinks 内存命名空间补丁
# ------------------------------------------------------------------------------
# [EN] Dynamic patch to override internal getBarcodeInfo & exported GDCquery_clinic
#      functions to bypass the GDC API schema changes causing "disease_response doesn't exist" errors.
# [ZH] 动态修补包内部 getBarcodeInfo 与 GDCquery_clinic 函数，绕过因 GDC API 字段变更（缺少 disease_response 列）导致的报错。
# ==============================================================================

ns <- asNamespace("TCGAbiolinks")
pkg_env <- as.environment("package:TCGAbiolinks")

# 1.1 Patch getBarcodeInfo (used internally by GDCprepare)
if (exists("getBarcodeInfo", envir = ns)) {
  unlockBinding("getBarcodeInfo", ns)
  
  patched_getBarcodeInfo <- function (barcode) 
  {
      baseURL <- "https://api.gdc.cancer.gov/cases/?"
      options.pretty <- "pretty=true"
      options.expand <- "expand=samples,follow_ups,project,diagnoses,diagnoses.treatments,annotations,family_histories,demographic,exposures"
      option.size <- paste0("size=", length(barcode))
      options.filter <- paste0("filters=", URLencode("{\"op\":\"or\",\"content\":[{\"op\":\"in\",\"content\":{\"field\":\"cases.submitter_id\",\"value\":["), 
          paste0("\"", paste(barcode, collapse = "\",\"")), URLencode("\"]}},"), 
          URLencode("{\"op\":\"in\",\"content\":{\"field\":\"submitter_sample_ids\",\"value\":["), 
          paste0("\"", paste(barcode, collapse = "\",\"")), URLencode("\"]}},"), 
          URLencode("{\"op\":\"in\",\"content\":{\"field\":\"submitter_aliquot_ids\",\"value\":["), 
          paste0("\"", paste(barcode, collapse = "\",\"")), URLencode("\"]}}"), 
          URLencode("]}"))
      url <- paste0(baseURL, paste(options.pretty, options.expand, 
          option.size, options.filter, sep = "&"))
      
      getURL <- TCGAbiolinks:::getURL
      json <- tryCatch(getURL(url, jsonlite::fromJSON, httr::timeout(600), simplifyDataFrame = TRUE), 
          error = function(e) {
              message(paste("Error: ", e, sep = " "))
              message("We will retry to access GDC again! URL:")
              jsonlite::fromJSON(httr::content(getURL(url, httr::GET, httr::timeout(600)), 
                  as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
          })
      results <- json$data$hits
      if (length(results) == 0) {
          return(data.frame(barcode, stringsAsFactors = FALSE))
      }
      submitter_id <- results$submitter_id
      submitter_aliquot_ids <- results$submitter_aliquot_ids
      
      df <- data.frame(submitter_id = submitter_id, stringsAsFactors = FALSE)
      
      if (!is.null(results$samples)) {
          samples <- data.table::rbindlist(results$samples, fill = TRUE)
          samples <- samples[match(barcode, samples$submitter_id), ]
          samples$sample_submitter_id <- stringr::str_extract_all(samples$submitter_id, 
              paste(barcode, collapse = "|")) %>% unlist %>% as.character
          tryCatch({
              samples$submitter_id <- stringr::str_extract_all(samples$submitter_id, 
                  paste(c(submitter_id, barcode), collapse = "|"), 
                  simplify = TRUE) %>% as.character
          }, error = function(e) {
              samples$submitter_id <- submitter_id
          })
          df <- samples[!is.na(samples$submitter_id), ]
          suppressWarnings({
              df[, c("updated_datetime", "created_datetime")] <- NULL
          })
      }
      if (!is.null(results$diagnoses)) {
          diagnoses <- data.table::rbindlist(lapply(results$diagnoses, function(x) if (is.null(x)) 
              data.frame(NA)
          else x), fill = TRUE)
          cols_to_remove <- intersect(c("updated_datetime", "created_datetime", "state", "days_to_last_follow_up"), colnames(diagnoses))
          for (col in cols_to_remove) {
              diagnoses[[col]] <- NULL
          }
          if (any(grepl("submitter_id", colnames(diagnoses)))) {
              diagnoses$submitter_id <- gsub("-diagnosis|_diagnosis.*|-DIAG|diag-", 
                  "", diagnoses$submitter_id)
          }
          else {
              diagnoses$submitter_id <- submitter_id
          }
          if (!any(df$submitter_id %in% diagnoses$submitter_id)) {
              diagnoses$submitter_id <- NULL
              df <- dplyr::bind_cols(df %>% as.data.frame, diagnoses %>% 
                  as.data.frame)
          }
          else {
              df <- dplyr::left_join(df, diagnoses, by = "submitter_id")
          }
      }
      if (!is.null(results$exposures)) {
          exposures <- data.table::rbindlist(lapply(results$exposures, function(x) if (is.null(x)) 
              data.frame(NA)
          else x), fill = TRUE)
          cols_to_remove <- intersect(c("updated_datetime", "created_datetime", "state"), colnames(exposures))
          for (col in cols_to_remove) {
              exposures[[col]] <- NULL
          }
          if (any(grepl("submitter_id", colnames(exposures)))) {
              exposures$submitter_id <- gsub("-exposure|_exposure.*|-EXP", 
                  "", exposures$submitter_id)
          }
          else {
              exposures$submitter_id <- submitter_id
          }
          if (!any(df$submitter_id %in% exposures$submitter_id)) {
              exposures$submitter_id <- NULL
              df <- dplyr::bind_cols(df, exposures)
          }
          else {
              df <- dplyr::left_join(df, exposures, by = "submitter_id")
          }
      }
      if (!is.null(results$follow_ups)) {
          follow_ups <- data.table::rbindlist(lapply(results$follow_ups, function(x) if (is.null(x)) 
              data.frame(NA)
          else x), fill = TRUE)
          cols_to_remove <- intersect(c("updated_datetime", "created_datetime", "state"), colnames(follow_ups))
          for (col in cols_to_remove) {
              follow_ups[[col]] <- NULL
          }
          if (any(grepl("submitter_id", colnames(follow_ups)))) {
              follow_ups$submitter_id <- gsub("_follow_up.*", "", 
                  follow_ups$submitter_id)
              
              # Safe selection & rename
              select_cols <- c("submitter_id", "days_to_follow_up")
              if ("disease_response" %in% colnames(follow_ups)) {
                  select_cols <- c(select_cols, "disease_response")
              }
              
              if (all(c("submitter_id", "days_to_follow_up") %in% colnames(follow_ups))) {
                  follow_ups_last <- follow_ups %>% dplyr::select(dplyr::all_of(select_cols)) %>% 
                      dplyr::filter(!is.na(submitter_id), !is.na(days_to_follow_up))
                  
                  if (nrow(follow_ups_last) > 0) {
                      follow_ups_last <- follow_ups_last %>% dplyr::group_by(submitter_id) %>% 
                          dplyr::filter(dplyr::row_number() == which.max(days_to_follow_up)) %>% 
                          dplyr::ungroup()
                      
                      if ("disease_response" %in% colnames(follow_ups_last)) {
                          follow_ups_last <- follow_ups_last %>% 
                              dplyr::rename_at(dplyr::vars(disease_response), .funs = function(x) paste0("follow_ups_", x))
                      }
                      
                      follow_ups_last <- follow_ups_last %>% 
                          dplyr::rename(days_to_last_follow_up = days_to_follow_up)
                      df <- dplyr::left_join(df, follow_ups_last, by = "submitter_id")
                  }
              }
          }
      }
      if (!is.null(results$demographic)) {
          demographic <- results$demographic
          cols_to_remove <- intersect(c("updated_datetime", "created_datetime", "state"), colnames(demographic))
          for (col in cols_to_remove) {
              demographic[[col]] <- NULL
          }
          if (any(grepl("submitter_id", colnames(demographic)))) {
              demographic$submitter_id <- gsub("-demographic|_demographic.*|-DEMO|demo-", 
                  "", results$demographic$submitter_id)
          }
          else {
              demographic$submitter_id <- submitter_id
          }
          if (!any(df$submitter_id %in% demographic$submitter_id)) {
              demographic$submitter_id <- NULL
              if (nrow(demographic) < nrow(df)) {
                  demographic <- plyr::ldply(1:length(results$submitter_sample_ids), 
                    .fun = function(x) {
                      demographic[x, ] %>% as.data.frame() %>% 
                        dplyr::slice(rep(dplyr::row_number(), sum(results$submitter_sample_ids[[x]] %in% 
                          barcode)))
                    })
              }
              df <- dplyr::bind_cols(df %>% as.data.frame, demographic)
          }
          else {
              df <- dplyr::left_join(df, demographic, by = "submitter_id")
          }
      }
      df$bcr_patient_barcode <- df$submitter_id %>% as.character()
      projects.info <- results$project
      projects.info <- results$project[, grep("state", colnames(projects.info), 
          invert = TRUE)]
      if (any(submitter_id %in% df$submitter_id)) {
          projects.info <- cbind(submitter_id = submitter_id, projects.info)
          suppressWarnings({
              df <- dplyr::left_join(df, projects.info, by = "submitter_id")
          })
      }
      else {
          if (nrow(projects.info) < nrow(df)) {
              projects.info <- plyr::ldply(1:length(results$submitter_sample_ids), 
                  .fun = function(x) {
                    projects.info[x, ] %>% as.data.frame() %>% 
                      dplyr::slice(rep(dplyr::row_number(), sum(results$submitter_sample_ids[[x]] %in% 
                        barcode)))
                  })
          }
          df <- dplyr::bind_cols(df, projects.info)
      }
      if (any(substr(barcode, 1, stringr::str_length(df$submitter_id)) %in% 
          df$submitter_id)) {
          df <- df[match(substr(barcode, 1, stringr::str_length(df$sample_submitter_id)), 
              df$sample_submitter_id), ]
          df <- df[!is.na(df$submitter_id), ]
      }
      else {
          idx <- sapply(substr(barcode, 1, stringr::str_length(df$submitter_aliquot_ids) %>% 
              max), FUN = function(x) {
              grep(x, df$submitter_aliquot_ids)
          })
          df <- df[idx, ]
      }
      df <- df %>% as.data.frame() %>% dplyr::select(which(colSums(is.na(df)) < 
          nrow(df)))
      return(df)
  }
  
  assign("getBarcodeInfo", patched_getBarcodeInfo, envir = ns)
  lockBinding("getBarcodeInfo", ns)
}

# 1.2 Patch GDCquery_clinic (used directly to fetch clinical data)
if (exists("GDCquery_clinic", envir = ns)) {
  fn_lines <- deparse(TCGAbiolinks::GDCquery_clinic)
  start_idx <- grep("follow_ups_last <- follow_ups %>% dplyr::select", fn_lines, fixed = TRUE)
  end_idx <- grep("df <- dplyr::full_join(df, follow_ups_last", fn_lines, fixed = TRUE) - 1
  
  if (length(start_idx) == 1 && length(end_idx) == 1) {
    patched_code <- c(
      "                select_cols <- c('submitter_id', 'days_to_follow_up')",
      "                if ('disease_response' %in% colnames(follow_ups)) {",
      "                    select_cols <- c(select_cols, 'disease_response')",
      "                }",
      "                follow_ups_last <- follow_ups %>% dplyr::select(dplyr::all_of(select_cols)) %>% ",
      "                  dplyr::filter(!is.na(submitter_id), !is.na(days_to_follow_up))",
      "                if (nrow(follow_ups_last) > 0) {",
      "                    follow_ups_last <- follow_ups_last %>% dplyr::group_by(submitter_id) %>% ",
      "                      dplyr::filter(dplyr::row_number() == which.max(days_to_follow_up)) %>% ",
      "                      dplyr::ungroup()",
      "                    if ('disease_response' %in% colnames(follow_ups_last)) {",
      "                        follow_ups_last <- follow_ups_last %>% dplyr::rename_at(dplyr::vars(disease_response), .funs = function(x) paste0('follow_ups_', x))",
      "                    }",
      "                    follow_ups_last <- follow_ups_last %>% dplyr::rename(days_to_last_follow_up = days_to_follow_up)",
      "                } else {",
      "                    follow_ups_last <- data.frame(submitter_id = character(), days_to_last_follow_up = numeric(), stringsAsFactors = FALSE)",
      "                }"
    )
    
    patched_fn_lines <- c(fn_lines[1:(start_idx-1)], patched_code, fn_lines[(end_idx+1):length(fn_lines)])
    patched_GDCquery_clinic <- eval(parse(text = paste(patched_fn_lines, collapse = "\n")))
    
    # Set environment to package namespace
    environment(patched_GDCquery_clinic) <- ns
    
    # Update Namespace
    unlockBinding("GDCquery_clinic", ns)
    assign("GDCquery_clinic", patched_GDCquery_clinic, envir = ns)
    lockBinding("GDCquery_clinic", ns)
    
    # Update Package Environment (Search Path)
    if (exists("GDCquery_clinic", envir = pkg_env)) {
      unlockBinding("GDCquery_clinic", pkg_env)
      assign("GDCquery_clinic", patched_GDCquery_clinic, envir = pkg_env)
      lockBinding("GDCquery_clinic", pkg_env)
    }
  }
}


# ==============================================================================
# SECTION 2: Create Output Folders | 第 2 部分：创建输出文件夹
# ------------------------------------------------------------------------------
# [EN] Sets up local directories for saving raw downloads and cleaned results.
# [ZH] 创建本地文件夹，分别用于存储下载的原始数据和最终清洗出的结果数据。
# ==============================================================================

dir.create("data_raw/tcga_laml", recursive = TRUE, showWarnings = FALSE)
dir.create("data_clean", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 3: Query TCGA-LAML Gene Expression Data | 第 3 部分：构建 TCGA-LAML 表达数据查询
# ------------------------------------------------------------------------------
# [EN] Configure query filters on the GDC portal to retrieve transcriptomic RNA-Seq (STAR - Counts) files for TCGA-LAML.
# [ZH] 配置 GDC 数据库查询过滤器，定位 TCGA-LAML 项目的转录组定量数据（STAR - Counts 类型）。
# ==============================================================================

query_expr <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)


# ==============================================================================
# SECTION 4: Download TCGA-LAML Raw Files | 第 4 部分：下载 TCGA-LAML 原始数据
# ------------------------------------------------------------------------------
# [EN] Downloads individual sample-level expression TSV files based on the query to local raw folder.
# [ZH] 执行数据下载，将查询到的每个样本的表达定量 TSV 散装文件下载至本地原始数据目录中。
# ==============================================================================

GDCdownload(
  query = query_expr,
  directory = "data_raw/tcga_laml"
)


# ==============================================================================
# SECTION 5: Merge Expression Matrix & Extract Gene Annotations | 第 5 部分：整理表达矩阵与基因注释
# ------------------------------------------------------------------------------
# [EN] Uses GDCprepare() to read and merge downloaded TSVs into a single TPM matrix and extract gene metadata.
# [ZH] 通过 GDCprepare() 批量读取下载的样本 TSV 并自动拼接合并为完整的 TPM 表达量矩阵，同时提取基因元数据。
# ==============================================================================

data_prep <- GDCprepare(
  query = query_expr,
  directory = "data_raw/tcga_laml"
)

cat("Available assays:\n")
print(assayNames(data_prep))

if (!"tpm_unstrand" %in% assayNames(data_prep)) {
  stop("The assay 'tpm_unstrand' was not found. Please check assayNames(data_prep).")
}

# Extract TPM matrix (STAR unstranded TPM)
tcga_tpm <- assay(data_prep, "tpm_unstrand")

# Check TPM dimensions and columns
cat("Dimension of TCGA-LAML TPM matrix:\n")
print(dim(tcga_tpm))

cat("Preview of TPM matrix:\n")
print(tcga_tpm[1:5, 1:min(5, ncol(tcga_tpm))])

# Save raw TPM expression matrix
write.csv(
  tcga_tpm,
  file = "data_raw/tcga_laml/tcga_laml_tpm_raw.csv"
)

saveRDS(
  tcga_tpm,
  file = "data_raw/tcga_laml/tcga_laml_tpm_raw.rds"
)

# Save gene annotation information
gene_info <- as.data.frame(rowRanges(data_prep))
# Clean gene info format slightly (matching original) and add clean Ensembl ID
gene_info <- data.frame(
  gene_id = gene_info$gene_id,
  ensembl_id_clean = str_remove(gene_info$gene_id, "\\..*$"),
  gene_name = gene_info$gene_name,
  gene_type = gene_info$gene_type,
  stringsAsFactors = FALSE
)

write.csv(
  gene_info,
  file = "data_raw/tcga_laml/tcga_laml_gene_info.csv",
  row.names = FALSE
)


# ==============================================================================
# SECTION 6: Download & Format Clinical Data | 第 6 部分：下载与整理临床数据
# ------------------------------------------------------------------------------
# [EN] Downloads patient clinical tables using GDCquery_clinic() and formats list columns for CSV compatibility.
# [ZH] 使用 GDCquery_clinic() 下载完整的患者临床信息，并将嵌套的 List 类型列进行字符化展平，以支持 CSV 导出。
# ==============================================================================

tcga_clinical <- GDCquery_clinic(
  project = "TCGA-LAML",
  type = "clinical"
)

# Convert list columns to character for CSV export (preserves data while allowing CSV writing)
tcga_clinical_csv <- tcga_clinical
list_cols <- sapply(tcga_clinical_csv, is.list)
for (col in names(list_cols)[list_cols]) {
  tcga_clinical_csv[[col]] <- sapply(tcga_clinical_csv[[col]], function(x) paste(x, collapse = ";"))
}

write.csv(
  tcga_clinical_csv,
  file = "data_raw/tcga_laml/tcga_laml_clinical.csv",
  row.names = FALSE
)

saveRDS(
  tcga_clinical,
  file = "data_raw/tcga_laml/tcga_laml_clinical.rds"
)

cat("Clinical data downloaded using TCGAbiolinks.\n")
cat("Dimension of clinical data:\n")
print(dim(tcga_clinical))

cat("Clinical data preview:\n")
print(head(tcga_clinical))


# ==============================================================================
# SECTION 7: Aggregate Duplicates by Unversioned Ensembl ID & Transform to log2(TPM + 1) | 第 7 部分：去版本号合并重复基因与对数转化
# ------------------------------------------------------------------------------
# [EN] Strips version suffix from Ensembl IDs, aggregates duplicate genes (e.g. PAR genes)
#      by mean TPM, and applies log2(TPM + 1) transformation.
# [ZH] 去除 Ensembl ID 的版本号后缀，对重复基因（如 PAR 拟常染色体基因区域）取 TPM 均值合并，并进行 log2(TPM + 1) 转化。
# ==============================================================================

# Convert TCGA matrix to data.table for fast aggregation
tcga_tpm_dt <- as.data.table(tcga_tpm, keep.rownames = "ensembl_id")
tcga_tpm_dt[, ensembl_id_clean := str_remove(ensembl_id, "\\..*$")]

# Retrieve expression columns
tcga_expr_cols <- setdiff(colnames(tcga_tpm_dt), c("ensembl_id", "ensembl_id_clean"))

# Aggregate by clean Ensembl ID (taking the mean of duplicates)
tcga_tpm_by_gene <- tcga_tpm_dt[, lapply(.SD, mean, na.rm = TRUE), by = ensembl_id_clean, .SDcols = tcga_expr_cols]

# Format as data.frame with clean Ensembl IDs as row names
tcga_expr_df <- as.data.frame(tcga_tpm_by_gene)
rownames(tcga_expr_df) <- tcga_expr_df$ensembl_id_clean
tcga_expr_df$ensembl_id_clean <- NULL

tcga_log2_tpm <- log2(tcga_expr_df + 1)

write.csv(
  tcga_log2_tpm,
  file = "data_clean/tcga_laml_expr_log2tpm.csv"
)

saveRDS(
  tcga_log2_tpm,
  file = "data_clean/tcga_laml_expr_log2tpm.rds"
)



# ==============================================================================
# SECTION 8: Save Group and Metadata Mapping | 第 8 部分：保存样本元数据映射
# ------------------------------------------------------------------------------
# [EN] Exports a mapping table of sample barcodes defining group (AML) and source (TCGA) groups.
# [ZH] 导出样本映射元数据表格，指定疾病分组（AML）和样本来源数据库（TCGA）。
# ==============================================================================

tcga_sample_info <- data.frame(
  sample_id = colnames(tcga_tpm),
  patient_id = substr(colnames(tcga_tpm), 1, 12),
  source = "TCGA",
  group = "AML",
  stringsAsFactors = FALSE
)

write.csv(
  tcga_sample_info,
  file = "data_clean/tcga_laml_sample_info.csv",
  row.names = FALSE
)

cat("TCGA-LAML download and preprocessing finished.\n")

