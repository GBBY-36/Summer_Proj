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


# ==============================================================================
# SECTION 3: Load Survival Data and Merge | 第 3 部分：加载生存数据并合并
# ------------------------------------------------------------------------------
# [EN] Load TCGA-LAML clinical database, calculate overall survival time and status,
#      and merge with the expression-based High/Low grouping results.
# [ZH] 加载 TCGA-LAML 临床生存数据，计算生存时间与生存状态，并与高/低表达分组结果合并。
# ==============================================================================

cat("\nLoading clinical survival data...\n")
clinical_file <- "data_raw/tcga_laml/tcga_laml_clinical.rds"
if (!file.exists(clinical_file)) {
  stop("Clinical data file 'data_raw/tcga_laml/tcga_laml_clinical.rds' not found. Please run scripts/01_download_tcga_laml.R first.")
}
tcga_clinical <- readRDS(clinical_file)

# Prepare clean survival metrics
# time: days_to_death if dead, days_to_last_follow_up if alive
# event: 1 if dead, 0 if alive
clinical_survival <- tcga_clinical %>%
  mutate(
    patient_id = submitter_id,
    time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    event = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time) & time > 0) %>%
  select(patient_id, time, event)

cat("Valid clinical records with survival info:", nrow(clinical_survival), "\n")

# Merge grouping results with survival metadata
# group_df has sample_id; we map to patient_id via tcga_sample_info
merged_data <- group_df %>%
  inner_join(tcga_sample_info[, c("sample_id", "patient_id")], by = "sample_id") %>%
  inner_join(clinical_survival, by = "patient_id")

cat("Samples aligned with expression grouping and survival data:", nrow(merged_data), "\n")


# ==============================================================================
# SECTION 4: Perform Kaplan-Meier & Log-rank Test | 第 4 部分：执行 Kaplan-Meier 和 Log-rank 检验
# ------------------------------------------------------------------------------
# [EN] For each of the 106 genes, run the Log-rank test (survdiff) to evaluate
#      whether there is a significant difference in survival between High and Low groups.
# [ZH] 针对 106 个基因，逐一进行 Log-rank 检验，分析高/低表达组之间是否存在显著的生存时间差异。
# ==============================================================================

library(survival)

cat("Running Log-rank tests for 106 safety-filtered genes...\n")
logrank_results <- list()

# Get the list of gene symbols (columns in merged_data, excluding sample_id, patient_id, time, event)
genes_list <- setdiff(colnames(merged_data), c("sample_id", "patient_id", "time", "event"))

for (gene in genes_list) {
  # Perform log-rank test
  # Formula: Surv(time, event) ~ gene_group
  # Wrap the gene in backticks to prevent syntax errors due to hyphens in gene symbols (e.g. NKX2-3)
  formula_str <- paste0("Surv(time, event) ~ ", "`", gene, "`")
  fit_diff <- tryCatch({
    survdiff(as.formula(formula_str), data = merged_data)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_diff)) next
  
  # Calculate p-value from Chi-square statistic
  df <- length(fit_diff$n) - 1
  p_val <- pchisq(fit_diff$chisq, df = df, lower.tail = FALSE)
  
  logrank_results[[gene]] <- data.frame(
    gene_symbol = gene,
    chisq = fit_diff$chisq,
    p_value = p_val,
    stringsAsFactors = FALSE
  )
}

# Bind list to dataframe and sort by p-value
logrank_df <- do.call(rbind, logrank_results) %>%
  arrange(p_value)

# Classify significance
logrank_df <- logrank_df %>%
  mutate(
    significance = ifelse(p_value < 0.05, "Significant", "Not Significant")
  )

cat("Log-rank analysis completed!\n")
cat("Total genes analyzed:", nrow(logrank_df), "\n")
cat("Significant prognostic genes (P < 0.05):", sum(logrank_df$p_value < 0.05), "\n")

cat("\nTop 10 most prognostic genes by Log-rank P-value:\n")
print(head(logrank_df, 10))


# ==============================================================================
# SECTION 5: Perform Univariate Cox Regression | 第 5 部分：执行单变量 Cox 回归分析
# ------------------------------------------------------------------------------
# [EN] Fit a univariate Cox proportional hazards regression model for each gene.
#      Evaluate the Hazard Ratio (HR) of High vs Low expression group.
# [ZH] 针对 106 个基因，逐一建立单变量 Cox 比例风险回归模型，评估高表达组相对于低表达组的风险比（HR）。
# ==============================================================================

cat("\nRunning Univariate Cox regression models...\n")
cox_results <- list()

for (gene in genes_list) {
  # Set the reference group to 'Low' so HR represents High vs Low
  merged_data[[gene]] <- factor(merged_data[[gene]], levels = c("Low", "High"))
  
  formula_str <- paste0("Surv(time, event) ~ ", "`", gene, "`")
  
  fit_cox <- tryCatch({
    coxph(as.formula(formula_str), data = merged_data)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_cox)) next
  
  sum_cox <- summary(fit_cox)
  coef_info <- sum_cox$coefficients
  conf_info <- sum_cox$conf.int
  
  hr_val <- coef_info[1, "exp(coef)"]
  p_val_cox <- coef_info[1, "Pr(>|z|)"]
  ci_lower <- conf_info[1, "lower .95"]
  ci_upper <- conf_info[1, "upper .95"]
  
  cox_results[[gene]] <- data.frame(
    gene_symbol = gene,
    coef = coef_info[1, "coef"],
    HR = hr_val,
    HR_CI_lower = ci_lower,
    HR_CI_upper = ci_upper,
    p_value = p_val_cox,
    stringsAsFactors = FALSE
  )
}

# Bind list to dataframe and sort by p-value
cox_df <- do.call(rbind, cox_results) %>%
  arrange(p_value)

# Apply user's threshold: HR > 1.5 and p_value < 0.05
cox_filtered_df <- cox_df %>%
  filter(HR > 1.5 & p_value < 0.05)

cat("Univariate Cox regression completed!\n")
cat("Total genes analyzed:", nrow(cox_df), "\n")
cat("Genes with HR > 1.5 & P < 0.05:", nrow(cox_filtered_df), "\n")

cat("\nSignificant Prognostic Hazard Genes (HR > 1.5 & P < 0.05):\n")
print(cox_filtered_df)


# ==============================================================================
# SECTION 6: Save All Results | 第 6 部分：保存分组与所有生存统计结果
# ==============================================================================

# Save grouping result
write.csv(group_df, file = "data_clean/tcga_laml_gene_groups.csv", row.names = FALSE)
saveRDS(group_df, file = "data_clean/tcga_laml_gene_groups.rds")

# Save Log-rank stats
write.csv(logrank_df, file = "data_clean/survival_logrank_results.csv", row.names = FALSE)
saveRDS(logrank_df, file = "data_clean/survival_logrank_results.rds")

# Save Cox regression stats
write.csv(cox_df, file = "data_clean/survival_cox_results.csv", row.names = FALSE)
saveRDS(cox_df, file = "data_clean/survival_cox_results.rds")

# Save filtered significant hazard genes (HR > 1.5 & P < 0.05)
write.csv(cox_filtered_df, file = "data_clean/survival_cox_sig_danger_genes.csv", row.names = FALSE)

cat("\nSaved output files in data_clean/:\n")
cat("  - tcga_laml_gene_groups.csv / .rds\n")
cat("  - survival_logrank_results.csv / .rds\n")
cat("  - survival_cox_results.csv / .rds\n")
cat("  - survival_cox_sig_danger_genes.csv\n")


