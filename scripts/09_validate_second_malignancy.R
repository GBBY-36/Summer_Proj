# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_bioc <- c("GEOquery", "Biobase")
required_cran <- c("data.table", "dplyr", "stringr", "survival", "ggplot2", "gridExtra")

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
library(survival)
library(ggplot2)
library(gridExtra)

# Create output directories
dir.create("Week6", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ==============================================================================
# SECTION 1: Load and Preprocess GSE10846 Dataset | 第 1 部分：加载并预处理 GSE10846
# ==============================================================================

eset_file <- "data_raw/gse10846/gse10846_eset.rds"
if (!file.exists(eset_file)) {
  dir.create("data_raw/gse10846", recursive = TRUE, showWarnings = FALSE)
  options(timeout = 100000)
  cat("Downloading GSE10846 Series Matrix from GEO...\n")
  gse_list <- getGEO(GEO = "GSE10846", GSEMatrix = TRUE, getGPL = FALSE, destdir = "data_raw/gse10846")
  if ("GPL570" %in% names(gse_list)) {
    gse10846_eset <- gse_list[["GPL570"]]
  } else {
    gse10846_eset <- gse_list[[1]]
  }
  saveRDS(gse10846_eset, file = eset_file)
} else {
  cat("Loading local GSE10846 ExpressionSet...\n")
  gse10846_eset <- readRDS(eset_file)
}

gse_expr_raw <- exprs(gse10846_eset)
gse_pheno <- pData(gse10846_eset)

# Map Probes to Gene Symbols
gpl_file <- "data_raw/gse13159/GPL570.annot.gz"
if (!file.exists(gpl_file)) {
  dir.create(dirname(gpl_file), recursive = TRUE, showWarnings = FALSE)
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    destfile = gpl_file,
    mode = "wb"
  )
}

gpl_data <- parseGEO(gpl_file)
gpl_table <- Table(gpl_data)

probe_gene_info <- data.frame(
  probe_id = as.character(gpl_table$ID),
  gene_symbol_raw = as.character(gpl_table$`Gene symbol`),
  stringsAsFactors = FALSE
) %>%
  mutate(
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1],
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---", gene_symbol != "NA")

expr_dt <- as.data.table(gse_expr_raw, keep.rownames = "probe_id")
probe_gene_dt <- as.data.table(probe_gene_info)
expr_gene_dt <- merge(expr_dt, probe_gene_dt[, .(probe_id, gene_symbol)], by = "probe_id")

sample_cols <- colnames(gse_expr_raw)
gse_expr_by_gene <- expr_gene_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol, .SDcols = sample_cols]

gse_expr_clean <- as.data.frame(gse_expr_by_gene)
rownames(gse_expr_clean) <- gse_expr_clean$gene_symbol
gse_expr_clean$gene_symbol <- NULL

# Truncate negative values and apply log2 transformation if not log scale
expr_range <- range(gse_expr_clean, na.rm = TRUE)
if (expr_range[2] > 100) {
  gse_expr_clean[gse_expr_clean < 0] <- 0
  gse_expr_final <- log2(gse_expr_clean + 1)
} else {
  gse_expr_final <- gse_expr_clean
}

saveRDS(gse_expr_final, file = "data_clean/gse10846_expr_clean.rds")

# ==============================================================================
# SECTION 2: Parse and Clean Phenotype Data | 第 2 部分：解析临床生存数据
# ==============================================================================

# Parse characteristics columns dynamically to avoid column offset errors
cat("Parsing DLBCL clinical data...\n")
parsed_clinical <- do.call(rbind, lapply(1:nrow(gse_pheno), function(i) {
  row_text <- paste(gse_pheno[i, ], collapse = "; ")
  status <- str_trim(str_match(row_text, "Follow up status: ([^;]+)")[, 2])
  years <- as.numeric(str_trim(str_match(row_text, "Follow up years: ([^;]+)")[, 2]))
  chemo <- str_trim(str_match(row_text, "Chemotherapy: ([^;]+)")[, 2])
  diag <- str_trim(str_match(row_text, "Final microarray diagnosis: ([^;]+)")[, 2])
  data.frame(
    sample_id = rownames(gse_pheno)[i],
    status = status,
    years = years,
    chemo = chemo,
    diag = diag,
    stringsAsFactors = FALSE
  )
}))

parsed_clinical <- parsed_clinical %>%
  mutate(
    event = ifelse(status == "DEAD", 1, ifelse(status == "ALIVE", 0, NA))
  )

# Filter samples with valid survival data and standard chemotherapy regimens
valid_idx <- which(!is.na(parsed_clinical$event) & !is.na(parsed_clinical$years) & 
                     parsed_clinical$chemo %in% c("CHOP-Like Regimen", "R-CHOP-Like Regimen"))

clean_expr <- gse_expr_final[, parsed_clinical$sample_id[valid_idx]]
clean_cli <- parsed_clinical[valid_idx, ]
clean_cli$chemo <- factor(clean_cli$chemo, levels = c("CHOP-Like Regimen", "R-CHOP-Like Regimen"))

saveRDS(clean_cli, file = "data_clean/gse10846_sample_info.rds")

# ==============================================================================
# SECTION 3: Perform Survival Analysis | 第 3 部分：执行生存分析
# ==============================================================================

targets <- c("IL1R2", "VNN2", "SLC15A3", "EGFL7", "CXCR1", "CXCR2", "CD14", "GZMB", "FPR2", "CMTM2")
dlbcl_survival_results <- list()

cat("Running survival analyses for top target candidates in DLBCL...\n")
for (g in targets) {
  val <- as.numeric(clean_expr[g, ])
  med <- median(val, na.rm = TRUE)
  grp <- ifelse(val >= med, "High", "Low")
  
  # Aligned clinical variables
  df_temp <- clean_cli
  df_temp$gene_expr <- val
  df_temp$gene_group <- factor(grp, levels = c("Low", "High"))
  
  # Log-rank test
  logrank_fit <- survdiff(Surv(years, event) ~ gene_group, data = df_temp)
  df <- length(logrank_fit$n) - 1
  logrank_p <- pchisq(logrank_fit$chisq, df = df, lower.tail = FALSE)
  
  # Univariate Cox
  uni_fit <- coxph(Surv(years, event) ~ gene_group, data = df_temp)
  uni_summary <- summary(uni_fit)
  uni_hr <- uni_summary$conf.int[1, 1]
  uni_lower <- uni_summary$conf.int[1, 3]
  uni_upper <- uni_summary$conf.int[1, 4]
  uni_p <- uni_summary$waldtest[3]
  
  # Multivariate adjusted Cox (adjust for chemo regimen)
  adj_fit <- coxph(Surv(years, event) ~ gene_group + chemo, data = df_temp)
  adj_summary <- summary(adj_fit)
  adj_hr <- adj_summary$conf.int[1, 1]
  adj_lower <- adj_summary$conf.int[1, 3]
  adj_upper <- adj_summary$conf.int[1, 4]
  adj_p <- adj_summary$coefficients[1, 5]
  
  dlbcl_survival_results[[g]] <- data.frame(
    gene_symbol = g,
    logrank_p = logrank_p,
    unadj_HR = uni_hr,
    unadj_HR_lower = uni_lower,
    unadj_HR_upper = uni_upper,
    unadj_cox_p = uni_p,
    adj_HR = adj_hr,
    adj_HR_lower = adj_lower,
    adj_HR_upper = adj_upper,
    adj_cox_p = adj_p,
    stringsAsFactors = FALSE
  )
}

dlbcl_res_df <- do.call(rbind, dlbcl_survival_results)
write.csv(dlbcl_res_df, file = "data_clean/gse10846_dlbcl_survival_results.csv", row.names = FALSE)
write.csv(dlbcl_res_df, file = "Week6/gse10846_dlbcl_survival_results.csv", row.names = FALSE)

print(dlbcl_res_df)

# ==============================================================================
# SECTION 4: Draw Kaplan-Meier Curves for Validated Targets | 第 4 部分：绘制 K-M 生存曲线
# ==============================================================================

# Custom ggplot2 KM plotting function to avoid survminer dependencies
draw_km_ggplot <- function(gene, data, clean_expr_data) {
  val <- as.numeric(clean_expr_data[gene, ])
  med <- median(val, na.rm = TRUE)
  grp <- ifelse(val >= med, "High", "Low")
  
  df_temp <- data
  df_temp$gene_group <- factor(grp, levels = c("Low", "High"))
  
  # Fit survival curve
  fit <- survfit(Surv(years, event) ~ gene_group, data = df_temp)
  
  # Compile plot dataframe
  fit_df <- data.frame(
    time = fit$time,
    surv = fit$surv,
    upper = fit$upper,
    lower = fit$lower,
    n.censor = fit$n.censor,
    group = rep(names(fit$strata), fit$strata)
  )
  fit_df$group <- gsub("gene_group=", "", fit_df$group)
  
  # Extract p-value
  p_val <- dlbcl_res_df[dlbcl_res_df$gene_symbol == gene, "logrank_p"]
  p_text <- paste0("Log-rank P = ", format.pval(p_val, digits = 3))
  
  # Plot
  p <- ggplot(fit_df, aes(x = time, y = surv, color = group)) +
    geom_step(size = 1.2) +
    labs(
      title = paste0("GSE10846 (DLBCL): ", gene),
      x = "Survival Time (Years)",
      y = "Overall Survival Probability",
      color = "Expression Group"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_color_manual(values = c("High" = "#D95F02", "Low" = "#1B9E77")) +
    annotate("text", x = 1.5, y = 0.15, label = p_text, fontface = "italic", size = 4.5, color = "black") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
    )
  return(p)
}

cat("Generating Kaplan-Meier plots for validated targets...\n")
p_vnn2 <- draw_km_ggplot("VNN2", clean_cli, clean_expr)
p_fpr2 <- draw_km_ggplot("FPR2", clean_cli, clean_expr)
p_cmtm2 <- draw_km_ggplot("CMTM2", clean_cli, clean_expr)

# Arrange and save as a panel
g_km <- grid.arrange(p_vnn2, p_fpr2, p_cmtm2, ncol = 3)
ggsave(filename = "Week6/km_curves_dlbcl_validated.png", plot = g_km, width = 12, height = 4.2, dpi = 300)
ggsave(filename = "results/km_curves_dlbcl_validated.png", plot = g_km, width = 12, height = 4.2, dpi = 300)

# ==============================================================================
# SECTION 5: Draw Cross-Dataset HR Forest Plot | 第 5 部分：绘制跨数据集 HR 森林对比图
# ==============================================================================

# Load AML survival results
aml_res <- read.csv("data_clean/survival_cox_results.csv") %>%
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, HR, HR_CI_lower, HR_CI_upper, p_value) %>%
  mutate(Cohort = "AML (TCGA LAML)")

# Clean columns and bind
dlbcl_forest_data <- dlbcl_res_df %>%
  select(gene_symbol, adj_HR, adj_HR_lower, adj_HR_upper, adj_cox_p) %>%
  rename(HR = adj_HR, HR_CI_lower = adj_HR_lower, HR_CI_upper = adj_HR_upper, p_value = adj_cox_p) %>%
  mutate(Cohort = "DLBCL (GSE10846, Adjusted)")

forest_all <- rbind(aml_res, dlbcl_forest_data)

# Reorder targets logically (e.g. consistently significant, AML-specific, lineage-opposing)
forest_all$gene_symbol <- factor(forest_all$gene_symbol, levels = rev(c("CMTM2", "FPR2", "GZMB", "IL1R2", "SLC15A3", "CXCR1", "CXCR2", "CD14", "EGFL7", "VNN2")))

p_forest <- ggplot(forest_all, aes(x = HR, y = gene_symbol, color = Cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = HR_CI_lower, xmax = HR_CI_upper), height = 0.3, size = 0.8, position = position_dodge(0.5)) +
  geom_point(size = 3.5, position = position_dodge(0.5)) +
  labs(
    title = "Cross-Dataset Prognostic Hazard Comparison (AML vs. DLBCL)",
    x = "Hazard Ratio (HR) with 95% CI",
    y = "Target Candidates",
    color = "Malignancy Cohort"
  ) +
  scale_color_manual(values = c("AML (TCGA LAML)" = "#E41A1C", "DLBCL (GSE10846, Adjusted)" = "#377EB8")) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

ggsave(filename = "Week6/forest_plot_cross_dataset.png", plot = p_forest, width = 8.5, height = 5.2, dpi = 300)
ggsave(filename = "results/forest_plot_cross_dataset.png", plot = p_forest, width = 8.5, height = 5.2, dpi = 300)

cat("All validations and visualizations completed successfully!\n")
