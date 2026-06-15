library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)
library(data.table)

# 1. Load the 131 robust biomarkers
degs_overlap <- read.csv("data_clean/sig_degs_overlap.csv", stringsAsFactors = FALSE)
cat("Loaded", nrow(degs_overlap), "robust biomarkers.\n")

# 2. Load GSE24759
cat("Loading GSE24759...\n")
gse <- getGEO(filename = "data_raw/gse24759/GSE24759_series_matrix.txt.gz", getGPL = FALSE)
expr_matrix <- exprs(gse)
pheno <- pData(gse)

# 3. Identify HSPC samples
hspc_idx <- which(grepl("Hematopoietic stem cell|Common myeloid progenitor|Granulocyte/monocyte progenitor|Megakaryocyte/ erythroid progenitor", 
                        pheno$title, ignore.case = TRUE))
hspc_samples <- pheno$geo_accession[hspc_idx]
hspc_titles <- pheno$title[hspc_idx]

cat("Found", length(hspc_samples), "HSPC samples:\n")
for (i in 1:length(hspc_samples)) {
  cat("  -", hspc_samples[i], ":", hspc_titles[i], "\n")
}

# 4. Map probes in GSE24759 to Gene Symbols using local GPL570 annot
cat("Mapping probes to Gene Symbols...\n")
gpl_file <- "data_raw/gse13159/GPL570.annot.gz"
gpl_data <- parseGEO(gpl_file)
gpl_table <- Table(gpl_data)

probe_gene <- data.frame(
  probe_id = as.character(gpl_table$ID),
  gene_symbol_raw = as.character(gpl_table$`Gene symbol`),
  stringsAsFactors = FALSE
) %>%
  mutate(
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1],
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---", gene_symbol != "NA")

# Create data.table for expression matrix
expr_dt <- as.data.table(expr_matrix, keep.rownames = "probe_id")
probe_gene_dt <- as.data.table(probe_gene)

# Merge
expr_gene <- merge(expr_dt, probe_gene_dt[, .(probe_id, gene_symbol)], by = "probe_id")

# Aggregate by gene symbol (mean)
expr_gene_mean <- expr_gene[, lapply(.SD, mean), by = gene_symbol, .SDcols = colnames(expr_matrix)]

cat("Aggregated matrix dimensions:", nrow(expr_gene_mean), "genes\n")

# Convert to data.frame
expr_gene_df <- as.data.frame(expr_gene_mean)
rownames(expr_gene_df) <- expr_gene_df$gene_symbol
expr_gene_df$gene_symbol <- NULL

# 5. Extract expression values for our 131 overlapping biomarkers in HSPC samples
overlap_genes_in_gse <- intersect(degs_overlap$gene_symbol, rownames(expr_gene_df))
cat("Found", length(overlap_genes_in_gse), "out of", nrow(degs_overlap), "overlapping biomarkers in GSE24759.\n")

hspc_expr <- expr_gene_df[overlap_genes_in_gse, hspc_samples]

# Calculate mean expression for each gene in HSPCs
gene_hspc_mean <- rowMeans(hspc_expr, na.rm = TRUE)

# Combine with gene metadata
result_df <- degs_overlap %>%
  filter(gene_symbol %in% overlap_genes_in_gse) %>%
  mutate(hspc_mean_expr = gene_hspc_mean[gene_symbol]) %>%
  arrange(desc(direction), desc(logFC_discovery))

cat("\nExpression of UP-regulated biomarkers in normal HSPCs:\n")
print(result_df %>% filter(direction == "Up") %>% select(gene_symbol, logFC_discovery, logFC_validation, hspc_mean_expr))

cat("\nSummary statistics of HSPC mean expression for UP-regulated genes:\n")
print(summary(result_df %>% filter(direction == "Up") %>% pull(hspc_mean_expr)))

cat("\nSummary statistics of HSPC mean expression for DOWN-regulated genes:\n")
print(summary(result_df %>% filter(direction == "Down") %>% pull(hspc_mean_expr)))
