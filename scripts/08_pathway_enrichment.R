# scripts/08_pathway_enrichment.R

# ==============================================================================
# SECTION 0: Setup and Package Loading | 第 0 部分：准备与加载包
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Install clusterProfiler if not present
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  cat("clusterProfiler not found. Installing from Bioconductor...\n")
  options(timeout = 100000)
  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("Week5", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 1: Load Input Genes | 第 1 部分：读取输入基因列表
# ==============================================================================

cat("Loading gene lists from Week 3 and Week 4...\n")
degs_filtered <- read.csv("data_clean/sig_degs_hspc_filtered.csv", stringsAsFactors = FALSE)
cox_filtered <- read.csv("data_clean/survival_cox_sig_danger_genes.csv", stringsAsFactors = FALSE)

up_genes <- degs_filtered %>% filter(direction == "Up") %>% pull(gene_symbol)
down_genes <- degs_filtered %>% filter(direction == "Down") %>% pull(gene_symbol)
prog_genes <- cox_filtered$gene_symbol

cat("Loaded:\n")
cat("  - Safety-filtered Up-regulated genes:", length(up_genes), "\n")
cat("  - Safety-filtered Down-regulated genes:", length(down_genes), "\n")
cat("  - Prognostic hazard genes:", length(prog_genes), "\n")

# ==============================================================================
# SECTION 2: Entrez ID Mapping | 第 2 部分：基因名转换为 Entrez ID
# ==============================================================================

cat("\nMapping Gene Symbols to Entrez IDs using org.Hs.eg.db...\n")
map_ids <- function(symbols) {
  if (length(symbols) == 0) return(NULL)
  tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    NULL
  })
}

ids_up <- map_ids(up_genes)
ids_down <- map_ids(down_genes)
ids_prog <- map_ids(prog_genes)

# ==============================================================================
# SECTION 3: Custom Enrichment Plotting Function | 第 3 部分：自定义富集分析绘图函数
# ==============================================================================

plot_enrichment_custom <- function(enrich_res, title, out_path_res, out_path_week5) {
  if (is.null(enrich_res)) return(FALSE)
  df <- as.data.frame(enrich_res)
  if (nrow(df) == 0) {
    cat("No significant terms enriched for:", title, "\n")
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No significant pathways enriched\n(", title, ")", sep=""), size = 4.5, fontface = "italic") + 
      theme_void()
    ggsave(out_path_res, plot = p, width = 7, height = 5, dpi = 300)
    file.copy(out_path_res, out_path_week5, overwrite = TRUE)
    return(FALSE)
  }
  
  df <- df %>% 
    arrange(p.adjust) %>% 
    head(15)
  
  df$GeneRatioVal <- sapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    parts[1] / parts[2]
  })
  
  df$Description <- sapply(df$Description, function(x) {
    if (nchar(x) > 55) {
      paste0(substr(x, 1, 52), "...")
    } else {
      x
    }
  })
  
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  p <- ggplot(df, aes(x = GeneRatioVal, y = Description, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "#C00000", high = "#2E75B6", name = "adj.P.Val") +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    theme_bw(base_size = 11) +
    labs(
      title = title,
      x = "Gene Ratio (Enriched / Total)",
      y = "Pathway / Term Description"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(color = "black", size = 9),
      axis.text.x = element_text(color = "black"),
      legend.title = element_text(face = "bold", size = 9)
    )
  
  ggsave(out_path_res, plot = p, width = 8.5, height = 6, dpi = 300)
  file.copy(out_path_res, out_path_week5, overwrite = TRUE)
  cat("Saved enrichment plot to:", out_path_res, "\n")
  return(TRUE)
}

# ==============================================================================
# SECTION 4: Run GO and KEGG Pathway Enrichment | 第 4 部分：执行 GO 和 KEGG 富集分析
# ==============================================================================

run_enrichment_analysis <- function(entrez_ids, prefix_title, file_suffix) {
  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    cat("Empty Entrez ID list for:", prefix_title, "\n")
    return(NULL)
  }
  
  cat("\n--- Running GO & KEGG for:", prefix_title, "---\n")
  
  # 1. GO Biological Process
  cat("Running GO Biological Process enrichment...\n")
  ego <- tryCatch({
    enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    NULL
  })
  
  plot_enrichment_custom(
    ego, 
    title = paste0(prefix_title, " GO Biological Process Enrichment"), 
    out_path_res = paste0("results/go_enrichment_bubble_", file_suffix, ".png"),
    out_path_week5 = paste0("Week5/go_enrichment_bubble_", file_suffix, ".png")
  )
  
  if (!is.null(ego)) {
    write.csv(as.data.frame(ego), file = paste0("data_clean/go_enrichment_results_", file_suffix, ".csv"), row.names = FALSE)
  }
  
  # 2. KEGG Pathways
  cat("Running KEGG Pathway enrichment (using tryCatch for network robustness)...\n")
  ekegg <- tryCatch({
    enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    cat("Warning: KEGG query failed due to network error. Skipping.\n")
    NULL
  })
  
  plot_enrichment_custom(
    ekegg, 
    title = paste0(prefix_title, " KEGG Pathway Enrichment"), 
    out_path_res = paste0("results/kegg_enrichment_bubble_", file_suffix, ".png"),
    out_path_week5 = paste0("Week5/kegg_enrichment_bubble_", file_suffix, ".png")
  )
  
  if (!is.null(ekegg)) {
    write.csv(as.data.frame(ekegg), file = paste0("data_clean/kegg_enrichment_results_", file_suffix, ".csv"), row.names = FALSE)
  }
}

# Run for up-regulated safety-filtered genes
run_enrichment_analysis(ids_up, "Up-regulated Safety Genes", "up")

# Run for prognostic hazard genes
run_enrichment_analysis(ids_prog, "Prognostic Hazard Genes", "prognostic")

# Run for down-regulated safety-filtered genes
run_enrichment_analysis(ids_down, "Down-regulated Safety Genes", "down")


# ==============================================================================
# SECTION 5: Druggability Assessment and Novel Target Screening | 第 5 部分：可靶性评估与新靶点筛选
# ==============================================================================

cat("\n======================================================================\n")
cat("Performing Druggability Comparison and Screening...\n")

# Combine all unique up-regulated candidates (from 19 safe up-regulated + 13 prognostic hazard)
all_candidates <- unique(c(up_genes, prog_genes))
cat("Total unique candidate genes for targeting (up-regulated in AML):", length(all_candidates), "\n")

# Define the Known Target List
known_targets <- c(
  "CD33", "IL3RA", "FLT3", "BCL2", "IDH1", "IDH2", "KIT", "CD47", 
  "HAVCR2", "WT1", "KMT2A", "TP53", "CD38", "NPM1"
)

# Manual subcellular localization and functional annotations for candidates
loc_map <- list(
  "IL1R2" = "Plasma membrane / Secreted",
  "SLC15A3" = "Lysosomal membrane / Plasma membrane",
  "NKX2-3" = "Nucleus",
  "NCF1C" = "Cytoplasm",
  "NCF1" = "Cytoplasm",
  "VNN2" = "Plasma membrane (GPI-anchored)",
  "CXCR1" = "Plasma membrane (GPCR)",
  "IL3RA" = "Plasma membrane",
  "CXCR2" = "Plasma membrane (GPCR)",
  "CD14" = "Plasma membrane / Secreted",
  "CMTM2" = "Plasma membrane / Transmembrane",
  "GZMB" = "Secreted / Extracellular space",
  "FPR2" = "Plasma membrane (GPCR)",
  "TPSAB1" = "Secreted (Extracellular)",
  "TPSB2" = "Secreted (Extracellular)",
  "CLIP2" = "Cytoplasm",
  "EGFL7" = "Secreted / Extracellular matrix",
  "CCNA1" = "Nucleus / Cytoplasm",
  "NRXN2" = "Plasma membrane",
  "PDE6G" = "Cytoplasm / Membrane",
  "CPXM1" = "Secreted (Extracellular)",
  "C1QTNF4" = "Secreted (Extracellular)",
  "COL24A1" = "Secreted / Extracellular matrix",
  "MEX3B" = "Cytoplasm / Nucleus",
  "NPW" = "Secreted (Extracellular)",
  "MAMDC2" = "Secreted / Extracellular matrix",
  "TRIM71" = "Cytoplasm",
  "UMODL1" = "Plasma membrane / Secreted",
  "IRX3" = "Nucleus"
)

# Compile comprehensive information for each candidate
target_info_list <- list()

for (gene in all_candidates) {
  # Subcellular localization lookup
  loc <- ifelse(!is.null(loc_map[[gene]]), loc_map[[gene]], "Unknown")
  
  # Novelty status
  status <- ifelse(gene %in% known_targets, "Known Target", "Novel Candidate")
  
  # Check if present in Week 3 safety-filtered DEG file
  match_deg <- degs_filtered %>% filter(gene_symbol == gene)
  logFC_disc <- ifelse(nrow(match_deg) > 0, match_deg$logFC_discovery, NA)
  logFC_valid <- ifelse(nrow(match_deg) > 0, match_deg$logFC_validation, NA)
  
  # Check if present in Week 4 prognostic file
  match_surv <- cox_filtered %>% filter(gene_symbol == gene)
  hr_val <- ifelse(nrow(match_surv) > 0, match_surv$HR, NA)
  p_val <- ifelse(nrow(match_surv) > 0, match_surv$p_value, NA)
  
  # Determine therapeutic modality
  modality <- "Small-molecule inhibitor"
  if (grepl("Plasma membrane|Secreted", loc)) {
    modality <- "Monoclonal antibody / ADC / CAR-T"
  }
  
  # Calculate priority score for screening:
  # Score starts with 0. 
  # Membrane/Secreted gets +2 (easier to target via immune therapy).
  # High risk in survival (HR > 1.5) gets +3.
  # Higher logFC (logFC > 2.0) gets +1.
  # Novel status gets +2 (academic novelty).
  score <- 0
  if (grepl("Plasma membrane|Secreted", loc)) score <- score + 2
  if (!is.na(hr_val) && hr_val > 1.5) score <- score + 3
  if (!is.na(logFC_disc) && logFC_disc > 2.0) score <- score + 1
  if (status == "Novel Candidate") score <- score + 2
  
  target_info_list[[gene]] <- data.frame(
    gene_symbol = gene,
    subcellular_localization = loc,
    novelty_status = status,
    modality_category = modality,
    logFC_discovery = logFC_disc,
    logFC_validation = logFC_valid,
    survival_HR = hr_val,
    survival_p_value = p_val,
    priority_score = score,
    stringsAsFactors = FALSE
  )
}

# Bind into dataframe
all_targets_df <- do.call(rbind, target_info_list) %>%
  arrange(desc(priority_score), desc(logFC_discovery))

# Filter for Novel Candidates
potential_novel_targets <- all_targets_df %>%
  filter(novelty_status == "Novel Candidate")

cat("Total novel candidate genes identified:", nrow(potential_novel_targets), "\n")
cat("Displaying prioritized novel target candidates:\n")
print(head(potential_novel_targets, 15))

# Save the full target table and the filtered novel target table
write.csv(all_targets_df, "data_clean/target_druggability_analysis.csv", row.names = FALSE)
write.csv(potential_novel_targets, "data_clean/potential_novel_targets_list.csv", row.names = FALSE)
write.csv(potential_novel_targets, "Week5/potential_novel_targets_list.csv", row.names = FALSE)

cat("\nOutputs successfully saved to data_clean/ and Week5/:\n")
cat("  - Week5/potential_novel_targets_list.csv\n")
cat("  - go/kegg plots in Week5/\n")

cat("\nPipeline script 08 completed successfully!\n")
