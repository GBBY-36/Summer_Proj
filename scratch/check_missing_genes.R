library(GEOquery)
library(Biobase)
library(dplyr)
library(stringr)

# Load robust biomarkers
degs_overlap <- read.csv("data_clean/sig_degs_overlap.csv", stringsAsFactors = FALSE)
cat("Robust biomarkers count:", nrow(degs_overlap), "\n")

# Load GSE24759
gse <- getGEO(filename = "data_raw/gse24759/GSE24759_series_matrix.txt.gz", getGPL = FALSE)
gse_probes <- rownames(gse)

# Load GPL570 annotation
gpl_file <- "data_raw/gse13159/GPL570.annot.gz"
gpl_data <- parseGEO(gpl_file)
gpl_table <- Table(gpl_data)

# Map all probes to gene symbols in GPL570
gpl_mapping <- gpl_table %>%
  select(ID, `Gene symbol`) %>%
  rename(probe_id = ID, gene_symbol_raw = `Gene symbol`) %>%
  mutate(
    gene_symbol = str_split(gene_symbol_raw, "///", simplify = TRUE)[, 1],
    gene_symbol = str_trim(gene_symbol)
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", gene_symbol != "---", gene_symbol != "NA")

cat("Total mapping rows in GPL570:", nrow(gpl_mapping), "\n")

# Probes in GSE24759
gse_mapping <- gpl_mapping %>%
  filter(probe_id %in% gse_probes)

cat("Probes in GSE24759 with valid gene symbol:", nrow(gse_mapping), "\n")
cat("Unique gene symbols in GSE24759 mapping:", length(unique(gse_mapping$gene_symbol)), "\n")

# Check which of our 131 robust biomarkers are not mapped
missing_genes <- setdiff(degs_overlap$gene_symbol, gse_mapping$gene_symbol)
cat("Missing genes (", length(missing_genes), "):\n", sep="")
print(missing_genes)

# For up-regulated genes
missing_up_genes <- setdiff(degs_overlap %>% filter(direction == "Up") %>% pull(gene_symbol), gse_mapping$gene_symbol)
cat("\nMissing UP-regulated genes (", length(missing_up_genes), "):\n", sep="")
print(missing_up_genes)

# Are these missing genes in the GPL570 annotation table at all?
gpl_all_symbols <- unique(gpl_mapping$gene_symbol)
missing_in_gpl <- setdiff(degs_overlap$gene_symbol, gpl_all_symbols)
cat("\nBiomarkers missing in GPL570 entirely (", length(missing_in_gpl), "):\n", sep="")
print(missing_in_gpl)
