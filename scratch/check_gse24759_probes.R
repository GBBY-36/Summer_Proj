library(GEOquery)
library(Biobase)

cat("Loading GSE24759...\n")
gse <- getGEO(filename = "data_raw/gse24759/GSE24759_series_matrix.txt.gz", getGPL = FALSE)

cat("\nPlatform annotation from ExpressionSet:\n")
print(annotation(gse))

cat("\nUnique platform IDs from phenotype data:\n")
print(table(pData(gse)$platform_id))

cat("\nFirst 10 probe IDs:\n")
print(head(rownames(gse), 10))

# Load our local GPL570 annotation
cat("\nLoading local GPL570 annotation...\n")
gpl_file <- "data_raw/gse13159/GPL570.annot.gz"
if (file.exists(gpl_file)) {
  gpl_data <- parseGEO(gpl_file)
  gpl_table <- Table(gpl_data)
  
  common_probes <- intersect(rownames(gse), gpl_table$ID)
  cat("Number of probes in GSE24759 that match GPL570 annotation:", length(common_probes), "\n")
  cat("Total probes in GSE24759:", nrow(gse), "\n")
} else {
  cat("GPL570 annotation file not found at", gpl_file, "\n")
}
