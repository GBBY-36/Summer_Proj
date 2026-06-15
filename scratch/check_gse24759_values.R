library(GEOquery)
library(Biobase)

cat("Loading GSE24759 local file with getGPL = FALSE...\n")
gse <- getGEO(filename = "data_raw/gse24759/GSE24759_series_matrix.txt.gz", getGPL = FALSE)
expr_matrix <- exprs(gse)

cat("Dimensions of expression matrix:\n")
print(dim(expr_matrix))

cat("Range of expression values:\n")
print(range(expr_matrix, na.rm = TRUE))

cat("First few rows and columns:\n")
print(expr_matrix[1:5, 1:5])

pheno <- pData(gse)
cat("Pheno columns:\n")
print(colnames(pheno))
