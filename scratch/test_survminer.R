if (!requireNamespace("survminer", quietly = TRUE)) {
  cat("Installing survminer...\n")
  install.packages("survminer", repos = "https://cloud.r-project.org")
}

library(survminer)
library(survival)

cat("survminer successfully loaded!\n")
