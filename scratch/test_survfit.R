library(survival)
library(dplyr)

cat("Loading data...\n")
group_df <- readRDS("data_clean/tcga_laml_gene_groups.rds")
tcga_sample_info <- read.csv("data_clean/tcga_laml_sample_info.csv", stringsAsFactors = FALSE)
clinical <- readRDS("data_raw/tcga_laml/tcga_laml_clinical.rds")

clinical_survival <- clinical %>%
  mutate(
    patient_id = submitter_id,
    time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    event = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time) & time > 0) %>%
  select(patient_id, time, event)

merged_data <- group_df %>%
  inner_join(tcga_sample_info[, c("sample_id", "patient_id")], by = "sample_id") %>%
  inner_join(clinical_survival, by = "patient_id")

cat("Testing survfit for NKX2-3...\n")
fit_diff <- survdiff(Surv(time, event) ~ `NKX2-3`, data = merged_data)
print(fit_diff)

fit_gene <- survfit(Surv(time, event) ~ `NKX2-3`, data = merged_data)
print(fit_gene)
