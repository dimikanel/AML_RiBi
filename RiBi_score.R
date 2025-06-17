# Load relevant packages 

# Load required libraries
library(tidyverse)
library(DESeq2)
library(limma)
library(umap)
library(glmnet)
library(survival)
library(finalfit)
library(survminer)
library(plotROC)
library(TCGAbiolinks)
library(SummarizedExperiment)



# Load Cohorts to be used in training the regression model


cs_counts <- read.csv('cs_counts.csv')
cs_meta <- read.csv('cs_meta.csv')

beat_counts <- read.csv('beat_counts.csv')
beat_meta <- read.csv('beat_meta.csv')

finn_counts <- read.csv('finn_counts.csv')
finn_meta <- read.csv('finn_meta.csv')


cs_meta <- cs_meta %>% filter(!Disease_Status == "Healthy")
cs_counts <- cs_counts %>% dplyr::rename(ensembl_ID = X) %>%
  distinct(ensembl_ID, .keep_all = TRUE) %>% select(ensembl_ID,cs_meta$patient)


beat_meta <- beat_meta %>% filter(!Disease_Status == "Healthy")
beat_counts <- beat_counts  %>% dplyr::rename(ensembl_ID = X) %>%
  distinct(ensembl_ID, .keep_all = TRUE) %>% select(ensembl_ID,beat_meta$patient)


finn_meta <- finn_meta %>% filter(!Disease_Status == "Healthy")
finn_counts <- finn_counts %>% dplyr::rename(ensembl_ID = X) %>%
  distinct(ensembl_ID, .keep_all = TRUE) %>% select(ensembl_ID,finn_meta$patient)


# Merge All Cohorts 
# Merge counts

all_counts <- list(cs_counts, beat_counts, finn_counts) %>%
  purrr::reduce(~ merge(.x, .y, by = "ensembl_ID")) %>%
  distinct(ensembl_ID, .keep_all = TRUE) %>%
  column_to_rownames("ensembl_ID")


# Clean duplicated columns

all_counts <- all_counts[, !duplicated(names(all_counts))]

# Merge all metadata

meta_all <- bind_rows(cs_meta, beat_meta, finn_meta) %>%
  distinct(patient, .keep_all = TRUE) %>%
  filter(patient %in% colnames(all_counts))


# Align columns and rows

all_counts <- all_counts[, meta_all$patient]
stopifnot(all.equal(colnames(all_counts), meta_all$patient))

table(meta_all$Disease_Status)
table(meta_all$cohort)

# DESeq2 Normalization 


dds <- DESeqDataSetFromMatrix(
  countData = all_counts,
  colData = meta_all %>% column_to_rownames("patient"),
  design = ~ 1
)

# Filter low-count genes

dds <- dds[rowSums(counts(dds)) >= 10, ]

# Apply variance-stabilizing transformation

vst_counts <- vst(dds)

# Remove batch effects

assay(vst_counts) <- limma::removeBatchEffect(assay(vst_counts), vst_counts$cohort)

# Extract normalized matrix

norm_counts <- assay(vst_counts)

write.csv(norm_counts, 'norm_counts.csv')

# Model Training 
# Load gene sets

ribi_terms <- read.csv("RiBi.only.terms.GSEA.AmiGO.csv")
de_names <- read.csv("DE_ClinSeq_BEAT_FINN_AML_Vs_Healthy_filtered_05.csv")

# Subset RiBi genes and DE genes

ribi_genes <- norm_counts %>% as.data.frame() %>%
  filter(rownames(.) %in% ribi_terms$ensembl_gene_id) %>%
  filter(rownames(.) %in% de_names$X)


# Filter metadata for survival-ready patients

surv_meta <- meta_all %>% filter(patient %in% colnames(ribi_genes), !is.na(time)) %>%
  select(patient, time, status) %>%
  column_to_rownames("patient")

# Transpose gene expression matrix to match glmnet input

x <- t(ribi_genes[, rownames(surv_meta)])
y <- surv_meta

# Fit Cox model

fit <- glmnet(as.matrix(x), as.matrix(y), family = "cox")

# Plot coefficient paths across lambda values

plot(fit, xvar = "lambda", label = TRUE)
title("Coefficient Paths")

set.seed(1)
cvfit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox")

# Partial likelihood deviance vs. log(lambda)

plot(cvfit)
title("Partial Likelihood Deviance (CV)", line = 2.5)


# Extract non-zero coefficients

coef_df <- coef(cvfit, s = cvfit$lambda.min)
features <- rownames(coef_df)[which(coef_df[, 1] != 0)]
weights <- coef_df[features, 1]


coef_table <- data.frame(gene = features, coefficient = as.numeric(weights))

# Save to CSV

write.csv(coef_table, "RiBi_model_coefficients.csv", row.names = FALSE)

# Compute RiBi Score 

# Subset and align matrix to selected features

score_matrix <- x[, features, drop = FALSE]
weighted_matrix <- t(score_matrix) * weights

# Score per patient

score_df <- data.frame(patient = rownames(score_matrix), RiBi_score = colSums(weighted_matrix, na.rm = TRUE))

# Classify score into quantile-based risk groups

score_df <- score_df %>%
  mutate(RiBi_status = case_when(
    RiBi_score > quantile(RiBi_score, 0.75) ~ "high",
    RiBi_score < quantile(RiBi_score, 0.25) ~ "low",
    TRUE ~ "moderate"
  ))


# Merge Score with Metadata

meta_with_score <- meta_all %>%
  left_join(score_df, by = "patient") %>%
  filter(!is.na(RiBi_score), time > 0)


# ROC Curve 

# Create ROC object

roc_obj <- ggplot(meta_with_score, aes(d = status, m = RiBi_score)) +
  geom_roc(color = "red", n.cuts = 0) +
  geom_abline(intercept = 0)


# Calculate AUC from the data

auc_val <- calc_auc(roc_obj)$AUC

# Add annotation and labels

roc_plot <- roc_obj + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_val, 3))) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  ggtitle("ROC - Training Cohort")


print(roc_plot)

# Survival Analysis 


dependent_os <- "Surv(time, status)"
explanatory <- "RiBi_status"

meta_with_score %>% surv_plot(dependent_os, explanatory, pval = TRUE)
subset(meta_with_score, time < 364*5) %>% surv_plot(dependent_os, explanatory, pval = TRUE)  # 5y survival

# Statistics (10y survival)

survdiff(Surv(time, status) ~ RiBi_status, data = subset(meta_with_score, time < 364*5))
pairwise_survdiff(Surv(time, status) ~ RiBi_status, data = subset(meta_with_score, time < 364*5))
coxph(Surv(time, status) ~ RiBi_status, data = subset(meta_with_score, time < 364*5))


# Multivariate COX proportional hazards ratio analysis

str(meta_with_score)
meta_with_score$age <- as.numeric(meta_with_score$age)

meta_with_score <- meta_with_score %>%
  filter(time < 364*5)

median_age <- median(meta_with_score$age, na.rm = TRUE)

meta_with_score <- meta_with_score %>%
  mutate(age_group = ifelse(age >= median_age, "High", "Low")) %>%
  mutate(age_group = factor(age_group, levels = c("Low", "High")))

write.csv(meta_with_score, 'meta_training.csv')


meta_with_score <- meta_with_score %>%
  mutate(across(c(RiBi_status, sex, NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53), as.factor))

# Fit Cox proportional hazards model

cox_multi <- coxph(Surv(time, status) ~ RiBi_status + sex + NPM1 + DNMT3A + IDH1 + IDH2 + FLT3 + TP53 + age_group, data = meta_with_score)

summary(cox_multi)
ggforest(cox_multi, data = meta_with_score, main = "Multivariable Cox Proportional Hazards Model")



# Load and Process TCGA Data 


# Expression data via TCGAbiolinks

query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga_raw <- GDCprepare(query)

# Extract and clean gene counts

tcga_counts <- assay(tcga_raw) %>% as.data.frame() %>%
  rownames_to_column("ensgene") %>% mutate(ensgene = gsub("\\..*", "", ensgene)) %>%
  distinct(ensgene, .keep_all = TRUE) %>% column_to_rownames("ensgene")

# Clean sample names

colnames(tcga_counts) <- sub("^([A-Z0-9-]+?)-\\d{2}[A-Z]-\\d{2}[A-Z]-\\d{4}-\\d{2}$", "\\1", colnames(tcga_counts))


# Load and Clean TCGA Metadata 

# Clinical info

tcga_meta <- read.delim("data_clinical_sample.txt", check.names = FALSE)
colnames(tcga_meta) <- tcga_meta[4, ]
tcga_meta <- tcga_meta[-(1:4), ]
colnames(tcga_meta)[2] <- "patient"

# Mutation file

tcga_mut <- read.delim("data_mutations.txt") %>%
  filter(Hugo_Symbol %in% c("NPM1", "DNMT3A", "IDH1", "IDH2", "FLT3", "TP53"))

# Build mutation status by gene

extract_status <- function(gene_symbol, colname) {
  patients <- tcga_mut %>%
    filter(Hugo_Symbol == gene_symbol) %>%
    pull(Matched_Norm_Sample_Barcode)
  mutate_column <- function(df) {
    df %>%
      mutate(!!colname := if_else(patient %in% patients, "mut", "wt"))
  }
  mutate_column
}

tcga_meta <- tcga_meta %>% dplyr::select(patient)
tcga_meta <- tcga_meta %>%
  extract_status("NPM1", "NPM1")() %>%
  extract_status("DNMT3A", "DNMT3A")() %>%
  extract_status("IDH1", "IDH1")() %>%
  extract_status("IDH2", "IDH2")() %>%
  extract_status("FLT3", "FLT3")() %>%
  extract_status("TP53", "TP53")()

# Load survival data


tcga_raw <- read.delim("data_clinical_patient.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(tcga_raw) 

colnames(tcga_raw) <- tcga_raw[5, ]
tcga_raw <- tcga_raw[-(1:5), ]
colnames(tcga_raw)
str(tcga_raw)

colnames(tcga_raw) 

# Drop APL and na samples

tcga_surv <- tcga_raw %>% filter(!CYTOGENETIC_CODE_OTHER %in% c('PML-RARA', 'N.D.')) %>%
  dplyr::select(patient = PATIENT_ID,
                OS_time = OS_MONTHS,
                OS_status = OS_STATUS,
                gender = SEX,
                age = AGE) %>%
  mutate(OS_status = as.numeric(str_detect(OS_status, "1:DECEASED")),
        sex = recode(gender, "Female" = "female", "Male" = "male"),
         age = as.numeric(age)) %>% dplyr::rename(time = OS_time) %>%
  dplyr::select(patient, time, status = OS_status, sex, age)



# Merge survival, mutation, cohort

tcga_meta <- inner_join(tcga_meta, tcga_surv, by = "patient") %>%
  mutate(cohort = "TCGA", Disease_Status = "Disease", stage = 'Initial Diagnosis') %>%
  dplyr::select(patient, time, status, cohort, Disease_Status, sex, age,
                NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53, stage)


tcga_meta <- tcga_meta %>% filter(patient %in% colnames(tcga_counts))

# Filter counts and match order

tcga_counts <- tcga_counts[, tcga_meta$patient]


write.csv(tcga_counts, 'tcga_counts.csv')
write.csv(tcga_meta, 'tcga_meta.csv')


#tcga_counts <- read.csv('tcga_counts.csv', row.names = 1, check.names = FALSE)
#tcga_meta <- read.csv('tcga_meta.csv')

# Normalize TCGA 

dds_tcga <- DESeqDataSetFromMatrix(
  countData = round(tcga_counts),
  colData = tcga_meta %>% column_to_rownames("patient"),
  design = ~ 1)

dds_tcga <- dds_tcga[rowSums(counts(dds_tcga)) >= 10, ]
vst_tcga <- vst(dds_tcga)
norm_tcga <- assay(vst_tcga)


tcga_score_matrix <- t(norm_tcga[features, , drop = FALSE])
tcga_weighted_matrix <- t(tcga_score_matrix) * weights

# Score per patient

tcga_score_df <- data.frame(patient = rownames(tcga_score_matrix), RiBi_score = colSums(tcga_weighted_matrix, na.rm = TRUE))

# Classify score into quantile-based risk groups

tcga_score_df <- tcga_score_df %>%
  mutate(RiBi_status = case_when(
    RiBi_score > quantile(RiBi_score, 0.75) ~ "high",
    RiBi_score < quantile(RiBi_score, 0.25) ~ "low",
    TRUE ~ "moderate"
  ))


# Merge Score with Metadata 

tcga_eval <- tcga_meta %>%
  left_join(tcga_score_df, by = "patient") %>%
  filter(time > 0)



# Build initial ROC ggplot object

roc_obj <- ggplot(tcga_eval, aes(d = status, m = RiBi_score)) +
  geom_roc(color = "darkblue", n.cuts = 0) +
  geom_abline(intercept = 0)

# Calculate AUC

auc_val <- calc_auc(roc_obj)$AUC

# Add annotation and labels

roc_plot <- roc_obj + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_val, 3))) +
  labs(title = "TCGA ROC", x = "1 - Specificity", y = "Sensitivity")

print(roc_plot)


# Survival


tcga_eval$time <- as.numeric(tcga_eval$time)
dependent_os <- "Surv(time, status)"
explanatory <- "RiBi_status"
tcga_eval %>% surv_plot(dependent_os, explanatory, pval = TRUE)
subset(tcga_eval, time < 12*5) %>% surv_plot(dependent_os, explanatory, pval = TRUE)


# Statistics (5y survval)

survdiff(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, time < 12*5))
pairwise_survdiff(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, time < 12*5))
coxph(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, time < 12*5))



## Multivariate COX analysis for TCGA

# Ensure age is numeric

tcga_eval$age <- as.numeric(tcga_eval$age)
tcga_eval <- tcga_eval %>% filter(time < 12*10)

# Create age group by median split

median_age_tcga <- median(tcga_eval$age, na.rm = TRUE)
tcga_eval <- tcga_eval %>%
  mutate(age_group = ifelse(age >= median_age_tcga, "High", "Low")) %>%
  mutate(age_group = factor(age_group, levels = c("Low", "High")))

# Convert relevant variables to factor (check that they exist in the dataset)

tcga_eval <- tcga_eval %>%
  mutate(across(c(RiBi_status, sex, NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53), ~ as.factor(.), .names = "factor_{.col}"))

# Rename factored columns if needed (or overwrite original columns)

tcga_eval <- tcga_eval %>%
  mutate(across(starts_with("factor_"), ~ ., .names = "{sub('factor_', '', .col)}"))


write.csv(tcga_eval, 'meta_validation.csv')

# Multivariate Cox model

cox_multi_tcga <- coxph(Surv(time, status) ~ RiBi_status + sex + NPM1 + DNMT3A + IDH1 + IDH2 + FLT3 + TP53 + age_group, data = tcga_eval)

# View summary

summary(cox_multi_tcga)

# Forest plot

ggforest(cox_multi_tcga, data = tcga_eval, main = "Multivariable Cox Proportional Hazards Model (TCGA)")






# Predictive Validation in the BEAT cohort - Keep treatment types with > 20 samples


beat_treat <- read.delim("beataml_wv1to4_clinical.txt") %>%
  select(dbgap_rnaseq_sample, cumulativeTreatmentTypes) %>%
  dplyr::rename(patient = dbgap_rnaseq_sample) %>% filter(!patient == '')

beat_pred_df <- left_join(beat_treat, meta_with_score, by = "patient") %>%
  filter(!is.na(cumulativeTreatmentTypes), cumulativeTreatmentTypes != "", !is.na(time)) %>%
  group_by(cumulativeTreatmentTypes) %>%
  filter(n() > 20) %>% 
  ungroup() %>% as.data.frame()
  

# Kaplan-Meier treatment per RiBi groups

beat_pred_df$time <- as.numeric(beat_pred_df$time)
beat_pred_df$status <- as.numeric(beat_pred_df$status)

dependent_os <- "Surv(time, status)"
explanatory <- c("RiBi_status", "cumulativeTreatmentTypes")
beat_pred_df %>% surv_plot(dependent_os, explanatory, pval = TRUE)
subset(beat_pred_df, time < 364*5) %>% surv_plot(dependent_os, explanatory, pval = TRUE) # 5y survival

write.csv(subset(beat_pred_df, time < 364*5), 'BEAT_predictive.csv')

# Statistics

fit_treat_1 <- surv_fit(Surv(time, status) ~ RiBi_status, data = subset(as.data.frame(beat_pred_df), time < 364*5), group.by = "cumulativeTreatmentTypes")
surv_pvalue(fit_treat_1)

fit_treat_2 <- surv_fit(Surv(time, status) ~ cumulativeTreatmentTypes, data = subset(as.data.frame(beat_pred_df), time < 364*5), group.by = "RiBi_status")
surv_pvalue(fit_treat_2)

# Multivariate

beat_pred_df %>% finalfit(dependent_os, explanatory)
beat_pred_df %>% hr_plot(dependent_os, explanatory)



