#setwd('C:/Users/Dimitrios Kanellis/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')
#setwd("C:/Users/Sofia/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready")
#setwd('C:/Users/dimka/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')

# Load required libraries

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
library(broom)
library(tidyverse)

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


finn_meta <- finn_meta %>% select(1:15) %>% filter(!Disease_Status == "Healthy")
finn_counts <- finn_counts %>% dplyr::rename(ensembl_ID = X) %>%
  distinct(ensembl_ID, .keep_all = TRUE) %>% select(ensembl_ID, finn_meta$patient)


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

# Requires glmnet 4.1-8 (see README). Fixed folds reproduce the published 31-gene signature.

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
coef_table
# Save to CSV

write.csv(coef_table, "RiBi_model_coefficients_rerun.csv", row.names = FALSE)

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
ggforest(cox_multi, data = meta_with_score, main = "Multivariable Cox Proportional Hazards Model", fontsize = 1)

# Save the hazards ratio analysis

cox_df <- tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value)

write.csv(cox_df, 'cox_multi_training.csv')


# Add the ELN2022 classification and run the HR analysis (training set only)

ELN2022_df <- read.csv('ELN2022_all.csv') %>% select(patient, eln2022)

meta_with_score_ELN <- merge(meta_with_score, ELN2022_df, by = 'patient')
cox_multi_ELN <- coxph(Surv(time, status) ~ RiBi_status + age_group + eln2022, data = meta_with_score_ELN)

summary(cox_multi_ELN)
ggforest(cox_multi_ELN, data = meta_with_score_ELN, main = "Multivariable Cox Proportional Hazards Model", fontsize = 1)


# Save the hazards ratio analysis

cox_df_eln <- tidy(cox_multi_ELN, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value)
write.csv(cox_df_eln, 'cox_multi_eln_training.csv')



# Load and Process TCGA Data 

library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

#files <- getResults(query)
#nrow(files)

GDCdownload(
  query,
  directory = "GDCdata",
  method = "api",
  files.per.chunk = 20
)

laml <- GDCprepare(query, directory = "GDCdata")

#assayNames(laml)

counts <- assay(laml, "unstranded")

# Retrieve metadata

sample_info <- as.data.frame(colData(laml))
gene_info <- as.data.frame(rowData(laml))

dim(counts)
counts[1:5, 1:5]


# Extract the raw counts
tcga_counts <- as.matrix(counts)
rownames(tcga_counts) <- sub("\\..*$", "", rownames(tcga_counts))

sum(duplicated(rownames(tcga_counts)))

if (anyDuplicated(rownames(tcga_counts))) {
  tcga_counts <- rowsum(
    tcga_counts,
    group = rownames(tcga_counts),
    reorder = FALSE
  )
}


tcga_metadata <- sample_info[match(colnames(tcga_counts), rownames(sample_info)),, drop = FALSE]

stopifnot(identical(colnames(tcga_counts), rownames(tcga_metadata)))

#dim(tcga_counts)
#tcga_counts[1:5, 1:5]
#colnames(tcga_metadata)



# Read clinical files

tcga_clin <- read.delim("data_clinical_sample.txt", comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
tcga_clin_surv <- read.delim("data_clinical_patient.txt", comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)


patient_id <- function(x) {
  substr(trimws(as.character(x)), 1, 12)
}

stopifnot("PATIENT_ID" %in% names(tcga_clin))

tcga_clin <- tcga_clin %>% transmute(patient = patient_id(PATIENT_ID)) %>% distinct()

# Read mutations

tcga_mut_all <- read.delim("data_mutations.txt", comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(all(c("Hugo_Symbol", "Tumor_Sample_Barcode") %in% names(tcga_mut_all)))

tcga_mut_all <- tcga_mut_all %>% mutate(patient = patient_id(Tumor_Sample_Barcode))

genes <- c("NPM1", "DNMT3A", "IDH1", "IDH2", "FLT3", "TP53")

tcga_mut <- tcga_mut_all %>% filter(Hugo_Symbol %in% genes)

patients_with_mutation_data <- unique(tcga_mut_all$patient)

for (gene in genes) {
  
  mutated_patients <- tcga_mut %>%
    filter(Hugo_Symbol == gene) %>%
    pull(patient) %>%
    unique()
  
  tcga_clin[[gene]] <- case_when(
    tcga_clin$patient %in% mutated_patients ~ "mut",
    tcga_clin$patient %in% patients_with_mutation_data ~ "wt",
    TRUE ~ NA_character_
  )
}


required_fields <- c("PATIENT_ID", "CYTOGENETIC_CODE_OTHER", "OS_MONTHS", "OS_STATUS", "SEX", "AGE")

stopifnot(all(required_fields %in% names(tcga_clin_surv)))

tcga_clin_surv <- tcga_clin_surv %>% filter(!CYTOGENETIC_CODE_OTHER %in% c("PML-RARA", "N.D.")) %>%
  transmute(patient = patient_id(PATIENT_ID), time = suppressWarnings(as.numeric(OS_MONTHS)),
    status = case_when(toupper(trimws(OS_STATUS)) %in% c("1:DECEASED", "DECEASED", "1") ~ 1,toupper(trimws(OS_STATUS)) %in%
        c("0:LIVING", "LIVING", "0") ~ 0,
      TRUE ~ NA_real_), sex = case_when(tolower(trimws(SEX)) == "female" ~ "female",
      tolower(trimws(SEX)) == "male" ~ "male", TRUE ~ NA_character_), age = suppressWarnings(as.numeric(AGE)))

stopifnot(!anyDuplicated(tcga_clin_surv$patient))


tcga_meta <- inner_join(tcga_clin, tcga_clin_surv, by = "patient") %>%
  mutate(cohort = "TCGA", Disease_Status = "Disease") %>%
  select(patient, time, status, cohort, Disease_Status, sex, age, all_of(genes))


stopifnot(identical(colnames(tcga_counts), rownames(tcga_metadata)))


if ("patient" %in% names(tcga_metadata)) {
  
  count_patients <- patient_id(tcga_metadata$patient)
  
} else if ("barcode" %in% names(tcga_metadata)) {
  
  count_patients <- patient_id(tcga_metadata$barcode)
  
} else {
  
  count_patients <- patient_id(colnames(tcga_counts))
}


stopifnot(all(!is.na(count_patients)), all(grepl("^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}$", count_patients)))


if (anyDuplicated(count_patients)) {
  stop("Multiple RNA-seq samples per patient. Select one sample per patient before continuing.")
}


meta_index <- match(count_patients, tcga_meta$patient)
keep <- !is.na(meta_index)

tcga_counts <- tcga_counts[, keep, drop = FALSE]
tcga_meta <- tcga_meta[meta_index[keep], , drop = FALSE]


tcga_meta$rna_sample_id <- colnames(tcga_counts)


colnames(tcga_counts) <- tcga_meta$patient
rownames(tcga_meta) <- tcga_meta$patient

stopifnot(identical(colnames(tcga_counts), tcga_meta$patient))


#dim(tcga_counts)
#dim(tcga_meta)
#table(tcga_meta$status, useNA = "ifany")

sapply(tcga_meta[, genes, drop = FALSE],
  function(x) table(factor(x, levels = c("wt", "mut")),
                    useNA = "ifany")
)

#write.csv(tcga_counts, 'tcga_counts.csv')
#write.csv(tcga_meta, 'tcga_meta.csv')


tcga_coldata <- as.data.frame(tcga_meta)
rownames(tcga_coldata) <- tcga_coldata$patient


stopifnot(all(colnames(tcga_counts) %in% rownames(tcga_coldata)))

tcga_coldata <- tcga_coldata[colnames(tcga_counts), , drop = FALSE]

stopifnot(identical(colnames(tcga_counts), rownames(tcga_coldata)))


dds_tcga <- DESeqDataSetFromMatrix(countData = round(as.matrix(tcga_counts)),
  colData = tcga_coldata, design = ~ 1)

dds_tcga <- dds_tcga[rowSums(counts(dds_tcga)) >= 10, ]
vst_tcga <- vst(dds_tcga)
norm_tcga <- assay(vst_tcga)


tcga_score_matrix <- t(norm_tcga[features, , drop = FALSE])
tcga_weighted_matrix <- t(tcga_score_matrix) * weights

# Score per patient

tcga_score_df <- data.frame(patient = rownames(tcga_score_matrix), RiBi_score = colSums(tcga_weighted_matrix, na.rm = TRUE))

tcga_score_df <- tcga_score_df %>% mutate(RiBi_status = case_when(RiBi_score > quantile(RiBi_score, 0.75) ~ "high",
    RiBi_score < quantile(RiBi_score, 0.25) ~ "low", TRUE ~ "moderate"))

tcga_eval <- tcga_meta %>% left_join(tcga_score_df, by = "patient") %>% filter(time > 0)

roc_obj <- ggplot(tcga_eval, aes(d = status, m = RiBi_score)) +
  geom_roc(color = "darkblue", n.cuts = 0) +
  geom_abline(intercept = 0)

auc_val <- calc_auc(roc_obj)$AUC

roc_plot <- roc_obj + annotate("text", x = 0.75, y = 0.25, label = paste("AUC =", round(auc_val, 3))) +
  labs(title = "TCGA ROC", x = "1 - Specificity", y = "Sensitivity")

print(roc_plot)


tcga_eval$time <- as.numeric(tcga_eval$time)
tcga_eval$year <- as.numeric(tcga_eval$time/12)
dependent_os <- "Surv(year, status)"
explanatory <- "RiBi_status"
tcga_eval %>% surv_plot(dependent_os, explanatory, pval = TRUE)
subset(tcga_eval, year < 5) %>% surv_plot(dependent_os, explanatory, pval = TRUE)


# Statistics (5y survval)

survdiff(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, year < 5))
pairwise_survdiff(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, year < 5))
coxph(Surv(time, status) ~ RiBi_status, data = subset(tcga_eval, year < 5))



# Multivariate COX analysis for TCGA

# Ensure age is numeric

tcga_eval$age <- as.numeric(tcga_eval$age)
tcga_eval <- tcga_eval %>% filter(time < 12*10)

# Create age group by median split

median_age_tcga <- median(tcga_eval$age, na.rm = TRUE)
tcga_eval <- tcga_eval %>%
  mutate(age_group = ifelse(age >= median_age_tcga, "High", "Low")) %>%
  mutate(age_group = factor(age_group, levels = c("Low", "High")))


write.csv(tcga_eval, 'meta_validation.csv')

# Multivariate Cox model

cox_multi_tcga <- coxph(Surv(time, status) ~ RiBi_status + sex + NPM1 + DNMT3A + IDH1 + IDH2 + FLT3 + TP53 + age_group, data = tcga_eval)

# View summary

summary(cox_multi_tcga)

# Forest plot

ggforest(cox_multi_tcga, data = tcga_eval, main = "Multivariable Cox Proportional Hazards Model (TCGA)")



# Merge all samples (training, validation) and run the multivariate COX

colnames(meta_with_score)
colnames(tcga_eval)

meta_with_score_1 <- meta_with_score %>% select(-X)

colnames(meta_with_score_1)
colnames(tcga_eval)

meta_with_score_1$year <- as.numeric(meta_with_score_1$time/365)
colnames(meta_with_score_1)
colnames(tcga_eval)

meta_with_score_2 <- meta_with_score_1 %>% select(1:16,18,17)
colnames(meta_with_score_2)
colnames(tcga_eval)

all_meta <- rbind(meta_with_score_2, tcga_eval)


# Fit Cox proportional hazards model

cox_all_meta <- coxph(Surv(time, status) ~ RiBi_status + sex + NPM1 + DNMT3A + IDH1 + IDH2 + FLT3 + TP53 + age_group, data = all_meta)

summary(cox_all_meta)
ggforest(cox_all_meta, data = meta_with_score, main = "Multivariable Cox Proportional Hazards Model", fontsize = 1)

# Save the hazards ratio analysis

cox_df_all_meta <- tidy(cox_all_meta, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value)
write.csv(cox_df_all_meta, 'cox_df_all_meta.csv')


# Repeat with ENL2022 

all_meta_ELN <- merge(all_meta, ELN2022_df, by = 'patient')

cox_all_meta_ELN <- coxph(Surv(year, status) ~ RiBi_status + age_group + eln2022, data = all_meta_ELN)

summary(cox_all_meta_ELN)
ggforest(cox_all_meta_ELN, data = all_meta_ELN, main = "Multivariable Cox Proportional Hazards Model", fontsize = 1)


# Save the hazards ratio analysis

cox_df_all_meta_eln <- tidy(cox_all_meta_ELN, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value)
write.csv(cox_df_all_meta_eln, 'cox_df_all_meta_eln.csv')



# Predictive Validation in the BEAT cohort - Keep treatment types with > 20 samples


beat_treat <- read.delim("beataml_wv1to4_clinical.txt") %>%
  select(dbgap_rnaseq_sample, cumulativeTreatmentTypes) %>%
  rename(patient = dbgap_rnaseq_sample)

table(beat_treat$cumulativeTreatmentTypes)

beat_pred_df <- left_join(beat_treat, meta_with_score, by = "patient") %>%
  filter(!is.na(cumulativeTreatmentTypes), cumulativeTreatmentTypes != "") %>% 
  filter(cumulativeTreatmentTypes %in% c('Standard Chemotherapy', 'Standard Chemotherapy|Bone Marrow Transplant')) %>%
  mutate(years = time / 365)
  

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


