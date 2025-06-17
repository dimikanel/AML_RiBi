# Load relevant packages

library(tidyverse)
library(DESeq2)
library(limma)
library(umap) 


# Load ClinSeq-AML cohort (consists of two batches that will be merged)

# Load and prepare ClinSeq.1 metadata

cs1_meta <-  read.csv('ClinSeq_meta.csv')


#Exclude APL (FAB M3) and NA aml_subtype samples

cs1_meta <-  cs1_meta[!(cs1_meta$aml_subtype == 'APL, FAB M3' | is.na(cs1_meta$aml_subtype)), ] 

# Rename first column

colnames(cs1_meta)[1] <- 'patient'

# Add new columns and convert days to years

cs1_meta <- cs1_meta %>% mutate(cohort = 'ClinSeq_1', Disease_Status = 'Disease') %>% dplyr::rename(time = time_to_status_days)

# Select relevant columns and rename 'time_to_status_days' -> 'time'

cs1_meta <- cs1_meta %>% select(patient, time, status, cohort, Disease_Status, sex, age,
         NPM1, DNMT3A, DNMT3A_R882, IDH1, IDH2_140, IDH2_172,
         FLT3_ITD, FLT3_HS, TP53) 

# Create combined mutation columns

cs1_meta <- cs1_meta %>% mutate(DNMT3A_all = case_when(DNMT3A == 'wt' ~ 'wt', DNMT3A == 'mut' ~ 'mut', DNMT3A_R882 == 'wt' ~ 'wt', DNMT3A_R882 == 'mut' ~ 'mut'),
                                IDH2_all = case_when(IDH2_140 == 'wt' ~ 'wt', IDH2_140 == 'mut' ~ 'mut', IDH2_172 == 'wt' ~ 'wt', IDH2_172 == 'mut' ~ 'mut'),
                                FLT3_all = case_when(FLT3_ITD == 'wt' ~ 'wt', FLT3_ITD == 'mut' ~ 'mut', FLT3_HS == 'wt' ~ 'wt', FLT3_HS == 'mut' ~ 'mut'))

# Final selection of columns and renaming combined columns

cs1_meta <- cs1_meta %>% select(patient, time, status, cohort, Disease_Status, sex, age,
         NPM1, DNMT3A_all, IDH1, IDH2_all, FLT3_all, TP53) %>%
  dplyr::rename(DNMT3A = DNMT3A_all, IDH2 = IDH2_all, FLT3 = FLT3_all)

# Load and subset counts matrix

cs1_counts <-  read.table('ClinSeq.counts.txt',header=TRUE) %>% rename_with(~ gsub('\\.T.*', '', .x)) %>% select(cs1_meta$patient)


### Load ClinSeq.2

## Load metadata, modify sample names to match countrs matrix colnames, discard APL samples


cs2_meta <- read.csv('KAW_CLINICAL_FORMAT.txt', header = T, sep = '\t') %>% mutate(ID = gsub('_', '', ID)) %>% 
  mutate(ID = sub('A$', '', ID)) %>% mutate(ID = sub('D$', '', ID)) %>% 
  filter(!a_snomedaml_Beskrivning %in% 'Akut promyelocytleukemi, FAB M3(98663)')

cs2_mut  <- read.delim('KAW_MUTATION_FORMAT.txt') %>% dplyr::rename(ID = sample_ID) %>% 
  mutate(ID = gsub('_', '', ID)) %>% mutate(ID = sub('A$', '', ID)) %>% 
  mutate(ID = sub('D$', '', ID)) %>% select(1, 'NPM1', 'DNMT3A', 'IDH1', 'IDH2', 'FLT3', 'TP53') 

cs2_meta <- merge(cs2_meta, cs2_mut, by = 'ID', all = T)

# Load counts matrix 

cs2_counts <- read.csv('AML_KAW_GENCODEv39.csv', sep = ';', fileEncoding = 'UTF-8-BOM') %>% dplyr::mutate(X = sub('\\..*', '', X))


# Some colnames are mislabeled. Rename some columns to match the sample IDS in the metadata df
#  Metadata IDs not found in colnames

diff1 <- setdiff(colnames(cs2_counts), cs2_meta$ID)

#  Colnames not found in Metadata IDs 

diff2 <- setdiff(cs2_meta$ID, colnames(cs2_counts))

# Rename the colnames that are found in metadata IDs

cs2_counts <- cs2_counts %>% dplyr::rename(ALG201519 = ALG1519,  ALG201533 = ALG1533,  ALG20157 = ALG157)  

# Update the metadata df to include all samples that have RNASeq data. Then subset counts matrix.
# Find missing IDs

missing_ids <- setdiff(colnames(cs2_counts), cs2_meta$ID)

# Exclude columns that contain normal samples (NBM) and the gene names (X)

missing_ids <- missing_ids[!str_detect(missing_ids, 'NBM|X')] 

#Create empty rows for missing IDs

new_rows <- data.frame(ID = missing_ids, 
                       matrix(NA, nrow = length(missing_ids), ncol = ncol(cs2_meta) - 1))
colnames(new_rows) <- colnames(cs2_meta)

# Bind original and new rows

cs2_meta <- bind_rows(cs2_meta, new_rows)

# Some of the APL samples were reintroduced in the metadata df so it needs to be cleaned again.

cs2_meta_APL <- read.csv('KAW_CLINICAL_FORMAT.txt', header = T, sep = '\t') %>% mutate(ID = gsub('_', '', ID)) %>% 
  mutate(ID = sub('A$', '', ID)) %>% mutate(ID = sub('D$', '', ID)) %>% 
  filter(a_snomedaml_Beskrivning %in% 'Akut promyelocytleukemi, FAB M3(98663)')

cs2_meta <- cs2_meta %>% filter(!ID %in% cs2_meta_APL$ID)

cs2_counts <- cs2_counts %>% select(X, any_of(cs2_meta$ID), matches('NBM')) %>% 
  distinct(X, .keep_all = T) %>% column_to_rownames('X')


# Reshape metadata

cs2_meta <- cs2_meta %>% mutate(patient = ID, age = a_alder, cohort = 'ClinSeq_2', Disease_Status = 'Disease',
                    status = VITALSTATUS, sex = recode(GENDER, 'M' = 'male', 'F' = 'female')) %>% dplyr::rename(time = nr_days_of_survival) %>%
                        mutate(across(c(NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53), ~ recode(as.character(.), `1` = 'mut', `0` = 'wt'))) %>%
                        dplyr::select(patient, time, status, cohort, Disease_Status, sex, age, NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53)


# Add healthy patients from another data frame

healthy_ids <- colnames(cs2_counts)[str_detect(colnames(cs2_counts), 'NBM')]

healthy_meta <- tibble(patient = healthy_ids, time = NA_real_, status = NA_real_,
                       cohort = 'ClinSeq_2', Disease_Status = 'Healthy', sex = NA_character_, age = NA_real_,
                       NPM1 = NA_character_, DNMT3A = NA_character_, IDH1 = NA_character_, IDH2 = NA_character_,
                       FLT3 = NA_character_, TP53 = NA_character_)


# Combine with main metadata

cs2_meta <- bind_rows(cs2_meta, healthy_meta)
cs2_meta <- cs2_meta %>% filter(patient %in% colnames(cs2_counts))

table(cs2_meta$Disease_Status,useNA = 'always')

# Ensure matching counts

cs2_counts <- cs2_counts %>% dplyr::select(cs2_meta$patient)

# The two ClinSeq batches share some common samples. 
# Find the common patients between ClinSeq,1 and ClinSe2. and merge all count matrices and metadata df excluding duplicates 

intersect(names(cs1_counts), names(cs2_counts))

cs1_counts <- cs1_counts %>% rownames_to_column('X')
cs2_counts <- cs2_counts %>% rownames_to_column('X')

cs_counts <- inner_join(cs1_counts, cs2_counts, by = 'X') %>%
  dplyr::select(-ends_with('.y')) %>%
  rename_with(~ sub('\\.x$', '', .x), ends_with('.x'))

cs_meta <- bind_rows(cs1_meta, cs2_meta) %>% distinct(patient, .keep_all = TRUE) #%>% mutate(cohort = 'ClinSeq')

# Keep only matched patients in counts and metadata

cs_meta <- cs_meta %>% filter(patient %in% colnames(cs_counts))
cs_counts <- cs_counts %>% dplyr::select(1, all_of(cs_meta$patient))

table(cs_meta$Disease_Status, useNA = 'always')
table(cs_meta$Disease_Status, cs_meta$status,useNA = 'always')

write.csv(cs_counts, 'cs_counts.csv')
write.csv(cs_meta, 'cs_meta.csv')



# Load BEAT-AML cohort

# Load counts matrix, filter only protein coding genes. 
# Since there might be a lack of consensus between gene name translation between all cohorts we will subset all 
# based on BEAT-AML data that are readily available

beat_counts <- read.delim('beataml_waves1to4_counts_dbgap.txt', header = TRUE, check.names = FALSE) %>%
  filter(biotype == 'protein_coding') %>% 
  dplyr::select(ensgene = 1, 5:711) %>% column_to_rownames('ensgene')


# Load metadata, discard samples without RNASeq data and keep only 'Healthy pooled CD34+' as healthy controls 
# to match similar controls from the other cohorts


beat_meta <- read.csv('beataml_waves1to4_sample_mapping.csv', header = T, sep = ';') %>% dplyr::select(4,6) %>% 
   filter(!rna_control %in% c('Healthy Individual CD34+','Healthy Individual BM MNC')) %>%
   mutate(Disease_Status = if_else(str_detect(rna_control, 'Healthy pooled CD34+'), 'Healthy', 'Disease')) %>%
   dplyr::select(1,3) %>% filter(!dbgap_rnaseq_sample == '' ) %>% arrange(dbgap_rnaseq_sample)

# Load clinical metadata. Discard APL and NA samples and those that do not refer to initial diagnosis (i.e. refractory etc)

beat_clin <- read.delim('beataml_wv1to4_clinical.txt') %>% filter(!dbgap_rnaseq_sample == '') %>% 
  filter(diseaseStageAtSpecimenCollection == 'Initial Diagnosis') %>% 
  filter(!specificDxAtInclusion %in% c('Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA', 'Unknown')) 

beat_meta <- merge(beat_clin, beat_meta, by = 'dbgap_rnaseq_sample', all = TRUE)
beat_meta <- beat_meta %>% dplyr::rename(sample = dbgap_rnaseq_sample)
beat_meta <- beat_meta %>% filter(sample %in% colnames(beat_counts))

# Mutation parsing from variantSummary

max_pipes <- max(str_count(beat_meta$variantSummary, '\\|'), na.rm = TRUE) + 1

beat_meta <- beat_meta %>%
  separate(variantSummary, into = paste0('mut', 1:max_pipes), sep = '\\|', fill = 'right') %>%
  rowwise() %>% mutate(DNMT3A = if_else(any(str_detect(c_across(starts_with('mut')), 'DNMT3A')), 'mut', 'wt'),
                       IDH1   = if_else(any(str_detect(c_across(starts_with('mut')), 'IDH1')), 'mut', 'wt'),
                       IDH2   = if_else(any(str_detect(c_across(starts_with('mut')), 'IDH2')), 'mut', 'wt'),
                       FLT3   = if_else(any(str_detect(c_across(starts_with('mut')), 'FLT3|positive')), 'mut', 'wt'),
                       TP53   = if_else(any(str_detect(c_across(starts_with('mut')), 'TP53')), 'mut', 'wt'),
                       NPM1   = if_else(str_detect(NPM1, 'positive|NPM1'), 'mut', 'wt'),
                       status = recode(vitalStatus, 'Dead' = 1, 'Alive' = 0, .default = NA_real_),
                       sex    = recode(consensus_sex, 'Female' = 'female', 'Male' = 'male'),
                       cohort = 'BEAT') %>% dplyr::rename(time = overallSurvival)%>% ungroup() %>%
                dplyr::select(sample, time, status, cohort, Disease_Status, sex, age = ageAtDiagnosis,
                NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53) %>%  dplyr::rename(patient = sample)


beat_meta <- beat_meta %>% filter(!(is.na(status) & Disease_Status != 'Healthy'))

# Align counts with metadata

beat_counts <- beat_counts[, beat_meta$patient]


table(beat_meta$Disease_Status)
table(beat_meta$Disease_Status, beat_meta$status,useNA = 'always')

write.csv(beat_counts, 'beat_counts.csv')
write.csv(beat_meta, 'beat_meta.csv')


# Load the FINF-AML cohort

# Load the counts matrix

finn_counts <- read.csv('File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv') %>% column_to_rownames('X')


# Load the metadata, rename, keep only samples with available RNASeq data and samples initially diagnosed as AML

finn_annot <- read.csv('File_0_Common_sample_annotation_252S.csv') %>%
  separate(Sample_ID, c('patient_ID_1', 'rep'), sep = '_(?=[^_]+$)', remove = FALSE) %>%
  filter(rep == '01', RNA.seq == 'available') %>% filter(Sample_ID %in% colnames(finn_counts)) %>% 
  dplyr::rename(stage = Disease.status) %>% mutate(stage = ifelse(stage == 'Diagnosis', 'Initial Diagnosis', stage)) %>% 
  filter(stage == 'Initial Diagnosis')


finn_counts <- finn_counts %>% dplyr::select(contains('Healthy'), finn_annot$Sample_ID)

# Load clinical metadata, discard APL and NA samples

finn_clin <- read.csv('File_1.1_Clinical_summary_186_Patients.csv', sep = ';') %>%
  dplyr::rename(patient_ID_1 = 1) %>%
  filter(ICD.O != '9866 Acute promyelocytic leuk.,t(15;17)(q22;q11-12)',!is.na(ICD.O),ICD.O != '')

# Merge the two meta df

finn_meta_pre <- merge(finn_annot,finn_clin, by = 'patient_ID_1')
colnames(finn_meta_pre)[3] <- 'patient'

# Load mutational metadata, reshape to match the rest of the cohorts


finn_mut <- read.csv('File_6_Binary_mutation_225S_57G.csv', sep = ';') %>%
  filter(Gene_Symbol %in% c('NPM1', 'DNMT3A', 'IDH1', 'IDH2', 'FLT3', 'TP53')) %>%
  dplyr::select(Gene_Symbol, starts_with('AML_')) %>%
  pivot_longer(cols = -Gene_Symbol, names_to = 'patient', values_to = 'mut') %>%
  pivot_wider(names_from = Gene_Symbol, values_from = mut) %>%
  separate(col = patient, into = c('patient_ID_1', 'rep'), sep = '_(?=[^_]+$)', remove = FALSE) %>%
  relocate(patient, patient_ID_1, rep)



# Merge mut with the rest of the metadata


finn_meta <- merge(finn_meta_pre, finn_mut, by = 'patient', all = T) %>% 
         filter(!is.na(Patient_ID)) %>% 
          mutate(cohort = 'FINN', Disease_Status = 'Disease', 
         status = as.numeric(os_event), sex = recode(Gender, 'Female' = 'female', 'Male' = 'male'),
         age = round(as.numeric(str_replace_all(Age.at.diagnosis, ',', '.')))) %>% dplyr::rename(time = os_time) %>%
         mutate(across(c(NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53),
                ~ recode(as.character(.), `1` = 'mut', `0` = 'wt'))) %>%
  dplyr::select(patient, time, status, cohort, Disease_Status, sex, age, NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53)


finn_annot <- finn_annot %>% dplyr::select(patient = Sample_ID, stage)
finn_meta <- merge(finn_meta, finn_annot, by = 'patient')


# Add 4 healthy patients 

healthy_ids <- colnames(finn_counts)[str_detect(colnames(finn_counts), 'Healthy')]

healthy_meta <- tibble(
  patient = healthy_ids, time = NA_real_, status = NA_real_, cohort = 'FINN',
  Disease_Status = 'Healthy', sex = NA_character_, age = NA_real_, NPM1 = NA_character_,
  DNMT3A = NA_character_, IDH1 = NA_character_, IDH2 = NA_character_,
  FLT3 = NA_character_, TP53 = NA_character_)


# Combine with main metadata

finn_meta <- bind_rows(finn_meta, healthy_meta)
finn_meta <- finn_meta %>% filter(patient %in% colnames(finn_counts), !(is.na(status) & Disease_Status == 'Disease'))

# Align FINN counts

finn_counts <- finn_counts %>% dplyr::select(finn_meta$patient)

table(finn_meta$Disease_Status)
table(finn_meta$Disease_Status, finn_meta$status,useNA = 'always')

write.csv(finn_counts, 'finn_counts.csv')
write.csv(finn_meta, 'finn_meta.csv')




# Merge All Cohorts 


beat_counts <- beat_counts %>% rownames_to_column('X')
finn_counts <- finn_counts %>% rownames_to_column('X')

common_genes <- Reduce(intersect, list(cs_counts$X, beat_counts$X, finn_counts$X))

cs_counts <- cs_counts %>% filter(X %in% common_genes)
beat_counts <- beat_counts  %>% filter(X %in% common_genes)
finn_counts <- finn_counts  %>% filter(X %in% common_genes)


# Step 3: Merge by 'X' using purrr::reduce

all_counts <- list(cs_counts, beat_counts, finn_counts) %>%
  purrr::reduce(~ merge(.x, .y, by = 'X')) %>%
  distinct(X, .keep_all = TRUE) %>%
  column_to_rownames('X')


# Clean duplicated columns

all_counts <- all_counts[, !duplicated(names(all_counts))]

# Merge all metadata

meta_all <- bind_rows(cs_meta, beat_meta, finn_meta) %>%
  distinct(patient, .keep_all = TRUE) %>%
  filter(patient %in% colnames(all_counts)) 


# Align columns and rows

all_counts <- all_counts[, meta_all$patient]
stopifnot(all.equal(colnames(all_counts), meta_all$patient))

write.csv(all_counts, 'all_counts_refined.csv')


# Check the composition of the merged metadata df

table(meta_all$Disease_Status, useNA = 'always')
table(meta_all$cohort, useNA = 'always')
table(meta_all$cohort, meta_all$Disease_Status, useNA = 'always')

write.csv(meta_all, 'meta_all_refined.csv')



# DESeq normalization and DE analysis

dds <- DESeqDataSetFromMatrix(countData = all_counts,
                             colData = meta_all,
                             design = ~ cohort + Disease_Status)


# Filter low-count genes

dds <- dds[rowSums(counts(dds)) >= 10, ]

# Apply variance-stabilizing transformation

vst_counts <- vst(dds)


# Remove batch effects and Run UMAP analysis

assay(vst_counts) <- removeBatchEffect(assay(vst_counts),vst_counts$cohort)

norm_counts <- assay(vst_counts) %>% t()
umap_results <- umap::umap(norm_counts)
umap_plot_df <- data.frame(umap_results$layout)

ggplot(umap_plot_df, aes(x = X1, y = X2, color = meta_all$Disease_Status)) + geom_point()


dds <- DESeq(dds) 

resultsNames(dds)
res <- results(dds, alpha=0.05, contrast = c('Disease_Status', 'Disease', 'Healthy')) 
sum(res$padj < 0.05, na.rm=TRUE) 

res_05 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0.5)
summary(res_05)

write.csv(as.data.frame(res_05),file = 'DE_ClinSeq_BEAT_FINN_AML_Vs_Healthy_filtered_05.csv')



