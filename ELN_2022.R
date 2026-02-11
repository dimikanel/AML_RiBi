setwd('C:/Users/Dimitrios Kanellis/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')
#setwd('C:/Users/Sofia/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')


library(tidyverse)
library(readxl)

# Load ClinSeq_1 data

CS_meta <-  read.csv('ClinSeq_meta.csv') %>% filter(!aml_subtype == 'APL, FAB M3') %>% 
  rename(patient = patient_id, karyotype = Curated.karyotype) %>% mutate(cohort = 'ClinSeq', cytogenetics = NA, CEBPA_bzip = NA, BCOR = NA, ZRSR2 = NA ) %>% 
  select(patient, cohort, cytogenetics, karyotype, NPM1, FLT3_ITD, CEBPA_bzip, TP53, ASXL1, BCOR, 
         EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2)

colnames(CS_meta)

# Load ClinSeq_2 (KAW) data


KAW_meta_clin <-  read.csv('KAW_CLINICAL_FORMAT.txt', header = T, sep = '\t') %>% filter(!a_snomedaml_Beskrivning %in% 'Akut promyelocytleukemi, FAB M3(98663)') %>%
  rename(patient = ID) %>% select(patient, a_karyotyp)
colnames(KAW_meta_clin)

KAW_meta_mut <- read.delim('KAW_MUTATION_FORMAT.txt') %>% rename(patient = sample_ID) %>%
  select(patient, TP53, NPM1, ASXL1, EZH2, RUNX1, SRSF2, STAG2, U2AF1, ZRSR2) %>%
  mutate(across(-patient, ~ case_when(.x == 1 ~ "mut", .x == 0 ~ "wt", TRUE ~ NA_character_)))

KAW_meta_FLT3 <- read.csv('KAW_Resultatsammanstallning_FLT3-ITD_ratio.csv', sep = ';') %>% rename(FLT3_ITD = Allelic_ratio, patient = Prov_id) %>%
  select(patient, FLT3_ITD) %>%  filter(grepl("ALBB|ALG16", patient)) %>% mutate(patient = sub("(.{3})$", "_\\1", patient)) %>%
  mutate(FLT3_ITD = case_when(is.na(FLT3_ITD) ~ NA_character_, FLT3_ITD == 0 ~ "wt", TRUE ~ "mut"))


KAW_meta_mut_alt <- read_excel('KAW_mutations.xlsx', sheet = 1) %>% rename(patient = CaseSamples) 

has_pos_280_340 <- function(x) {
  x <- as.character(x)
  nums <- str_extract_all(x, "\\d+")
  vapply(nums, function(v) {
    if (length(v) == 0) return(FALSE)
    any(as.integer(v) >= 280 & as.integer(v) <= 340)
  }, logical(1))
}

KAW_meta_mut_alt_1 <- KAW_meta_mut_alt %>%
  mutate(GeneSymbol = str_squish(str_replace_all(as.character(GeneSymbol), "[\u00A0]", " ")),
         ProteinVariant = as.character(ProteinVariant),
         CEBPA_present_row = str_detect(GeneSymbol, "\\bCEBPA\\b"),
         BCOR_present_row  = str_detect(GeneSymbol, "\\bBCOR\\b"),
         CEBPA_bzip_row = CEBPA_present_row & has_pos_280_340(ProteinVariant)) %>%
  group_by(patient) %>% summarise(CEBPA_bzip = if_else(any(CEBPA_bzip_row, na.rm = TRUE), "mut", "wt"),
                                  BCOR = if_else(any(BCOR_present_row,  na.rm = TRUE), "mut", "wt"), .groups = "drop")

KAW_meta_mut_alt_1


KAW_meta <- merge(KAW_meta_clin, KAW_meta_mut, by = 'patient')
KAW_meta_1 <- merge(KAW_meta, KAW_meta_FLT3, by = 'patient', all.x = T)
KAW_meta_2 <- merge(KAW_meta_1, KAW_meta_mut_alt_1, by = 'patient', all.x = T)


KAW_meta_3 <- KAW_meta_2 %>% rename(karyotype = a_karyotyp)  %>% 
             mutate(cohort = 'KAW', cytogenetics = NA, SF3B1 = NA ) %>%
            select(patient, cohort, cytogenetics, karyotype, NPM1, FLT3_ITD, CEBPA_bzip, TP53, ASXL1, BCOR, 
                   EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2)

# Merge ClinSeq and KAW data 

CS_KAW_meta <- rbind(CS_meta, KAW_meta_3)



# Classify based on ELN2022 (https://pubmed.ncbi.nlm.nih.gov/35797463/)


# Helper functions

has_pat <- function(x, pattern) {
  if_else(is.na(x), NA, str_detect(x, regex(pattern, ignore_case = TRUE)))
}

mr_genes <- c("ASXL1","BCOR","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")

is_mut <- function(x) {
  if_else(is.na(x), NA, as.character(x) %in% c("mut", "1", "TRUE", "T"))
}

CS_KAW_meta_eln <- CS_KAW_meta %>%
  mutate(
    FLT3_ITD_mut   = is_mut(FLT3_ITD),
    CEBPA_bzip_mut = is_mut(CEBPA_bzip),
    NPM1_mut       = is_mut(NPM1),
    TP53_mut       = is_mut(TP53),
    
    MR_gene_mut = case_when(
      if_all(all_of(mr_genes), ~ is.na(.x)) ~ NA,
      if_any(all_of(mr_genes), ~ is_mut(.x)) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    t_8_21 = has_pat(cytogenetics, "t\\(8;21\\)|RUNX1\\s*::\\s*RUNX1T1|RUNX1T1\\s*::\\s*RUNX1"),
    inv16  = has_pat(cytogenetics, "inv\\(16\\)|t\\(16;16\\)|CBFB\\s*::\\s*MYH11|MYH11\\s*::\\s*CBFB"),
    t_9_11 = has_pat(cytogenetics, "t\\(9;11\\)|MLLT3\\s*::\\s*KMT2A|KMT2A\\s*::\\s*MLLT3"),
    t_9_22 = has_pat(cytogenetics, "t\\(9;22\\)|BCR\\s*::\\s*ABL1|ABL1\\s*::\\s*BCR"),
    t_6_9  = has_pat(cytogenetics, "t\\(6;9\\)|DEK\\s*::\\s*NUP214|NUP214\\s*::\\s*DEK"),
    t_8_16 = has_pat(cytogenetics, "t\\(8;16\\)|KAT6A\\s*::\\s*CREBBP|CREBBP\\s*::\\s*KAT6A"),
    
    inv3_t3_3   = has_pat(cytogenetics, "inv\\(3\\)|t\\(3;3\\)"),
    mecom_rearr = has_pat(cytogenetics, "MECOM|EVI1"),
    
    mec1om_rearr = if_else(
      is.na(inv3_t3_3) & is.na(mecom_rearr),
      NA,
      coalesce(inv3_t3_3, FALSE) | coalesce(mecom_rearr, FALSE)
    ),
    
    kmt2a_any   = has_pat(cytogenetics, "KMT2A|\\bMLL\\b"),
    kmt2a_other = if_else(is.na(kmt2a_any), NA, kmt2a_any & !coalesce(t_9_11, FALSE)),
    
    minus5_del5q = has_pat(karyotype, "\\-5\\b|del\\(5q\\)|5q\\-"),
    minus7       = has_pat(karyotype, "\\-7\\b|del\\(7q\\)|7q\\-"),
    abn17p       = has_pat(karyotype, "abn\\(17p\\)|del\\(17p\\)|17p\\-|\\-17\\b"),
    
    plus8  = has_pat(karyotype, "\\+8\\b"),
    del9q  = has_pat(karyotype, "del\\(9q\\)"),
    del11q = has_pat(karyotype, "del\\(11q\\)"),
    del12p = has_pat(karyotype, "del\\(12p\\)"),
    del13q = has_pat(karyotype, "del\\(13q\\)"),
    del20q = has_pat(karyotype, "del\\(20q\\)"),
    
    monosomal_karyo = has_pat(karyotype, "monosomal|\\bmk\\b"),
    complex_karyo = case_when(
      is.na(karyotype) ~ NA,
      has_pat(karyotype, "complex") ~ TRUE,
      str_detect(str_squish(karyotype), "^(46,(XX|XY)|NK)$") ~ FALSE,
      TRUE ~ TRUE
    )
  )

CS_KAW_meta_eln_1 <- CS_KAW_meta_eln %>%
  mutate(
    eln2022 = case_when(
     
      coalesce(t_9_22, FALSE) | coalesce(t_6_9, FALSE) | coalesce(kmt2a_other, FALSE) |
        coalesce(t_8_16, FALSE) | coalesce(inv3_t3_3, FALSE) | coalesce(mec1om_rearr, FALSE) |
        coalesce(minus5_del5q, FALSE) | coalesce(minus7, FALSE) | coalesce(abn17p, FALSE) |
        coalesce(complex_karyo, FALSE) | coalesce(monosomal_karyo, FALSE) |
        coalesce(MR_gene_mut, FALSE) | coalesce(TP53_mut, FALSE) ~ "Adverse",
      
      coalesce(t_8_21, FALSE) | coalesce(inv16, FALSE) | coalesce(CEBPA_bzip_mut, FALSE) ~ "Favorable",
      
      coalesce(NPM1_mut, FALSE) & is.na(FLT3_ITD_mut) ~ NA_character_,   
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == FALSE) ~ "Favorable",
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == TRUE)  ~ "Intermediate",
      
    TRUE ~ "Intermediate"), eln2022 = factor(eln2022, levels = c("Favorable", "Intermediate", "Adverse")))


# Load BEAT AML data

BEAT_meta_1 <- read.csv('beataml_waves1to4_sample_mapping.csv', header = T, sep = ";") %>% filter(!grepl("Healthy", rna_control)) %>%
 filter(!dbgap_rnaseq_sample == '') %>% arrange(dbgap_rnaseq_sample)


BEAT_meta_2 <- read.delim('beataml_wv1to4_clinical.txt') %>% filter(!dbgap_rnaseq_sample == '') %>% 
                          arrange(dbgap_rnaseq_sample)

BEAT_meta <- merge(BEAT_meta_1, BEAT_meta_2, all=T)

BEAT_meta_x1 <- BEAT_meta %>% filter(specimenGroups == "Initial Acute Leukemia Diagnosis" & !specificDxAtInclusion %in% 
                              c("Acute promyelocytic leukaemia with t(15;17)(q22;q12); PML-RARA","Unknown")) %>% 
  mutate(cohort = 'BEAT') %>% rename(patient = dbgap_rnaseq_sample) %>% 
  select(patient, cohort, consensusAMLFusions, karyotype, otherCytogenetics, FLT3.ITD, NPM1, RUNX1, ASXL1, TP53, variantSummary ) %>% 
  separate(col = variantSummary, into = paste0("mut_", 1:11), sep = "\\|", fill = "right", extra = "drop") 


# Check if all ELN2022 relevant genes are screened in BEAT-AML

genes <- c("FLT3","DNMT3A","NPM1","TP53","ASXL1","BCOR","CEBPA","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")

mut_mat <- BEAT_meta_x1[, 7:13]

gene_present <- sapply(genes, function(g) {
  any(grepl(g, unlist(mut_mat)), na.rm = TRUE)
})

gene_present


# Add columns with mutational status

mutcols <- grep("^mut_\\d+$", names(BEAT_meta_x1), value = TRUE)

mut_mat <- as.matrix(BEAT_meta_x1[, mutcols])
row_mut_str <- apply(mut_mat, 1, function(x) paste(na.omit(x), collapse = " | "))

BEAT_meta_x2 <- BEAT_meta_x1

for (g in genes) {
  pat <- paste0("\\b", g, "\\b")
  BEAT_meta_x2[[g]] <- ifelse(grepl(pat, row_mut_str), "mut", "wt")
}

BEAT_meta_x2$row_mut_str <- row_mut_str

# For CEBPA try to find mutations in the bzip area (aa 280-340)

# Approximate CEBPA bZIP amino-acid window 

BZIP_MIN <- 280
BZIP_MAX <- 340

# Helper funtion: given one mutation string, return TRUE if it looks like a CEBPA bZIP in-frame event

is_cebpabzip_string <- function(s, bmin = BZIP_MIN, bmax = BZIP_MAX) {
  if (is.na(s) || !nzchar(s)) return(FALSE)
  if (!str_detect(s, "\\bCEBPA\\b")) return(FALSE)  
  p <- str_extract(s, "p\\.[^;\\)\\s]+") 
  if (is.na(p)) return(FALSE)
  if (str_detect(p, "fs|\\*|Ter|X")) return(FALSE) 
  pos <- suppressWarnings(as.integer(str_extract_all(p, "\\d+", simplify = TRUE))) 
  pos <- pos[!is.na(pos)]
  if (length(pos) == 0) return(FALSE)
  any(pos >= bmin & pos <= bmax) 
}

# Scan all columns with "mut" in their column names

mut_cols <- grep("mut", names(BEAT_meta_x2), value = TRUE)

BEAT_meta_x3 <- BEAT_meta_x2 %>% rowwise() %>% mutate(CEBPA_bzip = if_else(
      any(vapply(c_across(all_of(mut_cols)), is_cebpabzip_string, logical(1))), "mut", "wt")) %>% ungroup()

write.csv(BEAT_meta_x3, 'BEAT_meta_x3.csv')

# Compute ELN2022 for BEAT AML

#  helper functions (NA-safe) 

has_pat <- function(x, pattern) {
  if_else(
    is.na(x),
    NA,
    str_detect(x, regex(pattern, ignore_case = TRUE))
  )
}

# interpret mutation columns that are "mut"/"wt" (NA-safe)

is_mut <- function(x) {
  if_else(
    is.na(x),
    NA,
    as.character(x) %in% c("mut", "1", "TRUE", "T")
  )
}

mr_genes <- c("ASXL1","BCOR","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")

BEAT_meta_x3_eln <- BEAT_meta_x3 %>%
  mutate(
    cyto_text = str_squish(paste(consensusAMLFusions, otherCytogenetics, sep = " ; ")),
    
    NPM1_mut = is_mut(NPM1),
    TP53_mut = is_mut(TP53),
    CEBPA_bzip_mut = is_mut(CEBPA_bzip),
    
    FLT3_ITD_mut = case_when(
      is.na(`FLT3.ITD`) ~ NA,
      tolower(`FLT3.ITD`) == "positive" ~ TRUE,
      tolower(`FLT3.ITD`) == "negative" ~ FALSE,
      TRUE ~ NA
    ),
    
    MR_gene_mut = case_when(
      if_all(all_of(mr_genes), ~ is.na(.x)) ~ NA,
      if_any(all_of(mr_genes), ~ is_mut(.x)) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    t_8_21 = has_pat(cyto_text, "RUNX1\\s*[-:]{1,2}\\s*RUNX1T1|RUNX1T1\\s*[-:]{1,2}\\s*RUNX1|t\\(8;21\\)"),
    inv16  = has_pat(cyto_text, "CBFB\\s*[-:]{1,2}\\s*MYH11|MYH11\\s*[-:]{1,2}\\s*CBFB|inv\\(16\\)|t\\(16;16\\)"),
    t_9_11 = has_pat(cyto_text, "MLLT3\\s*[-:]{1,2}\\s*KMT2A|KMT2A\\s*[-:]{1,2}\\s*MLLT3|t\\(9;11\\)"),
    
    t_9_22 = has_pat(cyto_text, "BCR\\s*[-:]{1,2}\\s*ABL1|ABL1\\s*[-:]{1,2}\\s*BCR|t\\(9;22\\)"),
    t_6_9  = has_pat(cyto_text, "DEK\\s*[-:]{1,2}\\s*NUP214|NUP214\\s*[-:]{1,2}\\s*DEK|t\\(6;9\\)"),
    t_8_16 = has_pat(cyto_text, "KAT6A\\s*[-:]{1,2}\\s*CREBBP|CREBBP\\s*[-:]{1,2}\\s*KAT6A|t\\(8;16\\)"),
    
    inv3_t3_3   = has_pat(cyto_text, "inv\\(3\\)|t\\(3;3\\)"),
    mecom_rearr = has_pat(cyto_text, "MECOM|EVI1|GATA2\\s*[-:]{1,2}\\s*MECOM|RPN1\\s*[-:]{1,2}\\s*MECOM"),
    mec1om_rearr = if_else(
      is.na(inv3_t3_3) & is.na(mecom_rearr),
      NA,
      coalesce(inv3_t3_3, FALSE) | coalesce(mecom_rearr, FALSE)
    ),
    
    kmt2a_any   = has_pat(cyto_text, "KMT2A|\\bMLL\\b|KMT2A_re"),
    kmt2a_other = if_else(is.na(kmt2a_any), NA, kmt2a_any & !coalesce(t_9_11, FALSE)),
    
    minus5_del5q = has_pat(karyotype, "\\-5\\b|del\\(5q\\)|5q\\-"),
    minus7       = has_pat(karyotype, "\\-7\\b|del\\(7q\\)|7q\\-"),
    abn17p       = has_pat(karyotype, "abn\\(17p\\)|del\\(17p\\)|17p\\-|\\-17\\b"),
    
    monosomal_karyo = has_pat(karyotype, "monosomal|\\bmk\\b"),
    complex_karyo = case_when(
      is.na(karyotype) ~ NA,
      has_pat(karyotype, "complex") ~ TRUE,
      str_detect(str_to_upper(str_squish(karyotype)), "^(46,(XX|XY)(\\[\\d+\\])?|NK)$") ~ FALSE,
      TRUE ~ TRUE
    ),
    
    
    eln2022 = case_when(coalesce(t_9_22, FALSE) | coalesce(t_6_9, FALSE) | coalesce(kmt2a_other, FALSE) |
        coalesce(t_8_16, FALSE) | coalesce(inv3_t3_3, FALSE) | coalesce(mec1om_rearr, FALSE) |
        coalesce(minus5_del5q, FALSE) | coalesce(minus7, FALSE) | coalesce(abn17p, FALSE) |
        coalesce(complex_karyo, FALSE) | coalesce(monosomal_karyo, FALSE) |
        coalesce(MR_gene_mut, FALSE) | coalesce(TP53_mut, FALSE) ~ "Adverse",
      
      coalesce(t_8_21, FALSE) | coalesce(inv16, FALSE) | coalesce(CEBPA_bzip_mut, FALSE) ~ "Favorable",
      
      coalesce(NPM1_mut, FALSE) & is.na(FLT3_ITD_mut) ~ NA_character_,  
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == FALSE) ~ "Favorable",
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == TRUE)  ~ "Intermediate",
      
      TRUE ~ "Intermediate"
    ),
    
    eln2022 = factor(eln2022, levels = c("Favorable", "Intermediate", "Adverse"))) %>% select(-cyto_text)


# Load the Finnish cohort (keep initial registries for each patient)

FIN_meta_annot <- read.csv("File_0_Common_sample_annotation_252S.csv", sep = ",") %>%
  separate(Sample_ID, into = c("Sample_ID_1", "Sample_ID_2"), sep = "_(?=[^_]+$)", remove = FALSE) %>%
  filter(Sample_ID_2 == "01", RNA.seq == "available", !Sample_ID %in% c("AML_078_01", "AML_154_01"))


FIN_meta_clin <- read.csv("File_1.1_Clinical_summary_186_Patients.csv", header = TRUE, sep = ";") %>%
  mutate(ICD.O = na_if(ICD.O, ""), ICD.O = replace_na(ICD.O, "none")) %>% filter(!ICD.O %in% c(
      "Acute promyelocytic leuk.,t(15;17)(q22;q11-12)", "none")) %>% rename(patient = Patient_ID)

genes_keep <- c("FLT3","DNMT3A","NPM1","TP53","ASXL1","BCOR","CEBPA","EZH2","RUNX1",
                "SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")

FIN_meta_mut <- read.csv("File_6_Binary_mutation_225S_57G.csv", header = TRUE, sep = ";") %>%
  filter(Gene_Symbol %in% genes_keep) %>% select(Gene_Symbol, 2:227) %>%                 
  column_to_rownames("Gene_Symbol") %>% t() %>% as.data.frame() %>%
  rownames_to_column("patient") %>% filter(patient %in% FIN_meta_annot$Sample_ID) %>%
  mutate(patient = sub("^([^_]*_[^_]*).*", "\\1", patient))

FIN_meta <- merge(FIN_meta_clin, FIN_meta_mut, by = 'patient', all = T )
colnames(FIN_meta)

FIN_meta <- FIN_meta %>% mutate(cohort = 'FINN', CEBPA_bzip = NA, SF3B1 = NA, U2AF1 = NA, ZRSR2 = NA) %>%
  select(patient, cohort, Fusion_gene,karyotype_qp, complex_karyotype_qp,NPM1, FLT3_ITD, CEBPA_bzip,
         TP53, ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2)



# Calculate the ELN2022

#  NA-aware helper functions

has_pat <- function(x, pattern) {
  if_else(is.na(x), NA, str_detect(x, regex(pattern, ignore_case = TRUE)))
}

is_mut01 <- function(x) {
  case_when(
    is.na(x) ~ NA,
    x == 1   ~ TRUE,
    x == 0   ~ FALSE,
    TRUE     ~ NA
  )
}

mr_genes <- c("ASXL1","BCOR","EZH2","RUNX1","SF3B1","SRSF2","STAG2","U2AF1","ZRSR2")

FIN_meta_eln <- FIN_meta %>%
  mutate(
    
    NPM1_mut       = is_mut01(NPM1),
    FLT3_ITD_mut   = is_mut01(FLT3_ITD),
    CEBPA_bzip_mut = is_mut01(CEBPA_bzip),
    TP53_mut       = is_mut01(TP53),
    
    MR_gene_mut = case_when(
      if_all(all_of(mr_genes), ~ is.na(.x)) ~ NA,
      if_any(all_of(mr_genes), ~ is_mut01(.x)) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    
    minus5_del5q = has_pat(karyotype_qp, "\\-5\\b|del\\(5q\\)|5q\\-"),
    minus7       = has_pat(karyotype_qp, "\\-7\\b|del\\(7q\\)|7q\\-"),
    abn17p       = has_pat(karyotype_qp, "abn\\(17p\\)|del\\(17p\\)|17p\\-|\\-17\\b"),
    monosomal_karyo = has_pat(karyotype_qp, "monosomal|\\bmk\\b"),
    
    complex_karyo = case_when(
      !is.na(complex_karyotype_qp) ~ is_mut01(complex_karyotype_qp),
      is.na(karyotype_qp) ~ NA,
      has_pat(karyotype_qp, "complex") ~ TRUE,
      str_detect(str_to_upper(str_squish(karyotype_qp)),
                 "^(46,(XX|XY)(\\[(cp)?\\d+\\])?|NK)$") ~ FALSE,
      TRUE ~ TRUE
    ),
    
    
    t_8_21 = NA, inv16 = NA, t_9_11 = NA,
    t_9_22 = NA, t_6_9 = NA, t_8_16 = NA,
    inv3_t3_3 = NA, mec1om_rearr = NA, kmt2a_other = NA,
    
    eln2022 = case_when(coalesce(minus5_del5q, FALSE) | coalesce(minus7, FALSE) | coalesce(abn17p, FALSE) |
        coalesce(complex_karyo, FALSE) | coalesce(monosomal_karyo, FALSE) |
        coalesce(MR_gene_mut, FALSE) | coalesce(TP53_mut, FALSE) ~ "Adverse",
      
      coalesce(CEBPA_bzip_mut, FALSE) ~ "Favorable",
      
      coalesce(NPM1_mut, FALSE) & is.na(FLT3_ITD_mut) ~ NA_character_,   
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == FALSE) ~ "Favorable",
      coalesce(NPM1_mut, FALSE) & (FLT3_ITD_mut == TRUE)  ~ "Intermediate",
      
      TRUE ~ "Intermediate"), eln2022 = factor(eln2022, levels = c("Favorable", "Intermediate", "Adverse")))




# Merge all cohorts

CS_KAW_sub <- CS_KAW_meta_eln_1 %>% select(patient, cohort, eln2022)
BEAT_sub <-  BEAT_meta_x3_eln %>% select(patient,cohort, eln2022)
FINN_sub <- FIN_meta_eln %>% select(patient,cohort, eln2022)

ELN_all <- rbind(CS_KAW_sub, BEAT_sub, FINN_sub)

write.csv(ELN_all, 'ELN2022_all.csv')
