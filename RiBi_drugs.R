setwd('C:/Users/Dimitrios Kanellis/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')
#setwd('C:/Users/Sofia/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')


require(tidyverse)
require(matrixStats)
library(circlize)
library(ComplexHeatmap)

# Extract drug screen data from: https://aacrjournals.org/cancerdiscovery/article/12/2/388/678465/Implementing-a-Functional-Precision-Medicine-Tumor
# Load the FINN drug data

finn_drugs <- readxl::read_excel('cd-21-0410_supplementary_tables_suppst1.xlsx', sheet = 7, skip = 2)  

finn_drugs_ctr <- finn_drugs %>% filter(str_detect(Sample_ID,'Healthy')) %>% dplyr::select(1,3,18) %>%
                          pivot_wider(names_from = Chemical_compound, values_from = DSS) %>% column_to_rownames('Sample_ID')


finn_drugs_ctr_sum <- as.data.frame(rbind(colMeans(as.matrix(finn_drugs_ctr), na.rm = T),
                                                        colSds(as.matrix(finn_drugs_ctr), na.rm = T),
                                                        colMedians(as.matrix(finn_drugs_ctr), na.rm = T), 
                                                        colMads(as.matrix(finn_drugs_ctr), na.rm = T)))

rownames(finn_drugs_ctr_sum) <- c('mean', 'sd', 'median', 'mad')


finn_drugs_pat <- finn_drugs %>% filter(!str_detect(Sample_ID,'Healthy')) %>% dplyr::select(1,3,18) %>% 
   pivot_wider(names_from = Chemical_compound, values_from = DSS) %>% column_to_rownames('Sample_ID')


##Check DSS distribution between the healthy samples and the AML patients for Actinomycin D as an example

finn_drugs_pat_ActD <- data.frame(DSS = c(finn_drugs_ctr$Dactinomycin,finn_drugs_pat$Dactinomycin), 
                            status = rep(c("Healthy", "AML"), c(17,164)))

ggplot(finn_drugs_pat_ActD, aes(x = DSS, fill = status)) + geom_density(alpha = 0.5) + ggtitle('Dactinomycin')



# Plot the data in boxplot and compare DSS control Vs patient for Actinomycin D

dss_values <- c(finn_drugs_ctr$Dactinomycin, finn_drugs_pat$Dactinomycin)
groups <- c(rep("Control", nrow(finn_drugs_ctr)), rep("Patient", nrow(finn_drugs_pat)))

# Boxplot
boxplot(dss_values ~ groups,
        main = "Boxplot of ActD DSS: Control vs Patient",
        xlab = "Group",
        ylab = "DSS Value",
        col = c("lightblue", "lightgreen"))

# t-test
t_test_result <- t.test(finn_drugs_ctr$Dactinomycin, finn_drugs_pat$Dactinomycin)
print(t_test_result)



# Follow formulas of sDSS, zDDS and rDSS calculation from: https://www.nature.com/articles/s41596-023-00903-x
# Build a function to compute sDSS, zDSS, rDSS. 
# Check normality of DSS value distribution for each drug among healthy subjects
# If normal distribution: Z test; If non-normal: permutation test to assess the statisrical significance of a drug-patient combination (actionable drugs)



compute_dss_scores <- function(patient_dss, control_dss,
                               infer = TRUE,
                               fdr_threshold = 0.05,
                               n_perm = 10000,
                               normality_p_threshold = 0.05) {
  # Ensure column alignment
  
  stopifnot(all(colnames(patient_dss) == colnames(control_dss)))
  
  # Precompute control statistics
  
  mean_c <- colMeans(control_dss, na.rm = TRUE)
  sd_c <- apply(control_dss, 2, sd, na.rm = TRUE)
  median_c <- apply(control_dss, 2, median, na.rm = TRUE)
  mad_c <- apply(control_dss, 2, mad, constant = 1, na.rm = TRUE)
  
  # Avoid division by 0
  
  sd_c[is.na(sd_c) | sd_c == 0] <- 1e-6
  mad_c[is.na(mad_c) | mad_c == 0] <- 1e-6
  
  # Create matrices for broadcasting
  
  mean_mat <- matrix(rep(mean_c, each = nrow(patient_dss)), nrow = nrow(patient_dss))
  sd_mat <- matrix(rep(sd_c, each = nrow(patient_dss)), nrow = nrow(patient_dss))
  median_mat <- matrix(rep(median_c, each = nrow(patient_dss)), nrow = nrow(patient_dss))
  mad_mat <- matrix(rep(mad_c, each = nrow(patient_dss)), nrow = nrow(patient_dss))
  
  # Compute DSS scores
  
  sDSS <- patient_dss - mean_mat
  zDSS <- sDSS / (sd_mat + 1)
  rDSS <- (patient_dss - median_mat) / (mad_mat + 1)
  
  # Helper to convert matrices to long format
  
  make_long <- function(mat, label) {
    mat %>%
      as.data.frame() %>%
      rownames_to_column("Patient") %>%
      pivot_longer(cols = -Patient, names_to = "Drug", values_to = label)
  }
  
  # Create long format versions
  
  df_s <- make_long(sDSS, "sDSS")
  df_z <- make_long(zDSS, "zDSS")
  df_r <- make_long(rDSS, "rDSS")
  
  dss_long <- purrr::reduce(list(df_s, df_z, df_r), dplyr::full_join, by = c("Patient", "Drug"))
  
  # Return early if inference is not requested
  
  if (!infer) return(dss_long)
  
  # Begin inference
  
  infer_results <- list()
  
  for (drug in colnames(patient_dss)) {
    dss_pat <- patient_dss[, drug]
    dss_ctrl <- control_dss[, drug]
    
    if (all(is.na(dss_ctrl)) || all(is.na(dss_pat))) next
    
    # Shapiro-Wilk test for normality of controls
    
    is_normal <- FALSE
    if (length(unique(na.omit(dss_ctrl))) > 1) {
      norm_test <- tryCatch(shapiro.test(dss_ctrl), error = function(e) NULL)
      is_normal <- !is.null(norm_test) && norm_test$p.value > normality_p_threshold
    }
    
    for (i in seq_along(dss_pat)) {
      patient <- rownames(patient_dss)[i]
      z_val <- zDSS[i, drug]
      
      if (is.na(z_val)) {
        infer_results[[length(infer_results) + 1]] <- data.frame(
          Patient = patient, Drug = drug, p_value = NA, Test_Used = NA
        )
        next
      }
      
      if (is_normal) {
        # z-test
        
        p_val <- 1 - pnorm(z_val)
        test_used <- "z-test"
        
      } else {
        # Permutation test
        
        combined <- c(dss_pat[i], dss_ctrl)
        perm_z <- numeric(n_perm)
        
        for (k in seq_len(n_perm)) {
          perm <- sample(combined)
          perm_patient <- perm[1]
          perm_controls <- perm[-1]
          perm_mean <- mean(perm_controls, na.rm = TRUE)
          perm_sd <- sd(perm_controls, na.rm = TRUE)
          if (is.na(perm_sd) || perm_sd == 0) perm_sd <- 1e-6
          perm_z[k] <- (perm_patient - perm_mean) / (perm_sd + 1)
        }
        
        p_val <- (sum(perm_z >= z_val, na.rm = TRUE) + 1) / (n_perm + 1)
        test_used <- "permutation"
      }
      
      infer_results[[length(infer_results) + 1]] <- data.frame(
        Patient = patient, Drug = drug, p_value = p_val, Test_Used = test_used
      )
    }
  }
  
  # Compile results and adjust p-values
  
  inference_df <- bind_rows(infer_results)
  inference_df$FDR <- p.adjust(inference_df$p_value, method = "BH")
  inference_df$Significant <- inference_df$FDR < fdr_threshold
  
  # Final merge and return
  
  final_df <- dss_long %>%
    left_join(inference_df, by = c("Patient", "Drug"))
  
  return(final_df)
}



finn_drugs_pat <- finn_drugs_pat[, sort(names(finn_drugs_pat))] # alphabetical re-order
finn_drugs_ctr <- finn_drugs_ctr[, sort(names(finn_drugs_ctr))]


results <- compute_dss_scores(
  patient_dss = finn_drugs_pat,
  control_dss = finn_drugs_ctr, 
  fdr_threshold = 0.1,
  infer = TRUE,       # Set to FALSE for speed if you don’t need p-values
  n_perm = 10000, normality_p_threshold = 0.05     # Reduce for quick testing
)


results <- results %>% na.omit(Significant)
write.csv(results, 'sDSS.zDSS.rDSS.pval.csv')


results_sign <- results %>% filter(Significant == TRUE) 
write.csv(results_sign, 'sDSS.zDSS.rDSS.pval.sign.csv')

results <- read.csv('sDSS.zDSS.rDSS.pval.csv') %>% dplyr::select(2:10)
results_sign <- read.csv('sDSS.zDSS.rDSS.pval.sign.csv') %>% dplyr::select(2:10)


## Desing a rDSS distribution plot for all drugs


# Compute the threshold

threshold <- min(results_sign$rDSS, na.rm = TRUE)

# Compute the density manually

density_data <- density(results$rDSS, na.rm = TRUE)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Flag points beyond the threshold

density_df <- density_df |>
  mutate(highlight = x >= threshold)

# Plot

ggplot(density_df, aes(x = x, y = y)) +
  geom_area(aes(fill = highlight), alpha = 0.6, color = "black") +
  #geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "skyblue")) +
  xlab("rDSS") +
  ylab("Density") +
  xlim(-20, 20) +
  labs(fill = paste0("FDR <= 0.1")) +
  theme_minimal()



# Plot rDSS distribution for ActD

results_ActD <- results %>% filter(Drug == 'Dactinomycin') 
results_ActD_sign <- results_ActD %>% filter(Significant == TRUE)

ggplot(results_ActD, aes(x=rDSS)) +
  geom_histogram(binwidth=.5, colour="black", fill="white") + 
  geom_vline(aes(xintercept = min(results_ActD_sign$rDSS)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1) + xlab('rdss') + xlim(-20, 20) + theme_minimal()


# Draw a heatmap of rDSS values for the top 40 most variable drugs


results_hm <- results %>% dplyr::select(1,2,5) %>% 
  pivot_wider(names_from = Drug,values_from = rDSS) %>% column_to_rownames('Patient') %>% as.data.frame() %>% t()
 
# Calculate the variance of each row (drug)

variances <- apply(results_hm, 1, var)

# Find the 50 most variable rows

top_40_indices <- order(variances, decreasing = TRUE)[1:40]

# Subset the dataframe to include only the top 50 most variable rows

top_40_data <- results_hm[top_40_indices, ]


heatmap <- Heatmap(as.matrix(top_40_data), 
                   name = "expression", 
                   cluster_rows = TRUE, 
                   cluster_columns = TRUE,
                   show_row_names = TRUE, 
                   show_column_names = TRUE,
                   col = colorRamp2(c(min(as.matrix(top_40_data)), mean(as.matrix(top_40_data)), max(as.matrix(top_40_data))), c("navy", "white", "firebrick3")),
                   row_dend_reorder = TRUE, 
                   column_dend_reorder = TRUE,
                   heatmap_legend_param = list(title = "rDSS"),  row_names_gp = gpar(fontsize = 12), column_names_gp = gpar(fontsize = 4))

draw(heatmap)


# Add variable of culturing media used in the drug screen

drugs_media <- readxl::read_excel('cd-21-0410_supplementary_tables_suppst1.xlsx', sheet = 5, skip = 2) %>% 
             dplyr::select(2,3) %>% dplyr::rename(Patient = Sample_ID)

results_sign <- merge(results_sign, drugs_media, by = 'Patient', all.x = T)
head(results_sign,3)

# Add the RiBi status variable
# Keep the RiBi score derived from our analysis, subset the FINN data and re-calculate quantiles within the subset cohort 
# to tailor the "high," "moderate," and "low" definitions to the actual distribution of scores in the FINN cohort


finn_meta <- read.csv('RiBi_score_all_samples.csv') 
finn_meta <- finn_meta %>% filter(cohort == 'FINN') %>% select(2:16) %>%  
  mutate(RiBi_status = case_when(RiBi_score > quantile(finn_meta$RiBi_score, 0.75) ~ "high",
                                 RiBi_score > quantile(finn_meta$RiBi_score, 0.25) &
                                   RiBi_score < quantile(finn_meta$RiBi_score, 0.75) ~ "moderate",
                                 RiBi_score < quantile(finn_meta$RiBi_score, 0.25) ~ "low")) %>% dplyr::rename(Patient = patient)

write.csv(finn_meta, 'finn_meta.csv')


results_sign <- merge(results_sign, finn_meta, by = 'Patient')


# Add the drug class variable

drugs_class <- readxl::read_excel('cd-21-0410_supplementary_tables_suppst1.xlsx', sheet = 6, skip = 2) 
colnames(drugs_class)[2] <- 'Drug'
drugs_class <- drugs_class %>% dplyr::select(2,8,9) %>% as.data.frame()

results_sign <- merge(results_sign, as.data.frame(drugs_class), by = 'Drug') 


write.csv(results_sign, 'results_sign_all_meta.csv')

# Counting incidences per various groups

table(results_sign$Drug_Class, results_sign$RiBi_status) # Drug classes per RiBi group
table(results_sign$Drug_Class, results_sign$RiBi_status, results_sign$Medium) # Drug classes per RiBi group per culturing media

write.csv(table(results_sign$Drug_Class, results_sign$RiBi_status, results_sign$Medium), 'Drug_class_RiBi_status_Medium.csv')

# Focus on specific drug classes 

results_sign_kin_inh <- results_sign %>% filter(Drug_Class == 'Kinase inhibitor') # kinase inh per RiBi group
table(results_sign_kin_inh$Drug_Sub_Class, results_sign_kin_inh$RiBi_status)

results_sign_chemo <- results_sign %>% filter(Drug_Class == 'Chemotherapeutics') # Chemotherapeutics per RiBi group
table(results_sign_chemo$Drug_Sub_Class, results_sign_chemo$RiBi_status)

results_sign_diff <- results_sign %>% filter(Drug_Class == 'Differentiating/epigenetic modifier') # Differentiation inducer per RiBi group
table(results_sign_diff$Drug_Sub_Class, results_sign_diff$RiBi_status)

write.csv(table(results_sign_kin_inh$Drug_Sub_Class, results_sign_kin_inh$RiBi_status), 'Kin_inh_RiBi_status.csv')
write.csv(table(results_sign_chemo$Drug_Sub_Class, results_sign_chemo$RiBi_status), 'Chemotherapeutics_RiBi_status.csv')
write.csv(table(results_sign_diff$Drug_Sub_Class, results_sign_diff$RiBi_status), 'Diff_ind_RiBi_status.csv')


# Focus on RiBi inhibitors 


RiBi_inh <- read.csv('RiBi_inh_Curated.csv') %>% dplyr::rename(Drug = RiBi_inh)
results_sign_RiBi_inh <- merge(results_sign, RiBi_inh, by = 'Drug')

write.csv(results_sign_RiBi_inh,'results_sign_RiBi_inh.csv')

table(results_sign_RiBi_inh$Drug, results_sign_RiBi_inh$RiBi_status) # RiBi inhibitors per RiBi group

write.csv(table(results_sign_RiBi_inh$Drug, results_sign_RiBi_inh$RiBi_status), 'RiBi_inhibitors_RiBi.status.csv')


# Merge drugs with longitudinal RiBi score data

# Find samples with longitudinal data

results_sign_longi <- results_sign %>% filter(!is.na(stage)) %>% distinct(Patient, .keep_all = T)
results_sign_longi$prefix7 <- substr(results_sign_longi$Patient, 1, 7)
samples_with_shared_prefix <- results_sign_longi %>% group_by(prefix7) %>%
  filter(n() > 1) %>% ungroup() %>% arrange(Patient)
table(samples_with_shared_prefix$Patient)

results_sign_longi <- results_sign %>% filter(Patient %in% samples_with_shared_prefix$Patient)

write.csv(results_sign_longi, 'results_sign_longi.csv')


# count incidences

# Drug classes per stage

table(results_sign_longi$stage, results_sign_longi$Drug_Class)
table(results_sign_longi$stage, results_sign_longi$RiBi_status)
table(results_sign_longi$stage, results_sign_longi$Drug_Class, results_sign_longi$RiBi_status)
write.csv(table(results_sign_longi$stage, results_sign_longi$Drug_Class), 'longitudinal_drug_classes.csv')
write.csv(table(results_sign_longi$stage, results_sign_longi$Drug_Class, results_sign_longi$RiBi_status), 'longitudinal_drug_classes_RiBi_status.csv')


# Drug classes per stage for patient #AML_157 as en example

results_sign_longi_AML157 <- results_sign_longi %>% filter(str_detect(Patient, 'AML_157'))
write.csv(results_sign_longi_AML157, 'AML_157.csv')
table(results_sign_longi_AML157$stage, results_sign_longi_AML157$Drug_Class) 
write.csv(table(results_sign_longi_AML157$stage, results_sign_longi_AML157$Drug_Class), 'longitudinal_AML_157.csv')

# Chemotherapeutics per stage for patient #AML_157

results_sign_longi_AML157_chemo <- results_sign_longi_AML157 %>% filter(Drug_Class == 'Chemotherapeutics')
table(results_sign_longi_AML157_chemo$stage, results_sign_longi_AML157_chemo$Drug)
write.csv(table(results_sign_longi_AML157_chemo$stage, results_sign_longi_AML157_chemo$Drug), 'longitudinal_chemo_AML_157.csv')


