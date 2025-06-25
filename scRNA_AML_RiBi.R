#Load relevant packages

library(tidyverse)
library(Seurat)
library(scDblFinder)
library(RColorBrewer)
library(patchwork)
library(nortest)
library(Hmisc)
library(nortest)



# The scRNA data processing was done according to Lasry & colleagues (PMID: 36658429)
# Download the data from https://singlecell.broadinstitute.org/single_cell/study/SCP1987/an-inflammatory-state-remodels-the-immune-microenvironment-and-improves-risk-stratification-in-acute-myeloid-leukemia)

expression_data <- ReadMtx(mtx = 'RNA_soupX1.mtx.gz', 
                           features = 'features_RNA_soupX1.csv', cells = 'cells_RNA_soupX1.csv',feature.column = 1)


meta_data <- read.csv('metadata_clustering_w_header_upd.csv', row.names = 1)


seurat_obj <- CreateSeuratObject(counts = expression_data,
 project = 'Lasry.et.al', min.cells = 3, min.features = 100,  meta.data = meta_data)


# Normalize RNA data

seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)

# Identify variable features

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)


# Scale the data

seurat_obj <- ScaleData(seurat_obj)

# Perform PCA for dimensionality reduction
set.seed(123)


seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30)

# Step 7: UMAP and Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, k.param = 25)
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, n.neighbors = 25, min.dist = 0.3)

# Plot
DimPlot(seurat_obj, reduction = 'umap', label = TRUE) +
  ggtitle(paste0('UMAP - Resolution 2 - ', length(unique(Idents(seurat_obj))), ' clusters'))

DimPlot(seurat_obj, reduction = 'umap', label = F, group.by = 'ap_aml_age')

# Plot 10.000 randomly selected cells 

set.seed(123) 
cells_to_keep <- sample(Cells(seurat_obj), size = 10000)

# Subset the Seurat object

seurat_obj_subset <- subset(seurat_obj, cells = cells_to_keep)

# Plot the DimPlot

DimPlot(object = seurat_obj_subset, reduction = 'umap', group.by = 'ap_aml_age')


seurat_obj <-  readRDS('seurat_obj.rds')
# Calculate RiBi score for each cell (using our gene signature)

ribi_df <- read.csv('RiBi_model_coefficients_translated.csv', header = T, sep = ';') %>% dplyr::rename(gene_name = symbol)

# Penalise RiBi genes not expressed in cells 

gene_weights <- setNames(ribi_df$coefficient, ribi_df$gene_name)
present <- intersect(names(gene_weights), rownames(seurat_obj))
gene_weights <- gene_weights[present]

expr_matrix <- GetAssayData(seurat_obj, slot = 'data')[present, , drop = FALSE]
expressed <- expr_matrix > 0
penalty_value <- -1
score_matrix  <- expr_matrix * gene_weights + (1 - expressed) * (penalty_value * gene_weights)
penalized_score <- Matrix::colSums(score_matrix)
seurat_obj$RiBi_score <- penalized_score


# Ascribe RiBi status


low_threshold_RiBi <- quantile(seurat_obj@meta.data$RiBi_score, 0.25)
high_threshold_RiBi <- quantile(seurat_obj@meta.data$RiBi_score, 0.75)

seurat_obj$RiBi_status <- with(seurat_obj@meta.data, ifelse(RiBi_score <= low_threshold_RiBi, 'low',
         ifelse(RiBi_score >= high_threshold_RiBi, 'high', 'moderate')))



# Plot RiBi status

RiBi_palette <- c('low' = 'grey', 'moderate' = 'orange', 'high' = 'darkred')

# Reduce cell numbers

set.seed(123)  
cells_to_keep <- sample(Cells(seurat_obj), size = 10000)

# Subset the Seurat object

seurat_obj_subset <- subset(seurat_obj, cells = cells_to_keep)

# Plot the DimPlot

DimPlot(object = subset(seurat_obj_subset, ap_aml_age %in% c('adult_AML','control_39to53yrs')), 
        reduction = 'umap', group.by = 'RiBi_status', cols = RiBi_palette, , split.by = 'ap_aml_age')



# Calculate points per group

seurat_obj_sub <- subset(seurat_obj, subset = ap_aml_age %in% c('adult_AML','control_39to53yrs'))
table(seurat_obj_sub$RiBi_status,seurat_obj_sub$ap_aml_age)
write.csv(table(seurat_obj_sub$RiBi_status,seurat_obj_sub$ap_aml_age),'points_per_group.RiBi.AML.Healthy.Adults.csv')



# Plot RiBi status Vs cell types in adult AML


# Create the first DimPlot

cell_palette <- colorRampPalette(brewer.pal(9, 'Set3'))(30)


plot1a <- DimPlot(subset(seurat_obj_subset, ap_aml_age == 'adult_AML'), group.by = 'Broad_cell_identity',  
                  cols = cell_palette, label = T) +
  ggtitle('Plot colored by Variable 1')


# Create the second DimPlot

plot1b <- DimPlot(subset(seurat_obj_subset, ap_aml_age == 'adult_AML'), group.by = 'RiBi_status' , cols = RiBi_palette, label = F) +
  ggtitle('Plot colored by Variable 1')

# Arrange the plots side by side

combined_plot.1 <- plot1a | plot1b
print(combined_plot.1)


# Calculate the points per group

seurat_obj_sub <- subset(seurat_obj,  ap_aml_age == 'adult_AML')
table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status)
write.csv(table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status),'points_per_group.cell.types.adult_AML.csv')

seurat_obj_sub <- subset(seurat_obj,  ap_aml_age == 'control_39to53yrs')
table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status)
write.csv(table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status),'points_per_group.cell.types.adult.control.csv')


# Plot RiBi status in malignant Vs microenvironment cells

mal_palette <- c('malignant' = 'purple', 'microenvironment' = 'deepskyblue')

plot2a <- DimPlot(subset(seurat_obj_subset, ap_aml_age == 'adult_AML'), group.by = 'malignant',  cols = mal_palette, label = F) #+

plot2b <- DimPlot(subset(seurat_obj_subset, ap_aml_age == 'adult_AML'), group.by = 'RiBi_status', , cols = RiBi_palette, label = F) #+

combined_plot.2 <- plot2a | plot2b
print(combined_plot.2)


# Calculate the points per group

seurat_obj_sub <- subset(seurat_obj, ap_aml_age == 'adult_AML')
table(seurat_obj_sub$malignant,seurat_obj_sub$RiBi_status)
write.csv(table(seurat_obj_sub$malignant,seurat_obj_sub$RiBi_status),'points_per_group.RiBi.malignant.Vs.microenvironment.csv')


# calculate points per group

seurat_obj_sub <- subset(seurat_obj, ap_aml_age == 'adult_AML' & malignant == 'malignant')
table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status)
write.csv(table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status),'points_per_group.RiBi.malignant.csv')

seurat_obj_sub <- subset(seurat_obj, ap_aml_age == 'adult_AML' & malignant == 'microenvironment')
table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status)
write.csv(table(seurat_obj_sub$Broad_cell_identity,seurat_obj_sub$RiBi_status),'points_per_group.RiBi.microenvironment.csv')



# Calculate LSC17 score for each cell (procedure similar to the RiBi score calculation)
# Data acquired by K Ng & colleagues, 2016 (PMID: 27926740)

LSC17.df <- read.csv('LSC17.csv', header = T, sep = ',')


gene_weights_LSC17 <- setNames(LSC17.df$Weight, LSC17.df$gene_name)

present_LSC17 <- intersect(names(gene_weights_LSC17), rownames(seurat_obj))
gene_weights_LSC17 <- gene_weights_LSC17[present_LSC17]

expr_matrix_LSC17 <- GetAssayData(seurat_obj, slot = 'data')[present_LSC17, , drop = FALSE]
expressed_LSC17 <- expr_matrix_LSC17 > 0
penalty_value_LSC17 <- -1
score_matrix_LSC17  <- expr_matrix_LSC17 * gene_weights_LSC17 + (1 - expressed_LSC17) * (penalty_value_LSC17 * gene_weights_LSC17)
penalized_score_LSC17 <- Matrix::colSums(score_matrix_LSC17)
seurat_obj$LSC17_score <- penalized_score_LSC17

low_threshold_LSC17  <- quantile(seurat_obj@meta.data$LSC17_score, 0.2)
low_threshold_LSC17
high_threshold_LSC17  <- quantile(seurat_obj@meta.data$LSC17_score, 0.8)
high_threshold_LSC17

seurat_obj$LSC17_status <- with(seurat_obj@meta.data, ifelse(LSC17_score <= low_threshold_LSC17, 'low',
                                                            ifelse(LSC17_score >= high_threshold_LSC17, 'high', 'moderate')))



# Plot LSC17 status

LSC17_palette <- c('low' = 'cyan', 'moderate' = 'lightgreen', 'high' = 'gold4')


set.seed(123)  
cells_to_keep <- sample(Cells(seurat_obj), size = 10000)

# Subset the Seurat object

seurat_obj_subset <- subset(seurat_obj, cells = cells_to_keep)

# Plot the DimPlot

DimPlot(subset(seurat_obj_subset, ap_aml_age %in% c('adult_AML','control_39to53yrs')), 
        reduction = 'umap', group.by = 'LSC17_status', cols = LSC17_palette, , split.by = 'ap_aml_age', raster = F)


# Calculate points per group

seurat_obj_sub <- subset(seurat_obj, ap_aml_age %in% c('adult_AML','control_39to53yrs'))
table(seurat_obj_sub$LSC17_status,seurat_obj_sub$ap_aml_age)
write.csv(table(seurat_obj_sub$LSC17_status,seurat_obj_sub$ap_aml_age),'points_per_group.LSC17.adults.AML.Healthy.csv')


# Plot RiBi_status Vs LSC17_status

set.seed(123)  
cells_to_keep <- sample(Cells(seurat_obj), size = 10000)
seurat_obj_subset <- subset(seurat_obj, cells = cells_to_keep)

plot3a <- DimPlot(subset(seurat_obj_subset,  ap_aml_age == 'adult_AML'), group.by = 'RiBi_status',  cols = RiBi_palette, label = F) +
  ggtitle('Plot colored by RiBi.status')

# Create the second DimPlot

plot3b <- DimPlot(subset(seurat_obj_subset,  ap_aml_age == 'adult_AML'), group.by = 'LSC17_status',  cols = LSC17_palette, label = F) +
  ggtitle('Plot colored by LSC17.status')


combined_plot.3 <- plot3a | plot3b
print(combined_plot.3)


# Calculate points per group

seurat_obj_sub <- subset(seurat_obj, ap_aml_age == 'adult_AML')
table(seurat_obj_sub$RiBi_status, seurat_obj_sub$LSC17_status)
write.csv(table(seurat_obj_sub$RiBi_status, seurat_obj_sub$LSC17_status),'points_per_group.LSC17.RiBi.Adult.AML.csv')


seurat_obj_meta_df <- as.data.frame(seurat_obj@meta.data)
write.csv(seurat_obj_meta_df,'seurat_obj_meta_df.csv')
write.csv(subset(seurat_obj_meta_df,RiBi_status == 'low'),'seurat_obj_meta_df_LR.csv')
write.csv(subset(seurat_obj_meta_df,RiBi_status == 'moderate'),'seurat_obj_meta_df_MR.csv')
write.csv(subset(seurat_obj_meta_df,RiBi_status == 'high'),'seurat_obj_meta_df_HR.csv')


# Plot LSC17 values in different RiBi groups

summary_data <- seurat_obj_meta_df %>%
  group_by(RiBi_status) %>% dplyr::summarize(Mean = mean(LSC17_score, na.rm = TRUE), SD = sd(LSC17_score, na.rm = TRUE))


ggplot(summary_data, aes(x = RiBi_status, y = Mean, fill = RiBi_status)) +
  geom_bar(stat = 'identity', color = 'black', width = 0.6) +
  labs(title = 'Barplot with Standard Deviation no Error Bars',
       x = 'RiBi_status', y = 'Mean LSC17 score') +
  theme_minimal()

