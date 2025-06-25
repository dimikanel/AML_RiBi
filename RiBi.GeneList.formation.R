setwd('C:/Users/Dimitrios Kanellis/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')
#setwd('C:/Users/Sofia/Dropbox/Backup/Data_ppt/AMLSM/Latest/GitHub_Ready')

library(tidyverse)
library(biomaRt)


GSEA_df <- read.csv('RiBi.genes.GSEA.GO.0042254.csv') %>% distinct(Gene)
AmiGO_DF <- read.csv('RiBi.AmiGO.csv', sep = ';')

RiBi_df <- full_join(GSEA_df, AmiGO_DF, by = 'Gene') %>% dplyr::select(1) %>% dplyr::rename(external_gene_name = Gene)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",  version = 113)

genes <- RiBi_df$external_gene_name
gene_list <- getBM(filters = "external_gene_name",
                   attributes = c("ensembl_gene_id",
                                 "external_gene_name",'gene_biotype'),
                   values = genes,
                   mart = ensembl)

RiBi_df <- merge(RiBi_df, gene_list, by = "external_gene_name")
write.csv(RiBi_df, 'RiBi.only.terms.GSEA.AmiGO.csv')
