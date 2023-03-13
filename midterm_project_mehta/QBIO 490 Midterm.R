setwd("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee") 

if(!require(BiocManager)){
  install.packages("BiocManager")
}
library(BiocManager)

if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if(!require(maftools)) BiocManager::install("maftools")
library(maftools)

if(!require(ggplot2)) BiocManager::install("ggplot2")
library(ggplot2)

if(!require(SummarizedExperiment)) BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if(!require(survival)){
  install.packages("survival")
}
library(survival)

if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if(!require(DESeq)) BiocManager::install("DESeq2")
library(DESeq2)

if(!require(EnhancedVolcano)) BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#Query all data

#Clinical
clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
#GDCdownload(clin_query)
clinical <- GDCprepare_clinic(clin_query, clinical.info = "patient")

colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#MAF
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

#GDCdownload(maf_query)

maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

rna_query <- GDCquery(project ="TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

#GDCdownload(rna_query)

rna_se <- GDCprepare(rna_query)

#Create rna_clinical
na_mask <- !is.na(rna_se@colData$age_at_index)
rna_clinical <- as.data.frame(rna_se@colData[na_mask, ])
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F, T)
rna_clinical <- rna_clinical[, treatments_mask]
rna_clinical$age_category <- ifelse(rna_clinical$age_at_index <= 58, "young", "old")
def_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_clinical <- rna_clinical[def_mask, ]

#Create rna_genes
rna_genes <- as.data.frame(rna_se@rowRanges@elementMetadata)
rownames(rna_genes) <- rna_genes$gene_id

#Create rna_counts
rna_counts <- as.data.frame(rna_se@assays@data$unstranded[, na_mask])
rownames(rna_counts) <- rownames(rna_genes)
colnames(rna_counts) <- rownames(rna_clinical)
rna_counts <- rna_counts[ ,def_mask]

#Create outputs folder 
setwd("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/midterm_project_mehta")
dir.create("outputs")
setwd("outputs")

#Clinical data clean up

#Check for NAs in histological types
sum(is.na(clinical$histological_type))

#Remove empty values in histological types column
hist_type_mask <- ifelse(clinical$histological_type == "Infiltrating Ductal Carcinoma" | clinical$histological_type == "Infiltrating Lobular Carcinoma", T, F)
hist_type_cleaned_clinical <- clinical[hist_type_mask, ]
hist_type_cleaned_clinical$histological_type <- droplevels(hist_type_cleaned_clinical$histological_type)

#Create boxplot 
age_hist_type_boxplot <- boxplot(hist_type_cleaned_clinical$age_at_initial_pathologic_diagnosis~hist_type_cleaned_clinical$histological_type,
                                 xlab = "Histological Type", ylab = "Age at Diagnosis")

ggsave("age_hist_type_boxplot.png")

#making a survival time column for survival plots
hist_type_cleaned_clinical$survival_time <- ifelse(is.na(hist_type_cleaned_clinical$days_to_death),
                                                   hist_type_cleaned_clinical$survival_time <- hist_type_cleaned_clinical$days_to_last_followup,
                                                   hist_type_cleaned_clinical$survival_time <- hist_type_cleaned_clinical$days_to_death)

#remove any -Inf values in survival_time
inf_mask <- ifelse(hist_type_cleaned_clinical$survival_time == "-Inf", F, T)
hist_type_cleaned_clinical <- hist_type_cleaned_clinical[inf_mask, ]

#making a death event (T/F) column for survival plots
hist_type_cleaned_clinical$death_event <- ifelse(hist_type_cleaned_clinical$vital_status == "Alive", hist_type_cleaned_clinical$death_event <- FALSE,
                                           hist_type_cleaned_clinical$death_event <- TRUE)

#create survminer objects
survival_object <- Surv(time = hist_type_cleaned_clinical$survival_time, event = hist_type_cleaned_clinical$death_event)
fit_object <- survfit(survival_object ~ hist_type_cleaned_clinical$histological_type, data = hist_type_cleaned_clinical)

#create the plot
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot_hist_type <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

#show the plot
KM_plot_hist_type

ggsave("km_plot_hist_type.png")

#Create mask to identify patients within chemotherapy group
duct_car_mask <- ifelse(hist_type_cleaned_clinical$histological_type == 'Infiltrating Ductal Carcinoma', T, F)

#Apply chemotherapy mask to column of patient barcodes
duct_car_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[duct_car_mask]

#Subset maf_object to only contain chemotherapy group
duct_car_maf <- subsetMaf(maf = maf_object,
                              tsb = duct_car_patient_barcodes)

#Create mask to identify patients within hormone therapy group
duct_lob_mask <- ifelse(hist_type_cleaned_clinical$histological_type == 'Infiltrating Lobular Carcinoma', T, F)

#Apply hormone therapy mask to column of patient barcodes
duct_lob_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[duct_lob_mask]

#Subset maf_object to only contain hormone therapy group
duct_lob_maf <- subsetMaf(maf = maf_object,
                         tsb = duct_lob_patient_barcodes)
#Create co-oncoplot
coOncoplot(m1 = duct_car_maf, 
           m2 = duct_lob_maf, 
           m1Name = 'Infiltrating Ductal Carcinoma', 
           m2Name = 'Infiltrating Lobular Carcinoma',
           borderCol = NA)

ggsave("hist_type_coOncoplot.png")

#Create co-lollipop plot
lollipopPlot2(m1 = duct_car_maf, 
              m2 = duct_lob_maf, 
              m1_name = 'Infiltrating Ductal Carcinoma',
              m2_name = 'Infiltrating Lobular Carcinoma',
              gene = "TP53") ## pick any gene of your choosing to fill in here

ggsave("hist_type_co_lollipop_plot.png", plot = last_plot())


rna_clinical$primary_diagnosis <- as.factor(rna_clinical$primary_diagnosis)
summary(rna_clinical$primary_diagnosis)

prim_diag_mask <- ifelse(rna_clinical$primary_diagnosis == "Infiltrating duct carcinoma, NOS" | rna_clinical$primary_diagnosis == "Lobular carcinoma, NOS", T, F)
prim_diag_cleaned_rna_clinical <- rna_clinical[prim_diag_mask, ]
prim_diag_cleaned_rna_counts <- rna_counts[, prim_diag_mask]

prim_diag_cleaned_rna_clinical$primary_diagnosis <- droplevels(prim_diag_cleaned_rna_clinical$primary_diagnosis)
prim_diag_cleaned_rna_clinical$primary_diagnosis <- factor(prim_diag_cleaned_rna_clinical$primary_diagnosis)

prim_diag_cleaned_rna_clinical$race <- factor(prim_diag_cleaned_rna_clinical$race)
head(prim_diag_cleaned_rna_clinical$race)

sum(is.na(prim_diag_cleaned_rna_clinical$age_category))
sum(is.na(prim_diag_cleaned_rna_clinical$ajcc_pathologic_stage))
sum(is.na(prim_diag_cleaned_rna_clinical$race))

na_mask <- !is.na(prim_diag_cleaned_rna_clinical$ajcc_pathologic_stage)

prim_diag_cleaned_rna_clinical <-  prim_diag_cleaned_rna_clinical[na_mask,]

prim_diag_cleaned_rna_counts <-  prim_diag_cleaned_rna_counts[, na_mask]

row_sums <- rowSums(prim_diag_cleaned_rna_counts)

# create a boolean mask where genes with < 10 total counts are FALSE, and genes with >= 10 total counts are TRUE
low_counts_mask <- ifelse(row_sums >= 10, T, F)

# rewrite the rna_counts df, subsetting for only genes with >= 10 total counts
prim_diag_cleaned_rna_counts <- prim_diag_cleaned_rna_counts[low_counts_mask,]

#update rna_genes with the low_counts_mas
rna_genes <- rna_genes[low_counts_mask,]

colnames(prim_diag_cleaned_rna_counts) <- rownames(prim_diag_cleaned_rna_clinical)

dds <- DESeqDataSetFromMatrix(countData = prim_diag_cleaned_rna_counts,
                              colData = prim_diag_cleaned_rna_clinical,
                              design = ~primary_diagnosis)

dds_obj <- DESeq(dds) # note: this will likely take a long time (ie 45 minutes to 2 hours)

resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results <- results(dds_obj, format = "DataFrame", contrast = c("primary_diagnosis", "Infiltrating duct carcinoma, NOS", "Lobular carcinoma, NOS"))

results <- data.frame(rna_genes$gene_name, results@rownames, results@listData$log2FoldChange, results@listData$pvalue, results@listData$padj, -log10(results@listData$padj))
colnames(results) <- c("Gene Name", "Gene ID", "Log2 Fold Change", "P-value", "PAdj", "-log10(PAdj)")

sig_results <- results[results$PAdj < 0.05, ]
up_reg_results <- sig_results[order(sig_results$`Log2 Fold Change`, decreasing = TRUE), ]
up_reg_results <- up_reg_results[up_reg_results$`Log2 Fold Change` > 1, ]

down_reg_results <- sig_results[order(sig_results$`Log2 Fold Change`, decreasing = FALSE), ]
down_reg_results <- down_reg_results[down_reg_results$`Log2 Fold Change` < -1, ]

par(mar=c(1,1,1,1))
EnhancedVolcano(results, 
                lab = results$`Gene Name`, 
                x = "Log2 Fold Change", 
                y = "PAdj")

ggsave("hist_type_volcano_plot.png")

