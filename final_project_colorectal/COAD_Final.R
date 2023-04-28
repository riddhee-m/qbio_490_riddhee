#Set working directory
knitr::opts_knit$set(root.dir = normalizePath("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/final_project_colorectal/outputs"))

#Install and load all libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools")
if (!require("survival", quietly = TRUE))
  BiocManager::install("survival")
if (!require("survminer", quietly = TRUE))
  BiocManager::install("survminer")
if (!require("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")

library(BiocManager)
library(TCGAbiolinks)
library(maftools)
library(survival)
library(survminer)
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)

clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")

colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"
maf_query <- GDCquery(project = "TCGA-COAD", data.category = "Simple Nucleotide Variation", access = "open", data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf, clinicalData = clinic, isTCGA = TRUE)

rna_query <- GDCquery(project ="TCGA-COAD", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

# Determine if there are any NA values (there are none)
sum(is.na(clinic$age_at_initial_pathologic_diagnosis))


# Create a column in clinic with Young/Old 
young_mask <- ifelse((clinic$age_at_initial_pathologic_diagnosis < 50), T, F)
clinic$age_status <- ifelse(young_mask, "Early-Onset", "Late-Onset")


#making a survival time column for survival plots
clinic$survival_time <- ifelse(is.na(clinic$days_to_death),
                               clinic$survival_time <- clinic$days_to_last_followup,
                               clinic$survival_time <- clinic$days_to_death)

#remove any -Inf and NA values in survival_time
cleaned_clinic <- clinic
na_mask <- ifelse(is.na(clinic$survival_time), F, T)
cleaned_clinic <- cleaned_clinic[na_mask, ]
inf_mask <- ifelse(cleaned_clinic$survival_time == "-Inf", F, T)
cleaned_clinic <- cleaned_clinic[inf_mask, ]

#making a death event (T/F) column for survival plots
cleaned_clinic$death_event <- ifelse(cleaned_clinic$vital_status == "Alive",
                                     cleaned_clinic$death_event <- FALSE,
                                     cleaned_clinic$death_event <- TRUE)

#initializing a survival object
surv_object_age <- Surv(time = cleaned_clinic$survival_time, event = 
                          cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ cleaned_clinic$age_status, data = cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))

#show plot
KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/COAD_KM_plot_age.jpg")
KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                      axis.text = element_text(size=16),
                                                      legend.title = element_text(size=14),
                                                      legend.text = element_text(size=12))
KM_plot_age
dev.off()


################ Male Plot ####################
male_mask <- ifelse((cleaned_clinic$gender == 'MALE'), T, F)
male_cleaned_clinic <- cleaned_clinic[male_mask, ]

#initializing a survival object
surv_object_age <- Surv(time = male_cleaned_clinic$survival_time, event = 
                          male_cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ male_cleaned_clinic$age_status, data = male_cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
male_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                           axis.text = element_text(size=16),
                                                           legend.title = element_text(size=14),
                                                           legend.text = element_text(size=12))

#show plot
male_KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/male_COAD_KM_plot_age.jpg")
male_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                           axis.text = element_text(size=16),
                                                           legend.title = element_text(size=14),
                                                           legend.text = element_text(size=12))
male_KM_plot_age
dev.off()


################ Female Plot ####################
female_mask <- ifelse((cleaned_clinic$gender == 'FEMALE'), T, F)
female_cleaned_clinic <- cleaned_clinic[female_mask, ]

#initializing a survival object
surv_object_age <- Surv(time = female_cleaned_clinic$survival_time, event = 
                          female_cleaned_clinic$death_event)

#creating a fit object
age_fit <- survfit(surv_object_age ~ female_cleaned_clinic$age_status, data = female_cleaned_clinic)

#format and create KM plot
survplot_age <- ggsurvplot(age_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1),
                                                                                  "cm")),
                           legend = "right")

#save the plot to a variable
female_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                             axis.text = element_text(size=16),
                                                             legend.title = element_text(size=14),
                                                             legend.text = element_text(size=12))

#show plot
female_KM_plot_age

#save plot
jpeg("/Users/nataliefortunato/Documents/qbio_490_nataliefortunato/female_COAD_KM_plot_age.jpg")
female_KM_plot_age <- survplot_age$plot + theme_bw() + theme(axis.title = element_text(size=20),
                                                             axis.text = element_text(size=16),
                                                             legend.title = element_text(size=14),
                                                             legend.text = element_text(size=12))
female_KM_plot_age
dev.off()


################ Oncoplots ####################

#remove patients with no age information
age_na_mask <- ifelse(is.na(clinic$age_at_initial_pathologic_diagnosis), F, T)
age_cleaned_clinical <- clinic[age_na_mask, ]

#create a column in age_cleaned_clinical where patients are labeled "Young", "Middle", or "Old"
EO_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 50, T, F)
LO_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis > 50, T, F)
age_cleaned_clinical$diagnosis_status <- ifelse(EO_mask, "Early Onset", "Late Onset")

maf_object@clinical.data$diagnosis_status <- ifelse(maf_object@clinical.data$age_at_initial_pathologic_diagnosis > 50,'Late Onset', 'Early Onset')

#EO oncoplot 
EO_mask <- ifelse(maf_object@clinical.data$diagnosis_status == 'Early Onset', T, F)

EO_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[EO_mask]

EO_maf <- subsetMaf(maf = maf_object,
                    tsb = EO_patient_barcodes)

oncoplot(EO_maf, 10)

#LO oncoplot 
LO_mask <- ifelse(maf_object@clinical.data$diagnosis_status == 'Late Onset', T, F)

LO_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[LO_mask]

LO_maf <- subsetMaf(maf = maf_object,
                    tsb = LO_patient_barcodes)

oncoplot(LO_maf, 10)

################ Volcano Plot ####################

#Create rna_clinical
rna_clinical <- rna_se@colData[!is.na(rna_se@colData$age_at_index), ]

#OR

na_mask <- !is.na(rna_se@colData$age_at_index)
rna_clinical <- rna_se@colData[na_mask, ]

rna_clinical <- as.data.frame(rna_clinical)

#remove treatments
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F, T)

rna_clinical <- rna_clinical[, treatments_mask]

#Create rna_genes
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)

#Create rna_counts
rna_counts <- rna_se@assays@data$unstranded[, na_mask]
rna_counts <- as.data.frame(rna_counts)

#Create diagnosis status column in rna_clinical
rna_clinical$diagnosis_status <- ifelse(rna_clinical$age_at_index <= 50, "EO", "LO")
rna_clinical$diagnosis_status <- factor(rna_clinical$diagnosis_status)
summary(rna_clinical$diagnosis_status)

#Update row names and column names
rownames(rna_genes) <- rna_genes$gene_id
rownames(rna_counts) <- rownames(rna_genes)
colnames(rna_counts) <- rownames(rna_clinical)

#Remove Solid Tissue Normal patients
def_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_counts <- rna_counts[ ,def_mask]
rna_clinical <- rna_clinical[def_mask, ]
unique(rna_clinical$definition)

#Create rna_counts and rna_clinical for EO
EO_mask <- ifelse(rna_clinical$diagnosis_status == "EO", T, F)
EO_rna_counts <- rna_counts[ ,EO_mask]
EO_rna_clinical <- rna_clinical[EO_mask, ]

#Create rna_counts and rna_clinical for LO
LO_mask <- ifelse(rna_clinical$diagnosis_status == "LO", T, F)
LO_rna_counts <- rna_counts[ ,LO_mask]
LO_rna_clinical <- rna_clinical[LO_mask, ]

#check for NA's
sum(is.na(EO_rna_clinical$diagnosis_status))
sum(is.na(EO_rna_clinical$ajcc_pathologic_stage))
sum(is.na(EO_rna_clinical$gender))

sum(is.na(LO_rna_clinical$diagnosis_status))
sum(is.na(LO_rna_clinical$ajcc_pathologic_stage))
sum(is.na(LO_rna_clinical$gender))

#remove NA's
na_mask <- !is.na(rna_clinical$ajcc_pathologic_stage)

rna_clinical <-  rna_clinical[na_mask,]

rna_counts <-  rna_counts[, na_mask]

# use rowSums() to create a list with the total number of counts of each gene
row_sums <- rowSums(rna_counts)

# create a boolean mask where genes with < 10 total counts are FALSE, and genes with >= 10 total counts are TRUE
low_counts_mask <- ifelse(row_sums >= 10, T, F)

# rewrite the rna_counts df, subsetting for only genes with >= 10 total counts
rna_counts <- rna_counts[low_counts_mask,]

#update rna_genes with the low_counts_mas
rna_genes <- rna_genes[low_counts_mask,]

#Conduct DESeq2
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~diagnosis_status)

?DESeq
dds_obj <- DESeq(dds) # note: this will likely take a long time (ie 45 minutes to 2 hours)

?resultsNames
resultsNames(dds_obj)  # see what comparisons got run

# get the EO vs. LO comparison
?results

new_results <- results(dds_obj, format = "DataFrame", contrast = c("diagnosis_status", "EO", "LO")) # this is case sensitive so be careful to match it with your age_category factors closely!
  
new_results <- data.frame(rna_genes$gene_name, new_results@rownames, new_results@listData$log2FoldChange, new_results@listData$pvalue, new_results@listData$padj, -log10(new_results@listData$padj))

colnames(new_results) <- c("Gene Name", "Gene ID", "Log2 Fold Change", "P-value", "PAdj", "-log10(PAdj)")  
  
#Create enhanced volcano plot
par(mar=c(1,1,1,1))
EnhancedVolcano(new_results, 
                lab = results$`Gene Name`, 
                x = "Log2 Fold Change", 
                y = "PAdj")
jpeg("COAD_Volcano.jpg")

  