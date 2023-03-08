
#Setup

#Set working directory
setwd("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee")

#Read in clinical data csv
clin_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
clinical <- read.csv("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/analysis_data/brca_clinical_data.csv")

#Prepare drug and radiation data
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")

#Initialize  maf_object with clinical annotations
maf_query <- GDCquery( project = "TCGA-BRCA",
                       data.category = "Simple Nucleotide Variation",
                       access = "open",
                       data.type = "Masked Somatic Mutation",
                       workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

#GDCdownload(maf_query)
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf,
                       clinicalData = clinical,
                       isTCGA = TRUE)

#Coding Activity

#1. 

#Create mask to remove patients that don't fall under Chemotherapy or Hormone Therapy groups
therapy_types_mask <- ifelse(clinical_drug$therapy_types == "Chemotherapy" | clinical_drug$therapy_types == "Hormone Therapy", T, F)

#Apply mask
therapy_types_cleaned_clinical_drug <- clinical_drug[therapy_types_mask, ]

#Rewrite column as a factor
therapy_types_cleaned_clinical_drug$therapy_types <- as.factor(therapy_types_cleaned_clinical_drug$therapy_types)

str(therapy_types_cleaned_clinical_drug$therapy_types)

#2. 

#Create mask to identify patients within chemotherapy group
chemotherapy_mask <- ifelse(therapy_types_cleaned_clinical_drug$therapy_types == 'Chemotherapy', T, F)

#Apply chemotherapy mask to column of patient barcodes
chemotherapy_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[chemotherapy_mask]

#Subset maf_object to only contain chemotherapy group
chemotherapy_maf <- subsetMaf(maf = maf_object,
                          tsb = chemotherapy_patient_barcodes)

#Create mask to identify patients within hormone therapy group
hormone_mask <- ifelse(therapy_types_cleaned_clinical_drug$therapy_types == 'Hormone Therapy', T, F)

#Apply hormone therapy mask to column of patient barcodes
hormone_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[hormone_mask]

#Subset maf_object to only contain hormone therapy group
hormone_maf <- subsetMaf(maf = maf_object,
                          tsb = hormone_patient_barcodes)

#Create co-oncoplot
coOncoplot(m1 = chemotherapy_maf, 
           m2 = hormone_maf, 
           m1Name = 'Chemotherapy', 
           m2Name = 'Hormone Therapy', 
           borderCol = NA)

#PIK3CA has discrepancy in % mutated between chemotherapy and hormone therapy groups. The PIK3CA gene
#may allow the PI3K enzyme to become extremely active and allow cancer cells to grow. The chemotherapy group
#might have a lower % mutated than the hormone therapy group because it targets to kill fast growing cells in the body
#that have mutations in PIK3CA.

#3. 

#Subset maf_object to only contain patients that are positive for TP53
TP53_maf <- subsetMaf(maf = maf_object,
                      genes = 'TP53') ## fill in with your gene name as a string

#
TP53_patient_barcodes <- TP53_maf@clinical.data$Tumor_Sample_Barcode

#Find the number of TP53 positive patients
num_TP53_pos <- length(TP53_patient_barcodes)

#Find the number of patients in chemotherapy group
num_chemo <- length(chemotherapy_patient_barcodes)

#Find the number of patients in hormone therapy group
num_hormone <- length(hormone_patient_barcodes)

#Find the number of TP53 negative patients in chemotherapy group
num_TP53_neg_chemo <- length(chemotherapy_patient_barcodes) - length(TP53_patient_barcodes)

#Find the number of TP53 negative patients in hormone therapy group
num_TP53_neg_hormone <- length(hormone_patient_barcodes) - length(TP53_patient_barcodes)

#Find the number of TP53 positive patients in chemotherapy group
num_TP53_pos_chemo <- length(intersect(chemotherapy_patient_barcodes, TP53_patient_barcodes))

#Find the number of TP53 positive patients in hormone therapy group
num_TP53_pos_hormone <- length(intersect(hormone_patient_barcodes, TP53_patient_barcodes))

#Create contingency table with therapy type and TP53
contig <- matrix(c(num_TP53_pos_chemo,
                   num_TP53_neg_chemo,
                   num_TP53_pos_hormone,
                   num_TP53_neg_hormone)
                 , nrow = 2)

#Run Fisher's exact test
fisher.test(contig)

#Create mosaic plot
mosaicplot(contig)

#The output of the Fisher's exact test suggests that the relation between TP53 mutation and therapy type is statistically
#significant. The odds ratio shows that there is a 0.596% chance of being positive for TP53 mutation if therapy type is 
#chemotherapy compared to hormone therapy.
                 
#4.

#Create lollipop plot with descriptive names
lollipopPlot2(m1 = chemotherapy_maf, 
              m2 = hormone_maf, 
              m1_name = 'Chemotherapy',
              m2_name = 'Hormone Therapy',
              gene = "TP53") ## pick any gene of your choosing to fill in here

#There are no significant differences in mutations in terms of sites, types, and number between the two
#therapy type groups. 

#5. 

#Create overall survival status column
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$vital_status == 'Alive', T, F)

#Create mafSurvival KM plot based on mutations in TP53
mafSurvival(maf = maf_object,
            genes = "TP53", ## pick a gene of your choosing
            time = "days_to_last_followup", ## name of the column in maf_object@clinical.data containing survival time
            Status = "Overall_Survival_Status", ## name of the column that contains a boolean value for death events, you may need to recreate this... 
            isTCGA = TRUE)

#There is no significant difference between the survival probability of someone who is TP53 positive or negative, as the p-value is greater than 0.05. 
#This is because there is no difference in mutations between the two therapy groups.

