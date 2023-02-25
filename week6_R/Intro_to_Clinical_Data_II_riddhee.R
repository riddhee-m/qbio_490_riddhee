
#Necessary packages for KM plot
if(!require(survival)){
  install.packages("survival")
}
library(survival)

if(!require(survminer)){
  install.packages("survminer")
}
library(survminer)

if(!require(ggplot2)){
  install.packages("ggplot2")
}
library(ggplot2)

#Set working directory
setwd("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/analysis_data")

#Read clinical data csv
clinical <- read.csv("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/analysis_data/brca_clinical_data.csv")

#Prepare drug and radiation data
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug")
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation")


#1. 

#Check for NAs
sum(is.na(clinical$age_at_initial_pathologic_diagnosis))

#Age at diagnosis

#2.

#Age at diagnosis is a discrete variable.

#3.

#Check for NAs
sum(is.na(clinical_drug$therapy_types))

#Therapy types

#4. 

#Radiation type is a categorical variable. 

#5. 

#Hypothesis 1: The average age of diagnosis for patients treated with hormone therapy is 60 years.

#Hypothesis 2: Patients diagnosed under 50 are most likely to survive breast cancer.

#Hypothesis 3: Patients treated with hormone therapy are most likely to survive breast cancer.

#################################################################

#Coding Activity

#1. 

#Remove empty values in therapy types column
therapy_types_mask <- ifelse(clinical_drug$therapy_types == "", F, T)

therapy_types_cleaned_clinical <- clinical[therapy_types_mask, ]

therapy_types_cleaned_clinical_drug <- clinical_drug[therapy_types_mask, ]

#Create boxplot for comparing therapy types and diagnosis age
boxplot(therapy_types_cleaned_clinical$age_at_initial_pathologic_diagnosis~therapy_types_cleaned_clinical_drug$therapy_types)

#From the boxplot, I can determine the average age of diagnosis of patients being treated by different therapy types.
#I chose the boxplot because it can show the range of diagnosis ages of patients by therapy type rather than a scatterplot, which would show 
#a lot of points with no indication of the averages and quartiles.

#2. 

#remove patients with no age information
age_na_mask <- ifelse(is.na(clinical$age_at_initial_pathologic_diagnosis), F, T)
age_cleaned_clinical <- clinical[age_na_mask, ]

#create a column in age_cleaned_clinical where patients are labeled "Young", "Middle", or "Old"
young_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 35, T, F)
middle_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis >= 35 & age_cleaned_clinical$age_at_initial_pathologic_diagnosis <= 50, T, F)
old_mask <- ifelse(age_cleaned_clinical$age_at_initial_pathologic_diagnosis >= 35, T, F)
age_cleaned_clinical$age_status <- ifelse(young_mask, "Young", ifelse(middle_mask, "Middle", "Old"))

#making a survival time column for survival plots
age_cleaned_clinical$survival_time <- ifelse(is.na(age_cleaned_clinical$days_to_death),
                                             age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_last_followup,
                                             age_cleaned_clinical$survival_time <- age_cleaned_clinical$days_to_death)

#remove any -Inf values in survival_time
inf_mask <- ifelse(age_cleaned_clinical$survival_time == "-Inf", F, T)
age_cleaned_clinical <- age_cleaned_clinical[inf_mask, ]

#making a death event (T/F) column for survival plots
age_cleaned_clinical$death_event <- ifelse(age_cleaned_clinical$vital_status == "Alive", age_cleaned_clinical$death_event <- FALSE,
                                           age_cleaned_clinical$death_event <- TRUE)

#create survminer objects
survival_object <- Surv(time = age_cleaned_clinical$survival_time, event = age_cleaned_clinical$death_event)
fit_object <- survfit(survival_object ~ age_cleaned_clinical$age_status, data = age_cleaned_clinical)

#create the plot
survplot <- ggsurvplot(fit_object , pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
KM_plot_age <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

#show the plot
KM_plot_age

#5. (for age KM plot)

#Save KM plot as jpeg 
jpeg("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/week6_R/KM_plot_age.jpg")
KM_plot_age <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20),
                                                  axis.text = element_text(size = 16),
                                                  legend.title = element_text(size = 14),
                                                  legend.text = element_text(size = 12))
KM_plot_age
dev.off()

#3. 

#making a survival time column for survival plots
therapy_types_cleaned_clinical$survival_time <- ifelse(is.na(therapy_types_cleaned_clinical$days_to_death),
                                                            therapy_types_cleaned_clinical$survival_time <- therapy_types_cleaned_clinical$days_to_last_followup,
                                                            therapy_types_cleaned_clinical$survival_time <- therapy_types_cleaned_clinical$days_to_death)
#remove any -Inf values in survival_time
inf_mask <- ifelse(therapy_types_cleaned_clinical$survival_time == "-Inf", F, T)
therapy_types_inf_cleaned_clinical <- therapy_types_cleaned_clinical[inf_mask, ]
therapy_types_inf_cleaned_clinical_drug <- therapy_types_cleaned_clinical_drug[inf_mask, ]


#making a death event (T/F) column for survival plots
therapy_types_inf_cleaned_clinical$death_event <- ifelse(therapy_types_inf_cleaned_clinical$vital_status == "Alive", therapy_types_inf_cleaned_clinical$death_event <- FALSE, therapy_types_inf_cleaned_clinical$death_event <- TRUE)

#create survminer objects
survival_object <- Surv(time = therapy_types_inf_cleaned_clinical$survival_time, event = therapy_types_inf_cleaned_clinical$death_event)
fit_object <- survfit(survival_object ~ therapy_types_inf_cleaned_clinical_drug$therapy_types, data = therapy_types_inf_cleaned_clinical_drug)

#create the plot
survplot <- ggsurvplot(fit_object, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1,1,1,1), "cm")), legend = "right")
KM_plot_therapy_types <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16), legend.title = element_text(size = 14), legend.text = element_text(size = 12))

#show the plot
KM_plot_therapy_types

#4.Analyze your two KM plots. What do the KM plots suggest about the impact of your variables on breast cancer survival? What are the p-values? 
#Do the differences in survival appear to be significant? 

#The Age KM plot suggests that the "old" category of patients who were diagnosed at or above the age of 50, 
#they have a lower probability of surviving over time in comparison to the "young" and "middle" age categories.
#This supports my hypothesis. The p value is below 0.05, and the differences in survival appear to be significant.

#The Therapy Type KM plot suggests that patients treated with hormone therapy are least likely to survive breast cancer in
#comparison to other therapy types such as chemotherapy, ancillary, immunotherapy, targeted molecular therapy, and vaccine.
#This disproves my hypothesis. The p value is above 0.05, and the differences in survival appear to be significant.

#5. (for therapy types KM plot)

#Save KM plot as jpeg
jpeg("/Users/riddheemehta/Desktop/Code/QBIO/qbio_490_riddhee/week6_R/KM_plot_therapy_types.jpg")
KM_plot_therapy_types <- survplot$plot + theme_bw() + theme(axis.title = element_text(size = 20),
                                                  axis.text = element_text(size = 16),
                                                  legend.title = element_text(size = 14),
                                                  legend.text = element_text(size = 12))
KM_plot_therapy_types
dev.off()




