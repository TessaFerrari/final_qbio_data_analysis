# bar charts for cancer stage and mutation%

library(TCGAbiolinks)
library(maftools)

# Clinical
query_clin <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical",
                       file.type = "xml")
GDCdownload(query_clin)
clinic <- GDCprepare_clinic(query_clin, clinical.info = "patient")
# clinic <- GDCprepare_clinic(query_clin, clinical.info = "follow_up")
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"


# MAF
query_maf <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf)
maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep,
                       clinicalData = clinic,
                       isTCGA = TRUE)
clinical <- maf_object@clinical.data
NA_mask = is.na(clinical$race_list)
clinical_clean = clinical[!NA_mask]

white_patient_ids <- clinical_clean$Tumor_Sample_Barcode[clinical_clean$race_list == "WHITE"]
black_patient_ids <- clinical_clean$Tumor_Sample_Barcode[clinical_clean$race_list == "BLACK OR AFRICAN AMERICAN"]
asian_patient_ids <- clinical_clean$Tumor_Sample_Barcode[clinical_clean$race_list == "ASIAN"]
indian_patient_ids <- clinical_clean$Tumor_Sample_Barcode[clinical_clean$race_list == "AMERICAN INDIAN OR ALASKA NATIVE"]
white_maf = subsetMaf(maf = maf_object,
                      tsb = white_patient_ids)
black_maf = subsetMaf(maf = maf_object,
                      tsb = black_patient_ids)
asian_maf = subsetMaf(maf = maf_object,
                      tsb = white_patient_ids)
indian_maf = subsetMaf(maf = maf_object,
                       tsb = indian_patient_ids)


white_stages = white_maf@clinical.data$stage_event_pathologic_stage
white_stages_mask = is.na(white_stages)
white_stages_clean = white_stages[!white_stages_mask]
num_white_1 = sum(white_stages_clean == "Stage I" | white_stages_clean == "Stage IA" | white_stages_clean == "Stage IB" )
num_white_2 = sum(white_stages_clean == "Stage II" | white_stages_clean == "Stage IIA" | white_stages_clean == "Stage IIB")
num_white_3 = sum(white_stages_clean == "Stage III" | white_stages_clean == "Stage IIIA" | white_stages_clean == "Stage IIIB" |
                    white_stages_clean == "Stage IIIC")
num_white_4 = sum(white_stages_clean == "Stage IV")


black_stages = black_maf@clinical.data$stage_event_pathologic_stage
black_stages_mask = is.na(black_stages)
black_stages_clean = black_stages[!black_stages_mask]
num_black_1 = sum(black_stages_clean == "Stage I" | black_stages_clean == "Stage IA" | black_stages_clean == "Stage IB" )
num_black_2 = sum(black_stages_clean == "Stage II" | black_stages_clean == "Stage IIA" | black_stages_clean == "Stage IIB")
num_black_3 = sum(black_stages_clean == "Stage III" | black_stages_clean == "Stage IIIA" | black_stages_clean == "Stage IIIB" |
                    black_stages_clean == "Stage IIIC")
num_black_4 = sum(black_stages_clean == "Stage IV")

# asian_stages = asian_maf@clinical.data$stage_event_pathologic_stage

barplot(names=c("Stage I", "Stage II", "Stage III", "Stage IV"), c(num_white_1, num_white_2, num_white_3, num_white_4))
barplot(names=c("Stage I", "Stage II", "Stage III", "Stage IV"), c(num_black_1, num_black_2, num_black_3, num_black_4))

labelslist = c("Stage I", "Stage II", "Stage III", "Stage IV")
countswhite = c(num_white_1, num_white_2, num_white_3, num_white_4)
countsblack = c(num_black_1, num_black_2, num_black_3, num_black_4)

pie(countsblack, labels = labelslist, main="Cancer Stage For Black Patients")
pie(countswhite, labels = labelslist, main="Cancer Stage For Black Patients")


