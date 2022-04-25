# final project analysis

BiocManager::install(version='devel')

BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)

# Clinical
# query_clin <- GDCquery(project = "TCGA-LUAD",
#                        data.category = "Clinical",
#                        file.type = "xml")
# GDCdownload(query_clin)
# clinic <- GDCprepare_clinic(query_clin, clinical.info = "followup")
# colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"

# MAF
library(maftools)
query_maf <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      legacy = F)
GDCdownload(query_maf)
maf_prep <- GDCprepare(query_maf)
maf_object <- read.maf(maf = maf_prep,
                       clinicalData = clinic,
                       isTCGA = TRUE)

oncoplot(maf = maf_object,
         top = 10)

clinical <- maf_object@clinical.data
clinic <- GDCprepare_clinic(query_clin, clinical.info = "follow_up")
# unique(clinical$race_list)
# clean up NA's before splitting up
NA_mask = is.na(clinical$race_list)
clinical_clean = clinical[!NA_mask]

# create KM plot
# library(survival)
# library(survminer)
# clinical_clean$days_to_death <- ifelse(is.na(clinical_clean$days_to_death), 
#                                        clinical_clean$days_to_last_follow_up, 
#                                        clinical_clean$days_to_death)
# clinical_clean$death_event <- ifelse(clinical_clean$vital_status == "Alive", 0, 1)
# # We initialize a 'survival' object first, which contains the data we need.
# surv_object <- Surv(time = clinical_clean$days_to_death, 
#                     event = clinical_clean$death_event)
# 
# # We then create a fit object
# invasion_fit <- surv_fit( surv_object ~ clinical_clean$race_list, data = clinical_clean )
# 
# survplot = ggsurvplot(invasion_fit, 
#                       pval=TRUE, 
#                       ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
#                       legend = "right")
# 
# p = survplot$plot + 
#   theme_bw() +  # changes the appearance to be a bit prettier
#   theme(axis.title = element_text(size=20), # increase font sizes
#         axis.text = element_text(size=16),
#         legend.title = element_text(size=14),
#         legend.text = element_text(size=12))
# p


# now going to do some work with the maf data
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

# "/Users/btinsley/Desktop/qbio_student_group/final_qbio_data_analysis/white_black_oncoplot.png"
png("/Users/btinsley/Desktop/qbio_student_group/final_qbio_data_analysis/white_black_oncoplot.png", width=1800, height=700)
coOncoplot(m1 = white_maf, 
           m2 = black_maf, 
           m1Name = "WHITE", 
           m2Name = "BLACK",
           legendFontSize = 2.0,
           titleFontSize = 2,
           geneNamefont = 1.0,
           outer_mar = 4,
           gene_mar = 1.5)
dev.off()

lollipopPlot2(m1 = white_maf, 
              m2 = black_maf, 
              m1_name = "White",
              m2_name = "Black",
              gene = "TP53")
lollipopPlot2(m1 = white_maf, 
              m2 = black_maf, 
              m1_name = "White",
              m2_name = "Black",
              gene = "TTN")
lollipopPlot2(m1 = white_maf, 
              m2 = black_maf, 
              m1_name = "White",
              m2_name = "Black",
              gene = "MUC16")
lollipopPlot2(m1 = white_maf, 
              m2 = black_maf, 
              m1_name = "White",
              m2_name = "Black",
              gene = "CSMD3")



# Transcriptomics (RNASeq)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "STAR - Counts")
GDCdownload(query, method="api")
sum_exp = GDCprepare(query)
counts = assays(sum_exp)$unstranded
# remove counts for American Indians and for Not Reported
no_na = colData(sum_exp)[ colData(sum_exp)$race != "not reported" , ]
df_clean = no_na[no_na$race != "american indian or alaska native", ]
df_clean[df_clean$race == "black or african american", ] <- "black"
# also need to clean up counts
counts_clean = counts[, colData(sum_exp)$race != "not reported" & colData(sum_exp)$race != "american indian or alaska native"]
# counts_clean = counts[, colData(sum_exp)$race != "american indian or alaska native"]
# counts_clean = counts[, colData(sum_exp)$race != "not reported"]

colnames(colData(sum_exp))
colData(sum_exp)$race

# pick 3 genes from oncoplot, TP53, TTN, and MUC16
geneA_id_mask = (rowData(sum_exp)$gene_name == "TP53")
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$gene_id[geneA_id_mask]
ensembl_geneA

geneB_id_mask = (rowData(sum_exp)$gene_name == "TTN")
sum(geneB_id_mask) 
ensembl_geneB = rowData(sum_exp)$gene_id[geneB_id_mask]
ensembl_geneB

geneC_id_mask = (rowData(sum_exp)$gene_name == "MUC16")
sum(geneC_id_mask)
ensembl_geneC = rowData(sum_exp)$gene_id[geneC_id_mask] 
ensembl_geneC

geneD_id_mask = (rowData(sum_exp)$gene_name == "CSMD3")
sum(geneD_id_mask)
ensembl_geneD = rowData(sum_exp)$gene_id[geneD_id_mask] 
ensembl_geneD

# min(assays(sum_exp)$unstranded[ ensembl_geneA, ])
# max(assays(sum_exp)$unstranded[ ensembl_geneA, ])

# summary(assays(sum_exp)$unstranded[ ensembl_geneA, ])
# summary(assays(sum_exp)$unstranded[ ensembl_geneB, ])
# summary(assays(sum_exp)$unstranded[ ensembl_geneC, ])

# plot(assays(sum_exp)$unstranded[ ensembl_geneA, ],
#      assays(sum_exp)$unstranded[ ensembl_geneB, ],
#      xlab = "TP53", ylab = "TTN")

gene_countsA = counts_clean[geneA_id_mask, ]
gene_countsB = counts_clean[geneB_id_mask, ]
gene_countsC = counts_clean[geneC_id_mask, ]
gene_countsD = counts_clean[geneD_id_mask, ]
boxplot(gene_countsA ~ df_clean$race, xlab = "Race", ylab = "Counts", 
        main = "Counts for TP53 by Race", las=2)
boxplot(gene_countsB ~ df_clean$race, xlab = "Race", ylab = "Counts", 
        main = "Counts for TTN by Race", las=2)
boxplot(gene_countsC ~ df_clean$race, xlab = "Race", ylab = "Counts", 
        main = "Counts for MUC16 by Race", las=2)
boxplot(gene_countsD ~ df_clean$race, xlab = "Race", ylab = "Counts", 
        main = "Counts for CSMD3 by Race", las=2)
par(mar = c(3,3,3,3))


# pie chart
clinical$stage_event_pathologic_stage
# white_maf@clinical.data$stage_event_pathologic_stage
# black_maf@clinical.data$stage_event_pathologic_stage
# length(unique(white_maf@clinical.data$stage_event_pathologic_stage))
# length(unique(black_maf@clinical.data$stage_event_pathologic_stage))

tab <- table(white_maf@clinical.data$stage_event_pathologic_stage)
labelswhite = names(tab)
countswhite = c(tab[1], tab[2], tab[3], tab[4], tab[5], tab[6], tab[7],tab[8])
pie(countswhite, labels = labelswhite, main="Cancer Stage For White Patients")

tab2 <- table(black_maf@clinical.data$stage_event_pathologic_stage)
labelsblack = names(tab2)
labels = paste0(count, "%")
countsblack = c(tab2[1], tab2[2], tab2[3], tab2[4], tab2[5], tab2[6], tab2[7],tab2[8])
pie(countsblack, labels = labelsblack, main="Cancer Stage For Black Patients")


par(mar=c(0,0,0,0))     # Removes margins
par(mar=c(5,4,4,2)+0.1) # Default
par(mar=c(8,8,4,2))     # Larger bottom and left margins

