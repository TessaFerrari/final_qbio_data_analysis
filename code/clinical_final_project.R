#BiocManager::install(version='devel')
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)

query_clin <- GDCquery(project = "TCGA-LUAD", 
                       data.category = "Clinical",
                       file.type = "xml")
GDCdownload(query_clin)
clinic <- GDCprepare_clinic(query_clin, clinical.info = "patient")
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"

#SURVIVAL PLOT
library(survival)
library(survminer)
#adding info into the days_to_death column using days_to_last_known_alive and days_to_last_followup as substitutes
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_known_alive, clinic$days_to_death)
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_followup, clinic$days_to_death)
#create a boolean death_event column using vital_status
clinic$death_event = ifelse(clinic$vital_status == "Dead", T, F)

#removing the two patients that still do not have days_to_death info
death_mask=!is.na(clinic$days_to_death)
clinic=clinic[death_mask,]

#eliminating the clinic samples that have no race data and the 1 sample in american indian...
race_mask = (!is.na(clinic$race_list) & !clinic$race_list=="" & !clinic$race_list=="AMERICAN INDIAN OR ALASKA NATIVE" )
clinic = clinic[race_mask,]

# initialize a survival object for the time frame that we want
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# initialize a fit object for race 
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

# create the  survival plot, while inputting proper margins and putting the legend to the right
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# tweaking the formatting of the race plot
race_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# calling and showing the race survival plot
race_surv_plot



#getting rid of the empty string level in race_list
clinic$race_list = droplevels(clinic$race_list)

#graphing the boxplot
boxplot(clinic$age_at_initial_pathologic_diagnosis~clinic$race_list,
        xlab="Race",
        ylab='Age at Initial Diagnosis', main="Race vs Age at Initial Diagnosis", ylim=c(30,90))
