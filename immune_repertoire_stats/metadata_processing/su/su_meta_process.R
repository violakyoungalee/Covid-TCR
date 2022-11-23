library("data.table")
library("ggplot2")
library("umap")
library("Rtsne")
library("locStra")
library("Matrix")
library("scales")
library("lubridate")
library("stringr")
library("RColorBrewer")
library("dplyr")
library("plyr")
library("immunarch")


set.seed(20)



#################### 
# CD4 processing
#################### 
## from 254 COVID-19 blood draws (a draw near diagnosis (-BL) and a draw a few days later (-AC)) and 16 healthy donors.
## One blood draw was collected shortly after the initial clinical diagnosis (time = T1), and the second was collected a few days later (T2)
#samples_files <- substr(file.names,1,nchar(file.names)-4)
#samples_files <- samples_files[match(samplematch, samples_files)]

#meta <- fread("su_metadata.csv")
#metaHD <- fread("su_metadata.HD.csv")


meta <- fread("E-MTAB-9357.sdrf.txt")
clin_fused <- fread("clinical_meta_fused.csv",)
samplematch <- meta[,41]
samplematch <- samplematch$`Derived Array Data File`
samplematch <- substr(samplematch,1,nchar(samplematch)-7)
diseasestatus <- meta[,6]
diseasestatus <- diseasestatus$`Characteristics[disease]`
timepoint <-meta[,7]
timepoint <- timepoint$`Characteristics[sampling time point]`
timepoint[grep("BL", timepoint)] <- "T1"
timepoint[grep("AC", timepoint)] <- "T2"

individual <- meta[,5]
individual <- individual$`Characteristics[individual]`

#sex_healthy <- clin2$sex[match(individual, paste0("Healthy_",clin2$`Sample ID`))]
#age_healthy <- clin2$age[match(individual, paste0("Healthy_",clin2$`Sample ID`))]
#match(paste0(individual, ";",timepoint), paste0(individual, ";",timepoint))
newmeta <- clin_fused[match(paste0(individual, ";", timepoint), paste0(clin_fused$`Study Subject ID`,";", clin_fused$`Blood draw time point`)),]
temp <- clin_fused[match(individual, paste0("Healthy_", clin_fused$`Sample ID`)),]
newmeta[1:16,] <- temp[1:16,]
newmeta <- newmeta[,1:24]
newmeta <- cbind(data.frame(samplematch, diseasestatus, individual, timepoint), newmeta)
newmeta <- newmeta[,c(1,3:4,2,10:11, 8:9,12:28 )]
colnames(newmeta)[1:4] <- c("Sample", "Individual", "Timepoint", "Disease status")
newmeta$Sample <- paste0("su_cd4_VDJ_TCR_cd4_", substring(as.character(newmeta$Sample), 39))
write.table(newmeta, "su_meta_full_cd4.txt" ,quote = F,  row.names = F, sep = "\t")
write.table(newmeta[,1:6], "su_meta_sub_cd4.txt" ,quote = F,  row.names = F, sep = "\t")



#################### 
# CD8 processing
#################### 

## from 254 COVID-19 blood draws (a draw near diagnosis (-BL) and a draw a few days later (-AC)) and 16 healthy donors.
## One blood draw was collected shortly after the initial clinical diagnosis (time = T1), and the second was collected a few days later (T2)
#samples_files <- substr(file.names,1,nchar(file.names)-4)
#samples_files <- samples_files[match(samplematch, samples_files)]

meta <- fread("E-MTAB-9357.sdrf.txt")
clin_fused <- fread("clinical_meta_fused.csv",)
samplematch <- meta[,47]
samplematch <- samplematch$`Derived Array Data File`
samplematch <- substr(samplematch,1,nchar(samplematch)-7)
diseasestatus <- meta[,6]
diseasestatus <- diseasestatus$`Characteristics[disease]`
timepoint <-meta[,7]
timepoint <- timepoint$`Characteristics[sampling time point]`
timepoint[grep("BL", timepoint)] <- "T1"
timepoint[grep("AC", timepoint)] <- "T2"

individual <- meta[,5]
individual <- individual$`Characteristics[individual]`

#sex_healthy <- clin2$sex[match(individual, paste0("Healthy_",clin2$`Sample ID`))]
#age_healthy <- clin2$age[match(individual, paste0("Healthy_",clin2$`Sample ID`))]

#match(paste0(individual, ";",timepoint), paste0(individual, ";",timepoint))
newmeta <- clin_fused[match(paste0(individual, ";", timepoint), paste0(clin_fused$`Study Subject ID`,";", clin_fused$`Blood draw time point`)),]
temp <- clin_fused[match(individual, paste0("Healthy_", clin_fused$`Sample ID`)),]
newmeta[1:16,] <- temp[1:16,]
newmeta <- newmeta[,1:24]
newmeta <- cbind(data.frame(samplematch, diseasestatus, individual, timepoint), newmeta)
newmeta <- newmeta[,c(1,3:4,2,10:11, 8:9,12:28 )]
colnames(newmeta)[1:4] <- c("Sample", "Individual", "Timepoint", "Disease status")
newmeta$Sample <- paste0("su_cd8_VDJ_", substring(as.character(newmeta$Sample), 39))
write.table(newmeta, "su_meta_full_cd8.txt" ,quote = F,  row.names = F, sep = "\t")
write.table(newmeta[,1:6], "su_meta_sub_cd8.txt" ,quote = F,  row.names = F, sep = "\t")
