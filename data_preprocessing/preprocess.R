#### PREPROCESSING COVID TCR DATASETS 1.1.21
#### Zhang, Liao, Kusnadi, Su datasets 


#### Packages ----
library(data.table)
library(plyr)
library(dplyr)
library(strex)
library(tidyr)


#### LIAO ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses")

liao_meta <- fread("liao/liao_metadata.csv") ## metadata compiled manually from Supp Table 1

liao_TCR.files <- list.files(path="liao/liao_TCR", pattern="*.csv", full.names=TRUE, recursive=FALSE)
liao_TCR <- list()

for (i in 1:length(liao_TCR.files)){
  liao_TCR[[i]] <- fread(liao_TCR.files[i])
  liao_TCR[[i]] <- subset(liao_TCR[[i]], cdr3!="None") #omit rows with missing cdr3
  liao_TCR[[i]] <- subset(liao_TCR[[i]], chain!="Multi") #omit rows with "multi" instead of TRA or TRB
  liao_TCR[[i]] <- split(liao_TCR[[i]], liao_TCR[[i]]$barcode)
  liao_TCR[[i]] <- Filter(function(x){nrow(x) > 1}, liao_TCR[[i]]) # omit unpaired TRA/TRB
  liao_TCR[[i]] <- do.call(rbind, liao_TCR[[i]]) #join barcodes again
  liao_TCR[[i]] <- split(liao_TCR[[i]], liao_TCR[[i]]$chain) #split by TRA/TRB

  ## process TRA and TRB separately to keep chain with highest reads
  #TRA
  TRA <- liao_TCR[[i]]$TRA #only TRA
  TRA <- split(TRA, TRA$barcode) 
  for (j in 1:length(TRA)){
    if(nrow(TRA[[j]]) > 1) {
      TRA[[j]] <- TRA[[j]][which.max(TRA[[j]]$reads),]
    }
  } #keep TRA chain with higher reads 
  TRA <- do.call(rbind, TRA) 
  #TRB
  TRB <- liao_TCR[[i]]$TRB #only TRB
  TRB <- split(TRB, TRB$barcode) 
  for (j in 1:length(TRB)){
    if(nrow(TRB[[j]]) > 1) {
      TRB[[j]] <- TRB[[j]][which.max(TRB[[j]]$reads),]
    }
  } #keep TRB chain with higher reads
  TRB <- do.call(rbind, TRB) 
  
  ## combine TRA & TRB together 
  liao_TCR[[i]] <- rbind(TRA, TRB) 
  liao_TCR[[i]] <- dcast(setDT(liao_TCR[[i]]), barcode~chain, value.var=c('cdr3','cdr3_nt','v_gene', 'j_gene', 'c_gene', 'raw_clonotype_id'))
  
  ## omit rows with only TRA or TRB **** 1/9/21: KEEP ORPHAN CHAINS!
  #liao_TCR[[i]] <- liao_TCR[[i]] %>% subset(!is.na(liao_TCR[[i]]$cdr3_TRA)) 
  #liao_TCR[[i]] <- liao_TCR[[i]] %>% subset(!is.na(liao_TCR[[i]]$cdr3_TRB)) 

  ## add metadata
  liao_TCR[[i]]$patient_ID <- liao_TCR.files[i] %>% str_after_nth("/", 2) %>% str_before_nth("_", 2)  %>% str_after_nth("_", 1) #extract patient ID
  liao_TCR[[i]] <- left_join(liao_TCR[[i]], liao_meta, by="patient_ID")
  liao_TCR[[i]]$raw_clonotype_id_TRA <- paste0(liao_TCR[[i]]$patient_ID, "_", liao_TCR[[i]]$raw_clonotype_id_TRA) #add patient ID to clonotype (clonotypes are repeated across different patient samples)
  liao_TCR[[i]] <- data.frame(dataset = "Liao",
                              patient_ID = liao_TCR[[i]]$patient_ID, 
                              patient_subID = "NA",
                              cell_type = "T cell",
                              barcode = liao_TCR[[i]]$barcode,
                              TRA_cdr3 = liao_TCR[[i]]$cdr3_TRA,
                              TRA_cdr3_nt = liao_TCR[[i]]$cdr3_nt_TRA,
                              TRA_V_gene = liao_TCR[[i]]$v_gene_TRA, 
                              TRA_D_gene = "NA",
                              TRA_J_gene = liao_TCR[[i]]$j_gene_TRA, 
                              TRA_C_gene = liao_TCR[[i]]$c_gene_TRA, 
                              TRB_cdr3 = liao_TCR[[i]]$cdr3_TRB, 
                              TRB_cdr3_nt = liao_TCR[[i]]$cdr3_nt_TRB, 
                              TRB_V_gene = liao_TCR[[i]]$v_gene_TRB, 
                              TRB_D_gene = "NA",
                              TRB_J_gene = liao_TCR[[i]]$j_gene_TRB, 
                              TRB_C_gene = liao_TCR[[i]]$c_gene_TRB, 
                              clonotype = liao_TCR[[i]]$raw_clonotype_id_TRA, 
                              d_cond = liao_TCR[[i]]$disease_cond, 
                              sex = liao_TCR[[i]]$gender, 
                              age = liao_TCR[[i]]$age, 
                              race = "NA")
  ## save preprocessed files as separate patients first
  write.csv(liao_TCR[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/liao/liao_TCR_", unique(liao_TCR[[i]]$patient_ID), ".csv"))
}


## save aggregated patient file
liao_TCR <- rbindlist(liao_TCR)
liao_TCR <- as.data.frame(liao_TCR)
write.csv(liao_TCR, row.names=F, file = paste0("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr/norm_all features/liao_TCR_AllPatients.csv"))


#### ZHANG ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses")

zhang_TCR <- fread("zhang/zhang_TCR.csv")
zhang_meta <- fread("zhang/zhang_metadata.csv") 

zhang_TCR <- zhang_TCR[,-c(6,8)] 
zhang_TCR <- zhang_TCR %>% 
  mutate(CloneConsensusInfo = strsplit(CloneConsensusInfo, ";")) %>% 
  unnest(CloneConsensusInfo) # separates out the different consensus seq from one another
zhang_TCR <- separate(data = zhang_TCR, col = CloneConsensusInfo, into = c("consensus", "length", "tcr_chain", "v_gene", "d_gene", "j_gene", "c_gene", "AA_cdr3", "DNA_cdr3", "reads", "umis"), sep = "\\:")
zhang_TCR <- subset(zhang_TCR, !is.na(AA_cdr3)) #omit rows with missing cdr3

## Omit single TRA or TRB  **** 1/9/21: KEEP ORPHAN CHAINS!
#zhang_TCR <- split(zhang_TCR, zhang_TCR$CellName) 
#zhang_TCR <- Filter(function(x){nrow(x) > 1}, zhang_TCR) # omit unpaired TRA/TRB
#zhang_TCR <- do.call(rbind, zhang_TCR) #join barcodes again

## process TRA and TRB separately to keep chain with highest reads
zhang_TCR <- split(zhang_TCR, zhang_TCR$tcr_chain) #split by TRA/TRB

#TRA
TRA <- zhang_TCR$TRA #only TRA
TRA <- split(TRA, TRA$CellName) 
for (j in 1:length(TRA)){
  if(nrow(TRA[[j]]) > 1) {
    TRA[[j]] <- TRA[[j]][which.max(TRA[[j]]$reads),]
  }
} #keep TRA chain with higher reads 
TRA <- do.call(rbind, TRA) 
#TRB
TRB <- zhang_TCR$TRB #only TRB
TRB <- split(TRB, TRB$CellName) 
for (j in 1:length(TRB)){
  if(nrow(TRB[[j]]) > 1) {
    TRB[[j]] <- TRB[[j]][which.max(TRB[[j]]$reads),]
  }
} #keep TRB chain with higher reads
TRB <- do.call(rbind, TRB) 

## combine TRA & TRB together 
zhang_TCR <- rbind(TRA, TRB) 
rm(TRA,TRB)
zhang_TCR <- dcast(setDT(zhang_TCR), CellName~tcr_chain, value.var=c('FinalCellTypeTcell','AA_cdr3','DNA_cdr3','v_gene', 'd_gene', 'j_gene', 'c_gene', 'FinalCloneType'))

## omit rows with only TRA or TRB (if any, there should be none like this) **** 1/9/21: KEEP ORPHAN CHAINS!
#zhang_TCR <- zhang_TCR %>% subset(!is.na(zhang_TCR$AA_cdr3_TRA)) 
#zhang_TCR <- zhang_TCR %>% subset(!is.na(zhang_TCR$AA_cdr3_TRB)) 

## add metadata
zhang_TCR$patient_ID <- zhang_TCR$CellName %>% str_after_nth("_", 1) %>% str_before_nth("_", 1) #extract patient ID
zhang_TCR <- left_join(zhang_TCR, zhang_meta, by="patient_ID")
zhang_TCR <- data.frame(dataset = "Zhang",
                            patient_ID = zhang_TCR$patient_ID, 
                            patient_subID = zhang_TCR$patient_sub,
                            cell_type = zhang_TCR$FinalCellTypeTcell_TRA,
                            barcode = zhang_TCR$CellName,
                            TRA_cdr3 = zhang_TCR$AA_cdr3_TRA,
                            TRA_cdr3_nt = zhang_TCR$DNA_cdr3_TRA,
                            TRA_V_gene = zhang_TCR$v_gene_TRA, 
                            TRA_D_gene = zhang_TCR$d_gene_TRA, 
                            TRA_J_gene = zhang_TCR$j_gene_TRA, 
                            TRA_C_gene = zhang_TCR$c_gene_TRA, 
                            TRB_cdr3 = zhang_TCR$AA_cdr3_TRB, 
                            TRB_cdr3_nt = zhang_TCR$DNA_cdr3_TRB,
                            TRB_V_gene = zhang_TCR$v_gene_TRB, 
                            TRB_D_gene = zhang_TCR$d_gene_TRB, 
                            TRB_J_gene = zhang_TCR$j_gene_TRB, 
                            TRB_C_gene = zhang_TCR$c_gene_TRB, 
                            clonotype = zhang_TCR$FinalCloneType_TRA, 
                            d_cond = zhang_TCR$disease_cond, 
                            sex = zhang_TCR$gender, 
                            age = zhang_TCR$age, 
                            race = "NA")

#zhang_TCR[is.na(zhang_TCR)] <- "NA" ## replace all NAs with character NAs

## save preprocessed files as aggregated patients first
write.csv(zhang_TCR, row.names=F, file = paste0("preprocess_covidtcr/norm_all features/zhang_TCR_AllPatients.csv"))

## separate into individual files 
zhang_TCR <- as.data.frame(zhang_TCR)
zhang_TCR <- split(zhang_TCR, zhang_TCR$patient_ID)
for (i in 1:length(zhang_TCR)) {
  write.csv(zhang_TCR[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/zhang/zhang_TCR_", unique(zhang_TCR[[i]]$patient_ID), ".csv"))
}


#### KUSNADI ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses")

kusnadi_TCR <- fread("kusnadi/kusnadi_TCR.txt") 
kusnadi_TCR <- subset(kusnadi_TCR, orig.donor!="void") # omit patient w void patient_ID
#kusnadi_TCR <- kusnadi_TCR %>% subset(TRA.aa.chains.tag!="") %>% subset(TRB.aa.chains.tag!="")  **** 1/9/21: KEEP ORPHAN CHAINS!
kusnadi_TCR <- kusnadi_TCR[-which(kusnadi_TCR$TRA.aa.chains.tag == "" & kusnadi_TCR$TRB.aa.chains.tag == ""),]

## collect patient_ID and disease condition for editing (disease condition "void" are either HD or mild)
kusnadi_meta <- data.table(orig.donor = kusnadi_TCR[match(unique(kusnadi_TCR$orig.donor),kusnadi_TCR$orig.donor),]$orig.donor,
                           d_cond = kusnadi_TCR[match(unique(kusnadi_TCR$orig.donor),kusnadi_TCR$orig.donor),]$orig.severity_x)
kusnadi_meta$d_cond[25:32] <- "HD"
kusnadi_meta$d_cond[39:45] <- "mild"
kusnadi_meta$d_cond <- gsub("Mild","mild",kusnadi_meta$d_cond) 
kusnadi_meta$d_cond <- gsub("Moderate","moderate",kusnadi_meta$d_cond) 
kusnadi_meta$d_cond <- gsub("Severe","severe",kusnadi_meta$d_cond)

kusnadi_TCR <- left_join(kusnadi_TCR, kusnadi_meta, by="orig.donor")

kusnadi_TCR2 <- data.frame(dataset = "Kusnadi",
                             patient_ID = kusnadi_TCR$orig.donor, 
                             patient_subID = "NA",
                             cell_type = "SARS-Cov-2 reactive Memory CD8+",
                             barcode = kusnadi_TCR$barcode,
                             TRA_cdr3 = kusnadi_TCR$TRA.aa.chains.tag,
                             TRA_cdr3_nt = kusnadi_TCR$TRA.nt.chains.tag,
                             TRA_V_gene = "NA", 
                             TRA_D_gene = "NA",
                             TRA_J_gene = "NA", 
                             TRA_C_gene = "NA", 
                             TRB_cdr3 = kusnadi_TCR$TRB.aa.chains.tag, 
                             TRB_cdr3_nt = kusnadi_TCR$TRB.nt.chains.tag,
                             TRB_V_gene = "NA", 
                             TRB_D_gene = "NA",
                             TRB_J_gene = "NA", 
                             TRB_C_gene = "NA", 
                             clonotype = kusnadi_TCR$clonotype.tag, 
                             d_cond = kusnadi_TCR$d_cond,                           
                             sex = kusnadi_TCR$orig.sex, 
                             age = kusnadi_TCR$age, 
                             race = kusnadi_TCR$ethnicity)

#kusnadi_TCR2[is.na(kusnadi_TCR2)] <- "NA" ## replace all NAs with character NAs
rm(kusnadi_TCR)
## save preprocessed files as aggregated patients first
write.csv(kusnadi_TCR2, row.names=F, file = paste0("preprocess_covidtcr/norm_all features/kusnadi_TCR_AllPatients.csv"))

## separate into individual files 
kusnadi_TCR2 <- as.data.frame(kusnadi_TCR2)
kusnadi_TCR2 <- split(kusnadi_TCR2, kusnadi_TCR2$patient_ID)
for (i in 1:length(kusnadi_TCR2)) {
  write.csv(kusnadi_TCR2[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/kusnadi/kusnadi_TCR_", unique(kusnadi_TCR2[[i]]$patient_ID), ".csv"))
}


#### SU CD8+ ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses")

## metadata 1: diseased patients
su_meta <- fread("su/su_metadata.csv")
su_meta$ID <- gsub("INCOV0*","",su_meta$`Sample ID`) # add col with ID that matches the TCR data 
su_meta$ID <- gsub("-", "_", su_meta$ID)
su_meta$pat <- gsub("INCOV0*","Pt",su_meta$`Study Subject ID`) # add col with patient ID that includes both time point samples
WHO <- c("0", "1", "2", "1 or 2", "3", "4", "5", "6", "7")
d_cond <- c("HD", "mild", "mild", "mild", "moderate", "moderate", "severe", "severe", "severe")
su_meta$d_cond <- d_cond[match(su_meta$`Who Ordinal Scale`, WHO)]
su_meta <- data.frame(pat=su_meta$pat, patient=su_meta$ID, time_point=su_meta$`Blood draw time point`, d_cond=su_meta$d_cond, tenX=su_meta$`10X?`, 
                       sex=su_meta$Sex, age=su_meta$Age, ethnicity=su_meta$Ethnicity, race=su_meta$Race) 

## metadata 2: healthy donors
su_meta_HD <- fread("su/su_metadata.HD.csv")
su_meta_HD <- data.frame(pat=su_meta_HD$`Sample ID`, patient=su_meta_HD$`Sample ID`, time_point=NA, d_cond="HD", tenX=su_meta_HD$`10X?`, 
                          sex=su_meta_HD$sex, age=su_meta_HD$age, ethnicity=NA, race=NA) 

## combine metadata
su_meta <- rbind(su_meta, su_meta_HD); rm(su_meta_HD)
table(subset(su_meta, tenX == "yes")$d_cond) # patients with TCRs sequenced per d_cond

## load TCR files
tcr.files <- list.files(path="su/su_TCR_cd8", pattern="*.txt", full.names=TRUE, recursive=FALSE)
su_TCR <- list()
for (i in 1:length(tcr.files)){
  tcr <- fread(tcr.files[i])
  ID <- gsub(".*library_","",tcr.files[i]) # get patient ID
  ID <- gsub(".txt", "", ID)
  tcr <- data.frame(patient = ID, barcode = tcr$V1, 
                    TRA_cdr3 = tcr$TRA_1_cdr3, TRA_cdr3_nt = tcr$TRA_1_cdr3_nt, TRA_V_gene=tcr$TRA_1_v_gene, TRA_J_gene=tcr$TRA_1_j_gene, TRA_C_gene=tcr$TRA_1_c_gene, 
                    TRB_cdr3 = tcr$TRB_1_cdr3, TRB_cdr3_nt = tcr$TRB_1_cdr3_nt, TRB_V_gene=tcr$TRB_1_v_gene, TRB_J_gene=tcr$TRB_1_j_gene, TRB_C_gene=tcr$TRB_1_c_gene,
                   # TRA_cdr3_2 = tcr$TRA_2_cdr3, TRA_V_gene_2=tcr$TRA_2_v_gene, TRA_J_gene_2=tcr$TRA_2_j_gene, TRA_C_gene_2=tcr$TRA_2_c_gene, 
                   # TRB_cdr3_2 = tcr$TRB_2_cdr3, TRB_V_gene_2=tcr$TRB_2_v_gene, TRB_J_gene_2=tcr$TRB_2_j_gene, TRB_C_gene_2=tcr$TRB_2_c_gene,
                    clonotype=tcr$clonotype, clonotype_size=tcr$clonotype_size, chain_pairing=tcr$chain_pairing)
  su_TCR[[i]] <- tcr
  ## extracts only the first TRA and TRB
}
su_TCR <- rbindlist(su_TCR)

## combine TCR with metadata
su_TCR <- merge(su_TCR, su_meta, by="patient")

## process combined dataset 
#su_TCR <- subset(su_TCR, chain_pairing=="Single pair" | chain_pairing=="Extra alpha" | chain_pairing=="Extra beta" | chain_pairing=="Two full chains") **** 1/9/21: KEEP ORPHAN CHAINS!
su_TCR <- subset(su_TCR, chain_pairing=="Orphan alpha" | chain_pairing=="Orphan beta" | chain_pairing=="Single pair" | chain_pairing=="Extra alpha" | chain_pairing=="Extra beta" | chain_pairing=="Two full chains")
su_TCR <- data.frame(dataset = "Su",
                     patient_ID = su_TCR$patient, 
                     patient_subID = su_TCR$time_point,
                     cell_type = "CD8+",
                     barcode = su_TCR$barcode,
                     TRA_cdr3 = su_TCR$TRA_cdr3,
                     TRA_cdr3_nt = su_TCR$TRA_cdr3_nt,
                     TRA_V_gene = su_TCR$TRA_V_gene, 
                     TRA_D_gene = "NA",
                     TRA_J_gene = su_TCR$TRA_J_gene, 
                     TRA_C_gene = su_TCR$TRA_C_gene, 
                     TRB_cdr3 = su_TCR$TRB_cdr3, 
                     TRB_cdr3_nt = su_TCR$TRB_cdr3_nt,
                     TRB_V_gene = su_TCR$TRB_V_gene, 
                     TRB_D_gene = "NA",
                     TRB_J_gene = su_TCR$TRB_J_gene, 
                     TRB_C_gene = su_TCR$TRB_C_gene, 
                     clonotype = su_TCR$clonotype, 
                     d_cond = su_TCR$d_cond,                           
                     sex = su_TCR$sex, 
                     age = su_TCR$age, 
                     race = su_TCR$race)

#su_TCR[is.na(su_TCR)] <- "NA" ## replace all NAs with character NAs

## save preprocessed files as aggregated patients first
write.csv(su_TCR, row.names=F, file = paste0("preprocess_covidtcr/norm_all features/su_TCR_cd8_AllPatients.csv"))

## separate into individual files 
su_TCR <- as.data.frame(su_TCR)
su_TCR <- split(su_TCR, su_TCR$patient_ID)
for (i in 1:length(su_TCR)) {
  write.csv(su_TCR[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/su_cd8/su_TCR_cd8_", unique(su_TCR[[i]]$patient_ID), ".csv"))
}


#### SU CD4+ ----

## run code for su_meta (see above): obtain "su_meta"

## load TCR files
tcr.files <- list.files(path="su/su_TCR_cd4", pattern="*.txt", full.names=TRUE, recursive=FALSE)
su_TCR_cd4 <- list()
for (i in 1:length(tcr.files)){
  tcr <- fread(tcr.files[i])
  ID <- gsub(".*library_","",tcr.files[i]) # get patient ID
  ID <- gsub(".txt", "", ID)
  tcr <- data.frame(patient = ID, barcode = tcr$V1, 
                    TRA_cdr3 = tcr$TRA_1_cdr3, TRA_cdr3_nt = tcr$TRA_1_cdr3_nt, TRA_V_gene=tcr$TRA_1_v_gene, TRA_D_gene=tcr$TRA_1_d_gene, TRA_J_gene=tcr$TRA_1_j_gene, TRA_C_gene=tcr$TRA_1_c_gene, 
                    TRB_cdr3 = tcr$TRB_1_cdr3, TRB_cdr3_nt = tcr$TRB_1_cdr3_nt, TRB_V_gene=tcr$TRB_1_v_gene, TRB_D_gene=tcr$TRB_1_d_gene, TRB_J_gene=tcr$TRB_1_j_gene, TRB_C_gene=tcr$TRB_1_c_gene,
                    # TRA_cdr3_2 = tcr$TRA_2_cdr3, TRA_V_gene_2=tcr$TRA_2_v_gene, TRA_J_gene_2=tcr$TRA_2_j_gene, TRA_C_gene_2=tcr$TRA_2_c_gene, 
                    # TRB_cdr3_2 = tcr$TRB_2_cdr3, TRB_V_gene_2=tcr$TRB_2_v_gene, TRB_J_gene_2=tcr$TRB_2_j_gene, TRB_C_gene_2=tcr$TRB_2_c_gene,
                    clonotype=tcr$clonotype, clonotype_size=tcr$clonotype_size, chain_pairing=tcr$chain_pairing)
  su_TCR_cd4[[i]] <- tcr
  ## extracts only the first TRA and TRB
}
su_TCR_cd4 <- rbindlist(su_TCR_cd4)

## combine TCR with metadata
su_TCR_cd4 <- merge(su_TCR_cd4, su_meta, by="patient")

## process combined dataset 
#su_TCR_cd4 <- subset(su_TCR_cd4, chain_pairing=="Single pair" | chain_pairing=="Extra alpha" | chain_pairing=="Extra beta" | chain_pairing=="Two full chains") **** 1/9/21: KEEP ORPHAN CHAINS!
su_TCR_cd4 <- subset(su_TCR_cd4, chain_pairing=="Orphan alpha" | chain_pairing=="Orphan beta" | chain_pairing=="Single pair" | chain_pairing=="Extra alpha" | chain_pairing=="Extra beta" | chain_pairing=="Two full chains")
su_TCR_cd4 <- data.frame(dataset = "Su",
                     patient_ID = su_TCR_cd4$patient, 
                     patient_subID = su_TCR_cd4$time_point,
                     cell_type = "CD4+",
                     barcode = su_TCR_cd4$barcode,
                     TRA_cdr3 = su_TCR_cd4$TRA_cdr3,
                     TRA_cdr3_nt = su_TCR_cd4$TRA_cdr3_nt,
                     TRA_V_gene = su_TCR_cd4$TRA_V_gene, 
                     TRA_D_gene = su_TCR_cd4$TRA_D_gene,
                     TRA_J_gene = su_TCR_cd4$TRA_J_gene, 
                     TRA_C_gene = su_TCR_cd4$TRA_C_gene, 
                     TRB_cdr3 = su_TCR_cd4$TRB_cdr3, 
                     TRB_cdr3_nt = su_TCR_cd4$TRB_cdr3_nt,
                     TRB_V_gene = su_TCR_cd4$TRB_V_gene, 
                     TRB_D_gene = su_TCR_cd4$TRB_V_gene,
                     TRB_J_gene = su_TCR_cd4$TRB_J_gene, 
                     TRB_C_gene = su_TCR_cd4$TRB_C_gene, 
                     clonotype = su_TCR_cd4$clonotype, 
                     d_cond = su_TCR_cd4$d_cond,                           
                     sex = su_TCR_cd4$sex, 
                     age = su_TCR_cd4$age, 
                     race = su_TCR_cd4$race)

#su_TCR_cd4[is.na(su_TCR_cd4)] <- "NA" ## replace all NAs with character NAs

## save preprocessed files as aggregated patients first
write.csv(su_TCR_cd4, row.names=F, file = paste0("preprocess_covidtcr/norm_all features/su_TCR_cd4_AllPatients.csv"))

## separate into individual files 
su_TCR_cd4 <- as.data.frame(su_TCR_cd4)
su_TCR_cd4 <- split(su_TCR_cd4, su_TCR_cd4$patient_ID)
for (i in 1:length(su_TCR_cd4)) {
  write.csv(su_TCR_cd4[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/su_cd4/su_TCR_cd4_", unique(su_TCR_cd4[[i]]$patient_ID), ".csv"))
}

#### WEN ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses")

wen_meta <- fread("wen/wen_metadata.csv") 
wen_TCR <- fread("wen/wen_TCR.tsv")

## process TRA and TRB separately to keep chain with highest reads
wen_TCR <- split(wen_TCR, wen_TCR$locus) #split by TRA/TRB

#TRA
TRA <- wen_TCR$TRA #only TRA
TRA <- split(TRA, TRA$cell_id) 
for (j in 1:length(TRA)){
  if(nrow(TRA[[j]]) > 1) {
    TRA[[j]] <- TRA[[j]][which.max(TRA[[j]]$reads),]
  }
} #keep TRA chain with higher reads 
TRA <- do.call(rbind, TRA) 

#TRB
TRB <- wen_TCR$TRB #only TRB
TRB <- split(TRB, TRB$cell_id) 
for (j in 1:length(TRB)){
  if(nrow(TRB[[j]]) > 1) {
    TRB[[j]] <- TRB[[j]][which.max(TRB[[j]]$reads),]
  }
} #keep TRB chain with higher reads
TRB <- do.call(rbind, TRB) 

wen_TCR <- rbind(TRA, TRB) 
rm(TRA,TRB)
wen_TCR <- dcast(setDT(wen_TCR), cell_id~locus, value.var=c('junction_aa','junction','v_call', 'j_call', 'c_call', 'repertoire_id'))


## add metadata
wen_TCR$patient_ID_TRA <- wen_TCR$repertoire_id_TRA %>% str_after_nth("-", 1) %>% str_before_nth("-", 1) #extract patient ID
wen_TCR$patient_ID_TRB <- wen_TCR$repertoire_id_TRB %>% str_after_nth("-", 1) %>% str_before_nth("-", 1) #extract patient ID
wen_TCR <- wen_TCR %>% mutate(patient_ID = coalesce(patient_ID_TRA, patient_ID_TRB)) ##ID columns will have NAs, so coalesce here

wen_TCR <- left_join(wen_TCR, wen_meta, by="patient_ID")
wen_TCR$patient_ID <- gsub("Healthy_Control*","HD", wen_TCR$patient_ID)

wen_TCR <- data.frame(dataset = "Wen",
                        patient_ID = wen_TCR$patient_ID, 
                        patient_subID = "NA",
                        cell_type = "T cell",
                        barcode = wen_TCR$cell_id,
                        TRA_cdr3 = wen_TCR$junction_aa_TRA,
                        TRA_cdr3_nt = wen_TCR$junction_TRA,
                        TRA_V_gene = wen_TCR$v_call_TRA, 
                        TRA_D_gene = "NA", 
                        TRA_J_gene = wen_TCR$j_call_TRA, 
                        TRA_C_gene = wen_TCR$c_call_TRA, 
                        TRB_cdr3 = wen_TCR$junction_aa_TRB, 
                        TRB_cdr3_nt = wen_TCR$junction_TRB,
                        TRB_V_gene = wen_TCR$v_call_TRB, 
                        TRB_D_gene = "NA", 
                        TRB_J_gene = wen_TCR$j_call_TRB, 
                        TRB_C_gene = wen_TCR$c_call_TRB, 
                        clonotype = "NA", 
                        d_cond = wen_TCR$disease_cond, 
                        sex = wen_TCR$gender, 
                        age = wen_TCR$age, 
                        race = "NA")

#wen_TCR[is.na(wen_TCR)] <- "NA" ## replace all NAs with character NAs

## save preprocessed files as aggregated patients first
write.csv(wen_TCR, row.names=F, file = paste0("preprocess_covidtcr/norm_all features/wen_TCR_AllPatients.csv"))

## separate into individual files 
wen_TCR <- as.data.frame(wen_TCR)
wen_TCR <- split(wen_TCR, wen_TCR$patient_ID)
for (i in 1:length(wen_TCR)) {
  write.csv(wen_TCR[[i]], row.names=F, file = paste0("preprocess_covidtcr/norm_all features/wen/wen_TCR_", unique(wen_TCR[[i]]$patient_ID), ".csv"))
}


#### VDJ tools ----

### * VDJ: Liao ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

liao.files <- list.files(path="norm_all features/liao", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(liao.files)){
  liao_TCR <- fread(liao.files[i])
  liao_TCR <- subset(liao_TCR, !is.na(TRB_cdr3)) # 1/9/21: get rid of beta == NA, "", or "None"
  
  ID <- liao.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 3)

  liao_VDJ <- data.frame(CDR3nt = liao_TCR$TRB_cdr3_nt, 
                         CDR3aa = liao_TCR$TRB_cdr3,
                         V = liao_TCR$TRB_V_gene,
                         D = liao_TCR$TRB_D_gene,
                         J = liao_TCR$TRB_J_gene)
  liao_VDJ <- ddply(liao_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  liao_VDJ <- data.frame(liao_VDJ$V1, liao_VDJ$V1/sum(liao_VDJ$V1), liao_VDJ[,1:5], rep("NA", dim(liao_VDJ)[1]), rep("NA", dim(liao_VDJ)[1]), rep("NA", dim(liao_VDJ)[1]), rep("NA", dim(liao_VDJ)[1]))
  colnames(liao_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(liao_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/liao/liao_VDJ_", ID, ".txt"))
}

### * VDJ: Zhang ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

zhang.files <- list.files(path="norm_all features/zhang", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(zhang.files)){
  zhang_TCR <- fread(zhang.files[i])
  zhang_TCR <- subset(zhang_TCR, !is.na(TRB_cdr3)) # 1/9/21: get rid of beta == NA, "", or "None"

  ID <- zhang.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 3)
  
  zhang_VDJ <- data.frame(CDR3nt = zhang_TCR$TRB_cdr3_nt, 
                         CDR3aa = zhang_TCR$TRB_cdr3,
                         V = zhang_TCR$TRB_V_gene,
                         D = zhang_TCR$TRB_D_gene,
                         J = zhang_TCR$TRB_J_gene)
  zhang_VDJ <- ddply(zhang_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  zhang_VDJ <- data.frame(zhang_VDJ$V1, zhang_VDJ$V1/sum(zhang_VDJ$V1), zhang_VDJ[,1:5], rep("NA", dim(zhang_VDJ)[1]), rep("NA", dim(zhang_VDJ)[1]), rep("NA", dim(zhang_VDJ)[1]), rep("NA", dim(zhang_VDJ)[1]))
  colnames(zhang_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(zhang_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/zhang/zhang_VDJ_", ID, ".txt"))
}

### * VDJ: Kusnadi ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

kusnadi.files <- list.files(path="norm_all features/kusnadi", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(kusnadi.files)){
  kusnadi_TCR <- fread(kusnadi.files[i])
  kusnadi_TCR <- subset(kusnadi_TCR, TRB_cdr3 != "") # 1/9/21: get rid of beta == NA, "", or "None"
  
  ID <- kusnadi.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 3)
  
  kusnadi_VDJ <- data.frame(CDR3nt = kusnadi_TCR$TRB_cdr3_nt, 
                          CDR3aa = kusnadi_TCR$TRB_cdr3,
                          V = kusnadi_TCR$TRB_V_gene,
                          D = kusnadi_TCR$TRB_D_gene,
                          J = kusnadi_TCR$TRB_J_gene)
  kusnadi_VDJ <- ddply(kusnadi_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  kusnadi_VDJ <- data.frame(kusnadi_VDJ$V1, kusnadi_VDJ$V1/sum(kusnadi_VDJ$V1), kusnadi_VDJ[,1:5], rep("NA", dim(kusnadi_VDJ)[1]), rep("NA", dim(kusnadi_VDJ)[1]), rep("NA", dim(kusnadi_VDJ)[1]), rep("NA", dim(kusnadi_VDJ)[1]))
  colnames(kusnadi_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(kusnadi_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/kusnadi/kusnadi_VDJ_", ID, ".txt"))
}

### * VDJ: Su CD8 ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

su_cd8.files <- list.files(path="norm_all features/su_cd8", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(su_cd8.files)){
  su_cd8_TCR <- fread(su_cd8.files[i])
  su_cd8_TCR <- subset(su_cd8_TCR, TRB_cdr3 != "None") # 1/9/21: get rid of beta == NA, "", or "None"
  
  ID <- su_cd8.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 5)
  
  su_cd8_VDJ <- data.frame(CDR3nt = su_cd8_TCR$TRB_cdr3_nt, 
                            CDR3aa = su_cd8_TCR$TRB_cdr3,
                            V = su_cd8_TCR$TRB_V_gene,
                            D = su_cd8_TCR$TRB_D_gene,
                            J = su_cd8_TCR$TRB_J_gene)
  su_cd8_VDJ <- ddply(su_cd8_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  su_cd8_VDJ <- data.frame(su_cd8_VDJ$V1, su_cd8_VDJ$V1/sum(su_cd8_VDJ$V1), su_cd8_VDJ[,1:5], rep("NA", dim(su_cd8_VDJ)[1]), rep("NA", dim(su_cd8_VDJ)[1]), rep("NA", dim(su_cd8_VDJ)[1]), rep("NA", dim(su_cd8_VDJ)[1]))
  colnames(su_cd8_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(su_cd8_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/su_cd8/su_cd8_VDJ_", ID, ".txt"))
}

### * VDJ: Su CD4 ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

su_cd4.files <- list.files(path="norm_all features/su_cd4", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(su_cd4.files)){
  su_cd4_TCR <- fread(su_cd4.files[i])
  su_cd4_TCR <- subset(su_cd4_TCR, TRB_cdr3 != "None") # 1/9/21: get rid of beta == NA, "", or "None"
  
  ID <- su_cd4.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 3)
  
  su_cd4_VDJ <- data.frame(CDR3nt = su_cd4_TCR$TRB_cdr3_nt, 
                           CDR3aa = su_cd4_TCR$TRB_cdr3,
                           V = su_cd4_TCR$TRB_V_gene,
                           D = su_cd4_TCR$TRB_D_gene,
                           J = su_cd4_TCR$TRB_J_gene)
  su_cd4_VDJ <- ddply(su_cd4_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  su_cd4_VDJ <- data.frame(su_cd4_VDJ$V1, su_cd4_VDJ$V1/sum(su_cd4_VDJ$V1), su_cd4_VDJ[,1:5], rep("NA", dim(su_cd4_VDJ)[1]), rep("NA", dim(su_cd4_VDJ)[1]), rep("NA", dim(su_cd4_VDJ)[1]), rep("NA", dim(su_cd4_VDJ)[1]))
  colnames(su_cd4_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(su_cd4_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/su_cd4/su_cd4_VDJ_", ID, ".txt"))
}

### * VDJ: Wen ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

wen.files <- list.files(path="norm_all features/wen", pattern="*.csv", full.names=TRUE, recursive=FALSE)

for (i in 1:length(wen.files)){
  wen_TCR <- fread(wen.files[i])
  wen_TCR <- subset(wen_TCR, !is.na(TRB_cdr3)) # 1/9/21: get rid of beta == NA, "", or "None"
  
  ID <- wen.files[i] %>% str_before_nth(".csv", 1) %>% str_after_nth("_", 3)
  
  wen_VDJ <- data.frame(CDR3nt = wen_TCR$TRB_cdr3_nt, 
                         CDR3aa = wen_TCR$TRB_cdr3,
                         V = wen_TCR$TRB_V_gene,
                         D = wen_TCR$TRB_D_gene,
                         J = wen_TCR$TRB_J_gene)
  wen_VDJ <- ddply(wen_VDJ,.(CDR3nt, CDR3aa, V, D, J),nrow)
  wen_VDJ <- data.frame(wen_VDJ$V1, wen_VDJ$V1/sum(wen_VDJ$V1), wen_VDJ[,1:5], rep("NA", dim(wen_VDJ)[1]), rep("NA", dim(wen_VDJ)[1]), rep("NA", dim(wen_VDJ)[1]), rep("NA", dim(wen_VDJ)[1]))
  colnames(wen_VDJ) <- c("count", "frequency", "CDR3nt", "CDR3aa", "V", "D", "J", "VEnd", "Dstart", "Dend", "Jstart")
  
  ## export VDJ format 
  write.table(wen_VDJ, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/wen/wen_VDJ_", ID, ".txt"))
}

#### METADATA COMPILATION ----
setwd("C:/Users/stanl/Desktop/tcr data for analyses/preprocess_covidtcr")

## liao
liao_TCR <- fread("norm_all features/liao_TCR_AllPatients.csv")
liao_meta <- liao_TCR[!duplicated(liao_TCR$patient_ID)]
liao_meta <- data.frame(Sample = paste0("liao_VDJ_",liao_meta$patient_ID), 
                        d_cond = liao_meta$d_cond,
                        sex = liao_meta$sex, 
                        age = liao_meta$age, 
                        race = liao_meta$race)
write.table(liao_meta, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/liao_VDJ_metadata.txt"))

## zhang
zhang_TCR <- fread("norm_all features/zhang_TCR_AllPatients.csv")
zhang_meta <- zhang_TCR[!duplicated(zhang_TCR$patient_ID)]
zhang_meta <- data.frame(Sample = paste0("zhang_VDJ_",zhang_meta$patient_ID), 
                        d_cond = zhang_meta$d_cond,
                        sex = zhang_meta$sex, 
                        age = zhang_meta$age, 
                        race = zhang_meta$race)
write.table(zhang_meta, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/zhang_VDJ_metadata.txt"))

## kusnadi
kusnadi_TCR <- fread("norm_all features/kusnadi_TCR_AllPatients.csv")
kusnadi_meta <- kusnadi_TCR[!duplicated(kusnadi_TCR$patient_ID)]
kusnadi_meta <- data.frame(Sample = paste0("kusnadi_VDJ_",kusnadi_meta$patient_ID), 
                        d_cond = kusnadi_meta$d_cond,
                        sex = kusnadi_meta$sex, 
                        age = kusnadi_meta$age, 
                        race = kusnadi_meta$race)
write.table(kusnadi_meta, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/kusnadi_VDJ_metadata.txt"))

## wen
wen_TCR <- fread("norm_all features/wen_TCR_AllPatients.csv")
wen_meta <- wen_TCR[!duplicated(wen_TCR$patient_ID)]
wen_meta <- data.frame(Sample = paste0("wen_VDJ_",wen_meta$patient_ID), 
                        d_cond = wen_meta$d_cond,
                        sex = wen_meta$sex, 
                        age = wen_meta$age, 
                        race = wen_meta$race)
write.table(wen_meta, sep = "\t", row.names=F, file = paste0("norm_for VDJ tools/wen_VDJ_metadata.txt"))


## don't need to do Su et al 



