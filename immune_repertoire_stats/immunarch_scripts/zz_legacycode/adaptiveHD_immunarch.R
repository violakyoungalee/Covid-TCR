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
library("Seurat")


set.seed(20)


#################### 
# CD4 analysis
#################### 

immdata <- repLoad("../vdj_datasets/adaptiveHD")
meta <- fread("../vdj_datasets/adaptiveHD/metadata.txt")

exp_vol <- repExplore(immdata$data, .method = "volume")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats/clonotype_stats_total_adaptiveHD.csv" ,quote = F,  row.names = F)

#################### 
# CD4 analysis diversity
#################### 

div_chao <- repDiversity(immdata$data, "chao1")
div_chao <-  cbind(Sample = rownames(div_chao), div_chao, meta[match(rownames(div_chao), meta$Sample),])
write.csv(div_chao,"../immunarch_outputs/immunarch_stats/diversity_stats_chao1_adaptiveHD.csv" ,quote = F,  row.names = F)

div_ginisimp <- repDiversity(immdata$data, "gini.simp")
div_ginisimp <-  cbind(div_ginisimp, meta[match(div_ginisimp$Sample, meta$Sample),])
write.csv(div_ginisimp,"../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_adaptiveHD.csv" ,quote = F,  row.names = F)

div_hill <- repDiversity(immdata$data, "hill")
div_hill <-  cbind(div_hill, meta[match(div_hill$Sample, meta$Sample),])
write.csv(div_hill,"../immunarch_outputs/immunarch_stats/diversity_stats_hill_adaptiveHD.csv" ,quote = F,  row.names = F)


div_d50 <- repDiversity(immdata$data, "d50")
div_d50 <-  cbind(Sample = rownames(div_d50), div_d50, meta[match(rownames(div_d50), meta$Sample),])
write.csv(div_d50,"../immunarch_outputs/immunarch_stats/diversity_stats_d50_adaptiveHD.csv" ,quote = F,  row.names = F)


div_div <- repDiversity(immdata$data, "div")
div_div <-  cbind(div_div, meta[match(div_div$Sample, meta$Sample),])
write.csv(div_div,"../immunarch_outputs/immunarch_stats/diversity_stats_div_adaptiveHD.csv" ,quote = F,  row.names = F)



