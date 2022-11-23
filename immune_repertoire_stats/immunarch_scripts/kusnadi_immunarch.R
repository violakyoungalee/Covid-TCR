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

immdata <- repLoad("../vdj_datasets/kusnadi_reform")
meta <- fread("../vdj_datasets/kusnadi_reform/metadata.txt")

exp_vol <- repExplore(immdata$data, .method = "volume")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats/clonotype_stats_total_kusnadi.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats/clonotype_stats_clones_kusnadi.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats/clonotype_stats_len_kusnadi.csv" ,quote = F,  row.names = F)


exp_vol <- repExplore(immdata$data, .method = "volume", .col = "aa")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_kusnadi.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones", .col = "aa")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_clones_aa_kusnadi.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_kusnadi.csv" ,quote = F,  row.names = F)



#################### 
# CD4 analysis diversity
#################### 

div_chao <- repDiversity(immdata$data, "chao1")
div_chao <-  cbind(Sample = rownames(div_chao), div_chao, meta[match(rownames(div_chao), meta$Sample),])
write.csv(div_chao,"../immunarch_outputs/immunarch_stats/diversity_stats_chao1_kusnadi.csv" ,quote = F,  row.names = F)

div_ginisimp <- repDiversity(immdata$data, "gini.simp")
div_ginisimp <-  cbind(div_ginisimp, meta[match(div_ginisimp$Sample, meta$Sample),])
write.csv(div_ginisimp,"../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_kusnadi.csv" ,quote = F,  row.names = F)


div_hill <- repDiversity(immdata$data, "hill")
div_hill <-  cbind(div_hill, meta[match(div_hill$Sample, meta$Sample),])
write.csv(div_hill,"../immunarch_outputs/immunarch_stats/diversity_stats_hill_kusnadi.csv" ,quote = F,  row.names = F)


div_d50 <- repDiversity(immdata$data, "d50")
div_d50 <-  cbind(Sample = rownames(div_d50), div_d50, meta[match(rownames(div_d50), meta$Sample),])
write.csv(div_d50,"../immunarch_outputs/immunarch_stats/diversity_stats_d50_kusnadi.csv" ,quote = F,  row.names = F)


div_div <- repDiversity(immdata$data, "div")
div_div <-  cbind(div_div, meta[match(div_div$Sample, meta$Sample),])
write.csv(div_div,"../immunarch_outputs/immunarch_stats/diversity_stats_div_kusnadi.csv" ,quote = F,  row.names = F)


div_gini <- repDiversity(immdata$data, "gini")
div_gini <-  cbind(Sample = rownames(div_gini), div_gini, meta[match(rownames(div_gini), meta$Sample),])
write.csv(div_gini,"../immunarch_outputs/immunarch_stats/diversity_stats_gini_kusnadi.csv" ,quote = F,  row.names = F)

div_invsimp <- repDiversity(immdata$data, "inv.simp")
div_invsimp <-  cbind(div_invsimp, meta[match(div_invsimp$Sample, meta$Sample),])
write.csv(div_invsimp,"../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_kusnadi.csv" ,quote = F,  row.names = F)

