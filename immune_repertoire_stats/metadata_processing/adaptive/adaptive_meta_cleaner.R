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

metadata <- fread("ImmuneCODE-Repertoire-Tags-002.2.tsv")

d <- metadata[,c(3, 5, 7, 6, 8)]
colnames(d) <- c("Sample", "d_cond", "sex", "age", "race")
write.table(d, "adaptive_VDJ_metadata_reformatted.txt", quote = F, sep = "\t", row.names = F)




write.table(kusnadi, "kusnadi_VDJ_metadata_reformatted.txt", quote = F, sep = "\t", row.names = F)
lia <- fread("liao_VDJ_metadata.txt")
write.table(lia, "liao_VDJ_metadata_reformatted.txt", quote = F, sep = "\t", row.names = F)
zhang <- fread("zhang_VDJ_metadata.txt")
write.table(zhang, "zhang_VDJ_metadata_reformatted.txt", quote = F, sep = "\t", row.names = F)
wen <- fread("wen_VDJ_metadata.txt")
write.table(wen, "wen_VDJ_metadata_reformatted.txt", quote = F, sep = "\t", row.names = F)