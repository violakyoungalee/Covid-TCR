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

file.names <- dir("../vdj_datasets/wen", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/wen/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/wen_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}

file.names <- dir("../vdj_datasets/kusnadi", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/kusnadi/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/kusnadi_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}
file.names <- dir("../vdj_datasets/liao", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/liao/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/liao_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}
file.names <- dir("../vdj_datasets/zhang", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/zhang/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/zhang_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}

file.names <- dir("../vdj_datasets/su_cd4", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/su_cd4/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/su_cd4_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}
file.names <- dir("../vdj_datasets/su_cd8", pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("../vdj_datasets/su_cd8/", file.names[i]))
  write.table(file, paste0("../vdj_datasets/su_cd8_reform/", file.names[i]), quote = F, sep = "\t", row.names = F)
}