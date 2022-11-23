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

metadata <- fread("SampleOverview_06-03-2021_3-20-33_PM.tsv")
metadata <- metadata[grep("Healthy", metadata$sample_tags),]


racelist <- c("Asian or Pacific Islander", "Hispanic", "Black or African American", "Caucasian")
sexlist <- c("Male", "Female")
#smokelist <- c("Smoker","Non-smoker", "Unknown Smoking Status")


age <- unlist(lapply(strsplit(metadata$sample_tags, ", "), `[[`, 1))
age <- unlist(lapply(strsplit(age, " "), `[[`, 1))

sex <- NULL
for(i in 1:dim(metadata)[1]){
  tags <- unlist(strsplit(metadata$sample_tags, ", ")[i])
  sex <- c(sex, tags[which(tags %in% sexlist)])
}

race <- NULL
for(i in 1:dim(metadata)[1]){
  tags <- unlist(strsplit(metadata$sample_tags, ", ")[i])
  if(identical(tags[which(tags %in% racelist)], character(0))){
    race <- c(race, "NA")
  }else{
    race <- c(race, tags[which(tags %in% racelist)])
  }
}

d_cond <- rep("Healthy", dim(metadata)[1])


d <- data.frame(metadata$sample_name, d_cond, sex, age, race)
colnames(d) <- c("Sample", "d_cond", "sex", "age", "race")
write.table(d, "adaptiveHD_VDJ_metadata.txt", quote = F, sep = "\t", row.names = F)


