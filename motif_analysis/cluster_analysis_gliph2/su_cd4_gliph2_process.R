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
library("VennDiagram")
library("ggrepel")

set.seed(20)

immdata <- repLoad("../../vdj_datasets/su_cd4_reform/")
meta <- fread("../../vdj_datasets/su_cd4_reform/metadata.txt")

sampnames <- names(immdata$data)
sconds <- sampnames
sconds <- meta$`Who Ordinal Scale`[match(sconds, meta$Sample)]
sconds[which(sconds == "")] <- "HD"
sconds[union(union(which(sconds == "1"), which(sconds == "1 or 2")), which(sconds == "2"))] <- "mild"
sconds[union(which(sconds == "3"), which(sconds == "4"))] <- "moderate"
sconds[union(union(which(sconds == "5"), which(sconds == "6")), which(sconds == "7"))] <- "severe"

clonedf <- NULL
for(i in 1:length(sampnames)){
  print(i)
  isamp <- sampnames[i]
  icond <- sconds[i]
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(cdr3 = sampdata$CDR3.aa, 
                       clones = sampdata$Proportion, 
                       vgene = sampdata$V.name,
                       jgene = sampdata$J.name)
  tempdf <- aggregate(clones~.,data=tempdf,FUN=sum)
  tempdf <- cbind(tempdf[,1:3], 
                  "trav" = rep(NA, dim(tempdf)[1]), 
                  "cond" = paste0(rep(isamp, dim(tempdf)[1]), ":", rep(icond, dim(tempdf)[1])), 
                  "clones" = tempdf[,4])
  clonedf <- rbind(clonedf, tempdf)
}
write.table(clonedf, "su_cd4_gliph2_clones_all.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(clonedf[grep("HD", clonedf$cond),], "su_cd4_gliph2_clones_hd.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(clonedf[grep("mild", clonedf$cond),], "su_cd4_gliph2_clones_mild.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(clonedf[grep("moderate", clonedf$cond),], "su_cd4_gliph2_clones_moderate.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(clonedf[grep("severe", clonedf$cond),], "su_cd4_gliph2_clones_severe.txt", quote = F, row.names = F, col.names = F, sep = "\t")
