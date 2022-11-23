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
# CD4 analysis total clonotypes
#################### 

# TRIAL
immdata <- repLoad("../vdj_datasets/adaptive/")
meta <- fread("../vdj_datasets/adaptive/metadata.txt")

sampnames <- names(immdata$data)
allprop <- data.frame(cdr3 = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  sampdata <- sampdata[which(sampdata$Proportion > 0.0001),]
  tempdf <- data.frame(cdr3 = sampdata$CDR3.aa, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~cdr3,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  allprop <- merge(x = allprop, y = tempdf, by = "cdr3", all = TRUE)
}
allprop <- allprop[-which(is.na(allprop$cdr3)),]
allprop[is.na(allprop)] <- 0


test <- allprop[,-1]
test <- (test > 0)*1
write.csv(allprop,"allprop_cdr3_adaptive.csv" ,quote = F,  row.names = F)



# From HPC
metaHD <- fread("../vdj_datasets/adaptiveHD/metadata.txt")
allpropHD <- fread("allprop_cdr3_adaptiveHD.csv")


#matchedDzStat <- meta$d_cond[match(colnames(allprop), meta$Sample)]
#matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(allprop), meta$Sample)]

#allprop$meanprop <- rowMeans(allprop[,2:dim(allprop)[2]])


meanprop <- rowMeans(allpropHD[,-1])

allpropHD <- allpropHD[order(meanprop, decreasing = T),]


pltd <- allpropHD
pltd <- pltd[1:10,]
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value, group=variable)) + geom_jitter(size=3, pch=21, colour="black") # + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)

#pltd <- allprop_mod[, -which(colnames(allprop_mod) == "meanprop")]
pltd <- pltd[1:10,]
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value, group=variable)) + geom_jitter(size=3, pch=21, colour="black") + geom_violin()# + scale_fill_manual(values = c("#CD113B", "#053742"))




plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "D50")
pdf("box_d50.pdf", height = 10, width = 10)
plt
dev.off()

png("box_d50.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10d50, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "D50")
plt <- plt + coord_flip()
pdf("box_d50.fo.pdf", height = 10, width = 10)
plt
dev.off()



#imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
#imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)
#p1 <- vis(imm_ov1)
#p2 <- vis(imm_ov2, .text.size = 2)

#p1 + p2
#vis(imm_ov1, "heatmap2")

#repOverlapAnalysis(imm_ov1, "mds")

#repOverlapAnalysis(imm_ov1, "mds") %>% vis()




exp_vol <- repExplore(immdata$data, .method = "volume")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats/clonotype_stats_total_su_cd4.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats/clonotype_stats_clones_su_cd4.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats/clonotype_stats_len_su_cd4.csv" ,quote = F,  row.names = F)


exp_vol <- repExplore(immdata$data, .method = "volume", .col = "aa")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_su_cd4.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones", .col = "aa")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_clones_aa_su_cd4.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_su_cd4.csv" ,quote = F,  row.names = F)



#################### 
# CD4 analysis diversity
#################### 

div_chao <- repDiversity(immdata$data, "chao1")
div_chao <-  cbind(Sample = rownames(div_chao), div_chao, meta[match(rownames(div_chao), meta$Sample),])
write.csv(div_chao,"../immunarch_outputs/immunarch_stats/diversity_stats_chao1_su_cd4.csv" ,quote = F,  row.names = F)

div_ginisimp <- repDiversity(immdata$data, "gini.simp")
div_ginisimp <-  cbind(div_ginisimp, meta[match(div_ginisimp$Sample, meta$Sample),])
write.csv(div_ginisimp,"../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_su_cd4.csv" ,quote = F,  row.names = F)


div_hill <- repDiversity(immdata$data, "hill")
div_hill <-  cbind(div_hill, meta[match(div_hill$Sample, meta$Sample),])
write.csv(div_hill,"../immunarch_outputs/immunarch_stats/diversity_stats_hill_su_cd4.csv" ,quote = F,  row.names = F)


div_d50 <- repDiversity(immdata$data, "d50")
div_d50 <-  cbind(Sample = rownames(div_d50), div_d50, meta[match(rownames(div_d50), meta$Sample),])
write.csv(div_d50,"../immunarch_outputs/immunarch_stats/diversity_stats_d50_su_cd4.csv" ,quote = F,  row.names = F)


div_div <- repDiversity(immdata$data, "div")
div_div <-  cbind(div_div, meta[match(div_div$Sample, meta$Sample),])
write.csv(div_div,"../immunarch_outputs/immunarch_stats/diversity_stats_div_su_cd4.csv" ,quote = F,  row.names = F)


div_gini <- repDiversity(immdata$data, "gini")
div_gini <-  cbind(Sample = rownames(div_gini), div_gini, meta[match(rownames(div_gini), meta$Sample),])
write.csv(div_gini,"../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd4.csv" ,quote = F,  row.names = F)

div_invsimp <- repDiversity(immdata$data, "inv.simp")
div_invsimp <-  cbind(div_invsimp, meta[match(div_invsimp$Sample, meta$Sample),])
write.csv(div_invsimp,"../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_su_cd4.csv" ,quote = F,  row.names = F)



#################### 
# CD4 VDJ usage
#################### 

imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
pdf("../immunarch_outputs/immunarch_figs/geneUsage_immunarch_su_cd4_trbv.pdf", width = 20, height = 10)
print(vis(imm_gu, .by = "Disease status", .meta = immdata$meta))
print(vis(imm_gu, .by = "Disease status", .meta = immdata$meta, .plot = "box"))
dev.off()


imm_gu <- geneUsage(immdata$data, "hs.trbj", .norm = T)
pdf("../immunarch_outputs/immunarch_figs/geneUsage_immunarch_su_cd4_trbj.pdf", width = 20, height = 10)
print(vis(imm_gu, .by = "Disease status", .meta = immdata$meta))
print(vis(imm_gu, .by = "Disease status", .meta = immdata$meta, .plot = "box"))
dev.off()



ov <- repOverlap(immdata$data)
pdf("../immunarch_outputs/immunarch_figs/repOverlap_immunarch_su_cd4.pdf", width = 50, height = 50)
print(vis(ov, "circos"))
dev.off()



#################### 
# CD8 analysis
#################### 

immdata <- repLoad("../vdj_datasets/su_cd8_reform/")
meta <- fread("../vdj_datasets/su_cd8_reform/metadata.txt")

exp_vol <- repExplore(immdata$data, .method = "volume")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats/clonotype_stats_total_su_cd8.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats/clonotype_stats_clones_su_cd8.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats/clonotype_stats_len_su_cd8.csv" ,quote = F,  row.names = F)



exp_vol <- repExplore(immdata$data, .method = "volume", .col = "aa")
exp_vol <-  cbind(exp_vol, meta[match(exp_vol$Sample, meta$Sample),])
write.csv(exp_vol,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_su_cd8.csv" ,quote = F,  row.names = F)

exp_clones <- repExplore(immdata$data, .method = "clones", .col = "aa")
exp_clones <-  cbind(exp_clones, meta[match(exp_clones$Sample, meta$Sample),])
write.csv(exp_clones,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_clones_aa_su_cd8.csv" ,quote = F,  row.names = F)

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_len <-  cbind(exp_len, meta[match(exp_len$Sample, meta$Sample),])
write.csv(exp_len,"../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_su_cd8.csv" ,quote = F,  row.names = F)


#################### 
# CD8 analysis diversity
#################### 

div_chao <- repDiversity(immdata$data, "chao1")
div_chao <-  cbind(Sample = rownames(div_chao), div_chao, meta[match(rownames(div_chao), meta$Sample),])
write.csv(div_chao,"../immunarch_outputs/immunarch_stats/diversity_stats_chao1_su_cd8.csv" ,quote = F,  row.names = F)

div_ginisimp <- repDiversity(immdata$data, "gini.simp")
div_ginisimp <-  cbind(div_ginisimp, meta[match(div_ginisimp$Sample, meta$Sample),])
write.csv(div_ginisimp,"../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_su_cd8.csv" ,quote = F,  row.names = F)

div_hill <- repDiversity(immdata$data, "hill")
div_hill <-  cbind(div_hill, meta[match(div_hill$Sample, meta$Sample),])
write.csv(div_hill,"../immunarch_outputs/immunarch_stats/diversity_stats_hill_su_cd8.csv" ,quote = F,  row.names = F)


div_d50 <- repDiversity(immdata$data, "d50")
div_d50 <-  cbind(Sample = rownames(div_d50), div_d50, meta[match(rownames(div_d50), meta$Sample),])
write.csv(div_d50,"../immunarch_outputs/immunarch_stats/diversity_stats_d50_su_cd8.csv" ,quote = F,  row.names = F)


div_div <- repDiversity(immdata$data, "div")
div_div <-  cbind(div_div, meta[match(div_div$Sample, meta$Sample),])
write.csv(div_div,"../immunarch_outputs/immunarch_stats/diversity_stats_div_su_cd8.csv" ,quote = F,  row.names = F)


div_gini <- repDiversity(immdata$data, "gini")
div_gini <-  cbind(Sample = rownames(div_gini), div_gini, meta[match(rownames(div_gini), meta$Sample),])
write.csv(div_gini,"../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd8.csv" ,quote = F,  row.names = F)

div_invsimp <- repDiversity(immdata$data, "inv.simp")
div_invsimp <-  cbind(div_invsimp, meta[match(div_invsimp$Sample, meta$Sample),])
write.csv(div_invsimp,"../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_su_cd8.csv" ,quote = F,  row.names = F)



#################### 
# CD8 VDJ usage
#################### 

imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)
pdf("../immunarch_outputs/immunarch_figs/geneUsage_immunarch_su_cd8.pdf", width = 20, height = 10)
print(vis(imm_gu, .by = "Disease status", .meta = immdata$meta))
dev.off()


imm_gu <- geneUsage(immdata$data, "hs.trbj", .norm = T)






p1 <- vis(exp_vol, .by = c("Disease status"), .meta = immdata$meta)
p2 <- vis(exp_vol, .by = c("Disease status", "Sex"), .meta = immdata$meta)
pdf("CD4_stats_analyses/clonotype_stats_total_immunarch.pdf", width = 50, height = 20)
print(p1 + p2)
dev.off()

exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

pdf("CD4_stats_analyses/distro_lengths.pdf", width = 150, height = 20)
print(p1)
dev.off()
pdf("CD4_stats_analyses/distro_abundances.pdf", width = 150, height = 20)
print(p2)
dev.off()
pdf("CD4_stats_analyses/distro_clonotypes.pdf", width = 150, height = 20)
print(p3)
dev.off()

p4 <- vis(exp_len, .by = "Disease status", .meta = immdata$meta)
p5 <- vis(exp_cnt, .by = "Disease status", .meta = immdata$meta)
pdf("CD4_stats_analyses/distro_lengths_dz.pdf", width = 150, height = 20)
print(p4)
dev.off()
pdf("CD4_stats_analyses/distro_abundances_dz.pdf", width = 150, height = 20)
print(p5)
dev.off()



imm_pr <- repClonality(immdata$data, .method = "clonal.prop")
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_rare <- repClonality(immdata$data, .method = "rare")



imm_top <- repClonality(immdata$data[sample_ordered], .method = "top", .head = c(10, 100, 1000, 3000, 10000))
#imm_top <- repClonality(immdata$data[sample_ordered], .method = "top", .head = c(5, 10, 50, 100, 1000, 3000, 10000))

p4 <- vis(imm_top) #+ vis(imm_top, .by = "Disease status", .meta = immdata$meta)


imm_rare <- repClonality(immdata$data, .method = "rare")
imm_rare

sample_ordered <-subset(immdata$meta, timepoint == "not applicable")[,1]$Sample

sample_ordered <-c(subset(immdata$meta, timepoint == "not applicable")[,1]$Sample,
                   subset(immdata$meta, timepoint == "T1")[,1]$Sample,
                   subset(immdata$meta, timepoint == "T2")[,1]$Sample)

df <- imm_top[match(rownames(imm_top), sample_ordered),]
#df <- imm_top
dfdiff <- apply(df, 1, diff)
dfdiff <-cbind("10"= df[,1], t(dfdiff))
dfmelt <- melt(dfdiff)
colnames(dfmelt) <- c("sample", "indices", "prop")
dfmelt$indices <- as.character(dfmelt$indices)
dfmelt$sample <- factor(as.character(dfmelt$sample), levels = sample_ordered)
plt <- ggplot(data=dfmelt, aes(x=sample, y=prop, fill=indices)) + geom_bar(stat="identity") + scale_fill_brewer(palette = "Spectral")#+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
#plt <- plt + labs(x = "Sum variant freq. %", y = "Mean crRNA abundance in screen, log2 tpm")
#plt <- plt + ggtitle("Day 54 screen crRNA abundance vs cutting efficiency")
pdf("clon_prop.pdf", height = 20, width = 40)
plt
dev.off()

pdf("test.pdf", width = 50, height = 20)
#print(p1 + p2)
vis(imm_rare) + vis(imm_rare, .by = "Disease status", .meta = immdata$meta)
dev.off()


exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Timepoint"), .meta = immdata$meta)
p2 <- vis(exp_vol, .by = c("Disease status", "Sex"), .meta = immdata$meta)


exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
exp_cnt <- repExplore(immdata$data, .method = "count")
exp_vol <- repExplore(immdata$data, .method = "volume")

p1 <- vis(exp_len)
p2 <- vis(exp_cnt)
p3 <- vis(exp_vol)

p1
p4 <- vis(exp_len, .by = "Disease status", .meta = immdata$meta)
p5 <- vis(exp_cnt, .by = "Sex", .meta = immdata$meta)
p6 <- vis(exp_vol, .by = c("Status", "Sex"), .meta = immdata$meta)

p4
gu <- geneUsage(immdata$data)

imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)
p1 <- vis(imm_ov1, "heatmap2")



pdf("test.pdf", width = 50, height = 20)
#print(p1 + p2)
p4

dev.off()




#################### 
# CD8 processing
#################### 
path = "/Users/jonathan/Desktop/tcrcov/su_process/E-MTAB-9357.processed.4/"
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- fread(paste0("./E-MTAB-9357.processed.4/",file.names[i]))
  if(length(which(file$clonotype == "None")) > 0){
    test <- file[-which(file$clonotype == "None"),]
  }
  if(length(which(test$TRB_2_cdr3 != "None")) > 0){
    test <- test[-which(test$TRB_2_cdr3 != "None"),]
  }
  if(length(which(test$TRB_1_cdr3 == "None")) > 0){
    test <- test[-which(test$TRB_1_cdr3 == "None"),]
  }
  
  process <- data.frame(CDR3nt = test$TRB_1_cdr3_nt, 
                        CDR3aa = test$TRB_1_cdr3,
                        V = test$TRB_1_v_gene,
                        D = test$TRB_1_d_gene,
                        J = test$TRB_1_j_gene)
  temp  <- ddply(process,.(process$CDR3nt, process$CDR3aa, process$V, process$D, process$J),nrow)
  temp <- data.frame(temp$V1, temp$V1/sum(temp$V1), temp[,1:5], rep("NA", dim(temp)[1]), rep("NA", dim(temp)[1]), rep("NA", dim(temp)[1]), rep("NA", dim(temp)[1]))
  colnames(temp) <- c("count", "frequency", colnames(process), "VEnd", "Dstart", "Dend", "Jstart")
  write.table(temp, paste0("./CD8_vdjtools/",file.names[i]), quote = F, sep = "\t", row.names = F, col.names = T)
}



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
write.table(newmeta, "su_meta_full_cd8.txt" ,quote = F,  row.names = F, sep = "\t")
write.table(newmeta[,1:6], "su_meta_sub_cd8.txt" ,quote = F,  row.names = F, sep = "\t")

#################### 
# CD8 analysis
#################### 
immdata <- repLoad("./CD8_vdjtools/")

exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("Disease status"), .meta = immdata$meta)
p2 <- vis(exp_vol, .by = c("Disease status", "Sex"), .meta = immdata$meta)
p1 + p2
pdf("CD8_stats_prism/clonotype_stats_total_immunarch.pdf", width = 50, height = 20)
print(p1 + p2)
dev.off()
exp_vol <-  cbind(exp_vol, newmeta[match(exp_vol$Sample, newmeta$Sample),])
write.csv(exp_vol,"CD8_stats_prism/clonotype_stats_total.csv" ,quote = F,  row.names = F)






healthy <- subset(exp_vol, `Disease status` == "normal")$Volume
dz <- subset(exp_vol, `Disease status` == "COVID-19")$Volume
p <- wilcox.test(dz, healthy, alternative = "two.sided")
p$p.value # 0.0008214435







# VDJtools format
immdata <- repLoad("./test.txt")

ddply(process,.(process$count, process$CDR3nt,
                process$CDR3aa, process$V, process$D, process$J),nrow)
meta <- fread("E-MTAB-9357.sdrf.txt")
input <- input[-which(is.na(input$CloneConsensusInfo)),]



CDR3a <- gsub(".*TRAC:","",input$CloneConsensusInfo)
CDR3a <- gsub(":.*","",CDR3a)
CDR3b <- gsub(".*TRBC.:","",input$CloneConsensusInfo)
CDR3b <- gsub(":.*","",CDR3b)

TRBV <- gsub(".*TRBV","TRBV",input$CloneConsensusInfo)
TRBV <- gsub(":.*","",TRBV)
TRBJ <- gsub(".*TRBJ","TRBJ",input$CloneConsensusInfo)
TRBJ <- gsub(":.*","",TRBJ)

#CDR3a <- sapply(strsplit(input$CloneConsensusInfo, ":"), `[`, 8)    
#CDR3b <- sapply(strsplit(input$CloneConsensusInfo, ":"), `[`, 18)    
#TRBV <- sapply(strsplit(input$CloneConsensusInfo, ":"), `[`, 14) 
#TRBJ <- sapply(strsplit(input$CloneConsensusInfo, ":"), `[`, 16) 
count <- rep(1, length(CDR3b))
sub <- paste0(input$PatientID, ":", input$Condition)

out <- data.frame(CDR3b, TRBV, TRBJ, CDR3a, sub, count)
#out <- out[which(out$count > 1),]
fwrite(out, "zhangtcr.txt", col.names = F, sep = "\t")




####
play <- fread("zhang-et-al_combined_Single-cell-landscape-of-immuno.v3.1.csv")
play2 <- play

higher_count_indices <- which(unlist(lapply(strsplit(play$reads, ";"), order, decreasing = T)) == 1)
play$cell_subtype <- unlist(strsplit(play$cell_subtype, ";"))[higher_count_indices]
play$cell_barcode <- unlist(strsplit(play$cell_barcode, ";"))[higher_count_indices]
play$v_gene <- unlist(strsplit(play$v_gene, ";"))[higher_count_indices]
play$d_gene <- unlist(strsplit(play$d_gene, ";"))[higher_count_indices]
play$j_gene <- unlist(strsplit(play$j_gene, ";"))[higher_count_indices]
play$c_gene <- unlist(strsplit(play$c_gene, ";"))[higher_count_indices]
play$DNA_cdr3 <- unlist(strsplit(play$DNA_cdr3, ";"))[higher_count_indices]
play$AA_cdr3 <- unlist(strsplit(play$AA_cdr3, ";"))[higher_count_indices]
play$reads <- unlist(strsplit(play$reads, ";"))[higher_count_indices]
play$umis <- unlist(strsplit(play$umis, ";"))[higher_count_indices]

#temp <- play[grep(";",play$reads),]
#temp <- play$reads
#lapply(strsplit(play$reads, ";"), order)

play_beta <- subset(play, tcr_chain == "TRB")
play_alpha <- subset(play, tcr_chain == "TRA")



play_beta_conv <- subset(play_beta, disease_cond == "conv")
play_beta_hd <- subset(play_beta, disease_cond == "HD")
play_beta_mod <- subset(play_beta, disease_cond == "moderate")
play_beta_sev <- subset(play_beta, disease_cond == "severe")

play_alpha_conv <- subset(play_alpha, disease_cond == "conv")
play_alpha_hd <- subset(play_alpha, disease_cond == "HD")
play_alpha_mod <- subset(play_alpha, disease_cond == "moderate")
play_alpha_sev <- subset(play_alpha, disease_cond == "severe")



out <- play_beta_conv[,c(2,5,10:15,18)]
out2 <- play_alpha_conv
temp1 <- ddply(out,.(out$patient, out$disease_cond,
                      out$v_gene, out$d_gene, out$j_gene, out$c_gene,
                      out$DNA_cdr3, out$AA_cdr3, out$clonotype),nrow)
temp1 <- cbind(temp1, paste0(temp1$`out$patient`, ":", temp1$`out$disease_cond`))
temp1 <- cbind(temp1, out2$AA_cdr3[match(temp1$`out$clonotype`, out2$clonotype)])
temp3 <- NULL
for(i in unique(temp1$`out$patient`)){
  intermed <- subset(temp1, `out$patient` == i)
  intermed$freqprop <- intermed$V1/sum(intermed$V1)
  temp3 <- rbind(temp3, intermed)
}
#temp2 <- temp1[,c(8,3,5,12,11,10)]
temp2 <- temp3[,c(8,3,5,12,11,13)]
fwrite(temp2, "zhangtcr_v3_conv.txt", col.names = F, sep = "\t")


out <- play_beta_hd[,c(2,5,10:15,18)]
out2 <- play_alpha_hd
temp1 <- ddply(out,.(out$patient, out$disease_cond,
                     out$v_gene, out$d_gene, out$j_gene, out$c_gene,
                     out$DNA_cdr3, out$AA_cdr3, out$clonotype),nrow)
temp1 <- cbind(temp1, paste0(temp1$`out$patient`, ":", temp1$`out$disease_cond`))
temp1 <- cbind(temp1, out2$AA_cdr3[match(temp1$`out$clonotype`, out2$clonotype)])
temp3 <- NULL
for(i in unique(temp1$`out$patient`)){
  intermed <- subset(temp1, `out$patient` == i)
  intermed$freqprop <- intermed$V1/sum(intermed$V1)
  temp3 <- rbind(temp3, intermed)
}
#temp2 <- temp1[,c(8,3,5,12,11,10)]
temp2 <- temp3[,c(8,3,5,12,11,13)]
fwrite(temp2, "zhangtcr_v3_hd.txt", col.names = F, sep = "\t")

out <- play_beta_mod[,c(2,5,10:15,18)]
out2 <- play_alpha_mod
temp1 <- ddply(out,.(out$patient, out$disease_cond,
                     out$v_gene, out$d_gene, out$j_gene, out$c_gene,
                     out$DNA_cdr3, out$AA_cdr3, out$clonotype),nrow)
temp1 <- cbind(temp1, paste0(temp1$`out$patient`, ":", temp1$`out$disease_cond`))
temp1 <- cbind(temp1, out2$AA_cdr3[match(temp1$`out$clonotype`, out2$clonotype)])
temp3 <- NULL
for(i in unique(temp1$`out$patient`)){
  intermed <- subset(temp1, `out$patient` == i)
  intermed$freqprop <- intermed$V1/sum(intermed$V1)
  temp3 <- rbind(temp3, intermed)
}
#temp2 <- temp1[,c(8,3,5,12,11,10)]
temp2 <- temp3[,c(8,3,5,12,11,13)]
fwrite(temp2, "zhangtcr_v3_mod.txt", col.names = F, sep = "\t")

out <- play_beta_sev[,c(2,5,10:15,18)]
out2 <- play_alpha_sev
temp1 <- ddply(out,.(out$patient, out$disease_cond,
                     out$v_gene, out$d_gene, out$j_gene, out$c_gene,
                     out$DNA_cdr3, out$AA_cdr3, out$clonotype),nrow)
temp1 <- cbind(temp1, paste0(temp1$`out$patient`, ":", temp1$`out$disease_cond`))
temp1 <- cbind(temp1, out2$AA_cdr3[match(temp1$`out$clonotype`, out2$clonotype)])
temp3 <- NULL
for(i in unique(temp1$`out$patient`)){
  intermed <- subset(temp1, `out$patient` == i)
  intermed$freqprop <- intermed$V1/sum(intermed$V1)
  temp3 <- rbind(temp3, intermed)
}
#temp2 <- temp1[,c(8,3,5,12,11,10)]
temp2 <- temp3[,c(8,3,5,12,11,13)]
fwrite(temp2, "zhangtcr_v3_sev.txt", col.names = F, sep = "\t")


#agg <-  aggregate(.~id1+id2, out, mean)

#test <- fread("../gliph2/trial1/ref_CD48_v2.0.txt", header = F)

###

a = c(1, 1, 1, 2, 2, 3, 4, 4)
b = c(3.5, 3.5, 2.5, 2, 2, 1, 2.2, 7)
df <-data.frame(a,b)
library(plyr)
ddply(df,.df,nrow)

####

play <- fread("1040_TCRB.tsv")
play2 <- fread("../860011358_TCRB.tsv")
play3 <- fread("1041_TCRB.tsv")


play2 <- play2[-which(play2$productive_frequency == "na"),]




samples <- fread("samples.tsv")
meta <- fread("data/ImmuneCODE-Repertoire-Tags-002.tsv")

samples <- samples[order(samples$sample_name),]
meta <- meta[order(meta$sample_name),]
samples <- cbind(samples, meta[,c(2, 6:9)])




violin <- samples[-which(samples$`Biological Sex` == ""),]
violin <- violin[,-c(1,11:14,16,17)]
violin <- melt(violin,  id.vars = c("Biological Sex"))
plt <- ggplot(violin, aes(x=`Biological Sex`, y=value, fill=`Biological Sex`)) + geom_violin()  + geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#e4508f", "#556fb5"))
plt <- plt + facet_wrap(vars(variable), scales = "free")
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=30), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
pdf(paste0("violin_tcrrep_sex.pdf"), width = 12, height = 8)
print(plt)
dev.off()


#violin <- samples[-which(samples$`Tissue Source` == ""),]
violin <- samples[,-c(1,11:16)]
violin <- melt(violin,  id.vars = c("Tissue Source"))
plt <- ggplot(violin, aes(x=`Tissue Source`, y=value, fill=`Tissue Source`)) + geom_violin()  + geom_boxplot(width=0.1, fill="white")+scale_fill_manual(values=c("#e4508f", "#556fb5"))
plt <- plt + facet_wrap(vars(variable), scales = "free")
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=30), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
pdf(paste0("violin_tcrrep_tissue.pdf"), width = 12, height = 8)
print(plt)
dev.off()

violin <- samples[,-c(1,11:12,14:17)]
violin$Dataset <- substring(violin$Dataset, 10)
violin <- melt(violin,  id.vars = c("Dataset"))
plt <- ggplot(violin, aes(x=`Dataset`, y=value, fill=`Dataset`)) + geom_violin()  + geom_boxplot(width=0.1, fill="white")#+scale_fill_manual(values=c("#e4508f", "#556fb5"))
plt <- plt + facet_wrap(vars(variable), scales = "free")
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=10), axis.title=element_text(size=30), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
pdf(paste0("violin_tcrrep_dataset.pdf"), width = 14, height = 9)
print(plt)
dev.off()



violin <- samples[-which(samples$`Racial Group`== ""),]
violin <- violin[,-c(1,11:15,17)]
#violin$`Racial Group` <- substring(violin$`Racial Group`, 10)
violin <- melt(violin,  id.vars = c("Racial Group"))
plt <- ggplot(violin, aes(x=`Racial Group`, y=value, fill=`Racial Group`)) + geom_violin()  + geom_boxplot(width=0.1, fill="white")#+scale_fill_manual(values=c("#e4508f", "#556fb5"))
plt <- plt + facet_wrap(vars(variable), scales = "free")
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1), axis.text=element_text(size=10), axis.title=element_text(size=30), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
pdf(paste0("violin_tcrrep_race.pdf"), width = 16, height = 11)
print(plt)
dev.off()


violin <- samples[-which(samples$Age == ""),]
violin <- violin[-which(violin$Age == "08 Years"),]
violin <- violin[,-c(1,11:13,15:17)]
violin$Age <- substr(violin$Age, 1,2)
violin$Age <- as.numeric(violin$Age)
b <- seq(10,90,10)
#names <- c("Low", "Medium", "High")
violin$Age <- cut(violin$Age, breaks = b)
violin <- melt(violin,  id.vars = c("Age"))
plt <- ggplot(violin, aes(x=`Age`, y=value, fill=`Age`)) + geom_violin()  + geom_boxplot(width=0.1, fill="white")#+scale_fill_manual(values=c("#e4508f", "#556fb5"))
plt <- plt + facet_wrap(vars(variable), scales = "free")
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.line = element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=10), axis.title=element_text(size=30), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt #+ ylim(c(150,1150))# + ggtitle(paste0(i, " violin plot"))  + scale_y_continuous(limits=c(150, 1100))
pdf(paste0("violin_tcrrep_age.pdf"), width = 16, height = 8)
print(plt)
dev.off()
