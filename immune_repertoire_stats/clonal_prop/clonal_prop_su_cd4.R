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
library("ggsignif")
library("dichromat")
set.seed(20)



##### SELECTION ######

meta <- fread("../vdj_datasets/su_cd4_reform/metadata.txt")

meta_HD <- meta[which(meta$`Who Ordinal Scale` == ""),]
meta_mild <- meta[union(union(which(meta$`Who Ordinal Scale` == "1"), which(meta$`Who Ordinal Scale` == "1 or 2")), which(meta$`Who Ordinal Scale` == "2")),]
meta_mod <- meta[union(which(meta$`Who Ordinal Scale` == "3"), which(meta$`Who Ordinal Scale` == "4")),]
meta_severe <- meta[union(union(which(meta$`Who Ordinal Scale` == "5"), which(meta$`Who Ordinal Scale` == "6")), which(meta$`Who Ordinal Scale` == "7")),]
meta_HD$partition <- rep("HD", dim(meta_HD)[1])
meta_mild$partition <- rep("Mild", dim(meta_mild)[1])
meta_mod$partition <- rep("Moderate", dim(meta_mod)[1])
meta_severe$partition <- rep("Severe", dim(meta_severe)[1])


meta_HD <- meta_HD[sample(dim(meta_HD)[1], 16)]
meta_mild <- meta_mild[sample(dim(meta_mild)[1], 16)]
meta_mod <- meta_mod[sample(dim(meta_mod)[1], 16)]
meta_severe <- meta_severe[sample(dim(meta_severe)[1], 16)]

meta_pool <- rbind(meta_HD, meta_mild, meta_mod, meta_severe)
meta <- fwrite(meta_pool, "../vdj_datasets/subdata_clonal_prop/su_cd4/metadata.txt", quote = F, sep = "\t")
meta <- write.table(meta_pool$Sample, "../vdj_datasets/subdata_clonal_prop/su_cd4_sub.txt", quote = F, sep = "\t", row.names = F, col.names = F)




##### RUN ######


immdata <- repLoad("../vdj_datasets/subdata_clonal_prop/su_cd4/")
meta <- fread("../vdj_datasets/subdata_clonal_prop/su_cd4/metadata.txt")


clondata <- immdata$data
clondata <- clondata[match(meta$Sample, names(clondata))]
#imm_top <- repClonality(clondata, .method = "top")#, .head = c(10, 100, 1000, 3000, 30000, 1e+05))
#imm_top_raw <- repClonality(clondata, .method = "top")#, .head = c(10, 100, 1000, 3000, 30000, 1e+05))
imm_top_raw <- repClonality(clondata, .method = "top", .head = c(10, 50, 100, 500, 1000, 3000, 10000))
imm_top <- imm_top_raw
for(i in dim(imm_top)[2]:2){
  imm_top[,i] <- imm_top[,i] - imm_top[,i-1] 
}


#imm_top <- imm_top[match(filedf$Sample, rownames(imm_top)),]
#vis(imm_top)
#vis(imm_top) #+ vis(imm_top, .by = "Status", .meta = immdata$meta)
barstack <- melt(imm_top)
colnames(barstack) <- c("sample", "cloneindex", "prop")
barstack$cloneindex <- as.character(barstack$cloneindex)
barstack$cloneindex[which(barstack$cloneindex == "10")] <- "1:10"
barstack$cloneindex[which(barstack$cloneindex == "50")] <- "11:50"
barstack$cloneindex[which(barstack$cloneindex == "100")] <- "51:100"
barstack$cloneindex[which(barstack$cloneindex == "500")] <- "101:500"
barstack$cloneindex[which(barstack$cloneindex == "1000")] <- "501:1K"
barstack$cloneindex[which(barstack$cloneindex == "3000")] <- "1K:3K"
barstack$cloneindex[which(barstack$cloneindex == "10000")] <- "3K:10K"
barstack$cloneindex <- factor(barstack$cloneindex, levels = unique(barstack$cloneindex))
barstack$sample <- factor(as.character(barstack$sample), levels = meta$Sample)
barstack$sampleindexed <- paste0(meta$partition, seq(1,length(meta$partition)))[match(barstack$sample, meta$Sample)]
barstack$sampleindexed <- factor(as.character(barstack$sampleindexed), levels = unique(barstack$sampleindexed))


colorlist <- colorschemes$BluetoOrange.10
colorlist <- rev(colorlist[-c(1,5,6)])
#colorlist <- colorschemes$BluetoOrange.12
#colorlist <- rev(colorlist[c(1,3,5,7,9,11,12)])
plt <- ggplot(barstack, aes(x=sampleindexed, y=prop, fill=cloneindex)) + geom_bar(stat = "identity", position="stack", color="black") + scale_fill_manual(values =colorlist) #+ scale_fill_manual(values=c("#05668d", "#028090", "#00a896", "#02c39a", "#f0f3bd"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=26, colour = "black"),axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=30), axis.title.x=element_blank(), axis.line = element_blank(), plot.title=element_text(size=30),  legend.title=element_text(size=20), legend.text=element_text(size=15)) 
plt <- plt + ggtitle("Top clonal composition") + labs(x = "Sample", y = "Clonal proportion")
pdf("barstack_top_su_cd4.pdf", width = 14, height = 10)
plt
dev.off()



clondata <- immdata$data
clondata <- clondata[match(meta$Sample, names(clondata))]
imm_rare_raw <- repClonality(clondata, .method = "rare")#, .head = c(10, 100, 1000, 3000, 30000, 1e+05))
imm_rare <- imm_rare_raw
for(i in dim(imm_rare)[2]:2){
  imm_rare[,i] <- imm_rare[,i] - imm_rare[,i-1] 
}
imm_rare <- imm_rare*100
#vis(imm_rare)
barstack <- melt(imm_rare)
colnames(barstack) <- c("sample", "clonecount", "prop")
barstack$clonecount <- as.character(barstack$clonecount)
barstack$clonecount[which(barstack$clonecount == "1")] <- "1"
barstack$clonecount[which(barstack$clonecount == "3")] <- "2:3"
barstack$clonecount[which(barstack$clonecount == "10")] <- "4:10"
barstack$clonecount[which(barstack$clonecount == "30")] <- "11:30"
barstack$clonecount[which(barstack$clonecount == "100")] <- "31:100"
barstack$clonecount[which(barstack$clonecount == "MAX")] <- "101:MAX"
barstack$clonecount <- factor(barstack$clonecount, levels = unique(barstack$clonecount))
barstack$sample <- factor(as.character(barstack$sample), levels = meta$Sample)
barstack$sampleindexed <- paste0(meta$partition, seq(1,length(meta$partition)))[match(barstack$sample, meta$Sample)]
barstack$sampleindexed <- factor(as.character(barstack$sampleindexed), levels = unique(barstack$sampleindexed))

colorlist <- colorschemes$BluetoOrange.10
colorlist <- rev(colorlist[-c(1,5,6,7)])
#colorlist <- colorschemes$BluetoOrange.12
#colorlist <- rev(colorlist[c(1,3,5,7,9,11,12)])
plt <- ggplot(barstack, aes(x=sampleindexed, y=prop, fill=clonecount)) + geom_bar(stat = "identity", position="stack", color="black") + scale_fill_manual(values =colorlist) #+ scale_fill_manual(values=c("#05668d", "#028090", "#00a896", "#02c39a", "#f0f3bd"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=26, colour = "black"),axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=30), axis.title.x=element_blank(), axis.line = element_blank(), plot.title=element_text(size=30),  legend.title=element_text(size=20), legend.text=element_text(size=15)) 
plt <- plt + ggtitle("Rare clonal proportion") + labs(x = "Sample", y = "Occupied repertoire space (%)")
pdf("barstack_rare_su_cd4.pdf", width = 14, height = 10)
plt
dev.off()

