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


#################### 
# CD4 analysis total clonotypes
#################### 

immdata <- repLoad("../vdj_datasets/su_cd4_reform/")
meta <- fread("../vdj_datasets/su_cd4_reform/metadata.txt")


sampnames <- names(immdata$data)
vgeneprop <- data.frame(vgene = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(vgene = sampdata$V.name, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~vgene,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  vgeneprop <- merge(x = vgeneprop, y = tempdf, by = "vgene", all = TRUE)
}
vgeneprop <- vgeneprop[-which(is.na(vgeneprop$vgene)),]
vgeneprop[is.na(vgeneprop)] <- 0

write.csv(vgeneprop,"vgeneprop_vdj_su_cd4.csv" ,quote = F,  row.names = F)


matchedDzStat <- meta$`Disease status`[match(colnames(vgeneprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(vgeneprop), meta$Sample)]

# Su CD4 mean prop

#vgeneprop$meanprop <- rowMeans(vgeneprop[,2:dim(vgeneprop)[2]])

reindex <- substr(vgeneprop$vgene, 5, nchar(vgeneprop$vgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)

vgeneprop_vgene <- vgeneprop$vgene[reindex]
vgeneprop_HD <- vgeneprop[reindex,which(matchedWHO == "")]
vgeneprop_mild <- vgeneprop[reindex,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
vgeneprop_mod <- vgeneprop[reindex,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
vgeneprop_severe <- vgeneprop[reindex,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


pltd_HD_melt <- melt(cbind(vgeneprop_vgene, rep("HD",dim(vgeneprop_HD)[1]), vgeneprop_HD))
pltd_mild_melt <- melt(cbind(vgeneprop_vgene, rep("mild",dim(vgeneprop_mild)[1]), vgeneprop_mild))
pltd_mod_melt <- melt(cbind(vgeneprop_vgene, rep("moderate",dim(vgeneprop_mod)[1]), vgeneprop_mod))
pltd_severe_melt <- melt(cbind(vgeneprop_vgene, rep("severe",dim(vgeneprop_severe)[1]), vgeneprop_severe))
colnames(pltd_HD_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_mild_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_mod_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_severe_melt) <- c("vgene", "condition", "sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_mild_melt, pltd_mod_melt, pltd_severe_melt)
pltd$vgene <- factor(pltd$vgene, levels = unique(pltd$vgene))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white", outlier.colour=NA, position=position_dodge(width=0.9))
#plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + coord_flip()
pdf("su_cd4_vgene.pdf", height = 6, width = 20)
plt
dev.off()
png("su_cd4_vgene.png", height = 6, width = 20,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(pch=26, position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
pdf("su_cd4_vgene.fo.pdf", height = 6, width = 20)
plt
dev.off()





# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$vgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,49)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 49)[1]]
  part_mild <- parts[seq(i, length(parts), 49)[2]]
  part_mod <- parts[seq(i, length(parts), 49)[3]]
  part_sev <- parts[seq(i, length(parts), 49)[4]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                      pltd$prop[which(pltd$partition == part_mild)])
  
  wout2 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mod)])
  
  wout3 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_sev)])
  
  pvalouts <- c(wout1$p.value, wout2$p.value, wout3$p.value)
  print(pvalouts)
  pvalues <- cbind(pvalues, pvalouts)
}

dim(pvalues) <- c(1,147)
pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
dim(pvalues) <- c(3,49)

statvals <- as.data.frame(pvalues)
colnames(statvals) <- as.character(unique(pltd$vgene))
rownames(statvals) <- c("mild", "mod", "sev")
write.csv(statvals,"su_cd4_vgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"su_cd4_vgene_wilcoxon_pvalues_sig.csv" ,quote = F,  row.names = T)


vgeneprop_COV <- cbind(vgeneprop_mild, vgeneprop_mod, vgeneprop_severe)
vgeneset <- vgeneprop_vgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(vgeneset)){
  print(i)
  set1 <- as.numeric(vgeneprop_COV[i,-1])
  set2 <- as.numeric(vgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_vgene <- data.frame(vgene = vgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_vgene <- allprop_stats_vgene[order(allprop_stats_vgene$sig, decreasing = T),]
write.csv(allprop_stats_vgene,"su_cd4_allprop_stats_vgene.csv" ,quote = F,  row.names = F)

allprop_stats_vgene$hit <- rep("vgene", dim(allprop_stats_vgene)[1])
allprop_stats_vgene$hit[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc > 0))] <- "COVID-19"
allprop_stats_vgene$hit[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc < 0))] <- "HD"

allprop_stats_vgene$hitlabel <- rep("vgene", dim(allprop_stats_vgene)[1])
allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc > 0))] <- "COVID-19"
#allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$sig > 5), which(allprop_stats_vgene$logfc > 1))] <- "COVID-19"
allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc < 0))] <- "HD"

subdata <- subset(allprop_stats_vgene, hitlabel != "vgene")
plt <- ggplot(allprop_stats_vgene, aes(x=logfc, y=sig, color = hit)) + geom_point(size=4) + scale_colour_manual(values = c( "#11cbd7", "gray")) #"#A239EA", 
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) + labs(y = "beta value", x = "-log10( adj. p-value )")
plt <- plt + geom_hline(yintercept = -log10(0.05), linetype='dashed', color="#C6C8CA")
plt <- plt + geom_point(data=subset(allprop_stats_vgene, hit == "COVID-19"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Enriched"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt + geom_point(data=subset(allprop_stats_vgene, hit == "HD"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Depleted"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt  + geom_text_repel(data=subdata, mapping=aes(logfc, sig,label=vgene), size =6.5, segment.size = 0.2,color="black", box.padding = 0.1, point.padding = 0.1, show.legend = FALSE, seed = 11)
plt <- plt + labs(x = "logFC", y = "-log( adj. p-value )")
plt <- plt + scale_x_continuous(breaks = pretty_breaks(5))
pdf("su_cd4_allprop_stats_vgene_volcano.pdf", height = 8, width = 8.5)
plt
dev.off()




## Su CD4 jgene

sampnames <- names(immdata$data)
jgeneprop <- data.frame(jgene = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(jgene = sampdata$J.name, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~jgene,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  jgeneprop <- merge(x = jgeneprop, y = tempdf, by = "jgene", all = TRUE)
}
jgeneprop <- jgeneprop[-which(is.na(jgeneprop$jgene)),]
jgeneprop[is.na(jgeneprop)] <- 0

write.csv(jgeneprop,"jgeneprop_vdj_su_cd4.csv" ,quote = F,  row.names = F)


matchedDzStat <- meta$`Disease status`[match(colnames(jgeneprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(jgeneprop), meta$Sample)]

# Su CD4 mean prop

#jgeneprop$meanprop <- rowMeans(jgeneprop[,2:dim(jgeneprop)[2]])

reindex <- substr(jgeneprop$jgene, 5, nchar(jgeneprop$jgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)

jgeneprop_jgene <- jgeneprop$jgene[reindex]
jgeneprop_HD <- jgeneprop[reindex,which(matchedWHO == "")]
jgeneprop_mild <- jgeneprop[reindex,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
jgeneprop_mod <- jgeneprop[reindex,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
jgeneprop_severe <- jgeneprop[reindex,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


pltd_HD_melt <- melt(cbind(jgeneprop_jgene, rep("HD",dim(jgeneprop_HD)[1]), jgeneprop_HD))
pltd_mild_melt <- melt(cbind(jgeneprop_jgene, rep("mild",dim(jgeneprop_mild)[1]), jgeneprop_mild))
pltd_mod_melt <- melt(cbind(jgeneprop_jgene, rep("moderate",dim(jgeneprop_mod)[1]), jgeneprop_mod))
pltd_severe_melt <- melt(cbind(jgeneprop_jgene, rep("severe",dim(jgeneprop_severe)[1]), jgeneprop_severe))
colnames(pltd_HD_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_mild_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_mod_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_severe_melt) <- c("jgene", "condition", "sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_mild_melt, pltd_mod_melt, pltd_severe_melt)
pltd$jgene <- factor(pltd$jgene, levels = unique(pltd$jgene))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white", outlier.colour=NA, position=position_dodge(width=0.9))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits = c(0,35))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + ylim(0,35)
pdf("su_cd4_jgene.pdf", height = 6, width = 9)
plt
dev.off()
png("su_cd4_jgene.png", height = 6, width = 9,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(pch=26,position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits = c(0,35))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + ylim(0,35)
pdf("su_cd4_jgene.fo.pdf", height = 6, width = 9)
plt
dev.off()



# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$jgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,13)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 13)[1]]
  part_mild <- parts[seq(i, length(parts), 13)[2]]
  part_mod <- parts[seq(i, length(parts), 13)[3]]
  part_sev <- parts[seq(i, length(parts), 13)[4]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mild)])
  
  wout2 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mod)])
  
  wout3 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_sev)])
  
  pvalouts <- c(wout1$p.value, wout2$p.value, wout3$p.value)
  print(pvalouts)
  pvalues <- cbind(pvalues, pvalouts)
}

dim(pvalues) <- c(1,39)
pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
dim(pvalues) <- c(3,13)

statvals <- as.data.frame(pvalues)
colnames(statvals) <- as.character(unique(pltd$jgene))
rownames(statvals) <- c("mild", "mod", "sev")
write.csv(statvals,"su_cd4_jgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"su_cd4_jgene_wilcoxon_pvalues_sig.csv" ,quote = F,  row.names = T)



jgeneprop_COV <- cbind(jgeneprop_mild, jgeneprop_mod, jgeneprop_severe)
jgeneset <- jgeneprop_jgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(jgeneset)){
  print(i)
  set1 <- as.numeric(jgeneprop_COV[i,-1])
  set2 <- as.numeric(jgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_jgene <- data.frame(jgene = jgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_jgene <- allprop_stats_jgene[order(allprop_stats_jgene$sig, decreasing = T),]
write.csv(allprop_stats_jgene,"su_cd4_allprop_stats_jgene.csv" ,quote = F,  row.names = F)



allprop_stats_jgene$hit <- rep("jgene", dim(allprop_stats_jgene)[1])
allprop_stats_jgene$hit[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc > 0))] <- "COVID-19"
allprop_stats_jgene$hit[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc < 0))] <- "HD"

allprop_stats_jgene$hitlabel <- rep("jgene", dim(allprop_stats_jgene)[1])
allprop_stats_jgene$hitlabel[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc > 0))] <- "COVID-19"
allprop_stats_jgene$hitlabel[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc < 0))] <- "HD"

subdata <- subset(allprop_stats_jgene, hitlabel != "jgene")
plt <- ggplot(allprop_stats_jgene, aes(x=logfc, y=sig, color = hit)) + geom_point(size=4) + scale_colour_manual(values = c("gray"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) + labs(y = "beta value", x = "-log10( adj. p-value )")
plt <- plt + geom_hline(yintercept = -log10(0.05), linetype='dashed', color="#C6C8CA")
plt <- plt + geom_point(data=subset(allprop_stats_jgene, hit == "COVID-19"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Enriched"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt + geom_point(data=subset(allprop_stats_jgene, hit == "HD"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Depleted"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt  + geom_text_repel(data=subdata, mapping=aes(logfc, sig,label=jgene), size =6.5, segment.size = 0.2,color="black", box.padding = 0.1, point.padding = 0.1, show.legend = FALSE, seed = 11)
plt <- plt + labs(x = "logFC", y = "-log( adj. p-value )")
plt <- plt + scale_x_continuous(breaks = pretty_breaks(5))
pdf("su_cd4_allprop_stats_jgene_volcano.pdf", height = 8, width = 8.5)
plt
dev.off()



#################### 
# cd8 analysis total clonotypes
#################### 

immdata <- repLoad("../vdj_datasets/su_cd8_reform/")
meta <- fread("../vdj_datasets/su_cd8_reform/metadata.txt")


sampnames <- names(immdata$data)
vgeneprop <- data.frame(vgene = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(vgene = sampdata$V.name, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~vgene,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  vgeneprop <- merge(x = vgeneprop, y = tempdf, by = "vgene", all = TRUE)
}
vgeneprop <- vgeneprop[-which(is.na(vgeneprop$vgene)),]
vgeneprop[is.na(vgeneprop)] <- 0

write.csv(vgeneprop,"vgeneprop_vdj_su_cd8.csv" ,quote = F,  row.names = F)


matchedDzStat <- meta$`Disease status`[match(colnames(vgeneprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(vgeneprop), meta$Sample)]

# Su cd8 mean prop

#vgeneprop$meanprop <- rowMeans(vgeneprop[,2:dim(vgeneprop)[2]])

reindex <- substr(vgeneprop$vgene, 5, nchar(vgeneprop$vgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)

vgeneprop_vgene <- vgeneprop$vgene[reindex]
vgeneprop_HD <- vgeneprop[reindex,which(matchedWHO == "")]
vgeneprop_mild <- vgeneprop[reindex,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
vgeneprop_mod <- vgeneprop[reindex,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
vgeneprop_severe <- vgeneprop[reindex,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


pltd_HD_melt <- melt(cbind(vgeneprop_vgene, rep("HD",dim(vgeneprop_HD)[1]), vgeneprop_HD))
pltd_mild_melt <- melt(cbind(vgeneprop_vgene, rep("mild",dim(vgeneprop_mild)[1]), vgeneprop_mild))
pltd_mod_melt <- melt(cbind(vgeneprop_vgene, rep("moderate",dim(vgeneprop_mod)[1]), vgeneprop_mod))
pltd_severe_melt <- melt(cbind(vgeneprop_vgene, rep("severe",dim(vgeneprop_severe)[1]), vgeneprop_severe))
colnames(pltd_HD_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_mild_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_mod_melt) <- c("vgene", "condition", "sample", "prop")
colnames(pltd_severe_melt) <- c("vgene", "condition", "sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_mild_melt, pltd_mod_melt, pltd_severe_melt)
pltd$vgene <- factor(pltd$vgene, levels = unique(pltd$vgene))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white", outlier.colour=NA, position=position_dodge(width=0.9))
#plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + coord_flip()
pdf("su_cd8_vgene.pdf", height = 6, width = 20)
plt
dev.off()
png("su_cd8_vgene.png", height = 6, width = 20,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(pch=26, position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
pdf("su_cd8_vgene.fo.pdf", height = 6, width = 20)
plt
dev.off()



# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$vgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,49)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 49)[1]]
  part_mild <- parts[seq(i, length(parts), 49)[2]]
  part_mod <- parts[seq(i, length(parts), 49)[3]]
  part_sev <- parts[seq(i, length(parts), 49)[4]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mild)])
  
  wout2 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mod)])
  
  wout3 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_sev)])
  
  pvalouts <- c(wout1$p.value, wout2$p.value, wout3$p.value)
  print(pvalouts)
  pvalues <- cbind(pvalues, pvalouts)
}

dim(pvalues) <- c(1,147)
pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
dim(pvalues) <- c(3,49)

statvals <- as.data.frame(pvalues)
colnames(statvals) <- as.character(unique(pltd$vgene))
rownames(statvals) <- c("mild", "mod", "sev")
write.csv(statvals,"su_cd8_vgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"su_cd8_vgene_wilcoxon_pvalues_sig.csv" ,quote = F,  row.names = T)


vgeneprop_COV <- cbind(vgeneprop_mild, vgeneprop_mod, vgeneprop_severe)
vgeneset <- vgeneprop_vgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(vgeneset)){
  print(i)
  set1 <- as.numeric(vgeneprop_COV[i,-1])
  set2 <- as.numeric(vgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_vgene <- data.frame(vgene = vgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_vgene <- allprop_stats_vgene[order(allprop_stats_vgene$sig, decreasing = T),]
write.csv(allprop_stats_vgene,"su_cd8_allprop_stats_vgene.csv" ,quote = F,  row.names = F)

allprop_stats_vgene$hit <- rep("vgene", dim(allprop_stats_vgene)[1])
allprop_stats_vgene$hit[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc > 0))] <- "COVID-19"
allprop_stats_vgene$hit[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc < 0))] <- "HD"

allprop_stats_vgene$hitlabel <- rep("vgene", dim(allprop_stats_vgene)[1])
allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc > 0))] <- "COVID-19"
#allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$sig > 5), which(allprop_stats_vgene$logfc > 1))] <- "COVID-19"
allprop_stats_vgene$hitlabel[intersect(which(allprop_stats_vgene$padj < 0.05), which(allprop_stats_vgene$logfc < 0))] <- "HD"

subdata <- subset(allprop_stats_vgene, hitlabel != "vgene")
plt <- ggplot(allprop_stats_vgene, aes(x=logfc, y=sig, color = hit)) + geom_point(size=4) + scale_colour_manual(values = c( "#11cbd7", "gray")) #"#A239EA", 
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) + labs(y = "beta value", x = "-log10( adj. p-value )")
plt <- plt + geom_hline(yintercept = -log10(0.05), linetype='dashed', color="#C6C8CA")
plt <- plt + geom_point(data=subset(allprop_stats_vgene, hit == "COVID-19"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Enriched"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt + geom_point(data=subset(allprop_stats_vgene, hit == "HD"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Depleted"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt  + geom_text_repel(data=subdata, mapping=aes(logfc, sig,label=vgene), size =6.5, segment.size = 0.2,color="black", box.padding = 0.1, point.padding = 0.1, show.legend = FALSE, seed = 11)
plt <- plt + labs(x = "logFC", y = "-log( adj. p-value )")
plt <- plt + scale_x_continuous(breaks = pretty_breaks(5))
pdf("su_cd8_allprop_stats_vgene_volcano.pdf", height = 8, width = 8.5)
plt
dev.off()



## Su cd8 jgene

sampnames <- names(immdata$data)
jgeneprop <- data.frame(jgene = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(jgene = sampdata$J.name, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~jgene,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  jgeneprop <- merge(x = jgeneprop, y = tempdf, by = "jgene", all = TRUE)
}
jgeneprop <- jgeneprop[-which(is.na(jgeneprop$jgene)),]
jgeneprop[is.na(jgeneprop)] <- 0

write.csv(jgeneprop,"jgeneprop_vdj_su_cd8.csv" ,quote = F,  row.names = F)


matchedDzStat <- meta$`Disease status`[match(colnames(jgeneprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(jgeneprop), meta$Sample)]

# Su cd8 mean prop

#jgeneprop$meanprop <- rowMeans(jgeneprop[,2:dim(jgeneprop)[2]])

reindex <- substr(jgeneprop$jgene, 5, nchar(jgeneprop$jgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)

jgeneprop_jgene <- jgeneprop$jgene[reindex]
jgeneprop_HD <- jgeneprop[reindex,which(matchedWHO == "")]
jgeneprop_mild <- jgeneprop[reindex,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
jgeneprop_mod <- jgeneprop[reindex,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
jgeneprop_severe <- jgeneprop[reindex,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


pltd_HD_melt <- melt(cbind(jgeneprop_jgene, rep("HD",dim(jgeneprop_HD)[1]), jgeneprop_HD))
pltd_mild_melt <- melt(cbind(jgeneprop_jgene, rep("mild",dim(jgeneprop_mild)[1]), jgeneprop_mild))
pltd_mod_melt <- melt(cbind(jgeneprop_jgene, rep("moderate",dim(jgeneprop_mod)[1]), jgeneprop_mod))
pltd_severe_melt <- melt(cbind(jgeneprop_jgene, rep("severe",dim(jgeneprop_severe)[1]), jgeneprop_severe))
colnames(pltd_HD_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_mild_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_mod_melt) <- c("jgene", "condition", "sample", "prop")
colnames(pltd_severe_melt) <- c("jgene", "condition", "sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_mild_melt, pltd_mod_melt, pltd_severe_melt)
pltd$jgene <- factor(pltd$jgene, levels = unique(pltd$jgene))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white", outlier.colour=NA, position=position_dodge(width=0.9))
#plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + ylim(0,35)
pdf("su_cd8_jgene.pdf", height = 6, width = 9)
plt
dev.off()
png("su_cd8_jgene.png", height = 6, width = 9,units = "in", res = 300)
plt
dev.off()
plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(pch=26,position=position_jitterdodge(dodge.width=0.9), size=1)#+ scale_color_manual(values = c("#787A91", "#FFC947", "#FF7600", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Percentage")
#plt <- plt + ylim(0,35)
pdf("su_cd8_jgene.fo.pdf", height = 6, width = 9)
plt
dev.off()



# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$jgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,13)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 13)[1]]
  part_mild <- parts[seq(i, length(parts), 13)[2]]
  part_mod <- parts[seq(i, length(parts), 13)[3]]
  part_sev <- parts[seq(i, length(parts), 13)[4]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mild)])
  
  wout2 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_mod)])
  
  wout3 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_sev)])
  
  pvalouts <- c(wout1$p.value, wout2$p.value, wout3$p.value)
  print(pvalouts)
  pvalues <- cbind(pvalues, pvalouts)
}

dim(pvalues) <- c(1,39)
pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
dim(pvalues) <- c(3,13)

statvals <- as.data.frame(pvalues)
colnames(statvals) <- as.character(unique(pltd$jgene))
rownames(statvals) <- c("mild", "mod", "sev")
write.csv(statvals,"su_cd8_jgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"su_cd8_jgene_wilcoxon_pvalues_sig.csv" ,quote = F,  row.names = T)




jgeneprop_COV <- cbind(jgeneprop_mild, jgeneprop_mod, jgeneprop_severe)
jgeneset <- jgeneprop_jgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(jgeneset)){
  print(i)
  set1 <- as.numeric(jgeneprop_COV[i,-1])
  set2 <- as.numeric(jgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_jgene <- data.frame(jgene = jgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_jgene <- allprop_stats_jgene[order(allprop_stats_jgene$sig, decreasing = T),]
write.csv(allprop_stats_jgene,"su_cd8_allprop_stats_jgene.csv" ,quote = F,  row.names = F)


allprop_stats_jgene$hit <- rep("jgene", dim(allprop_stats_jgene)[1])
allprop_stats_jgene$hit[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc > 0))] <- "COVID-19"
allprop_stats_jgene$hit[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc < 0))] <- "HD"

allprop_stats_jgene$hitlabel <- rep("jgene", dim(allprop_stats_jgene)[1])
allprop_stats_jgene$hitlabel[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc > 0))] <- "COVID-19"
allprop_stats_jgene$hitlabel[intersect(which(allprop_stats_jgene$padj < 0.05), which(allprop_stats_jgene$logfc < 0))] <- "HD"

subdata <- subset(allprop_stats_jgene, hitlabel != "jgene")
plt <- ggplot(allprop_stats_jgene, aes(x=logfc, y=sig, color = hit)) + geom_point(size=4) + scale_colour_manual(values = c("gray"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) + labs(y = "beta value", x = "-log10( adj. p-value )")
plt <- plt + geom_hline(yintercept = -log10(0.05), linetype='dashed', color="#C6C8CA")
plt <- plt + geom_point(data=subset(allprop_stats_jgene, hit == "COVID-19"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Enriched"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt + geom_point(data=subset(allprop_stats_jgene, hit == "HD"), size = 2, show.legend = F)# + geom_text_repel(data=subset(volcano, hit == "Depleted"), mapping=aes(beta,lq,label=gene), size = 6, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt  + geom_text_repel(data=subdata, mapping=aes(logfc, sig,label=jgene), size =6.5, segment.size = 0.2,color="black", box.padding = 0.1, point.padding = 0.1, show.legend = FALSE, seed = 11)
plt <- plt + labs(x = "logFC", y = "-log( adj. p-value )")
plt <- plt + scale_x_continuous(breaks = pretty_breaks(5))
pdf("su_cd8_allprop_stats_jgene_volcano.pdf", height = 8, width = 8.5)
plt
dev.off()
