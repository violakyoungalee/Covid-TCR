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
allprop <- data.frame(cdr3 = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(cdr3 = sampdata$CDR3.aa, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~cdr3,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  allprop <- merge(x = allprop, y = tempdf, by = "cdr3", all = TRUE)
}
allprop <- allprop[-which(is.na(allprop$cdr3)),]
allprop[is.na(allprop)] <- 0
write.csv(allprop,"allprop_cdr3_su_cd4.csv" ,quote = F,  row.names = F)




# Su CD4 mean prop

# Scatter
allprop <- fread("allprop_cdr3_su_cd4.csv")
allprop <- as.data.frame(allprop)

matchedDzStat <- meta$`Disease status`[match(colnames(allprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(allprop), meta$Sample)]
#temp1 <- apply(allprop[,-1],1,var)
#temp2 <- rowMeans(allprop[,-1])
allprop <- allprop[which(rowSums(allprop[,-1]) > 0.005),]

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$meanprop <- rowMeans(allprop_HD)
allprop_mild$meanprop <- rowMeans(allprop_mild)
allprop_mod$meanprop <- rowMeans(allprop_mod)
allprop_severe$meanprop <- rowMeans(allprop_severe)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3


q <- data.frame(cdr3 = allprop_mild$cdr3 ,
                diff = (allprop_mild$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Mild"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point( size = 2) + scale_colour_manual(values = c("gray","#11cbd7","#F0A500"))# + scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Mild"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Mild"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5))
pdf("waterfall_cd4_mildvshd.pdf", height = 6, width = 10)
plt
dev.off()


q <- data.frame(cdr3 = allprop_mod$cdr3 ,
                diff = (allprop_mod$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Moderate"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point(size = 2) + scale_colour_manual(values = c("gray","#11cbd7","#FF5200"))#+ scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Moderate"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Moderate"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5))
pdf("waterfall_cd4_modvshd.pdf", height = 6, width = 10)
plt
dev.off()

q <- data.frame(cdr3 = allprop_severe$cdr3 ,
                diff = (allprop_severe$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Severe"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point(size = 2)+ scale_colour_manual(values = c("gray","#11cbd7","#810000")) #+ scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Severe"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Severe"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5))
#plt <- plt +coord_cartesian(clip = "off")
pdf("waterfall_cd4_severevshd.pdf", height = 6, width = 10)
plt
dev.off()


#allprop$meanprop <- rowMeans(allprop[,2:dim(allprop)[2]])

allprop <- fread("allprop_cdr3_su_cd4.csv")
allprop <- as.data.frame(allprop)

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$meanprop <- rowMeans(allprop_HD)
allprop_mild$meanprop <- rowMeans(allprop_mild)
allprop_mod$meanprop <- rowMeans(allprop_mod)
allprop_severe$meanprop <- rowMeans(allprop_severe)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3

allprop_HD <- allprop_HD[order(allprop_HD$meanprop, decreasing = T),]
allprop_mild <- allprop_mild[order(allprop_mild$meanprop, decreasing = T),]
allprop_mod <- allprop_mod[order(allprop_mod$meanprop, decreasing = T),]
allprop_severe <- allprop_severe[order(allprop_severe$meanprop, decreasing = T),]

#plot(allprop_severe$meanprop, allprop_mod$meanprop)



write.csv(allprop_HD,"allprop_cdr3_su_cd4_processed_HD.csv" ,quote = F,  row.names = F)
write.csv(allprop_mild,"allprop_cdr3_su_cd4_processed_mild.csv" ,quote = F,  row.names = F)
write.csv(allprop_mod,"allprop_cdr3_su_cd4_processed_mod.csv" ,quote = F,  row.names = F)
write.csv(allprop_severe,"allprop_cdr3_su_cd4_processed_severe.csv" ,quote = F,  row.names = F)




ns <- 15
pltd_HD <- data.frame(cdr3 = allprop_HD$cdr3, meanprop = allprop_HD[, which(colnames(allprop_HD) == "meanprop")])
pltd_HD$meanprop <- pltd_HD$meanprop*100
pltd_HD <- pltd_HD[1:ns,]
pltd_HD$group <- rep("HD", ns)

pltd_mild <- data.frame(cdr3 = allprop_mild$cdr3, meanprop = allprop_mild[, which(colnames(allprop_mild) == "meanprop")])
pltd_mild$meanprop <- pltd_mild$meanprop*100
pltd_mild <- pltd_mild[1:ns,]
pltd_mild$group <- rep("mild", ns)

pltd_mod <- data.frame(cdr3 = allprop_mod$cdr3, meanprop = allprop_mod[, which(colnames(allprop_mod) == "meanprop")])
pltd_mod$meanprop <- pltd_mod$meanprop*100
pltd_mod <- pltd_mod[1:ns,]
pltd_mod$group <- rep("moderate", ns)

pltd_severe <- data.frame(cdr3 = allprop_severe$cdr3, meanprop = allprop_severe[, which(colnames(allprop_severe) == "meanprop")])
pltd_severe$meanprop <- pltd_severe$meanprop*100
pltd_severe <- pltd_severe[1:ns,]
pltd_severe$group <- rep("severe", ns)

pltd <- pltd_HD
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#B2B1B9") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Healthy Donor")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_su_cd4_mean_hd.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mild
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FFF9B0") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Mild")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_su_cd4_mean_mild.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mod
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FB9300") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Moderate")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_su_cd4_mean_moderate.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_severe
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#CD113B") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Severe")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_su_cd4_mean_severe.pdf", height = 8, width = 6)
plt
dev.off()




# Venn diagrams
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0)]))


pdf("venn_cdr3_su_cd4_meanprop_0threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()


pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.00001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.00001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.00001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.00001)]))

pdf("venn_cdr3_su_cd4_meanprop_0.00001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.0001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.0001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.0001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.0001)]))


pdf("venn_cdr3_su_cd4_meanprop_0.0001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

write.table(setdiff(intersect(intersect(pattern_mild, pattern_mod), pattern_sev), pattern_hd), "venn_cdr3_su_cd4_meanprop_0.0001threshold_setdiffHD-intersectDZ.txt", row.names = F, col.names = F, quote = F)




# Su CD4 max prop

#allprop$maxprop <- rowmaxs(allprop[,2:dim(allprop)[2]])

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$maxprop <- apply(allprop_HD, 1, FUN=max)
allprop_mild$maxprop <- apply(allprop_mild, 1, FUN=max)
allprop_mod$maxprop <- apply(allprop_mod, 1, FUN=max)
allprop_severe$maxprop <- apply(allprop_severe, 1, FUN=max)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3

allprop_HD <- allprop_HD[order(allprop_HD$maxprop, decreasing = T),]
allprop_mild <- allprop_mild[order(allprop_mild$maxprop, decreasing = T),]
allprop_mod <- allprop_mod[order(allprop_mod$maxprop, decreasing = T),]
allprop_severe <- allprop_severe[order(allprop_severe$maxprop, decreasing = T),]


ns <- 15
ymaxval <- 25
pltd_HD <- data.frame(cdr3 = allprop_HD$cdr3, maxprop = allprop_HD[, which(colnames(allprop_HD) == "maxprop")])
pltd_HD$maxprop <- pltd_HD$maxprop*100
pltd_HD <- pltd_HD[1:ns,]
pltd_HD$group <- rep("HD", ns)

pltd_mild <- data.frame(cdr3 = allprop_mild$cdr3, maxprop = allprop_mild[, which(colnames(allprop_mild) == "maxprop")])
pltd_mild$maxprop <- pltd_mild$maxprop*100
pltd_mild <- pltd_mild[1:ns,]
pltd_mild$group <- rep("mild", ns)

pltd_mod <- data.frame(cdr3 = allprop_mod$cdr3, maxprop = allprop_mod[, which(colnames(allprop_mod) == "maxprop")])
pltd_mod$maxprop <- pltd_mod$maxprop*100
pltd_mod <- pltd_mod[1:ns,]
pltd_mod$group <- rep("moderate", ns)

pltd_severe <- data.frame(cdr3 = allprop_severe$cdr3, maxprop = allprop_severe[, which(colnames(allprop_severe) == "maxprop")])
pltd_severe$maxprop <- pltd_severe$maxprop*100
pltd_severe <- pltd_severe[1:ns,]
pltd_severe$group <- rep("severe", ns)

pltd <- pltd_HD
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#B2B1B9") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Healthy Donor")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/2)
pdf("cdr3_su_cd4_max_hd.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mild
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FFF9B0") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Mild")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/2)
pdf("cdr3_su_cd4_max_mild.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mod
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FB9300") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Moderate")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/2)
pdf("cdr3_su_cd4_max_moderate.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_severe
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#CD113B") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Severe")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/2)
pdf("cdr3_su_cd4_max_severe.pdf", height = 8, width = 6)
plt
dev.off()

# Venn diagrams max
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0)]))


pdf("venn_cdr3_su_cd4_maxprop_0threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()


pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0.001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0.001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0.001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0.001)]))

pdf("venn_cdr3_su_cd4_maxprop_0.001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0.01)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0.01)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0.01)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0.01)]))


pdf("venn_cdr3_su_cd4_maxprop_0.01threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

write.table(setdiff(intersect(intersect(pattern_mild, pattern_mod), pattern_sev), pattern_hd), "venn_cdr3_su_cd4_maxprop_0.01threshold_setdiffHD-intersectDZ.txt", row.names = F, col.names = F, quote = F)







#################### 
# CD8 analysis total clonotypes
#################### 

immdata <- repLoad("../vdj_datasets/su_cd8_reform/")
meta <- fread("../vdj_datasets/su_cd8_reform/metadata.txt")

sampnames <- names(immdata$data)
allprop <- data.frame(cdr3 = c(NA))
for(i in 1:length(sampnames)){
  print(i)
  sampdata <- immdata$data[[i]]
  tempdf <- data.frame(cdr3 = sampdata$CDR3.aa, prop = sampdata$Proportion)
  tempdf <- aggregate(prop~cdr3,data=tempdf,FUN=sum)
  colnames(tempdf)[2] <- sampnames[i]
  allprop <- merge(x = allprop, y = tempdf, by = "cdr3", all = TRUE)
}
allprop <- allprop[-which(is.na(allprop$cdr3)),]
allprop[is.na(allprop)] <- 0

write.csv(allprop,"allprop_cdr3_su_cd8.csv" ,quote = F,  row.names = F)



# Su CD8 mean prop



# Scatter
allprop <- fread("allprop_cdr3_su_cd8.csv")
allprop <- as.data.frame(allprop)

matchedDzStat <- meta$`Disease status`[match(colnames(allprop), meta$Sample)]
matchedWHO <- meta$`Who Ordinal Scale`[match(colnames(allprop), meta$Sample)]
#temp1 <- apply(allprop[,-1],1,var)
#temp2 <- rowMeans(allprop[,-1])
allprop <- allprop[which(rowSums(allprop[,-1]) > 0.005),]

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$meanprop <- rowMeans(allprop_HD)
allprop_mild$meanprop <- rowMeans(allprop_mild)
allprop_mod$meanprop <- rowMeans(allprop_mod)
allprop_severe$meanprop <- rowMeans(allprop_severe)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3


q <- data.frame(cdr3 = allprop_mild$cdr3 ,
                diff = (allprop_mild$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Mild"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point( size = 2) + scale_colour_manual(values = c("gray","#11cbd7","#F0A500"))# + scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Mild"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE,nudge_y = -0.001, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Mild"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_y = 0.007,nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5), limits = c(-0.027, 0.02))
pdf("waterfall_cd8_mildvshd.pdf", height = 6, width = 10)
plt
dev.off()


q <- data.frame(cdr3 = allprop_mod$cdr3 ,
                diff = (allprop_mod$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Moderate"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point(size = 2) + scale_colour_manual(values = c("gray","#11cbd7","#FF5200"))#+ scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Moderate"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_y = -0.001,nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Moderate"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_y = 0.007,nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5), limits = c(-0.027, 0.02))
pdf("waterfall_cd8_modvshd.pdf", height = 6, width = 10)
plt
dev.off()

q <- data.frame(cdr3 = allprop_severe$cdr3 ,
                diff = (allprop_severe$meanprop - allprop_HD$meanprop))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "Severe"
plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point(size = 2)+ scale_colour_manual(values = c("gray","#11cbd7","#810000")) #+ scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "Severe"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_y = -0.002, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "Severe"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE,nudge_y = 0.005,  nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5), limits = c(-0.027, 0.027))
#plt <- plt +coord_cartesian(clip = "off")
pdf("waterfall_cd8_severevshd.pdf", height = 6, width = 10)
plt
dev.off()

#allprop$meanprop <- rowMeans(allprop[,2:dim(allprop)[2]])

allprop <- fread("allprop_cdr3_su_cd8.csv")
allprop <- as.data.frame(allprop)

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$meanprop <- rowMeans(allprop_HD)
allprop_mild$meanprop <- rowMeans(allprop_mild)
allprop_mod$meanprop <- rowMeans(allprop_mod)
allprop_severe$meanprop <- rowMeans(allprop_severe)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3


allprop_HD <- allprop_HD[order(allprop_HD$meanprop, decreasing = T),]
allprop_mild <- allprop_mild[order(allprop_mild$meanprop, decreasing = T),]
allprop_mod <- allprop_mod[order(allprop_mod$meanprop, decreasing = T),]
allprop_severe <- allprop_severe[order(allprop_severe$meanprop, decreasing = T),]

write.csv(allprop_HD,"allprop_cdr3_su_cd8_processed_HD.csv" ,quote = F,  row.names = F)
write.csv(allprop_mild,"allprop_cdr3_su_cd8_processed_mild.csv" ,quote = F,  row.names = F)
write.csv(allprop_mod,"allprop_cdr3_su_cd8_processed_mod.csv" ,quote = F,  row.names = F)
write.csv(allprop_severe,"allprop_cdr3_su_cd8_processed_severe.csv" ,quote = F,  row.names = F)



ns <- 15
ymaxval <- 3
pltd_HD <- data.frame(cdr3 = allprop_HD$cdr3, meanprop = allprop_HD[, which(colnames(allprop_HD) == "meanprop")])
pltd_HD$meanprop <- pltd_HD$meanprop*100
pltd_HD <- pltd_HD[1:ns,]
pltd_HD$group <- rep("HD", ns)

pltd_mild <- data.frame(cdr3 = allprop_mild$cdr3, meanprop = allprop_mild[, which(colnames(allprop_mild) == "meanprop")])
pltd_mild$meanprop <- pltd_mild$meanprop*100
pltd_mild <- pltd_mild[1:ns,]
pltd_mild$group <- rep("mild", ns)

pltd_mod <- data.frame(cdr3 = allprop_mod$cdr3, meanprop = allprop_mod[, which(colnames(allprop_mod) == "meanprop")])
pltd_mod$meanprop <- pltd_mod$meanprop*100
pltd_mod <- pltd_mod[1:ns,]
pltd_mod$group <- rep("moderate", ns)

pltd_severe <- data.frame(cdr3 = allprop_severe$cdr3, meanprop = allprop_severe[, which(colnames(allprop_severe) == "meanprop")])
pltd_severe$meanprop <- pltd_severe$meanprop*100
pltd_severe <- pltd_severe[1:ns,]
pltd_severe$group <- rep("severe", ns)

pltd <- pltd_HD
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#B2B1B9") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Healthy Donor")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 4)
pdf("cdr3_su_cd8_mean_hd.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mild
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FFF9B0") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Mild")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 4)
pdf("cdr3_su_cd8_mean_mild.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mod
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FB9300") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Moderate")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 4)
pdf("cdr3_su_cd8_mean_moderate.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_severe
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#CD113B") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Severe")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 4)
pdf("cdr3_su_cd8_mean_severe.pdf", height = 8, width = 6)
plt
dev.off()

# Venn diagrams
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0)]))


pdf("venn_cdr3_su_cd8_meanprop_0threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()


pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.00001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.00001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.00001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.00001)]))

pdf("venn_cdr3_su_cd8_meanprop_0.00001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.0001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.0001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.0001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.0001)]))


pdf("venn_cdr3_su_cd8_meanprop_0.0001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

write.table(setdiff(intersect(intersect(pattern_mild, pattern_mod), pattern_sev), pattern_hd), "venn_cdr3_su_cd8_meanprop_0.0001threshold_setdiffHD-intersectDZ.txt", row.names = F, col.names = F, quote = F)




# Su CD8 max prop

#allprop$maxprop <- rowmaxs(allprop[,2:dim(allprop)[2]])

allprop_HD <- allprop[,which(matchedWHO == "")]
allprop_mild <- allprop[,union(union(which(matchedWHO == "1"), which(matchedWHO == "1 or 2")), which(matchedWHO == "2"))]
allprop_mod <- allprop[,union(which(matchedWHO == "3"), which(matchedWHO == "4"))]
allprop_severe <- allprop[,union(union(which(matchedWHO == "5"), which(matchedWHO == "6")), which(matchedWHO == "7"))]


allprop_HD$maxprop <- apply(allprop_HD, 1, FUN=max)
allprop_mild$maxprop <- apply(allprop_mild, 1, FUN=max)
allprop_mod$maxprop <- apply(allprop_mod, 1, FUN=max)
allprop_severe$maxprop <- apply(allprop_severe, 1, FUN=max)

allprop_HD$cdr3 <- allprop$cdr3
allprop_mild$cdr3 <- allprop$cdr3
allprop_mod$cdr3 <- allprop$cdr3
allprop_severe$cdr3 <- allprop$cdr3

allprop_HD <- allprop_HD[order(allprop_HD$maxprop, decreasing = T),]
allprop_mild <- allprop_mild[order(allprop_mild$maxprop, decreasing = T),]
allprop_mod <- allprop_mod[order(allprop_mod$maxprop, decreasing = T),]
allprop_severe <- allprop_severe[order(allprop_severe$maxprop, decreasing = T),]


ns <- 15
ymaxval <- 50
pltd_HD <- data.frame(cdr3 = allprop_HD$cdr3, maxprop = allprop_HD[, which(colnames(allprop_HD) == "maxprop")])
pltd_HD$maxprop <- pltd_HD$maxprop*100
pltd_HD <- pltd_HD[1:ns,]
pltd_HD$group <- rep("HD", ns)

pltd_mild <- data.frame(cdr3 = allprop_mild$cdr3, maxprop = allprop_mild[, which(colnames(allprop_mild) == "maxprop")])
pltd_mild$maxprop <- pltd_mild$maxprop*100
pltd_mild <- pltd_mild[1:ns,]
pltd_mild$group <- rep("mild", ns)

pltd_mod <- data.frame(cdr3 = allprop_mod$cdr3, maxprop = allprop_mod[, which(colnames(allprop_mod) == "maxprop")])
pltd_mod$maxprop <- pltd_mod$maxprop*100
pltd_mod <- pltd_mod[1:ns,]
pltd_mod$group <- rep("moderate", ns)

pltd_severe <- data.frame(cdr3 = allprop_severe$cdr3, maxprop = allprop_severe[, which(colnames(allprop_severe) == "maxprop")])
pltd_severe$maxprop <- pltd_severe$maxprop*100
pltd_severe <- pltd_severe[1:ns,]
pltd_severe$group <- rep("severe", ns)

pltd <- pltd_HD
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#B2B1B9") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Healthy Donor")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/5)
pdf("cdr3_su_cd8_max_hd.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mild
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FFF9B0") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Mild")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/5)
pdf("cdr3_su_cd8_max_mild.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_mod
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#FB9300") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Moderate")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/5)
pdf("cdr3_su_cd8_max_moderate.pdf", height = 8, width = 6)
plt
dev.off()

pltd <- pltd_severe
cdr3list <- pltd$cdr3
pltd <- melt(pltd)
pltd$cdr3 <- factor(pltd$cdr3, levels = cdr3list)
plt <- ggplot(pltd, aes(x=cdr3, y=value)) + geom_bar(stat = "identity", color="black", fill="#CD113B") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Max CDR3 usage (%)")
plt <- plt + ggtitle("Severe")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,ymaxval) + coord_fixed(ratio = 1/5)
pdf("cdr3_su_cd8_max_severe.pdf", height = 8, width = 6)
plt
dev.off()

# Venn diagrams max
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0)]))


pdf("venn_cdr3_su_cd8_maxprop_0threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()


pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0.001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0.001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0.001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0.001)]))

pdf("venn_cdr3_su_cd8_maxprop_0.001threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$maxprop > 0.01)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$maxprop > 0.01)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$maxprop > 0.01)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$maxprop > 0.01)]))


pdf("venn_cdr3_su_cd8_maxprop_0.01threshold.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_sev),
  area3 = length(pattern_mild),
  area4 = length(pattern_mod),
  n12 = length(intersect(pattern_hd, pattern_sev)),
  n13 = length(intersect(pattern_hd, pattern_mild)),
  n14 = length(intersect(pattern_hd, pattern_mod)),
  n23 = length(intersect(pattern_sev, pattern_mild)),
  n24 = length(intersect(pattern_sev, pattern_mod)),
  n34 = length(intersect(pattern_mild, pattern_mod)),
  n123 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mild)),
  n124 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_mod)),
  n134 = length(intersect(intersect(pattern_hd, pattern_mild), pattern_mod)),
  n234 = length(intersect(intersect(pattern_sev, pattern_mild), pattern_mod)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_sev), pattern_mild), pattern_mod)),
  category = c("HD","Severe", "Mild", "Moderate"),
  fill = c("white", "#ffcfdf", "#fefdca", "#FFC074"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()

write.table(setdiff(intersect(intersect(pattern_mild, pattern_mod), pattern_sev), pattern_hd), "venn_cdr3_su_cd8_maxprop_0.01threshold_setdiffHD-intersectDZ.txt", row.names = F, col.names = F, quote = F)








##### LEGACY
# Stats analysis

allprop_comp_mild <- cbind(allprop_HD, allprop_mild)
allprop_comp_mild <- allprop_comp_mild[which(rowSums(allprop_comp_mild) > 0.002),]
cdr3count <- dim(allprop_comp_mild)[1]
w <- NULL
for(i in 1:cdr3count){
  print(i)
  wout <- wilcox.test(as.numeric(allprop_comp_mild[i,1:16]),
                      as.numeric(allprop_comp_mild[i,17:124]))
  w <- c(w, wout$p.value)
}
w_mildvshd <- p.adjust(w, "BH")
#w_mildvshad[which(is.na(w_mildvshad))] <- 1

allprop_comp_mod <- cbind(allprop_HD, allprop_mod)
allprop_comp_mod <- allprop_comp_mod[which(rowSums(allprop_comp_mod) > 0.002),]
cdr3count <- dim(allprop_comp_mod)[1]
w <- NULL
for(i in 1:cdr3count){
  print(i)
  wout <- wilcox.test(as.numeric(allprop_comp_mod[i,1:16]),
                      as.numeric(allprop_comp_mod[i,17:109]))
  w <- c(w, wout$p.value)
}
w_modvshd <- p.adjust(w, "BH")


allprop_comp_severe <- cbind(allprop_HD, allprop_severe)
#allprop_comp_severe <- allprop_comp_severe[which(rowSums(allprop_comp_severe) > 0.002),]
#allprop_comp_severe <- log10(allprop_comp_severe + 1)
cdr3count <- dim(allprop_comp_severe)[1]
w <- NULL
for(i in 1:cdr3count){
  print(i)
  wout <- wilcox.test(as.numeric(allprop_comp_severe[i,1:16]),
                      as.numeric(allprop_comp_severe[i,17:65]))
  w <- c(w, wout$p.value)
}
w_severevshd <- w
w_severevshd[which(is.na(w_severevshd))] <- 1
w_severevshd_adj <- p.adjust(w_severevshd, "BH")

