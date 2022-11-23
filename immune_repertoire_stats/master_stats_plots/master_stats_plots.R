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
library("ggpubr")


set.seed(20)


#################### 
# Chao1
#################### 


chao1_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_chao1_adaptiveHD.csv")
chao1_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_chao1_adaptive.csv")
chao1_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_chao1_su_cd4.csv")
chao1_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_chao1_su_cd8.csv")
chao1_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_chao1_zhang.csv")
chao1_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_chao1_liao.csv")
chao1_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_chao1_wen.csv")
#chao1_adaptive <- chao1_adaptive[-which(is.na(chao1_adaptive[,6])),]

chao1_pool_estimates <- c(chao1_adaptiveHD$Estimator,
                          chao1_adaptive$Estimator,
                          chao1_su_cd4$Estimator,
                          chao1_su_cd8$Estimator,
                          chao1_zhang$Estimator,
                          #chao1_liao$Estimator,
                          chao1_wen$Estimator)

chao1_pool_dataset <- c(rep("AdaptiveHD",dim(chao1_adaptiveHD)[1]),
                        rep("Adaptive",dim(chao1_adaptive)[1]),
                        rep("Su_CD4",dim(chao1_su_cd4)[1]),
                        rep("Su_CD8",dim(chao1_su_cd8)[1]),
                        rep("Zhang",dim(chao1_zhang)[1]),
                        #rep("Liao",dim(chao1_liao)[1]),
                        rep("Wen",dim(chao1_wen)[1]))

chao1_pool_disease <- c(chao1_adaptiveHD$d_cond,
                        chao1_adaptive$d_cond,
                        chao1_su_cd4$`Disease status`,
                        chao1_su_cd8$`Disease status`,
                        chao1_zhang$d_cond,
                        #chao1_liao$d_cond,
                        chao1_wen$d_cond)
chao1_pool_disease[which(chao1_pool_disease == "normal")] <- "Healthy"
chao1_pool_disease[which(chao1_pool_disease == "HD")] <- "Healthy"
chao1_pool_disease[which(chao1_pool_disease == "conv")] <- "COVID-19"
chao1_pool_disease[which(chao1_pool_disease == "moderate")] <- "COVID-19"
chao1_pool_disease[which(chao1_pool_disease == "severe")] <- "COVID-19"
chao1_pool_disease[which(chao1_pool_disease == "")] <- "COVID-19"
chao1_pool_disease[which(chao1_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#chao1_pool_dataset[which(is.na(chao1_pool_disease))] <- "COVID-19"

pltd <- data.frame(chao1 = chao1_pool_estimates, 
                   dataset = chao1_pool_dataset, 
                   disease = chao1_pool_disease,
                   partition = paste0(chao1_pool_dataset, "_", chao1_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10chao1 <- log10(pltd$chao1)


# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$log10chao1[which(pltd$partition == parts[i])],
                      pltd$log10chao1[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8",  "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_chao1.csv" ,quote = F,  row.names = F)




plt <- ggplot(pltd, aes(x=partition, y=log10chao1, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Chao1 Estimator )")
plt <- plt + coord_flip()
pdf("box_chao1.pdf", height = 10, width = 10)
plt
dev.off()

png("box_chao1.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10chao1, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Chao1 Estimator )")
plt <- plt + coord_flip()
pdf("box_chao1.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# Ginisimp
#################### 

ginisimp_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_ginisimp_adaptiveHD.csv")
ginisimp_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_ginisimp_adaptive.csv")
ginisimp_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_su_cd4.csv")
ginisimp_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_su_cd8.csv")
ginisimp_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_zhang.csv")
ginisimp_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_liao.csv")
ginisimp_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_ginisimp_wen.csv")
#ginisimp_adaptive <- ginisimp_adaptive[-which(is.na(ginisimp_adaptive[,6])),]

ginisimp_pool_estimates <- c(ginisimp_adaptiveHD$Value,
                          ginisimp_adaptive$Value,
                          ginisimp_su_cd4$Value,
                          ginisimp_su_cd8$Value,
                          ginisimp_zhang$Value,
                          #ginisimp_liao$Value,
                          ginisimp_wen$Value)

ginisimp_pool_dataset <- c(rep("AdaptiveHD",dim(ginisimp_adaptiveHD)[1]),
                        rep("Adaptive",dim(ginisimp_adaptive)[1]),
                        rep("Su_CD4",dim(ginisimp_su_cd4)[1]),
                        rep("Su_CD8",dim(ginisimp_su_cd8)[1]),
                        rep("Zhang",dim(ginisimp_zhang)[1]),
                        #rep("Liao",dim(ginisimp_liao)[1]),
                        rep("Wen",dim(ginisimp_wen)[1]))

ginisimp_pool_disease <- c(ginisimp_adaptiveHD$d_cond,
                        ginisimp_adaptive$d_cond,
                        ginisimp_su_cd4$`Disease status`,
                        ginisimp_su_cd8$`Disease status`,
                        ginisimp_zhang$d_cond,
                        #ginisimp_liao$d_cond,
                        ginisimp_wen$d_cond)
ginisimp_pool_disease[which(ginisimp_pool_disease == "normal")] <- "Healthy"
ginisimp_pool_disease[which(ginisimp_pool_disease == "HD")] <- "Healthy"
ginisimp_pool_disease[which(ginisimp_pool_disease == "conv")] <- "COVID-19"
ginisimp_pool_disease[which(ginisimp_pool_disease == "moderate")] <- "COVID-19"
ginisimp_pool_disease[which(ginisimp_pool_disease == "severe")] <- "COVID-19"
ginisimp_pool_disease[which(ginisimp_pool_disease == "")] <- "COVID-19"
ginisimp_pool_disease[which(ginisimp_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#ginisimp_pool_dataset[which(is.na(ginisimp_pool_disease))] <- "COVID-19"

pltd <- data.frame(ginisimp = ginisimp_pool_estimates, 
                   dataset = ginisimp_pool_dataset, 
                   disease = ginisimp_pool_disease,
                   partition = paste0(ginisimp_pool_dataset, "_", ginisimp_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

#pltd$log10ginisimp <- log10(pltd$ginisimp)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$ginisimp[which(pltd$partition == parts[i])],
                      pltd$ginisimp[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_ginisimp.csv" ,quote = F,  row.names = F)



plt <- ggplot(pltd, aes(x=partition, y=ginisimp, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Gini-Simpson index")
plt <- plt + coord_flip()
pdf("box_ginisimp.pdf", height = 10, width = 10)
plt
dev.off()

png("box_ginisimp.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=ginisimp, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Gini-Simpson index")
plt <- plt + coord_flip()
pdf("box_ginisimp.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# Total clonotypes
#################### 

total_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_total_adaptiveHD.csv")
total_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_total_adaptive.csv")
total_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_su_cd4.csv")
total_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_su_cd8.csv")
total_zhang <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_zhang.csv")
total_liao <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_liao.csv")
total_wen <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_wen.csv")
#total_adaptive <- total_adaptive[-which(is.na(total_adaptive[,6])),]

total_pool_estimates <- c(total_adaptiveHD$Volume,
                          total_adaptive$Volume,
                          total_su_cd4$Volume,
                          total_su_cd8$Volume,
                          total_zhang$Volume,
                          #total_liao$Volume,
                          total_wen$Volume)

total_pool_dataset <- c(rep("AdaptiveHD",dim(total_adaptiveHD)[1]),
                        rep("Adaptive",dim(total_adaptive)[1]),
                        rep("Su_CD4",dim(total_su_cd4)[1]),
                        rep("Su_CD8",dim(total_su_cd8)[1]),
                        rep("Zhang",dim(total_zhang)[1]),
                        #rep("Liao",dim(total_liao)[1]),
                        rep("Wen",dim(total_wen)[1]))

total_pool_disease <- c(total_adaptiveHD$d_cond,
                        total_adaptive$d_cond,
                        total_su_cd4$`Disease status`,
                        total_su_cd8$`Disease status`,
                        total_zhang$d_cond,
                        #total_liao$d_cond,
                        total_wen$d_cond)
total_pool_disease[which(total_pool_disease == "normal")] <- "Healthy"
total_pool_disease[which(total_pool_disease == "HD")] <- "Healthy"
total_pool_disease[which(total_pool_disease == "conv")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "moderate")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "severe")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#total_pool_dataset[which(is.na(total_pool_disease))] <- "COVID-19"

pltd <- data.frame(total = total_pool_estimates, 
                   dataset = total_pool_dataset, 
                   disease = total_pool_disease,
                   partition = paste0(total_pool_dataset, "_", total_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10total <- log10(pltd$total)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$log10total[which(pltd$partition == parts[i])],
                      pltd$log10total[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_total.csv" ,quote = F,  row.names = F)


plt <- ggplot(pltd, aes(x=partition, y=log10total, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Total unique clonotypes )")
plt <- plt + coord_flip()
pdf("box_total.pdf", height = 10, width = 10)
plt
dev.off()

png("box_total.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10total, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Total unique clonotypes )")
plt <- plt + coord_flip()
pdf("box_total.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# D50
#################### 


clones_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_clones_adaptiveHD.csv")
clones_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_clones_adaptive.csv")
clones_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_clones_su_cd4.csv")
clones_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_clones_su_cd8.csv")
clones_zhang <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_clones_zhang.csv")
clones_liao <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_clones_liao.csv")
clones_wen <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_clones_wen.csv")
#clones_adaptive <- clones_adaptive[-which(is.na(clones_adaptive[,6])),]

clones_pool_estimates <- c(clones_adaptiveHD$Clones,
                          clones_adaptive$Clones,
                          clones_su_cd4$Clones,
                          clones_su_cd8$Clones,
                          clones_zhang$Clones,
                          #clones_liao$Clones,
                          clones_wen$Clones)


d50_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_d50_adaptiveHD.csv")
d50_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_d50_adaptive.csv")
d50_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_d50_su_cd4.csv")
d50_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_d50_su_cd8.csv")
d50_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_d50_zhang.csv")
d50_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_d50_liao.csv")
d50_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_d50_wen.csv")
#d50_adaptive <- d50_adaptive[-which(is.na(d50_adaptive[,6])),]

d50_pool_estimates <- c(d50_adaptiveHD$Clones,
                          d50_adaptive$Clones,
                          d50_su_cd4$Clones,
                          d50_su_cd8$Clones,
                          d50_zhang$Clones,
                          #d50_liao$Clones,
                          d50_wen$Clones)
d50_pool_estimates <- d50_pool_estimates/clones_pool_estimates

d50_pool_dataset <- c(rep("AdaptiveHD",dim(d50_adaptiveHD)[1]),
                        rep("Adaptive",dim(d50_adaptive)[1]),
                        rep("Su_CD4",dim(d50_su_cd4)[1]),
                        rep("Su_CD8",dim(d50_su_cd8)[1]),
                        rep("Zhang",dim(d50_zhang)[1]),
                        #rep("Liao",dim(d50_liao)[1]),
                        rep("Wen",dim(d50_wen)[1]))

d50_pool_disease <- c(d50_adaptiveHD$d_cond,
                        d50_adaptive$d_cond,
                        d50_su_cd4$`Disease status`,
                        d50_su_cd8$`Disease status`,
                        d50_zhang$d_cond,
                        #d50_liao$d_cond,
                        d50_wen$d_cond)
d50_pool_disease[which(d50_pool_disease == "normal")] <- "Healthy"
d50_pool_disease[which(d50_pool_disease == "HD")] <- "Healthy"
d50_pool_disease[which(d50_pool_disease == "conv")] <- "COVID-19"
d50_pool_disease[which(d50_pool_disease == "moderate")] <- "COVID-19"
d50_pool_disease[which(d50_pool_disease == "severe")] <- "COVID-19"
d50_pool_disease[which(d50_pool_disease == "")] <- "COVID-19"
d50_pool_disease[which(d50_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#d50_pool_dataset[which(is.na(d50_pool_disease))] <- "COVID-19"

pltd <- data.frame(d50 = d50_pool_estimates, 
                   dataset = d50_pool_dataset, 
                   disease = d50_pool_disease,
                   partition = paste0(d50_pool_dataset, "_", d50_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10d50 <- log10(pltd$d50)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$d50[which(pltd$partition == parts[i])],
                      pltd$d50[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_d50.csv" ,quote = F,  row.names = F)


plt <- ggplot(pltd, aes(x=partition, y=d50, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "D50")
plt <- plt + coord_flip()
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


#################### 
# div
#################### 

div_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_div_adaptiveHD.csv")
div_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_div_adaptive.csv")
div_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_div_su_cd4.csv")
div_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_div_su_cd8.csv")
div_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_div_zhang.csv")
div_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_div_liao.csv")
div_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_div_wen.csv")
#div_adaptive <- div_adaptive[-which(is.na(div_adaptive[,6])),]

div_pool_estimates <- c(div_adaptiveHD$Value,
                          div_adaptive$Value,
                          div_su_cd4$Value,
                          div_su_cd8$Value,
                          div_zhang$Value,
                          #div_liao$Value,
                          div_wen$Value)

div_pool_dataset <- c(rep("AdaptiveHD",dim(div_adaptiveHD)[1]),
                        rep("Adaptive",dim(div_adaptive)[1]),
                        rep("Su_CD4",dim(div_su_cd4)[1]),
                        rep("Su_CD8",dim(div_su_cd8)[1]),
                        rep("Zhang",dim(div_zhang)[1]),
                        #rep("Liao",dim(div_liao)[1]),
                        rep("Wen",dim(div_wen)[1]))

div_pool_disease <- c(div_adaptiveHD$d_cond,
                        div_adaptive$d_cond,
                        div_su_cd4$`Disease status`,
                        div_su_cd8$`Disease status`,
                        div_zhang$d_cond,
                        #div_liao$d_cond,
                        div_wen$d_cond)
div_pool_disease[which(div_pool_disease == "normal")] <- "Healthy"
div_pool_disease[which(div_pool_disease == "HD")] <- "Healthy"
div_pool_disease[which(div_pool_disease == "conv")] <- "COVID-19"
div_pool_disease[which(div_pool_disease == "moderate")] <- "COVID-19"
div_pool_disease[which(div_pool_disease == "severe")] <- "COVID-19"
div_pool_disease[which(div_pool_disease == "")] <- "COVID-19"
div_pool_disease[which(div_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#div_pool_dataset[which(is.na(div_pool_disease))] <- "COVID-19"

pltd <- data.frame(div = div_pool_estimates, 
                   dataset = div_pool_dataset, 
                   disease = div_pool_disease,
                   partition = paste0(div_pool_dataset, "_", div_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10div <- log10(pltd$div)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$log10div[which(pltd$partition == parts[i])],
                      pltd$log10div[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_div.csv" ,quote = F,  row.names = F)


plt <- ggplot(pltd, aes(x=partition, y=log10div, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( div )")
plt <- plt + coord_flip()
pdf("box_div.pdf", height = 10, width = 10)
plt
dev.off()

png("box_div.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10div, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( div )")
plt <- plt + coord_flip()
pdf("box_div.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# gini
#################### 
gini_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_gini_adaptiveHD.csv")
gini_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_gini_adaptive.csv")
gini_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd4.csv")
gini_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd8.csv")
gini_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_zhang.csv")
gini_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_liao.csv")
gini_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_wen.csv")
#gini_adaptive <- gini_adaptive[-which(is.na(gini_adaptive[,6])),]

gini_pool_estimates <- c(gini_adaptiveHD$V1,
                          gini_adaptive$V1,
                          gini_su_cd4$V1,
                          gini_su_cd8$V1,
                          gini_zhang$V1,
                          #gini_liao$V1,
                          gini_wen$V1)

gini_pool_dataset <- c(rep("AdaptiveHD",dim(gini_adaptiveHD)[1]),
                        rep("Adaptive",dim(gini_adaptive)[1]),
                        rep("Su_CD4",dim(gini_su_cd4)[1]),
                        rep("Su_CD8",dim(gini_su_cd8)[1]),
                        rep("Zhang",dim(gini_zhang)[1]),
                        #rep("Liao",dim(gini_liao)[1]),
                        rep("Wen",dim(gini_wen)[1]))

gini_pool_disease <- c(gini_adaptiveHD$d_cond,
                        gini_adaptive$d_cond,
                        gini_su_cd4$`Disease status`,
                        gini_su_cd8$`Disease status`,
                        gini_zhang$d_cond,
                        #gini_liao$d_cond,
                        gini_wen$d_cond)
gini_pool_disease[which(gini_pool_disease == "normal")] <- "Healthy"
gini_pool_disease[which(gini_pool_disease == "HD")] <- "Healthy"
gini_pool_disease[which(gini_pool_disease == "conv")] <- "COVID-19"
gini_pool_disease[which(gini_pool_disease == "moderate")] <- "COVID-19"
gini_pool_disease[which(gini_pool_disease == "severe")] <- "COVID-19"
gini_pool_disease[which(gini_pool_disease == "")] <- "COVID-19"
gini_pool_disease[which(gini_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#gini_pool_dataset[which(is.na(gini_pool_disease))] <- "COVID-19"

pltd <- data.frame(gini = gini_pool_estimates, 
                   dataset = gini_pool_dataset, 
                   disease = gini_pool_disease,
                   partition = paste0(gini_pool_dataset, "_", gini_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10gini <- log10(pltd$gini)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$gini[which(pltd$partition == parts[i])],
                      pltd$gini[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_gini.csv" ,quote = F,  row.names = F)


plt <- ggplot(pltd, aes(x=partition, y=gini, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Gini coefficient")
plt <- plt + coord_flip()
pdf("box_gini.pdf", height = 10, width = 10)
plt
dev.off()

png("box_gini.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10gini, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "Gini coefficient")
plt <- plt + coord_flip()
pdf("box_gini.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# invsimp
#################### 
invsimp_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_invsimp_adaptiveHD.csv")
invsimp_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_invsimp_adaptive.csv")
invsimp_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_su_cd4.csv")
invsimp_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_su_cd8.csv")
invsimp_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_zhang.csv")
invsimp_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_liao.csv")
invsimp_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_wen.csv")
#invsimp_adaptive <- invsimp_adaptive[-which(is.na(invsimp_adaptive[,6])),]

invsimp_pool_estimates <- c(invsimp_adaptiveHD$Value,
                         invsimp_adaptive$Value,
                         invsimp_su_cd4$Value,
                         invsimp_su_cd8$Value,
                         invsimp_zhang$Value,
                         #invsimp_liao$Value,
                         invsimp_wen$Value)

invsimp_pool_dataset <- c(rep("AdaptiveHD",dim(invsimp_adaptiveHD)[1]),
                       rep("Adaptive",dim(invsimp_adaptive)[1]),
                       rep("Su_CD4",dim(invsimp_su_cd4)[1]),
                       rep("Su_CD8",dim(invsimp_su_cd8)[1]),
                       rep("Zhang",dim(invsimp_zhang)[1]),
                       #rep("Liao",dim(invsimp_liao)[1]),
                       rep("Wen",dim(invsimp_wen)[1]))

invsimp_pool_disease <- c(invsimp_adaptiveHD$d_cond,
                       invsimp_adaptive$d_cond,
                       invsimp_su_cd4$`Disease status`,
                       invsimp_su_cd8$`Disease status`,
                       invsimp_zhang$d_cond,
                       #invsimp_liao$d_cond,
                       invsimp_wen$d_cond)
invsimp_pool_disease[which(invsimp_pool_disease == "normal")] <- "Healthy"
invsimp_pool_disease[which(invsimp_pool_disease == "HD")] <- "Healthy"
invsimp_pool_disease[which(invsimp_pool_disease == "conv")] <- "COVID-19"
invsimp_pool_disease[which(invsimp_pool_disease == "moderate")] <- "COVID-19"
invsimp_pool_disease[which(invsimp_pool_disease == "severe")] <- "COVID-19"
invsimp_pool_disease[which(invsimp_pool_disease == "")] <- "COVID-19"
invsimp_pool_disease[which(invsimp_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#invsimp_pool_dataset[which(is.na(invsimp_pool_disease))] <- "COVID-19"

pltd <- data.frame(invsimp = invsimp_pool_estimates, 
                   dataset = invsimp_pool_dataset, 
                   disease = invsimp_pool_disease,
                   partition = paste0(invsimp_pool_dataset, "_", invsimp_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10invsimp <- log10(pltd$invsimp)

# stats
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
for(i in seq(1,length(parts), 2)){
  print(i)
  wout <- wilcox.test(pltd$log10invsimp[which(pltd$partition == parts[i])],
                      pltd$log10invsimp[which(pltd$partition == parts[i+1])])
  pvalues <- c(pvalues, wout$p.value)
}
statvals <- data.frame(partition = c("Adaptive", "Su_CD4", "Su_CD8", "Wen", "Zhang"), pvalues)
statvals <- statvals[c(1:3,5,4),]
write.csv(statvals,"wilcoxon_pvalues_invsimp.csv" ,quote = F,  row.names = F)


plt <- ggplot(pltd, aes(x=partition, y=log10invsimp, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Inverse Simpson index )")
plt <- plt + coord_flip()
pdf("box_invsimp.pdf", height = 10, width = 10)
plt
dev.off()

png("box_invsimp.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10invsimp, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Inverse Simpson index )")
plt <- plt + coord_flip()
pdf("box_invsimp.fo.pdf", height = 10, width = 10)
plt
dev.off()






#################### 
# LEN ANALYSIS
#################### 
len_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_len_adaptiveHD.csv")
len_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/clonotype_stats_len_adaptive.csv")
len_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_su_cd4.csv")
len_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_su_cd8.csv")
len_zhang <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_zhang.csv")
len_liao <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_liao.csv")
len_wen <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_wen.csv")
#len_adaptive <- len_adaptive[-which(is.na(len_adaptive[,6])),]



len_pool_estimates <- c(len_adaptiveHD$Count,
                            len_adaptive$Count,
                            len_su_cd4$Count,
                            len_su_cd8$Count,
                            len_zhang$Count,
                            #len_liao$Count,
                            len_wen$Count)

len_pool_lengths <- c(len_adaptiveHD$Length,
                        len_adaptive$Length,
                        len_su_cd4$Length,
                        len_su_cd8$Length,
                        len_zhang$Length,
                        #len_liao$Length,
                        len_wen$Length)


len_pool_dataset <- c(rep("AdaptiveHD",dim(len_adaptiveHD)[1]),
                          rep("Adaptive",dim(len_adaptive)[1]),
                          rep("Su_CD4",dim(len_su_cd4)[1]),
                          rep("Su_CD8",dim(len_su_cd8)[1]),
                          rep("Zhang",dim(len_zhang)[1]),
                          #rep("Liao",dim(len_liao)[1]),
                          rep("Wen",dim(len_wen)[1]))

len_pool_disease <- c(len_adaptiveHD$d_cond,
                          len_adaptive$d_cond,
                          len_su_cd4$`Disease status`,
                          len_su_cd8$`Disease status`,
                          len_zhang$d_cond,
                          #len_liao$d_cond,
                          len_wen$d_cond)

len_pool_disease[which(len_pool_disease == "normal")] <- "Healthy"
len_pool_disease[which(len_pool_disease == "HD")] <- "Healthy"
len_pool_disease[which(len_pool_disease == "conv")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "moderate")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "severe")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#len_pool_dataset[which(is.na(len_pool_disease))] <- "COVID-19"

len_sample <- c(len_adaptiveHD$Sample,
                len_adaptive$Sample,
                len_su_cd4$Sample,
                len_su_cd8$Sample,
                len_zhang$Sample,
                #len_liao$Sample,
                len_wen$Sample)

total_sample <- c(total_adaptiveHD$Sample,
                  total_adaptive$Sample,
                  total_su_cd4$Sample,
                  total_su_cd8$Sample,
                  total_zhang$Sample,
                  #total_liao$Sample,
                  total_wen$Sample)

len_pool_estimates_nl_factor <- total_pool_estimates[match(len_sample, total_sample)]
len_pool_estimates_nl <- len_pool_estimates / len_pool_estimates_nl_factor
### BUILD DF
pltd <- data.frame(lenlength = len_pool_lengths, 
                   lencount = len_pool_estimates_nl,
                   dataset = len_pool_dataset, 
                   disease = len_pool_disease,
                   partition = paste0(len_pool_dataset, "_", len_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)

#pltd <- pltd[intersect(which(pltd$lenlength >= 15), which(pltd$lenlength <= 81)),]
pltd <- pltd[which(pltd$lenlength %in% seq(15, 81, 3)),]
  
#pltd$partition <- factor(pltd$partition , levels = partition_levels)
pltd$dataset[which(pltd$dataset == "AdaptiveHD")] <- "Adaptive"
#pltd$dataset <- factor(pltd$dataset , levels = unique(pltd$dataset))
pltd$lenlength <- factor(pltd$lenlength, levels = sort(unique(pltd$lenlength), decreasing = F))

#pltd$log10len <- log10(pltd$len)
#pltd <- melt(pltd,id.vars=c('lenlength'), measure.vars=c('lencount'))
  
pltd$dataset <- factor(pltd$dataset, levels = c("Adaptive", "Su_CD4", "Su_CD8", "Zhang", "Wen"))

plt <- ggplot(pltd, aes(x=lenlength, y=lencount, fill = disease)) 
plt <- plt + geom_point(position=position_jitterdodge(),size=2.5, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(alpha=1, lwd=0.5, color="#FFC947",  outlier.shape = NA)#+ scale_fill_manual(values = c("white"))
plt <- plt+ facet_grid(dataset ~ .)
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=16), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "CDR3 length (NT)", y = "Clonotype frequency")
pdf("box_len.pdf", height = 14, width = 14)
plt
dev.off()

png("box_len.png", height = 14, width = 14,units = "in", res = 300)
plt
dev.off()


plt <- ggplot(pltd, aes(x=lenlength, y=lencount, fill = disease)) 
plt <- plt + geom_point(position=position_jitterdodge(),size=2.5, pch=26, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt+ facet_grid(dataset ~ .)
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=16), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "CDR3 length (NT)", y = "Clonotype frequency")
pdf("box_len.fo.pdf", height = 14, width = 14)
plt
dev.off()



# stats are IMPERFECT
pvalues <- NULL
parts <- as.character(unique(pltd$partition))
partslen <- as.numeric(unique(pltd$lenlength))

for(i in seq(1,length(parts), 2)){
  print(i)
  for(j in partslen){
    print(j)
    if(length(intersect(which(pltd$lenlength == j),which(pltd$partition == parts[i]))) > 0){
      if(length(intersect(which(pltd$lenlength == j),which(pltd$partition == parts[i+1]))) > 0){
        wout <- wilcox.test(pltd$lencount[intersect(which(pltd$lenlength == j),which(pltd$partition == parts[i]))],
                            pltd$lencount[intersect(which(pltd$lenlength == j),which(pltd$partition == parts[i+1]))],
                            exact=F)
        pvalues <- c(pvalues, wout$p.value) 
      }
    }
  }
}
# stats are IMPERFECT


#################### 
# AMINO ACID
####################
# AMINO ACID
####################
# AMINO ACID
#################### 

#################### 
# Total clonotypes AA
#################### 

total_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive_aa/clonotype_stats_total_aa_adaptiveHD.csv")
total_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive_aa/clonotype_stats_total_aa_adaptive.csv")
total_su_cd4 <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_su_cd4.csv")
total_su_cd8 <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_su_cd8.csv")
total_zhang <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_zhang.csv")
total_liao <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_liao.csv")
total_wen <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_wen.csv")
#total_adaptive <- total_adaptive[-which(is.na(total_adaptive[,6])),]

total_pool_estimates <- c(total_adaptiveHD$Volume,
                          total_adaptive$Volume,
                          total_su_cd4$Volume,
                          total_su_cd8$Volume,
                          total_zhang$Volume,
                          #total_liao$Volume,
                          total_wen$Volume)

total_pool_dataset <- c(rep("AdaptiveHD",dim(total_adaptiveHD)[1]),
                        rep("Adaptive",dim(total_adaptive)[1]),
                        rep("Su_CD4",dim(total_su_cd4)[1]),
                        rep("Su_CD8",dim(total_su_cd8)[1]),
                        rep("Zhang",dim(total_zhang)[1]),
                        #rep("Liao",dim(total_liao)[1]),
                        rep("Wen",dim(total_wen)[1]))

total_pool_disease <- c(total_adaptiveHD$d_cond,
                        total_adaptive$d_cond,
                        total_su_cd4$`Disease status`,
                        total_su_cd8$`Disease status`,
                        total_zhang$d_cond,
                        #total_liao$d_cond,
                        total_wen$d_cond)
total_pool_disease[which(total_pool_disease == "normal")] <- "Healthy"
total_pool_disease[which(total_pool_disease == "HD")] <- "Healthy"
total_pool_disease[which(total_pool_disease == "conv")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "moderate")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "severe")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "")] <- "COVID-19"
total_pool_disease[which(total_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#total_pool_dataset[which(is.na(total_pool_disease))] <- "COVID-19"

pltd <- data.frame(total = total_pool_estimates, 
                   dataset = total_pool_dataset, 
                   disease = total_pool_disease,
                   partition = paste0(total_pool_dataset, "_", total_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10total <- log10(pltd$total)

plt <- ggplot(pltd, aes(x=partition, y=log10total, fill=disease)) + geom_jitter(size=3, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(fill="white", alpha=0.5, lwd=1, color="black",  outlier.shape = NA)
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Total unique clonotypes )")
plt <- plt + coord_flip()
pdf("box_total_aa.pdf", height = 10, width = 10)
plt
dev.off()

png("box_total_aa.png", height = 10, width = 10,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=partition, y=log10total, color=disease)) + geom_jitter(size=3, pch=26) + scale_color_manual(values = c("#CD113B", "#053742"))
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "", y = "log10( Total unique clonotypes )")
plt <- plt + coord_flip()
pdf("box_total_aa.fo.pdf", height = 10, width = 10)
plt
dev.off()



#################### 
# LEN ANALYSIS AMINO ACID
#################### 
len_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive_aa/clonotype_stats_len_aa_adaptiveHD.csv")
len_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive_aa/clonotype_stats_len_aa_adaptive.csv")
len_su_cd4 <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_su_cd4.csv")
len_su_cd8 <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_su_cd8.csv")
len_zhang <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_zhang.csv")
len_liao <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_liao.csv")
len_wen <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_wen.csv")
#len_adaptive <- len_adaptive[-which(is.na(len_adaptive[,6])),]



len_pool_estimates <- c(len_adaptiveHD$Count,
                        len_adaptive$Count,
                        len_su_cd4$Count,
                        len_su_cd8$Count,
                        len_zhang$Count,
                        #len_liao$Count,
                        len_wen$Count)

len_pool_lengths <- c(len_adaptiveHD$Length,
                      len_adaptive$Length,
                      len_su_cd4$Length,
                      len_su_cd8$Length,
                      len_zhang$Length,
                      #len_liao$Length,
                      len_wen$Length)


len_pool_dataset <- c(rep("AdaptiveHD",dim(len_adaptiveHD)[1]),
                      rep("Adaptive",dim(len_adaptive)[1]),
                      rep("Su_CD4",dim(len_su_cd4)[1]),
                      rep("Su_CD8",dim(len_su_cd8)[1]),
                      rep("Zhang",dim(len_zhang)[1]),
                      #rep("Liao",dim(len_liao)[1]),
                      rep("Wen",dim(len_wen)[1]))

len_pool_disease <- c(len_adaptiveHD$d_cond,
                      len_adaptive$d_cond,
                      len_su_cd4$`Disease status`,
                      len_su_cd8$`Disease status`,
                      len_zhang$d_cond,
                      #len_liao$d_cond,
                      len_wen$d_cond)

len_pool_disease[which(len_pool_disease == "normal")] <- "Healthy"
len_pool_disease[which(len_pool_disease == "HD")] <- "Healthy"
len_pool_disease[which(len_pool_disease == "conv")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "moderate")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "severe")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "")] <- "COVID-19"
len_pool_disease[which(len_pool_disease == "COVID-19 Positive")] <- "COVID-19"
#len_pool_dataset[which(is.na(len_pool_disease))] <- "COVID-19"

len_sample <- c(len_adaptiveHD$Sample,
                len_adaptive$Sample,
                len_su_cd4$Sample,
                len_su_cd8$Sample,
                len_zhang$Sample,
                #len_liao$Sample,
                len_wen$Sample)

total_sample <- c(total_adaptiveHD$Sample,
                  total_adaptive$Sample,
                  total_su_cd4$Sample,
                  total_su_cd8$Sample,
                  total_zhang$Sample,
                  #total_liao$Sample,
                  total_wen$Sample)

len_pool_estimates_nl_factor <- total_pool_estimates[match(len_sample, total_sample)]
len_pool_estimates_nl <- len_pool_estimates / len_pool_estimates_nl_factor
### BUILD DF
pltd <- data.frame(lenlength = len_pool_lengths, 
                   lencount = len_pool_estimates_nl,
                   lencountclones = len_pool_estimates,
                   dataset = len_pool_dataset, 
                   disease = len_pool_disease,
                   partition = paste0(len_pool_dataset, "_", len_pool_disease))
#pltd <- pltd[-which(pltd$country == ""),]
pltd <- pltd[order(pltd$partition),]
partition_levels <- c("AdaptiveHD_Healthy", "Adaptive_COVID-19",
                      "Su_CD4_Healthy", "Su_CD4_COVID-19",
                      "Su_CD8_Healthy", "Su_CD8_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19",
                      "Wen_Healthy", "Wen_COVID-19")
partition_levels <- rev(partition_levels)

#pltd <- pltd[intersect(which(pltd$lenlength >= 15), which(pltd$lenlength <= 81)),]
potentialengths <- seq(15, 81, 3)/3
pltd <- pltd[which(pltd$lenlength %in% potentialengths),]

#pltd$partition <- factor(pltd$partition , levels = partition_levels)
pltd$dataset[which(pltd$dataset == "AdaptiveHD")] <- "Adaptive"
#pltd$dataset <- factor(pltd$dataset , levels = unique(pltd$dataset))
pltd$lenlength <- factor(pltd$lenlength, levels = sort(unique(pltd$lenlength), decreasing = F))

#pltd$log10len <- log10(pltd$len)
#pltd <- melt(pltd,id.vars=c('lenlength'), measure.vars=c('lencount'))
pltd$dataset <- factor(pltd$dataset, levels = c("Adaptive", "Su_CD4", "Su_CD8", "Zhang", "Wen"))


plt <- ggplot(pltd, aes(x=lenlength, y=lencount, fill = disease)) 
plt <- plt + geom_point(position=position_jitterdodge(),size=2.5, pch=21, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt + geom_boxplot(alpha=1, lwd=0.5, color="#FFC947",  outlier.shape = NA)#+ scale_fill_manual(values = c("white"))
plt <- plt+ facet_grid(dataset ~ .)
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=16), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "CDR3 length (AA)", y = "Clonotype frequency")
pdf("box_len_aa.pdf", height = 14, width = 14)
plt
dev.off()

png("box_len_aa.png", height = 14, width = 14,units = "in", res = 300)
plt
dev.off()


plt <- ggplot(pltd, aes(x=lenlength, y=lencount, fill = disease)) 
plt <- plt + geom_point(position=position_jitterdodge(),size=2.5, pch=26, colour="black") + scale_fill_manual(values = c("#CD113B", "#053742"))
plt <- plt+ facet_grid(dataset ~ .)
plt <- plt + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=16), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6))
plt <- plt + labs(x= "CDR3 length (AA)", y = "Clonotype frequency")
pdf("box_len_aa.fo.pdf", height = 14, width = 14)
plt
dev.off()



