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
# Ginisimp
#################### 


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



#################### 
# Ginisimp
#################### 

ginisimp_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_ginisimp_adaptiveHD.csv")
ginisimp_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_ginisimp_adaptive.csv")
ginisimp_adaptive <- ginisimp_adaptive[-which(is.na(ginisimp_adaptive[,6])),]
ginisimp_pooled <- rbind(ginisimp_adaptiveHD, ginisimp_adaptive)
  

# Age
ginisimp_pooled_age <- ginisimp_pooled[-which(ginisimp_pooled$age == ""),]



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
total_wen <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_wen.csv")
total_liao <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_liao.csv")
total_zhang <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_total_zhang.csv")

total_adaptive <- total_adaptive[-which(is.na(total_adaptive[,6])),]

total_pool_estimates <- c(total_adaptiveHD$Volume,
                          total_adaptive$Volume,
                          total_su_cd4$Volume,
                          total_su_cd8$Volume,
                          total_wen$Volume,
                          #total_liao$Volume,
                          total_zhang$Volume)

total_pool_dataset <- c(rep("AdaptiveHD",dim(total_adaptiveHD)[1]),
                        rep("Adaptive",dim(total_adaptive)[1]),
                        rep("Su_CD4",dim(total_su_cd4)[1]),
                        rep("Su_CD8",dim(total_su_cd8)[1]),
                        rep("Wen",dim(total_wen)[1]),
                        #rep("Liao",dim(total_liao)[1]),
                        rep("Zhang",dim(total_zhang)[1]))

total_pool_disease <- c(total_adaptiveHD$d_cond,
                        total_adaptive$d_cond,
                        total_su_cd4$`Disease status`,
                        total_su_cd8$`Disease status`,
                        total_wen$d_cond,
                        #total_liao$d_cond,
                        total_zhang$d_cond)
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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10total <- log10(pltd$total)

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
# gini
#################### 
gini_adaptiveHD <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_gini_adaptiveHD.csv")
gini_adaptive <- fread("../immunarch_outputs/immunarch_stats_adaptive/diversity_stats_gini_adaptive.csv")
gini_su_cd4 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd4.csv")
gini_su_cd8 <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_su_cd8.csv")
gini_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_wen.csv")
gini_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_liao.csv")
gini_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_gini_zhang.csv")

gini_adaptive <- gini_adaptive[-which(is.na(gini_adaptive[,6])),]

gini_pool_estimates <- c(gini_adaptiveHD$V1,
                          gini_adaptive$V1,
                          gini_su_cd4$V1,
                          gini_su_cd8$V1,
                          gini_wen$V1,
                          #gini_liao$V1,
                          gini_zhang$V1)

gini_pool_dataset <- c(rep("AdaptiveHD",dim(gini_adaptiveHD)[1]),
                        rep("Adaptive",dim(gini_adaptive)[1]),
                        rep("Su_CD4",dim(gini_su_cd4)[1]),
                        rep("Su_CD8",dim(gini_su_cd8)[1]),
                        rep("Wen",dim(gini_wen)[1]),
                        #rep("Liao",dim(gini_liao)[1]),
                        rep("Zhang",dim(gini_zhang)[1]))

gini_pool_disease <- c(gini_adaptiveHD$d_cond,
                        gini_adaptive$d_cond,
                        gini_su_cd4$`Disease status`,
                        gini_su_cd8$`Disease status`,
                        gini_wen$d_cond,
                        #gini_liao$d_cond,
                        gini_zhang$d_cond)
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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10gini <- log10(pltd$gini)

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
invsimp_wen <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_wen.csv")
invsimp_liao <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_liao.csv")
invsimp_zhang <- fread("../immunarch_outputs/immunarch_stats/diversity_stats_invsimp_zhang.csv")

invsimp_adaptive <- invsimp_adaptive[-which(is.na(invsimp_adaptive[,6])),]

invsimp_pool_estimates <- c(invsimp_adaptiveHD$Value,
                         invsimp_adaptive$Value,
                         invsimp_su_cd4$Value,
                         invsimp_su_cd8$Value,
                         invsimp_wen$Value,
                         #invsimp_liao$Value,
                         invsimp_zhang$Value)

invsimp_pool_dataset <- c(rep("AdaptiveHD",dim(invsimp_adaptiveHD)[1]),
                       rep("Adaptive",dim(invsimp_adaptive)[1]),
                       rep("Su_CD4",dim(invsimp_su_cd4)[1]),
                       rep("Su_CD8",dim(invsimp_su_cd8)[1]),
                       rep("Wen",dim(invsimp_wen)[1]),
                       #rep("Liao",dim(invsimp_liao)[1]),
                       rep("Zhang",dim(invsimp_zhang)[1]))

invsimp_pool_disease <- c(invsimp_adaptiveHD$d_cond,
                       invsimp_adaptive$d_cond,
                       invsimp_su_cd4$`Disease status`,
                       invsimp_su_cd8$`Disease status`,
                       invsimp_wen$d_cond,
                       #invsimp_liao$d_cond,
                       invsimp_zhang$d_cond)
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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
partition_levels <- rev(partition_levels)
pltd$partition <- factor(pltd$partition , levels = partition_levels)

pltd$log10invsimp <- log10(pltd$invsimp)

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
len_wen <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_wen.csv")
len_liao <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_liao.csv")
len_zhang <- fread("../immunarch_outputs/immunarch_stats/clonotype_stats_len_zhang.csv")

len_adaptive <- len_adaptive[-which(is.na(len_adaptive[,6])),]



len_pool_estimates <- c(len_adaptiveHD$Count,
                            len_adaptive$Count,
                            len_su_cd4$Count,
                            len_su_cd8$Count,
                            len_wen$Count,
                            #len_liao$Count,
                            len_zhang$Count)

len_pool_lengths <- c(len_adaptiveHD$Length,
                        len_adaptive$Length,
                        len_su_cd4$Length,
                        len_su_cd8$Length,
                        len_wen$Length,
                        #len_liao$Length,
                        len_zhang$Length)


len_pool_dataset <- c(rep("AdaptiveHD",dim(len_adaptiveHD)[1]),
                          rep("Adaptive",dim(len_adaptive)[1]),
                          rep("Su_CD4",dim(len_su_cd4)[1]),
                          rep("Su_CD8",dim(len_su_cd8)[1]),
                          rep("Wen",dim(len_wen)[1]),
                          #rep("Liao",dim(len_liao)[1]),
                          rep("Zhang",dim(len_zhang)[1]))

len_pool_disease <- c(len_adaptiveHD$d_cond,
                          len_adaptive$d_cond,
                          len_su_cd4$`Disease status`,
                          len_su_cd8$`Disease status`,
                          len_wen$d_cond,
                          #len_liao$d_cond,
                          len_zhang$d_cond)

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
                len_wen$Sample,
                #len_liao$Sample,
                len_zhang$Sample)

total_sample <- c(total_adaptiveHD$Sample,
                  total_adaptive$Sample,
                  total_su_cd4$Sample,
                  total_su_cd8$Sample,
                  total_wen$Sample,
                  #total_liao$Sample,
                  total_zhang$Sample)

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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
partition_levels <- rev(partition_levels)

#pltd <- pltd[intersect(which(pltd$lenlength >= 15), which(pltd$lenlength <= 81)),]
pltd <- pltd[which(pltd$lenlength %in% seq(15, 81, 3)),]
  
#pltd$partition <- factor(pltd$partition , levels = partition_levels)
pltd$dataset[which(pltd$dataset == "AdaptiveHD")] <- "Adaptive"
#pltd$dataset <- factor(pltd$dataset , levels = unique(pltd$dataset))
pltd$lenlength <- factor(pltd$lenlength, levels = sort(unique(pltd$lenlength), decreasing = F))

#pltd$log10len <- log10(pltd$len)
#pltd <- melt(pltd,id.vars=c('lenlength'), measure.vars=c('lencount'))
  
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
total_wen <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_wen.csv")
total_liao <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_liao.csv")
total_zhang <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_total_aa_zhang.csv")

total_adaptive <- total_adaptive[-which(is.na(total_adaptive[,6])),]

total_pool_estimates <- c(total_adaptiveHD$Volume,
                          total_adaptive$Volume,
                          total_su_cd4$Volume,
                          total_su_cd8$Volume,
                          total_wen$Volume,
                          #total_liao$Volume,
                          total_zhang$Volume)

total_pool_dataset <- c(rep("AdaptiveHD",dim(total_adaptiveHD)[1]),
                        rep("Adaptive",dim(total_adaptive)[1]),
                        rep("Su_CD4",dim(total_su_cd4)[1]),
                        rep("Su_CD8",dim(total_su_cd8)[1]),
                        rep("Wen",dim(total_wen)[1]),
                        #rep("Liao",dim(total_liao)[1]),
                        rep("Zhang",dim(total_zhang)[1]))

total_pool_disease <- c(total_adaptiveHD$d_cond,
                        total_adaptive$d_cond,
                        total_su_cd4$`Disease status`,
                        total_su_cd8$`Disease status`,
                        total_wen$d_cond,
                        #total_liao$d_cond,
                        total_zhang$d_cond)
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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
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
len_wen <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_wen.csv")
len_liao <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_liao.csv")
len_zhang <- fread("../immunarch_outputs/immunarch_stats_aa/clonotype_stats_len_aa_zhang.csv")

len_adaptive <- len_adaptive[-which(is.na(len_adaptive[,6])),]



len_pool_estimates <- c(len_adaptiveHD$Count,
                        len_adaptive$Count,
                        len_su_cd4$Count,
                        len_su_cd8$Count,
                        len_wen$Count,
                        #len_liao$Count,
                        len_zhang$Count)

len_pool_lengths <- c(len_adaptiveHD$Length,
                      len_adaptive$Length,
                      len_su_cd4$Length,
                      len_su_cd8$Length,
                      len_wen$Length,
                      #len_liao$Length,
                      len_zhang$Length)


len_pool_dataset <- c(rep("AdaptiveHD",dim(len_adaptiveHD)[1]),
                      rep("Adaptive",dim(len_adaptive)[1]),
                      rep("Su_CD4",dim(len_su_cd4)[1]),
                      rep("Su_CD8",dim(len_su_cd8)[1]),
                      rep("Wen",dim(len_wen)[1]),
                      #rep("Liao",dim(len_liao)[1]),
                      rep("Zhang",dim(len_zhang)[1]))

len_pool_disease <- c(len_adaptiveHD$d_cond,
                      len_adaptive$d_cond,
                      len_su_cd4$`Disease status`,
                      len_su_cd8$`Disease status`,
                      len_wen$d_cond,
                      #len_liao$d_cond,
                      len_zhang$d_cond)

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
                len_wen$Sample,
                #len_liao$Sample,
                len_zhang$Sample)

total_sample <- c(total_adaptiveHD$Sample,
                  total_adaptive$Sample,
                  total_su_cd4$Sample,
                  total_su_cd8$Sample,
                  total_wen$Sample,
                  #total_liao$Sample,
                  total_zhang$Sample)

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
                      "Wen_Healthy", "Wen_COVID-19",
                      "Zhang_Healthy", "Zhang_COVID-19")
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



