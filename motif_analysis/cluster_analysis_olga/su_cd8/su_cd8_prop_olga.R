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

for(i in c("all", "hd", "mild", "moderate", "severe")){
  input <- fread(paste0("su_cd8_gliph2_clones_",i,"_cluster.csv"), fill = T)
  beta <- input$TcRb
  v <- input$V
  j <- input$J
  df <- data.frame(beta, v, j)
  df <- df[-which(df$beta == ""), ]
  write.table(df, paste0("su_cd8_clones_",i,"_cluster_olgainput.tsv"), quote = F, col.names = F, row.names = F, sep = "\t")

}

# olga-compute_pgen -i out.tsv --humanTRB -o out_pgens.tsv --v_in 1 --j_in 2
# olga-compute_pgen -i out.tsv --humanTRB -o out_pgens.tsv
#input <- input[-which(input$TcRb == ""), ]




i<-"hd"
input <- fread(paste0("su_cd8_gliph2_clones_",i,"_cluster.csv"), fill = T)
input <- input[-which(input$TcRb == ""), ]
pgenscores <- fread(paste0("su_cd8_clones_",i,"_cluster_pgens.tsv"))
pooled <- cbind(input, pgen = pgenscores$V2, logpgen = -log2(pgenscores$V2))
#pooled <- pooled[which(pooled$Fisher_score < 0.05),]



agg_pgen <-  aggregate(pooled$pgen, by = list(pooled$index), FUN = median)
#agg_logpgen <-  aggregate(pooled$logpgen, by = list(pooled$index), FUN = median)
agg_freq <-  aggregate(pooled$Freq, by = list(pooled$index), FUN = median)
agg_nsub <-  aggregate(pooled$number_subject, by = list(pooled$index), FUN = median)

aggdf <- data.frame(pattern = pooled$pattern[-which(duplicated(pooled$index))], 
                    freq = agg_freq$x, 
                    nsub = as.integer(agg_nsub$x),
                    pgen = agg_pgen$x,
                    logpgen = -log2(agg_pgen$x))

pooled_hd <- pooled
pooled_mild <- pooled
pooled_mod <- pooled
pooled_severe <- pooled

aggdf_hd <- aggdf
aggdf_mild <- aggdf
aggdf_mod <- aggdf
aggdf_severe <- aggdf

aggdf_hd <- aggdf_hd[-which(aggdf_hd$pattern == "single"),]
aggdf_mild <- aggdf_mild[-which(aggdf_mild$pattern == "single"),]
aggdf_mod <- aggdf_mod[-which(aggdf_mod$pattern == "single"),]
aggdf_severe <- aggdf_severe[-which(aggdf_severe$pattern == "single"),]


pattern_hd <- unique(aggdf_hd$pattern)
pattern_mild <- unique(aggdf_mild$pattern)
pattern_mod <- unique(aggdf_mod$pattern)
pattern_severe <- unique(aggdf_severe$pattern)

i <- "severe"
dfplt <- aggdf_severe
hits <- rep("clone", dim(dfplt)[1])
#hits[which(dfplt$pattern %in% setdiff(pattern_mod, pattern_hd))] <- "hit"
hits[which(dfplt$logpgen > 45)] <- "hit"
hits[which(dfplt$freq > 0.06)] <- "hit"
dfplt$hits <- hits
plt <- ggplot(dfplt, aes(x=freq , y=logpgen, size=nsub)) + geom_point(shape = 1,stroke = 1, color="black", alpha=0.3) #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),  axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_blank(), plot.title=element_text(size=30), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt  + labs(x = "Median cluster frequency", y = "-log2( Median cluster pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(3,15))
#plt <- plt + scale_x_continuous(breaks = pretty_breaks(n=7), limits = c(0,0.2))
plt <- plt + scale_y_continuous(breaks = pretty_breaks(n=5), limits = c(20,60))
plt <- plt + geom_point(data=subset(dfplt, hits == "hit"), aes(x=freq , y=logpgen, size=nsub), shape = 21,stroke = 1,  fill="#8E2657", alpha=0.5)
plt <- plt + geom_point(data=subset(dfplt, hits == "hit"), aes(x=freq , y=logpgen, size=nsub), shape = 1,stroke = 1,  color="#8E2657", alpha=1)
#plt <- plt + scale_colour_manual(values = c("#f6830f", "#2d6187", "#5aa469", "#9d0191"))
plt <- plt + geom_text_repel(data=subset(dfplt, hits == "hit"), mapping=aes(x=freq , y=logpgen, size=nsub, label =pattern), color="#8E2657", size = 12, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
plt <- plt + guides(color = guide_legend(override.aes = list(color = c("black") ) ) )
pdf(paste0("clusterplots/su_cd8_clones_",i,"_clusterplot.pdf"), height = 8, width = 10)
plt
dev.off()


pattern_dz <- intersect(intersect(pattern_mild, pattern_mod), pattern_severe)
pattern_dz_specific <- setdiff(pattern_dz, pattern_hd)

pooled_hd_temp <- pooled_hd
pooled_hd_temp$condition <- rep("HD", dim(pooled_hd_temp)[1])
pooled_mild_temp <- pooled_mild
pooled_mild_temp$condition <- rep("mild", dim(pooled_mild_temp)[1])
pooled_mod_temp <- pooled_mod
pooled_mod_temp$condition <- rep("moderate", dim(pooled_mod_temp)[1])
pooled_severe_temp <- pooled_severe
pooled_severe_temp$condition <- rep("severe", dim(pooled_severe_temp)[1])

pooled_master <- rbind(pooled_hd_temp, pooled_mild_temp, pooled_mod_temp, pooled_severe_temp)
fwrite(pooled_master, "clusterplots/su_cd8_clones_clonotypes_all.csv")
#pooled_master <- pooled_master[which(pooled_master$pattern %in% pattern_dz),]
pooled_master <- pooled_master[-which(pooled_master$pattern == "single"),]
pooled_master <- pooled_master[order(pooled_master$Freq, decreasing = T),]
#pooled_master <- pooled_master[order(pooled_master$logpgen, decreasing = T),]
pooled_master <- pooled_master[which(pooled_master$pattern %in% unique(pooled_master$pattern)[1:35]),]
pooled_master$pattern <- factor(pooled_master$pattern, levels = as.character(unique(pooled_master$pattern)))

plt <- ggplot(pooled_master, aes(x=Freq , y=logpgen, fill=condition)) + geom_point(shape = 21, size =4, alpha=1) #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  axis.text=element_text(size=15), axis.title=element_text(size=30), 
                                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.line = element_blank(), plot.title=element_text(size=30), 
                                strip.text.x = element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt + facet_wrap(~ pattern, ncol=7)
plt <- plt  + labs(x = "Clonotype frequency", y = "-log2( Clonotype pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(3,15))
plt <- plt + scale_x_continuous(breaks = pretty_breaks(n=5))#, limits = c(0,0.2))
pdf(paste0("clusterplots/su_cd8_clones_clonotypes_facet.pdf"), height = 8, width = 12)
plt
dev.off()

fwrite(pooled_master, "clusterplots/su_cd8_clones_clonotypes_dzspecific.csv")





pattern_all <- union(union(union(pattern_mild, pattern_mod), pattern_severe),pattern_hd)
pattern_dz <- intersect(intersect(pattern_mild, pattern_mod), pattern_severe)
pattern_dz_specific <- setdiff(pattern_dz, pattern_hd)


pdf("clusterplots/su_cd8_clones_clonotypes_dz_venn.pdf", height = 6, width = 6)
draw.triple.venn(area1 = length(pattern_mild), area2 = length(pattern_mod), area3 = length(pattern_severe), 
                 n12 = length(intersect(pattern_mild, pattern_mod)), 
                 n23 = length(intersect(pattern_mod, pattern_severe)), 
                 n13 = length(intersect(pattern_mild, pattern_severe)), 
                 n123 = length(intersect(intersect(pattern_mild, pattern_mod), pattern_severe)), 
                 cex = 2,
                 alpha     = 0.2,
                 #cat.cex = 2,
                 category = c("Mild", "Moderate", "Severe"),lwd = 3, col = "black",
                 fill = c("#FFC93C", "#FC5404", "#CD113B"));
dev.off()


pdf("clusterplots/su_cd8_clones_clonotypes_dz_specific_venn.pdf", height = 6, width = 6)
grid.newpage()
venn.plot <- draw.pairwise.venn(area1      = length(pattern_hd),
                                area2      = length(pattern_dz),
                                cross.area = length(intersect(pattern_hd, pattern_dz)),
                                category   = c("HD", "COVID-19"),
                                fill = c("#356895", "#ED6663"),
                                cex = 2,
                                alpha     = 0.2,
                                scaled     = FALSE)
dev.off()
write.table(pattern_dz_specific, "clusterplots/su_cd8_clones_clonotypes_dz_specific_venn_setdiff.csv", quote = F, row.names = F, col.names = F)




a<-unique(pattern_severe)
b<-unique(pattern_hd)
c<-unique(union(pattern_severe, pattern_hd))
q <- length(intersect(a,b))
m <- length(a)
n <- length(c) - length(a)
k <- length(b)
phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 


