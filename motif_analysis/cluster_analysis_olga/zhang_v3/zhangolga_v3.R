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
library("VennDiagram")
library("ggrepel")


set.seed(20)


# PROCESSING
input <- fread("zhangtcr_v3_conv_cluster.csv", fill = T)
beta <- input$TcRb
v <- input$V
j <- input$J
df <- data.frame(beta, v, j)
df <- df[-which(df$beta == ""), ]
write.table(df, "zhangtcr_v3_conv_cluster_olgainput.tsv", quote = F, col.names = F, row.names = F, sep = "\t")

# olga-compute_pgen -i out.tsv --humanTRB -o out_pgens.tsv --v_in 1 --j_in 2
# olga-compute_pgen -i out.tsv --humanTRB -o out_pgens.tsv
#input <- input[-which(input$TcRb == ""), ]



##sort(c(df_conv$logpgen, df_hd$logpgen, df_sev$logpgen, df_mod$logpgen), decreasing = T)
### LOAD DATA

input <- fread("zhangtcr_v3_hd_cluster.csv", fill = T)
input <- input[-which(input$TcRb == ""), ]
pgenscores <- fread("zhangtcr_v3_hd_cluster_olgainput_pgens.tsv")
pooled <- cbind(input, pgen = pgenscores$V2, logpgen = -log2(pgenscores$V2))
clust_hd <- pooled
pooled <- pooled[which(pooled$Fisher_score < 0.05),]
agg_logpgen <-  aggregate(pooled$logpgen, by = list(pooled$index), FUN = median)
agg_freq <-  aggregate(pooled$Freq, by = list(pooled$index), FUN = median)
agg_nsub <-  aggregate(pooled$number_subject, by = list(pooled$index), FUN = median)

#agg_pattern <-  aggregate(pooled$pattern, by = list(pooled$index), FUN = [[])
#plot(agg_freq$x, agg_logpgen$x)
df_hd <- data.frame(pattern = pooled$pattern[-which(duplicated(pooled$index))], 
                    freq = agg_freq$x, 
                    logpgen = agg_logpgen$x,
                    nsub = as.integer(agg_nsub$x))



input <- fread("zhangtcr_v3_conv_cluster.csv", fill = T)
input <- input[-which(input$TcRb == ""), ]
pgenscores <- fread("zhangtcr_v3_conv_cluster_olgainput_pgens.tsv")
pooled <- cbind(input, pgen = pgenscores$V2, logpgen = -log2(pgenscores$V2))
clust_conv <- pooled
pooled <- pooled[which(pooled$Fisher_score < 0.05),]
agg_logpgen <-  aggregate(pooled$logpgen, by = list(pooled$index), FUN = median)
agg_freq <-  aggregate(pooled$Freq, by = list(pooled$index), FUN = median)
agg_nsub <-  aggregate(pooled$number_subject, by = list(pooled$index), FUN = median)

#agg_pattern <-  aggregate(pooled$pattern, by = list(pooled$index), FUN = [[])
#plot(agg_freq$x, agg_logpgen$x)
df_conv <- data.frame(pattern = pooled$pattern[-which(duplicated(pooled$index))], 
                      freq = agg_freq$x, 
                      logpgen = agg_logpgen$x,
                      nsub = as.integer(agg_nsub$x))



input <- fread("zhangtcr_v3_mod_cluster.csv", fill = T)
input <- input[-which(input$TcRb == ""), ]
pgenscores <- fread("zhangtcr_v3_mod_cluster_olgainput_pgens.tsv")
pooled <- cbind(input, pgen = pgenscores$V2, logpgen = -log2(pgenscores$V2))
clust_mod <- pooled
pooled <- pooled[which(pooled$Fisher_score < 0.05),]
agg_logpgen <-  aggregate(pooled$logpgen, by = list(pooled$index), FUN = median)
agg_freq <-  aggregate(pooled$Freq, by = list(pooled$index), FUN = median)
agg_nsub <-  aggregate(pooled$number_subject, by = list(pooled$index), FUN = median)

#agg_pattern <-  aggregate(pooled$pattern, by = list(pooled$index), FUN = [[])
#plot(agg_freq$x, agg_logpgen$x)
df_mod <- data.frame(pattern = pooled$pattern[-which(duplicated(pooled$index))], 
                     freq = agg_freq$x, 
                     logpgen = agg_logpgen$x,
                     nsub = as.integer(agg_nsub$x))



input <- fread("zhangtcr_v3_sev_cluster.csv", fill = T)
input <- input[-which(input$TcRb == ""), ]
pgenscores <- fread("zhangtcr_v3_sev_cluster_olgainput_pgens.tsv")
pooled <- cbind(input, pgen = pgenscores$V2, logpgen = -log2(pgenscores$V2))
clust_sev <- pooled
pooled <- pooled[which(pooled$Fisher_score < 0.05),]
agg_logpgen <-  aggregate(pooled$logpgen, by = list(pooled$index), FUN = median)
agg_freq <-  aggregate(pooled$Freq, by = list(pooled$index), FUN = median)
agg_nsub <-  aggregate(pooled$number_subject, by = list(pooled$index), FUN = median)

#agg_pattern <-  aggregate(pooled$pattern, by = list(pooled$index), FUN = [[])
#plot(agg_freq$x, agg_logpgen$x)
df_sev <- data.frame(pattern = pooled$pattern[-which(duplicated(pooled$index))], 
                     freq = agg_freq$x, 
                     logpgen = agg_logpgen$x,
                     nsub = as.integer(agg_nsub$x))


#### Venn diagram


pattern_hd <- unique(as.character(df_hd$pattern))
pattern_mod <- unique(as.character(df_mod$pattern))
pattern_sev <- unique(as.character(df_sev$pattern))
pattern_conv <- unique(as.character(df_conv$pattern))


pdf("venn_cluster_pattern_full.pdf", height = 10, width = 10)
venn.plot <- draw.quad.venn(
  area1 = length(pattern_hd),
  area2 = length(pattern_mod),
  area3 = length(pattern_sev),
  area4 = length(pattern_conv),
  n12 = length(intersect(pattern_hd, pattern_mod)),
  n13 = length(intersect(pattern_hd, pattern_sev)),
  n14 = length(intersect(pattern_hd, pattern_conv)),
  n23 = length(intersect(pattern_mod, pattern_sev)),
  n24 = length(intersect(pattern_mod, pattern_conv)),
  n34 = length(intersect(pattern_sev, pattern_conv)),
  n123 = length(intersect(intersect(pattern_hd, pattern_mod), pattern_sev)),
  n124 = length(intersect(intersect(pattern_hd, pattern_mod), pattern_conv)),
  n134 = length(intersect(intersect(pattern_hd, pattern_sev), pattern_conv)),
  n234 = length(intersect(intersect(pattern_mod, pattern_sev), pattern_conv)),
  n1234 = length(intersect(intersect(intersect(pattern_hd, pattern_mod), pattern_sev), pattern_conv)),
  category = c("HD", "Moderate", "Severe", "Convalescent"),
  fill = c("white", "#fefdca", "#ffcfdf", "#a5dee5"),
  lty = "solid",
  cex = 2,
  cat.cex = 2,
  #cat.col = c("orange", "red", "green", "blue"),
  margin = 0.1
);
dev.off()



require(gridExtra)
g1 <- draw.triple.venn(area1 = length(pattern_mod), area2 = length(pattern_sev), area3 = length(pattern_conv), 
                       n12 = length(intersect(pattern_mod, pattern_sev)), n23 = length(intersect(pattern_sev, pattern_conv)), n13 = length(intersect(pattern_mod, pattern_conv)), 
                       n123 = length(intersect(intersect(pattern_mod, pattern_sev), pattern_conv)), 
                       category = c("Moderate", "Severe", "Convalescent"),lwd = 3, col = "black",
                       lty = "solid",
                       cex = 2,
                       cat.cex = 2,
                       #cat.col = c("orange", "red", "green", "blue"),
                       margin = 0.1,
                       fill = c("#fefdca", "#ffcfdf", "#a5dee5"))
pdf("venn_cluster_pattern_covid.pdf", height = 10, width = 10)
grid.arrange(gTree(children=g1))
dev.off()


unique_pattern_mod <-  setdiff(setdiff(setdiff(pattern_mod,pattern_sev), pattern_conv), pattern_hd)
unique_pattern_sev<-  setdiff(setdiff(setdiff(pattern_sev,pattern_mod), pattern_conv), pattern_hd)
unique_pattern_conv <-  setdiff(setdiff(setdiff(pattern_conv, pattern_sev), pattern_mod), pattern_hd)
unique_pattern_hd <-  setdiff(setdiff(setdiff(pattern_hd,pattern_sev), pattern_mod), pattern_conv)



shared_pattern <-  intersect(pattern_mod,pattern_hd)
#shared_pattern_sev<-  intersect(intersect(pattern_sev,pattern_mod), pattern_conv)
#shared_pattern_conv <-  intersect(intersect(pattern_conv, pattern_sev), pattern_mod)
#shared_pattern_hd <-  intersect(intersect(intersect(pattern_hd,pattern_sev), pattern_mod), pattern_conv)


####

# firts arg is df input
# second arg is default cluster set
# next three args are set remove clusters
# next three args are set remove clusters labels
setlabelsub <- function(dfinput, A, B, C, D, Blabel, Clabel, Dlabel){
  AB <- intersect(A,B)
  AC <- intersect(A,C)
  AD <- intersect(A,D)
  Amultiple <- union(union(intersect(AB,AC), intersect(AB,AD)), intersect(AC,AD))
  
  dfAB <- subset(dfinput, pattern %in% AB)
  dfAC <- subset(dfinput, pattern %in% AC)
  dfAD <- subset(dfinput, pattern %in% AD)
  dfAmultiple <- subset(dfinput, pattern %in% Amultiple)
  dfrebind <- rbind(dfAB, dfAC, dfAD, dfAmultiple)
  dfrebind$type <- c(rep(Blabel, dim(dfAB)[1]), 
                     rep(Clabel, dim(dfAC)[1]), 
                     rep(Dlabel, dim(dfAD)[1]), 
                     rep("Multiple", dim(dfAmultiple)[1]))
  return(dfrebind)
}


####



dfplt <- df_hd
dfsubunique <- setlabelsub(df_hd, pattern_hd, pattern_sev, pattern_mod, pattern_conv, "Severe", "Moderate", "Convalescent")
dfsubunique$type <- factor(dfsubunique$type, levels = c("Severe", "Moderate", "Convalescent", "Multiple"))
dfsubunique <- dfsubunique[rev(1:dim(dfsubunique)[1]),]
plt <- ggplot(dfplt, aes(x=freq , y=logpgen, size=nsub)) + geom_point(shape = 1,stroke = 1, color="darkgray") #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt  + labs(x = "Median cluster frequency", y = "-log2( Median cluster pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(6,12))
plt <- plt + scale_x_continuous(breaks = pretty_breaks(), limits = c(0,0.015))
plt <- plt + scale_y_continuous(breaks = pretty_breaks(), limits = c(15,50))
plt <- plt + geom_point(data=dfsubunique, aes(x=freq , y=logpgen, size=nsub, color = type), shape = 1,stroke = 1)
plt <- plt + scale_colour_manual(values = c("#bb2205", "#f6830f", "#2d6187", "#9d0191"))
plt <- plt + geom_text_repel(data=subset(dfsubunique, logpgen > 30), mapping=aes(x=freq , y=logpgen, size=nsub, color=type, label =pattern), size = 12, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
pdf("clusterplots/clusterplots_hd.pdf", height = 8, width = 10)
plt
dev.off()


dfplt <- df_conv
dfsubunique <- setlabelsub(df_conv, pattern_conv, pattern_sev, pattern_mod, pattern_hd, "Severe", "Moderate", "HD")
dfsubunique$type <- factor(dfsubunique$type, levels = c("Severe", "Moderate", "HD", "Multiple"))
dfsubunique <- dfsubunique[rev(1:dim(dfsubunique)[1]),]
plt <- ggplot(dfplt, aes(x=freq , y=logpgen, size=nsub)) + geom_point(shape = 1,stroke = 1, color="darkgray") #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt  + labs(x = "Median cluster frequency", y = "-log2( Median cluster pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(6,12))
plt <- plt + scale_x_continuous(breaks = pretty_breaks(), limits = c(0,0.015))
plt <- plt + scale_y_continuous(breaks = pretty_breaks(), limits = c(15,50))
plt <- plt + geom_point(data=dfsubunique, aes(x=freq , y=logpgen, size=nsub, color = type), shape = 1,stroke = 1)
plt <- plt + scale_colour_manual(values = c("#bb2205", "#f6830f", "#5aa469", "#9d0191"))
plt <- plt + geom_text_repel(data=subset(dfsubunique, logpgen > 30), mapping=aes(x=freq , y=logpgen, size=nsub, color=type, label =pattern), size = 12, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
pdf("clusterplots/clusterplots_conv.pdf", height = 8, width = 10)
plt
dev.off()




dfplt <- df_mod
dfsubunique <- setlabelsub(df_mod, pattern_mod, pattern_sev, pattern_conv, pattern_hd, "Severe", "Convalescent", "HD")
dfsubunique$type <- factor(dfsubunique$type, levels = c("Severe", "Convalescent", "HD", "Multiple"))
dfsubunique <- dfsubunique[rev(1:dim(dfsubunique)[1]),]
plt <- ggplot(dfplt, aes(x=freq , y=logpgen, size=nsub)) + geom_point(shape = 1,stroke = 1, color="darkgray") #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt  + labs(x = "Median cluster frequency", y = "-log2( Median cluster pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(6,12))
plt <- plt + scale_x_continuous(breaks = pretty_breaks(), limits = c(0,0.015))
plt <- plt + scale_y_continuous(breaks = pretty_breaks(), limits = c(15,50))
plt <- plt + geom_point(data=dfsubunique, aes(x=freq , y=logpgen, size=nsub, color = type), shape = 1,stroke = 1)
plt <- plt + scale_colour_manual(values = c("#bb2205", "#2d6187", "#5aa469", "#9d0191"))
plt <- plt + geom_text_repel(data=subset(dfsubunique, logpgen > 30), mapping=aes(x=freq , y=logpgen, size=nsub, color=type, label =pattern), size = 12, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
pdf("clusterplots/clusterplots_mod.pdf", height = 8, width = 10)
plt
dev.off()


dfplt <- df_sev
dfsubunique <- setlabelsub(df_sev, pattern_sev,pattern_mod, pattern_conv, pattern_hd, "Moderate", "Convalescent", "HD")
dfsubunique$type <- factor(dfsubunique$type, levels = c("Moderate", "Convalescent", "HD", "Multiple"))
dfsubunique <- dfsubunique[rev(1:dim(dfsubunique)[1]),]
plt <- ggplot(dfplt, aes(x=freq , y=logpgen, size=nsub)) + geom_point(shape = 1,stroke = 1, color="darkgray") #+ scale_colour_manual(values = c("red", "gray", "gray90"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=25), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=15), legend.text=element_text(size=15))
plt <- plt  + labs(x = "Median cluster frequency", y = "-log2( Median cluster pGen score )")
plt <- plt + scale_size(breaks = 1:5, range = c(6,12))
plt <- plt + scale_x_continuous(breaks = pretty_breaks(), limits = c(0,0.015))
plt <- plt + scale_y_continuous(breaks = pretty_breaks(), limits = c(15,50))
plt <- plt + geom_point(data=dfsubunique, aes(x=freq , y=logpgen, size=nsub, color = type), shape = 1,stroke = 1)
plt <- plt + scale_colour_manual(values = c("#f6830f", "#2d6187", "#5aa469", "#9d0191"))
plt <- plt + geom_text_repel(data=subset(dfsubunique, logpgen > 30), mapping=aes(x=freq , y=logpgen, size=nsub, color=type, label =pattern), size = 12, box.padding = 0.2, point.padding = 0.6, show.legend = FALSE, seed = 17)
pdf("clusterplots/clusterplots_sev.pdf", height = 8, width = 10)
plt
dev.off()





#write.table(intersect(intersect(intersect(downgene1$gene, downgene2$gene), downgene3$gene), downgene4$gene), "venn_TINKvSPNK_sig100_down_intersect.txt", row.names = F, col.names = F, quote = F)
#write.table(union(union(union(downgene1$gene, downgene2$gene), downgene3$gene), downgene4$gene), "venn_TINKvSPNK_sig100_down_union.txt", row.names = F, col.names = F, quote = F)



#dfconv <- df
#dfmod <- df
#dfsev <- df
#dfconv <- dfconv[order(dfconv$logpgen, decreasing = T),]
#dfmod <- dfmod[order(dfmod$logpgen, decreasing = T),]
#dfsev <- dfsev[order(dfsev$logpgen, decreasing = T),]


### plotting

### LEGACY CODE
#shared_pattern_sev <-  intersect(pattern_mod,pattern_sev)
#shared_pattern_conv <-  intersect(pattern_mod,pattern_conv)
#shared_pattern_hd <-  intersect(pattern_mod,pattern_hd)
#shared_pattern_multiple <-  intersect(intersect(pattern_mod,pattern_conv), shared_pattern_sev)
#dfsubunique_sev <- subset(dfplt, pattern %in% shared_pattern_sev)
#dfsubunique_conv <- subset(dfplt, pattern %in% shared_pattern_conv)
#dfsubunique_all <- subset(dfplt, pattern %in% shared_pattern_all)
#dfsubunique <- rbind(dfsubunique_sev, dfsubunique_conv, dfsubunique_all)
#dfsubunique$type <- c(rep("Severe", dim(dfsubunique_sev)[1]), 
#                      rep("Convalescent", dim(dfsubunique_conv)[1]), 
#                      rep("All", dim(dfsubunique_all)[1]))

