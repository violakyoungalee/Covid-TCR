library(data.table)
setwd("~/Dropbox/MLGROUP/tcrcov/JP/kmers_cd4")
set.seed(20)

kmersize <- 6
kmermaster <- fread(paste0("su_cd4_",kmersize,"mers.csv"))
varexp <- apply(as.data.frame(kmermaster[,-1]), 1, var)
kmermaster <- kmermaster[order(varexp, decreasing = T),]
kmermaster <- kmermaster[c(1:20000),]
kmermaster <- na.omit(kmermaster)

kmernames <- kmermaster$Kmer

library(janitor)
su <- as.data.frame(t(kmermaster),index=TRUE)
su <- su %>% row_to_names(row_number = 1)
head(su)

su <- cbind(Patient_ID = rownames(su), su)
rownames(su) <- 1:nrow(su)
head(su)

setwd("~/Dropbox/MLGROUP/tcrcov/JP/kmers_cd4")
metadata <- fread("metadata.txt")
wtsamples <- metadata$Sample[which(metadata$`Who Ordinal Scale` == "")]
expsamples1 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "1"), which(metadata$`Who Ordinal Scale` == "1 or 2")), which(metadata$`Who Ordinal Scale` == "2"))]
expsamples2 <- metadata$Sample[union(which(metadata$`Who Ordinal Scale` == "3"), which(metadata$`Who Ordinal Scale` == "4"))]
expsamples3 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "5"), which(metadata$`Who Ordinal Scale` == "6")), which(metadata$`Who Ordinal Scale` == "7"))]
expsamples4 <- union(union(expsamples1, expsamples2), expsamples3)
expid1 <- "mild"
expid2 <- "moderate"
expid3 <- "severe"
expid4 <- "pooled"

su$d_cond <- ifelse(su$Patient_ID %in% expsamples3,expid3,ifelse(su$Patient_ID %in% expsamples2, expid2, ifelse(su$Patient_ID %in% expsamples1, expid1, "HD")))

write.csv(su,"su_cd4_6mers_cleaned.csv", row.names = FALSE)


# normalized
library(data.table)
setwd("~/Dropbox/MLGROUP/tcrcov/JP/kmers_cd8")
set.seed(20)

kmersize <- 7
kmermaster <- fread(paste0("su_cd8_",kmersize,"mers.csv"))
varexp <- apply(as.data.frame(kmermaster[,-1]), 1, var)
kmermaster <- kmermaster[order(varexp, decreasing = T),]
kmermaster <- kmermaster[c(1:20000),]
kmermaster <- na.omit(kmermaster)

kmernames <- kmermaster$Kmer

pcinput <- kmermaster[,-1]
pcinput <- t(t(pcinput)/colSums(pcinput))
pcinput <- as.data.frame(t(pcinput))

# pcinput is the normalized kmer dataset
colnames(pcinput) <- kmernames
head(pcinput)

su <- cbind(Patient_ID = rownames(pcinput), pcinput)
rownames(su) <- 1:nrow(su)
head(su)

setwd("~/Dropbox/MLGROUP/tcrcov/JP/kmers_cd8")
metadata <- fread("metadata.txt")
wtsamples <- metadata$Sample[which(metadata$`Who Ordinal Scale` == "")]
expsamples1 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "1"), which(metadata$`Who Ordinal Scale` == "1 or 2")), which(metadata$`Who Ordinal Scale` == "2"))]
expsamples2 <- metadata$Sample[union(which(metadata$`Who Ordinal Scale` == "3"), which(metadata$`Who Ordinal Scale` == "4"))]
expsamples3 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "5"), which(metadata$`Who Ordinal Scale` == "6")), which(metadata$`Who Ordinal Scale` == "7"))]
expsamples4 <- union(union(expsamples1, expsamples2), expsamples3)
expid1 <- "mild"
expid2 <- "moderate"
expid3 <- "severe"
expid4 <- "pooled"

su$d_cond <- ifelse(su$Patient_ID %in% expsamples3,expid3,ifelse(su$Patient_ID %in% expsamples2, expid2, ifelse(su$Patient_ID %in% expsamples1, expid1, "HD")))

write.csv(su,"su_cd8_7mers_normalized.csv", row.names = FALSE)
