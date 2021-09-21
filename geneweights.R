#Remove NAs from ensembl translation
df <- df[!is.na(df$ensembl_gene),]

#in a for loop, find which gene weights have an associated ensembl name and build a list
gene.w <- list()
for (j in 1:length(gene.set.collection)){
  gwi <- c()
  for (i in 1:length(names(gene.set.collection[[j]]))){
    gwi <- c(gwi, which(df$ensembl_gene==names(gene.set.collection[[j]][i])))
  }
  gwi <- intersect(gwi, which(df$V1 == names(gene.set.collection)[j]))
  gene.w[[j]] <- (as.double(df[gwi,3]))
}

gene.w
