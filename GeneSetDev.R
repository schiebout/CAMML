# load data
celldata <- celldex::HumanPrimaryCellAtlasData()

#labels
labs <- unique(celldata$label.main)
labs <- sort(labs)

#cell counts
counts=celldata@assays@data$logcounts
colcount <- celldata@assays@data$logcounts

colnames(counts) <- celldata$label.main
colcount <- colcount[,which(colnames(counts) %in% labs)]
counts <- counts[,which(colnames(counts) %in% labs)]

#DE analysis
d <- DGEList(counts=exp(colcount), group = ifelse(colnames(counts) == labs[1],1,0))
d <- calcNormFactors(d)
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
et12 <- exactTest(d1, pair = c(1,2))
gp <- gp[order(gp[,1], decreasing = T),]

#keep cells > logfc 5
r <- rownames(gp[gp[,1]>5,])
gw <- (gp[gp[,1]>5,1])
v <- cbind(rep(labs[1], length(r)), r,gw)

#repeat for every cell typ
for (i in 2:length(labs)){
  d <- DGEList(counts=exp(colcount), group = ifelse(colnames(counts) == labs[i],1,0))
  d <- calcNormFactors(d)
  d1 <- estimateCommonDisp(d, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  et12 <- exactTest(d1, pair = c(1,2))
  
  gp <- gp[order(gp[,1], decreasing = T),]
  r <- rownames(gp[gp[,1]>5,])
  
  gw <- (gp[gp[,1]>5,1])
  
  v <- rbind(v,cbind(rep(labs[i], length(r)), r,gw))
}

#if intersecting with MSigDB C8, simple intersection calls were used

df <- data.frame(v)

#ensembl tranforming from VAM vignette

# Get the entrez gene IDs that are mapped to an Ensembl ID
entrez2ensembl = mappedkeys(org.Hs.egSYMBOL2EG)
# Convert to a list
entrez2ensembl = as.list(org.Hs.egSYMBOL2EG[entrez2ensembl])

msigdb.entrez.ids = (df$r)
num.ids = length(msigdb.entrez.ids)
msigdb.ensembl.ids = rep(NA, num.ids)

for (i in 1:num.ids) {
  entrez.id = msigdb.entrez.ids[i]
  id.index = (which(names(entrez2ensembl) == entrez.id))
  if (length(id.index > 0)) {
    # only use the first mapped ensembl id
    msigdb.ensembl.ids[i] =(entrez2ensembl[[id.index]][1])
  }
}

entrez2ensembl = mappedkeys(org.Mm.egENSEMBL)
# Convert to a list
entrez2ensembl = as.list(org.Mm.egENSEMBL[entrez2ensembl])
msigdb.entrez.ids <- msigdb.ensembl.ids
num.ids = length(msigdb.entrez.ids)
msigdb.ensembl.ids = rep(NA, num.ids)

for (i in 1:num.ids) {
  entrez.id = msigdb.entrez.ids[i] 
  id.index = (which(names(entrez2ensembl) == entrez.id))
  if (length(id.index > 0)) {
    # only use the first mapped ensembl id
    msigdb.ensembl.ids[i] =(entrez2ensembl[[id.index]][1])
  }
}

df$ensembl_gene = msigdb.ensembl.ids

#set up gene set collection from VAM

(gene.set.names = unique(df$V1))

gene.set.collection = list()
for (i in 1:num.sets) {
  gene.set.name = gene.set.names[i]
  gene.set.rows = which(df$V1 == gene.set.name)
  gene.set.ensembl.ids = df$ensembl_gene[gene.set.rows]
  gene.set.collection[[i]] = gene.set.ensembl.ids
}

(names(gene.set.collection) = gene.set.names)

#find relevant ensembl IDs
feature.data = read.delim("yourdata/features/genes.tsv.gz",
                          header = FALSE, stringsAsFactors = FALSE)
ensembl.ids = feature.data[,1]
gene.names = feature.data[,2]
genes.after.QC = rownames(seurat@assays$RNA@counts)
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]









