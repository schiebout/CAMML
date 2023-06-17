hpca.se <- celldex::ImmGenData()

labs <- unique(hpca.se$label.main)
labs <- sort(labs)
labs <- labs[c(1,4,5,10,15,19)]
labs <- sort(labs)


counts=hpca.se@assays@data$logcounts
colcount <- hpca.se@assays@data$logcounts
#colnames(colcount) <- 1:1561
colnames(counts) <- hpca.se$label.main
colcount <- colcount[,which(colnames(counts) %in% labs)]
counts <- counts[,which(colnames(counts) %in% labs)]
length(colnames(colcount))

i <- 1
d <- DGEList(counts=exp(colcount), group = ifelse(colnames(counts) == labs[1],1,0))
d <- calcNormFactors(d)
d1 <- estimateCommonDisp(d, verbose=TRUE)
d1 <- estimateTagwiseDisp(d1)

et12 <- exactTest(d1, pair = c(1,2))
#p <- et12$table$PValue
gp <- et12$table#[which(p.adjust(p, method =  "bonferroni", n = length(p)) < .05),]
gp <- gp[order(gp[,1], decreasing = TRUE),]
r <- rownames(gp[gp[,1]>5,])
#r <- rownames(gp[c(1:50),])
gw <- (gp[gp[,1]>5,1])
#gw <- (gp[c(1:50),1])
v <- cbind(rep(labs[1], length(r)), r,gw)
gene.weight <- list()
gene.weight[[1]] <- gw


for (i in 2:length(labs)){
  d <- DGEList(counts=exp(colcount), group = ifelse(colnames(counts) == labs[i],1,0))
  d <- calcNormFactors(d)
  d1 <- estimateCommonDisp(d, verbose=TRUE)
  d1 <- estimateTagwiseDisp(d1)
  print(labs[i])
  et12 <- exactTest(d1, pair = c(1,2))
  #p <- et12$table$PValue
  gp <- et12$table#[which(p.adjust(p, method =  "bonferroni", n = length(p)) < .05),]
  gp <- gp[order(gp[,1], decreasing = TRUE),]
  r <- rownames(gp[gp[,1]>5,])
  #r <- rownames(gp[c(1:50),])
  gw <- (gp[gp[,1]>5,1])
  #gw <- (gp[c(1:50),1])
  print(length(r))
  v <- rbind(v,cbind(rep(labs[i], length(r)), r,gw))
  gene.weight[[i]] <- gw
}

df <- data.frame(v)
#saveRDS(df,"df.RDS")

pbmc <- readRDS("CITESeq_B16_IL21_10_100_30_3000_TRUE_FALSE_FALSE.RDS")
DefaultAssay(pbmc) <- "RNA"
pbmc[['SCT']] <- NULL
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
# pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = .25)
pbmc <- RunUMAP(pbmc, dims = 1:30)
UMAPPlot(pbmc)

# Get the entrez gene IDs that are mapped to an Ensembl ID
entrez2ensembl = mappedkeys(org.Mm.egSYMBOL2EG)
# Convert to a list
entrez2ensembl = as.list(org.Mm.egSYMBOL2EG[entrez2ensembl])
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


#df <- na.omit(df)
#df <- df[-which(df$V1 == "Pre-B_cell_CD34-"),]
table(df$V1)
# Create a gene.set.collection list of Ensembl IDs
(gene.set.names = unique(df$V1))

num.sets = length(gene.set.names)

gene.set.collection = list()
for (i in 1:num.sets) {
  gene.set.name = gene.set.names[i]
  gene.set.rows = which(df$V1 == gene.set.name)
  gene.set.ensembl.ids = df$ensembl_gene[gene.set.rows]
  gene.set.collection[[i]] = gene.set.ensembl.ids
}

(names(gene.set.collection) = gene.set.names)

feature.data = read.delim("../../../../Desktop/10x_CITEseq_9_15_2020/B16_outs/filtered_feature_bc_matrix/features.tsv.gz",
                          header = FALSE, stringsAsFactors = FALSE)
(ensembl.ids = (feature.data[,1]))
gene.names =(feature.data[,2])
#rownames(pbmc@assays$RNA@counts) <- gsub("\\..*","",rownames(pbmc@assays$RNA@counts))
genes.after.QC = (rownames(pbmc@assays$RNA))
indices.to.keep = unlist(sapply(genes.after.QC, function(x) {which(gene.names == x)[1]}))
ensembl.ids = ensembl.ids[indices.to.keep]

gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids,
                                              gene.set.collection=gene.set.collection)

df <- df[!is.na(df$ensembl_gene),]

gene.w <- list()
for (j in 1:length(gene.set.collection)){
  gwi <- c()
  for (i in 1:length(names(gene.set.collection[[j]]))){
    gwi <- c(gwi, which(df$ensembl_gene==names(gene.set.collection[[j]][i])))
  }
  gwi <- intersect(gwi, which(df$V1 == names(gene.set.collection)[j]))
  gene.w[[j]] <- (as.double(df[gwi,3]))
}

pbmc = vamForSeurat(seurat.data=pbmc, gene.weights = gene.w,
                    gene.set.collection=gene.set.collection,
                    center=F, gamma=T, sample.cov=F, return.dist=T)
#set.seed(22)
pbmc.markers = FindAllMarkers(pbmc, assay="VAMcdf", only.pos = TRUE, logfc.threshold = .01)
#set.seed(22)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
#set.seed(22)
DefaultAssay(object = pbmc) = "VAMcdf"
top.pathways <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#set.seed(22)
DoHeatmap(pbmc, slot="data", features = top.pathways$gene, size=2, label=T,disp.min=0 ,disp.max = 1, angle = 90) +
  scale_fill_gradient(low = "black", high = "lightskyblue1") +
  #theme(axis.text.y = element_text(angle = 135, hjust=.25)) +
  labs(fill='Score\n') + NoLegend()

###ChIMP
cd4 <- kmeans(pbmc@assays$ADT@data[1,], centers = 2)
if (median(pbmc@assays$ADT@data[1,][cd4$cluster==1]) > median(pbmc@assays$ADT@data[1,][cd4$cluster==2])){
  cd4$cluster[cd4$cluster == 2] <- 0
} else{cd4$cluster <- cd4$cluster-1}

cd8 <- kmeans(pbmc@assays$ADT@data[2,], centers = 2)
if (median(pbmc@assays$ADT@data[2,][cd8$cluster==1]) > median(pbmc@assays$ADT@data[2,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[6,] <- pbmc@assays$VAMcdf@data[6,]* pmax(cd4$cluster, cd8$cluster)

cd8 <- kmeans(pbmc@assays$ADT@data[9,], centers = 2)
if (median(pbmc@assays$ADT@data[9,][cd8$cluster==1]) > median(pbmc@assays$ADT@data[9,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[1,] <- pbmc@assays$VAMcdf@data[1,]*cd8$cluster


cd8 <- kmeans(pbmc@assays$ADT@data[10,], centers = 2)
if (median(pbmc@assays$ADT@data[10,][cd8$cluster==1]) > median(pbmc@assays$ADT@data[10,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[2,] <- pbmc@assays$VAMcdf@data[2,]*cd8$cluster


cd8 <- kmeans(pbmc@assays$ADT@data[4,], centers = 2)
if (median(pbmc@assays$ADT@data[4,][cd8$cluster==1]) < median(pbmc@assays$ADT@data[4,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[3,] <- pbmc@assays$VAMcdf@data[3,]*cd8$cluster

cd8 <- kmeans(pbmc@assays$ADT@data[8,], centers = 2)
if (median(pbmc@assays$ADT@data[8,][cd8$cluster==1]) > median(pbmc@assays$ADT@data[8,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[4,] <- pbmc@assays$VAMcdf@data[4,]*cd8$cluster

cd8 <- kmeans(pbmc@assays$ADT@data[5,], centers = 2)
if (median(pbmc@assays$ADT@data[5,][cd8$cluster==1]) > median(pbmc@assays$ADT@data[5,][cd8$cluster==2])){
  cd8$cluster[cd8$cluster == 2] <- 0
} else{cd8$cluster <- cd8$cluster-1}

pbmc@assays$VAMcdf@data[5,] <- pbmc@assays$VAMcdf@data[5,]*cd8$cluster

DoHeatmap(pbmc, slot="data", features = top.pathways$gene, size=2, label=T,disp.min=0 ,disp.max = 1, angle = 90) +
  scale_fill_gradient(low = "black", high = "lightskyblue1") +
  #theme(axis.text.y = element_text(angle = 135, hjust=.25)) +
  labs(fill='Score\n') + NoLegend()


doubletpbmc <- pbmc


#isolate T and M
fer <- c()
fert <- c()
for (i in 1:length(colnames(doubletpbmc))){
  fer[i] <- sum(doubletpbmc@assays$VAMcdf[4,i] > .75 & doubletpbmc@assays$VAMcdf[6,i] > .75)
  if (fer[i] > 0){
    fert <- c(fert, i)
  }
}

doubletpbmc <- doubletpbmc[,fert]
doubletpbmc

#make fakers
fer <- c()
fert <- c()
for (i in 1:length(colnames(pbmc))){
  fer[i] <- sum(pbmc@assays$VAMcdf[,i] > 0.5)
  fert <- c(fert, rownames(pbmc@assays$VAMcdf[which(pbmc@assays$VAMcdf[,i] > 0.5),i]))
}

plot(fer, pch = 20)
table(fer)

tees <- c()

for (i in 1:length(colnames(pbmc))){
  if (fer[i] == 1){
    if (as.numeric(pbmc@assays$VAMcdf[6,i]) > .75){
      tees <- c(tees,i)
    }
  }
}

mac <- c()

for (i in 1:length(colnames(pbmc))){
  if (fer[i] == 1){
    if (as.numeric(pbmc@assays$VAMcdf[4,i]) > .75){
      mac <- c(mac,i)
    }
  }
}


set.seed(22)
tees <- sample(tees, length(tees), replace = F)

set.seed(22)
mac <- sample(mac, length(mac), replace = F)

doublet <- rowSums(cbind(pbmc@assays$RNA@counts[,tees[1]],pbmc@assays$RNA@counts[,mac[1]]))

for (i in 2:min(length(tees),length(mac))){
 doublet <- cbind(doublet,rowSums(cbind(pbmc@assays$RNA@counts[,tees[i]],pbmc@assays$RNA@counts[,mac[i]])))
}


colnames(doublet) <- 1:(min(length(tees),length(mac)))


true.data <- doubletpbmc@assays$RNA@counts
synth.data <- doublet

true <- CreateSeuratObject(counts = true.data, 
                           project = "True Doublets", min.cells = (ncol(true.data)*.5), min.features = 100)

synth <- CreateSeuratObject(counts = synth.data, 
                            project = "Synthetic Doublets", min.cells = (ncol(synth.data)*.5), min.features = 100)
combo = merge(x = true, y = c(synth))



ribosomal.pattern = "^Rp[sl]"
gene.names = rownames(combo@assays$RNA)
ribosomal.genes = grep(pattern = ribosomal.pattern, x = gene.names, value = TRUE)
message("Number of ribosomal genes: ", length(ribosomal.genes))
not.ribosomal.genes = which(!(gene.names %in% ribosomal.genes))
message("Number of non-ribosomal genes: ", length(not.ribosomal.genes))
combo = combo[not.ribosomal.genes,]

mt.pattern = "^mt-"
gene.names = rownames(combo@assays$RNA)
mt.genes = grep(pattern = mt.pattern, x = gene.names, value = TRUE)
message("Number of mt genes: ", length(mt.genes))
not.mt.genes = which(!(gene.names %in% mt.genes))
message("Number of non-mt genes: ", length(not.mt.genes))
combo = combo[not.mt.genes,]


combo <- NormalizeData(combo)
combo <- FindVariableFeatures(combo, selection.method = "vst", nfeatures = 2000)
combo <- ScaleData(combo)

combo <- RunPCA(combo)
combo <- FindNeighbors(combo, dims = 1:30)
combo <- FindClusters(combo, resolution = .25)

combo <- RunUMAP(combo, dims = 1:30)

UMAPPlot(combo, group.by = "orig.ident") + ggtitle("Identities")+
  theme(plot.title = element_text(hjust = 0)) 
UMAPPlot(combo) + ggtitle("Clusters")

cl <- which(data.frame(table(combo$orig.ident,Idents(combo))["Synthetic Doublets",]) == 0)-1

Idents(combo) <- ifelse(Idents(combo) %in% cl,"True","Synthetic")

comp.markers <- FindAllMarkers(combo, only.pos = T)

top.c.markers <- comp.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
DotPlot(combo, features = unique(top.c.markers$gene))+ coord_flip()
