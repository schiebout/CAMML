CAMML <- function(seurat, gene.weights, gene.set.collection,){
  #run VAM with ensembl gene set 
  #gene.set.collection = createGeneSetCollection(gene.ids=ensembl.ids, gene.set.collection=gene.set.collection)
  
  seurat = vamForSeurat(seurat.data=seurat, gene.weights = gene.weights, gene.set.collection=gene.set.collection, center=F, gamma=T, sample.cov=F, return.dist=T)
  
  DefaultAssay(object = pbmc) = "VAMcdf"
  
  #initalize all label options
  le <- length(seurat$orig.ident)
  celz <- list(c(1:le))
  celztop2 <- list(c(1:le))
  celzfold <- list(c(1:le))
  celz2 <- list(c(1:le))
  celz3 <- list(c(1:le))
  celz10p <- list(c(1:le))
  
  #find labels
  for (i in 1:le){
    #take data and label
    dif <- data.frame(seurat@assays$VAMcdf[,i])
    dif <- data.frame(dif[order(-dif[,1]), , drop = FALSE])
    #for cells with no scores, skip
    if (max(dif) == 0){
      celztop2[[i]]<- NA
      celzfold[[i]]<- NA
      celz10p[[i]] <- NA
      celz[[i]] <- "none"
      next
    }
    #top cells
    celz[i]<- rownames(dif)[1]
    dif2 <- data.frame(dif[-1,], row.names=rownames(dif)[-1])
    celz2[i] <- rownames(dif2)[1]
    celz3[i] <- rownames(dif2)[2]
    
    #top2
    celztop2[[i]] <- data.frame(dif[c(1:2),], row.names=rownames(dif)[c(1:2)])
    
    #fold change
    me <- mean(dif[,1])
    fold <- which(dif[,1] > me*2)
    celzfold[[i]] <- data.frame(dif[fold,],row.names=rownames(dif)[fold])
    
    #percent dif
    inde <- which(dif > max(dif*.9))
    celz10p[[i]] <- data.frame(dif[inde,],row.names=rownames(dif)[inde])
  }
  return(seurat, celz, celztop2, celzfold, celz10p)
}
