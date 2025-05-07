#DE detection#########
#library(RaceID)
library(OmnipathR)
library(magrittr)
library(BiocParallel)
#library(MAST)
library(reshape2)
library(trend)
library(Seurat)
library(SeuratData)
#library(SeuratDisk)
library(ggrepel)
library(ggplot2)
library(gtable)
library(ggpubr)
library(gridExtra)
library(gridGraphics)
library(grid)
library(presto)
library(fgsea)
library(readxl)
suppressMessages(library(tidyverse))
#library(destiny)
library(fitdistrplus)


library(bayNorm)
library(uwot)
library(Rtsne)
library(irlba)
library(igraph)

library(readbitmap)
library(grid)
library(tibble)
library(readr)

library(rdist)
library(dplyr)
library(plyr)

library(cowplot)
library(crayon)
#library(scatterpie)
library(viridis)
library(RColorBrewer)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)

library(data.table)
library(abind)

library(stringr)

library(future)
library(pbapply)
#library(limma)
library(future.apply)
library(BioQC)

library(purrr)

library(scuttle)
#library(scran)
#library(scater)
#font types#
#library(extrafont)
#useful when use other font types
# library(extrafont) 
# font_import()
# loadfonts(device = "win")

#library(OneR)
#library(scalop) #scalop::sigScores

#library(ggalluvial)

##ternary####
#library(ggtern)
library(vcd)


library(doRNG)
#trajectory
#library(slingshot)
library(SingleCellExperiment)
library(mclust, quietly = TRUE)
library(rhdf5)
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

#library(infercnv)

#rotate coordinates
#library(spdep)
#cdata<-Rotation(cdata,270*pi/180)

#load for stRNA-seq
library(rjson)
library(rhdf5)

#load plot utils#########
#source('D:/utils_plot.R', echo=FALSE)
#source('/home/clustor2/ma/w/wt215/BTC_project/BTC_R/utils_plot.R', echo=FALSE)
source('/home/clustor2/ma/w/wt215/RFILES/utils_plot.R', echo=FALSE)

read_excel_allsheets <- function(filename, tibble = FALSE,...) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X,...))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#fast PCA######

fastPCA_wt<-function(object,weight.by.var=TRUE,npcs=30,...){
  #object: rows: cells columns: genes
  pca.results <- irlba(A = (object), nv = npcs,...)
  #Weight the cell embeddings by the variance of each PC
  
  if (weight.by.var) {
    cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  }else{
    cell.embeddings <- pca.results$u
  }
  sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
  
  rownames(cell.embeddings)<-rownames(object)
  return(list(
    cell.embeddings=cell.embeddings,
    pca.results=pca.results
  ))
}



#seurat norm#######
Seurat_norm<-function(Data,scale.factor = 10000,scale.max = 10){
  ff=LogNormalize(Data,scale.factor = scale.factor)
  ff<-apply(ff,1,function(x){(x-mean(x))/sd(x)})
  ff<-t(ff)
  ff[ff>scale.max]<-scale.max
  return(ff)
}


####dimension reduction DIY############
DimReduce_self<-function(
  object,
  npcs=50,
  approx=TRUE,
  n_neighbors = 15,
  n_components = 2,
  numHVG=2000,
  method='umap',
  use_seurat=FALSE,
  norm_seurat=FALSE,
  regress_dat=NULL,
  regress_var=NULL,
  resolution_num=0.5,
  resolution_vec=NULL,
  weight.by.var=TRUE,
  markers_seurat=FALSE,
  dim.embed=2,
  pca_reduced=20,
  pca_neighbour=10,...){
  #object: rows: cells columns: genes
  message(
    cyan$bold(
      '
method: umap, tsne, tsne_seurat.
set resolution when using tsne_seurat
  '))
  
  if(norm_seurat){
    object=t(Seurat_norm(t(object)))
  }
  
  
  npcs <- min(npcs, nrow(x = object) - 1)  

  if(!use_seurat){
    if(approx){
      pca.results <- irlba(A = (object), nv = npcs,...)
      #Weight the cell embeddings by the variance of each PC
      
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
      }else{
        cell.embeddings <- pca.results$u
      }
      sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
      
      rownames(cell.embeddings)<-rownames(object)
      
    } else{
      pca.results <- prcomp(x = (object), rank. = npcs, ...)
      
      if (weight.by.var) {
        cell.embeddings <-pca.results$x %*% diag(pca.results$sdev[1:npcs]^2)
      }else{
        cell.embeddings <- pca.results$x
      }

      sdev <- pca.results$sdev
      
      rownames(cell.embeddings)<-rownames(object)
    }
    
    
    if(method=='umap'){
      
      uout<-uwot::umap(cell.embeddings,n_neighbors = n_neighbors,n_components = n_components,...)
      
      if(class(uout)[1]!='list'){
        rownames(uout)<-rownames(object)
      }
      return(list(
        pca.results=pca.results,
        pca.cell.embeddings=cell.embeddings,
        uout=uout)
      )
    } 
    
    if(method=='tsne'){
      #set.seed(1)
      #For tsne, output dims should be either 1, 2 or 3
      tsne_out<-Rtsne(cell.embeddings,initial_dims=npcs,dims = n_components,...)

      return(list(
        pca.results=pca.results,
        pca.cell.embeddings=cell.embeddings,
        tsne_out=tsne_out)
      )
    }
  }

  if(use_seurat){
    #For tsne, output dims should be either 1, 2 or 3
    #row gene, column cell
    temp <- CreateSeuratObject(counts= t(object))
    
    #qq<-t(t(qq)/colSums(qq))*10000 normalization
    temp <- NormalizeData(object = temp, normalization.method = "LogNormalize", scale.factor = 10000)
    temp <- FindVariableFeatures(object = temp, mean.function = ExpMean, dispersion.function = LogVMR,nfeatures=numHVG)
    
    if(!is.null(regress_dat) & !is.null(regress_var)){
      
      for(reind in 1:dim(regress_dat)[2]){
        temp@meta.data[colnames(regress_dat)[reind]]<-regress_dat[,reind]
      }
      
      #temp@meta.data['percent.mito']<-mitper
      temp <- ScaleData(object = temp, vars.to.regress = regress_var)
    } else{
      #qq<-apply(qq,1,function(x){(x-mean(x))/sd(x)})
      temp <- ScaleData(object = temp, do.scale = TRUE,do.center = TRUE)
    }
    
    temp <- RunPCA(object = temp, pc.genes =temp@assays$RNA@var.features, do.print = TRUE, pcs.print = 1:5, genes.print = 5,npcs = npcs)
  
    if(method=='tsne' | method=='tsne_seurat'){
      temp <- RunTSNE(object = temp, dims = 1:pca_reduced,reduction='pca',dim.embed=dim.embed,check_duplicates = FALSE)
      temp <- FindNeighbors(temp, dims = 1:pca_neighbour)
      temp <- FindClusters(object = temp,resolution = resolution_num, verbose= 0)
      Embeddings<-temp@reductions$tsne@cell.embeddings
      Reduction_out<-temp@reductions$tsne
      
    }
    if(method=='umap' | method=='umap_seurat'){
      temp <- RunUMAP(object = temp, dims = 1:pca_reduced,reduction='pca',n.components=dim.embed)
      temp <- FindNeighbors(temp, dims = 1:pca_neighbour)
      temp <- FindClusters(object = temp, resolution = resolution_num, verbose= 0)
      Embeddings<-temp@reductions$umap@cell.embeddings
      Reduction_out<-temp@reductions$umap
    }
    
    scale.data<-temp@assays$RNA@scale.data
    pca_out<-temp@reductions$pca
    itentities<-Idents(temp)
    
    identities_list<-list()
    #try different resolutions
    if(!is.null(resolution_vec)){
      for(iden_ind in 1:length(resolution_vec)){
        print(paste('Resolution:',resolution_vec[iden_ind],sep=''))
        temp <- FindClusters(object = temp, resolution = resolution_vec[iden_ind], verbose= 0)
        identities_list[[iden_ind]]<-Idents(temp)
      }
      names(identities_list)<-resolution_vec
    }

    
    ##Final return
  
    if(markers_seurat){
      Allmarkers <- FindAllMarkers(temp,...) 
      return(list(
        scale.data=scale.data,
        Embeddings=Embeddings,
        itentities=itentities,
        HVG_seurat=temp@assays$RNA@var.features,
        Reduction_out=Reduction_out,
        pca_out=pca_out,
        Allmarkers=Allmarkers,
        identities_list=identities_list
      )
      )
      
    } else{
      return(list(
        scale.data=scale.data,
        Embeddings=Embeddings,
        itentities=itentities,
        HVG_seurat=temp@assays$RNA@var.features,
        Reduction_out=Reduction_out,
        pca_out=pca_out,
        identities_list=identities_list
      )
      )
    }


    
  }
  

}

##rerun clustering########
# x.seurat <- CreateSeuratObject(counts =bay_allall$Bay_out,assay = 'RNA')
# x.seurat@reductions$tsne<-dimenout_allall_tsne$Reduction_out
# x.seurat@reductions$pca<-dimenout_allall_tsne$pca_out
# 
# x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
# x.seurat <- FindClusters(object = x.seurat,resolution = 0.6, verbose= 0)


####clustering_DIY###########



clustering_DIY<-function(data,PCAforNN=NULL,assay_name='bayNorm',uout=NULL,dims_range=1:dim(uout)[2],resolution=0.1,annoy.metric="euclidean",algorithm=1,...){
  
  message(
    cyan$bold(
      '
annoy.metric: euclidean, cosine, manhattan, and hamming
algorithm: Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
  '))


  
  x.seurat_self <- CreateSeuratObject(counts= data,assay = assay_name)
  
  if(!is.null(uout)){
    rownames(uout)<-colnames(data)
    x.seurat_self@reductions['umap']<-CreateDimReducObject(
      embeddings = uout,
      key = 'UMAP',
      assay = assay_name,
      global = TRUE
    )
  }

  
  if(is.null(PCAforNN)){
    x.seurat_self <- FindNeighbors(x.seurat_self, dims = dims_range,annoy.metric = annoy.metric,reduction ='umap',...)
  } else{
    rownames(PCAforNN)<-colnames(data)
    x.seurat_self@reductions['pca']<-CreateDimReducObject(
      embeddings = PCAforNN,
      key = 'PCA',
      assay = assay_name,
      global = TRUE
    )
    x.seurat_self <- FindNeighbors(x.seurat_self, dims = dims_range,annoy.metric = annoy.metric,reduction ='pca',...)
  }

  
  x.seurat_self <- FindClusters(x.seurat_self, resolution = resolution,algorithm=algorithm,...)
  
  cell_identity<-as.factor(Idents(x.seurat_self))
  return(list(cell_identity=cell_identity,x.seurat_self=x.seurat_self))
}


FindAllMarkers_DIY<-function(x.seurat_self,assay_name='bayNorm',...){
  
#,only.pos=TRUE,min.pct=0.1,logfc.threshold = 0.25,test.use='wilcox',max.cells.per.ident=Inf,
  markers_out <- FindAllMarkers(x.seurat_self,assay=assay_name,...)
  
  
  markers_list<-list()
  
  for(i in 1:length(unique(markers_out$cluster))){
    markers_list[[i]]<-markers_out[which(markers_out$cluster==unique(markers_out$cluster)[i]),]
  }
  names(markers_list)<-paste('Markers:',(unique(markers_out$cluster)),sep='')
  
  message(
    cyan$bold(
      '
#markers_out %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  '))
  
  
  
  return(markers_list)
}


FindMarkers_DIY<-function(x.seurat_self,assay_name='bayNorm',only.pos=TRUE,min.pct=0.1,logfc.threshold = 0.25,test.use='wilcox',max.cells.per.ident=Inf,...){
  
  markers_out <- FindMarkers(x.seurat_self,assay=assay_name,...)
  
  
  return(markers_out)
}












# assay_name<-'bayNorm'
# 

#pca
# npcs = 50
# npcs <- min(npcs, nrow(x = object) - 1)
# pca.results <- irlba(A = t(bayout$Bay_out), nv = npcs)
# pca_embeddings<-pca.results$u %*% diag(pca.results$d)
# rownames(pca_embeddings)<-colnames(object)
# feature.loadings <- pca.results$v
# sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
# total.variance<-sum(apply(t(x = bayout$Bay_out),1,var))


# x.seurat_self <- CreateSeuratObject(counts= bayout$Bay_out,assay = assay_name)
# # x.seurat_self@reductions['pca']<-CreateDimReducObject(
# #   embeddings = pca_embeddings,
# #   loadings = feature.loadings,
# #   assay = assay_name,
# #   stdev = sdev,
# #   key = 'pca',
# #   misc = list(total.variance = total.variance)
# # )
# x.seurat_self@reductions['umap']<-CreateDimReducObject(
#   embeddings = uout,
#   key = 'UMAP',
#   assay = assay_name,
#   global = TRUE
# )
# x.seurat_self <- FindNeighbors(x.seurat_self, dims = 1:10,annoy.metric = "euclidean",reduction ='umap')
# x.seurat_self <- FindClusters(x.seurat_self, resolution = 0.1)
# 
# cell_identity<-as.factor(Idents(x.seurat_self))
# num_classes<-length(unique(cell_identity))







########SCnorm_runMAST3##########
# SCnorm_runMAST3 <- function(Data, NumCells) {
SCnorm_runMAST3 <- function(G1,G2) {
  
  
  if(length(dim(G1))==2){
    #resultss<-SCnorm_runMAST(Data, NumCells)
    resultss<-SCnorm_runMAST(G1,G2)
  } else if(length(dim(G1))==3){
    
    library(foreach)
    resultss<- foreach(sampleind=1:dim(G1)[3],.combine=cbind)%do%{
      print(sampleind)
      
      #qq<-SCnorm_runMAST(Data[,,sampleind], NumCells)
      qq<-SCnorm_runMAST(G1[,,sampleind], G2[,,sampleind])
      return(qq$adjpval)
    }
  }
  message(cyan$bold('Note: logFC: G2/G1'))
  
  return(resultss) 
}
# SCnorm_runMAST <- function(Data, NumCells) {
#   NA_ind<-which(rowSums(Data)==0)
#   Data = as.matrix(log2(Data+1))
#   G1<-Data[,seq(1,NumCells[1])]
#   G2<-Data[,-seq(1,NumCells[1])]
SCnorm_runMAST <- function(G1,G2) {
  G1 = as.matrix(log2(G1+1))
  G2 = as.matrix(log2(G2+1))

  qq_temp<- rowMeans(G2)-rowMeans(G1)
  NumCells<-c(dim(G1)[2],dim(G2)[2])
  
  Data<-cbind(G1,G2)
  if(sum(duplicated(colnames(Data)))>0){
    colnames(G1)<-paste('G1_',colnames(G1),sep='')
    Data<-cbind(G1,G2)
  }
  

  #qq_temp[NA_ind]=NA
  numGenes = dim(G1)[1]
  datalong = melt(Data)
  Cond = c(rep("c1", NumCells[1]*numGenes), rep("c2", NumCells[2]*numGenes))
  dataL = cbind(datalong, Cond)
  colnames(dataL) = c("gene","cell","value","Cond")
  dataL$gene = factor(dataL$gene)
  dataL$cell = factor(dataL$cell)
  vdata = FromFlatDF(dataframe = dataL, idvars = "cell", primerid = "gene", measurement = "value",  id = numeric(0), cellvars = "Cond", featurevars = NULL,  phenovars = NULL)
  
  
  zlm.output = zlm(~ Cond, vdata, method='glm', ebayes=TRUE)
  zlm.lr = lrTest(zlm.output, 'Cond')
  gpval = zlm.lr[,,'Pr(>Chisq)']
  adjpval = p.adjust(gpval[,1],"BH") ## Use only pvalues from the continuous part.
  adjpval = adjpval[rownames(Data)]
  return(list(adjpval=adjpval, logFC_re=qq_temp,pval_ori=gpval[rownames(Data),1]))
}

Find_Marker_self<-function(Data,cluster_labels,Ncores=4,samplee=c(1,1),seed=12300,path_save=NULL){
  
  unique_labels<-unique(cluster_labels)
  
  if(length(unique_labels)==2){
    ll=1
  } else{
    ll=length(unique_labels)
  }
  
  #if(Ncores==0){
    #mc.cores = Ncores
    options(mc.cores=Ncores)
    
    markers_list<-list()
    for(i in 1:ll){
      print(paste('Cluster:', i,sep=''))
      g1=which(cluster_labels==unique_labels[i])
      g2=which(cluster_labels!=unique_labels[i])
      
      #markers_list[[i]]<- SCnorm_runMAST3(Data[,c(g1,g2)], NumCells=c(length(g1),length(g2)))
      gg1<-Data[,g1] 
      gg2<-Data[,g2]
      
      
      if(sum(samplee)==2){
        markers_list[[i]]<- SCnorm_runMAST3(gg1,gg2)
      } else{
        set.seed(seed)
        k1<-sample(seq(1,dim(gg1)[2]),round(dim(gg1)[2]*samplee[1]))
        k2<-sample(seq(1,dim(gg2)[2]),round(dim(gg2)[2]*samplee[2]))
        markers_list[[i]]<- SCnorm_runMAST3(gg1[,k1],gg2[,k2])
        
      }
      
      if(!is.null(path_save)){
        ff=markers_list[[i]]
        save(ff,file=paste(path_save,paste('Markers_',unique_labels[i],sep=''),'.RData',sep=''))
      }
      
    }
    
    #parallel MAST
  # } else{
  #   
  #   cluster <- makeCluster(NCores, type = "SOCK")
  #   registerDoSNOW(cluster)
  #   getDoParWorkers()
  #   
  #   iterations <- ll
  #   pb <- txtProgressBar(max = iterations, style = 3)
  #   progress <- function(n) setTxtProgressBar(pb, n)
  #   opts <- list(progress = progress)
  #   
  #   markers_list=foreach(i = 1:ll,.options.snow = opts,.packages = c('MAST','reshape2','crayon'),.export = c('SCnorm_runMAST3','SCnorm_runMAST'))%dopar%{
  #     
  #     g1=which(cluster_labels==unique_labels[i])
  #     g2=which(cluster_labels!=unique_labels[i])
  #     
  #     
  #     ff<- SCnorm_runMAST3(Data[,g1], Data[,g2])
  #     return(ff)
  #     
  #   }
  #   
  #   close(pb)
  #   stopCluster(cluster)
  # }
  
    if(length(unique_labels)==2){
      names(markers_list)<-paste('Markers:',unique_labels[1],sep='')
    } else{
      names(markers_list)<-paste('Markers:',unique_labels,sep='')
    }
  
  
  
  
  message(cyan$bold('
Note: logFC: G2/G1 \n
Markers of the i^th cluster: 
i=1
which(Markers_list[[i]]$adjpval<0.05 & Markers_list[[i]]$logFC_re<0)
                    
                    '))
  
  
  return(markers_list)
}




#wilcox from Seurat#####
WilcoxDETest_wt<-function(
  data.use,
  cells.1,
  cells.2,
  verbose = TRUE,
  using_method='BioQC',
  NCores=4,
  ...
) {
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  j <- seq_len(length.out = length(x = cells.1))
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  
  if(using_method=='Seurat'){
    p_val <- my.sapply(
      X = 1:nrow(x = data.use),
      FUN = function(x) {
        return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, statistics = data.use[x, ])), 1))
      }
    )
  }else if(using_method=='foreach'){
    
    
    cluster <- makeCluster(NCores, type = "SOCK")
    registerDoSNOW(cluster)
    getDoParWorkers()
    
    iterations <- dim(data.use)[1]
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    p_val <- foreach(
      Geneind = 1:iterations,
      .combine = c,
      .options.snow = opts,
      .packages = c('Matrix')) %dopar% {
        
        wout<-wilcox.test(data.use[Geneind, cells.1],data.use[Geneind, cells.2])
        return(wout$p.value)
      }
    close(pb)
    stopCluster(cluster)
    
  } else if(using_method=='BioQC'){
    p_val <-wmwTest(t(as.matrix(data.use)), c(rep(TRUE,length(cells.1)),rep(FALSE,length(cells.2))), valType="p.two.sided")
  }
  
  
  
  return(data.frame(p_val, row.names = rownames(x = data.use)))
}

WilcoxDETest_markers_wt<-function(data.use,clusters){
  classes<-levels(clusters)
  
  
  wout_list<-list()
  for(i in 1:length(classes)){
    print(i)
    
    cells.1<-which(clusters==classes[i])
    cells.2=seq(1,dim(data.use)[2])[-cells.1]
    wout<-WilcoxDETest_wt(
      data.use=data.use,
      cells.1=cells.1,
      cells.2=cells.2,
      verbose = TRUE,
      using_method='BioQC',
      NCores=4
    )
    G1=data.use[,cells.1]
    G2=data.use[,cells.2]
    G1 = as.matrix(log2(G1+1))
    G2 = as.matrix(log2(G2+1))
    qq_temp<- rowMeans(G2)-rowMeans(G1)
    
    wout['padj']<-p.adjust(wout$p_val,method = 'BH')
    wout['logFC']<-qq_temp
    wout<-wout[order(wout$logFC,decreasing=FALSE),]
    wout_list[[i]]<-wout
  }
  names(wout_list)<-classes
  return(wout_list)
  
}








some_genes_fun<-function(wout_H,thres_pval=0.05,thres_logFC=(-0.005),thres_logFC_low=NULL,limitnum=20,returnlist=FALSE){
  wout_TOP_H_l<-lapply(seq(1,length(wout_H)),function(i){
    x=wout_H[[i]]
    genes=rownames(x)[which(x$padj<thres_pval & x$logFC<thres_logFC)]
    if(!is.null(thres_logFC_low)){
      genes=rownames(x)[union(which(x$padj<thres_pval & x$logFC<thres_logFC),which(x$padj<thres_pval & x$logFC>thres_logFC_low))]
    }
    
    if(length(genes)>limitnum){
      genes=genes[order(x[genes,'logFC'],decreasing=FALSE)[1:limitnum]]
    }
    clusters<-rep(names(wout_H)[i],length(genes))
    return(cbind(genes,clusters))
  })
  wout_TOP_H<-do.call(rbind,wout_TOP_H_l)
  #wout_TOP_H<-wout_TOP_H[!duplicated(wout_TOP_H[,1]),]
  some_genes=wout_TOP_H[,1]
  #some_genes[duplicated(some_genes)]<-paste(some_genes[duplicated(some_genes)],'.duplicated',sep='')
  
  names(some_genes)=wout_TOP_H[,2]
  outlist<-lapply(wout_TOP_H_l,function(x){return(x[,1])})
  names(outlist)<-names(wout_H)
  
  if(returnlist){
    return(outlist)
  }else{
    return(some_genes)
  }

}



#noisy genes#########


noiseBaseFit_modi<-function (x, step = 0.01, thr = 0.05)
{
  m <- apply(x, 1, mean)
  v <- apply(x, 1, var)
  f <- m > 0
  lm <- log2(m)[f]
  lv <- log2(v)[f]
  lmN <- c()
  lvN <- c()
  for (i in 1:round(1/step, 0)) {
    f <- lm > quantile(lm, (i - 1) * step) & lm <= quantile(lm,
                                                            i * step)
    vi <- lv[f]
    mi <- lm[f]
    q <- quantile(vi, thr)
    f <- vi < q
    lmN <- append(lmN, mi[f])
    lvN <- append(lvN, vi[f])
  }
  nfit <- lm(lvN ~ lmN + I(lmN^2))
  list(nfit = nfit, lm_used = lmN, lv_used = lvN)
}

noise_raceid<-function(input_dat){
  muu<-rowMeans(input_dat)
  sdd<-apply(input_dat,1,sd)
  oot<-noiseBaseFit_modi(input_dat,step=0.01,thr=0.05)
  lm_out<-oot$nfit
  myVar<-2^(coef(lm_out)[1]+coef(lm_out)[2]*log2(muu)+coef(lm_out)[3]*(log2(muu))^2)
  tt2<-(input_dat-muu)^2
  vij<-rowMeans(tt2)/(myVar)
  vijlog2<-rowMeans(tt2)/log2(myVar)
  #return(list(vij=vij,vijlog2=vijlog2,noise_fit=oot))
  return(list(vij=vij))
}

#object=CNV
noise_seurat<-function(object,...){
  #FindVariableFeatures in seurat is applied on the input count data, if not exist then applied on the data in data (normalized) slot.
  vij_dat<-FindVariableFeatures(object, mean.function = ExpMean, dispersion.function = LogVMR,... )
  vij<-vij_dat$variance.standardized
  names(vij)<-rownames(object)
  return(list(vij_dat=vij_dat,vij=vij))
}






#####dropout functions###########
dropout_fun<-function(data){
  xx<-rowMeans(data)
  yy<-apply(data,1,function(x){length(which(x==0))/length(x)})
  message("
plot(drout$mean,drout$dropouts,log='x',pch=16,xlab='Mean',ylab='Dropout rates')
lines(sort(drout$mean),exp(-sort(drout$mean)),col=3)")
  
  return(list(mean=xx,dropouts=yy,nCells=(1-yy)*dim(data)[2]))
}

## list available HCL color palettes########
# hcl.pals("qualitative")
# hcl.pals("sequential")
# hcl.pals("diverging")
# hcl.pals("divergingx")
#pie(rep(1, 12), col = hcl.colors(12, "Tropic"), main = "HCL")

library(RColorBrewer)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
# pie(rep(1, 12), col = SpatialColors(12), main = "HCL")



#bayNorm relevant########
library(Matrix)


bayLL_tempfun<-function(rawdata,BETA,labels,labels_order){
  bout_list<-list()
  for(i in 1:length(labels_order)){
    print(i)
    usedind<-names(labels)[which(labels==labels_order[i])]
    bout<-bayNorm(Data=rawdata[,usedind],BETA_vec = BETA[usedind],mean_version = TRUE)
    bout_list[[i]]<-bout$Bay_out
  }
  names(bout_list)<-labels_order
  return(bout_list)
}


Beta_default<-function(Data,mean_beta=0.06){
  xx<-Matrix::colSums(Data)
  BETA_vec <- xx/mean(xx)*mean_beta
  BETA_vec[BETA_vec<=0]<-min(BETA_vec[BETA_vec>0])
  BETA_vec[BETA_vec>=1]<-max(BETA_vec[BETA_vec<1])
  return(BETA_vec)
}


######map vector of numbers to colors###########
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


###utils from seurat ###########
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}


#Inspired by a need for a reciprocal %||%, %iff% is a binary operator that returns the value of the right-hand side if the left-hand side is not NULL
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}


#frequent used colors##########
# message(
#   cyan$bold(
#     '
# Run plot_grid(qq2) to plot \n
# DIY color:                  \n
#   Discrete: scale_color_manual(values=SpatialColors(length(unique(cell_identity)))) \n
#   Continuous: scale_color_continuous(type = "viridis")
#   '))



#trajectory slingshot#########

Find_trajectory<-function(LowerDim,ClusteringLabels,start.clus = '0',...){

  qq=slingshot(LowerDim, ClusteringLabels, start.clus = start.clus,...)
  sce <- qq@curves
  
  traje_list<-list()
  for(i in 1:length(sce)){
    tempp<-as.data.frame(sce[[i]]$s)
    tempp<-tempp[sce[[i]]$ord,]
    colnames(tempp)<-c('X','Y')
    tempp['curve']<-rep(i,dim(tempp)[1])
    traje_list[[i]]<-tempp
  }
  traje_dat<-do.call(rbind,traje_list)
  traje_dat$curve<-as.factor(traje_dat$curve)
  return(traje_dat)
}
# plot(LowerDim,col=ClusteringLabels,pch=16)
# lines(qq)
# 
# lines(sce$curve2)
# 
# 



######Supervised learning######
##Random Forest#######
library(randomForest)
library("foreach")
library("doSNOW")

RF_selfun<-function(
  X,Y,
  ncores=4,
  ntree = rep(250, 4),
  mtry=100,
  importance=TRUE,
  show_info=FALSE,
  seed=NULL){
  
  
  if(show_info){
    message(
      cyan$bold(     '
#X: Rows: samples; columns: features
#Y: Labels
#top100_fea<-rownames(rf_imp)[1:100]
#yhat<-predict(rf_out,newdata=t(test))
  '))
  }

  cl<-makeCluster(ncores, type="SOCK")
  registerDoSNOW(cl)
  
  if(!is.null(seed)){
    system.time(rf_out <- foreach(ntree = ntree, .combine = randomForest::combine, .packages = "randomForest") %dorng%  randomForest(x=X, y=Y, ntree = ntree,mtry=mtry,importance=importance))
  }else{
    system.time(rf_out <- foreach(ntree = ntree, .combine = randomForest::combine, .packages = "randomForest") %dopar%  randomForest(x=X, y=Y, ntree = ntree,mtry=mtry,importance=importance))
  }
  
  stopCluster(cl)
  
  #yhat<-predict(rf_out,newdata=t(test))
  
  rf_imp<-rf_out$importance[order(rf_out$importance[,1],decreasing=T),]
  #top100_fea<-rownames(rf_imp)[1:100]
  return(list(rf_out=rf_out,rf_imp=rf_imp))
}



####load velocity data under windows environment (spliced unspliced, taken from velocity.R and Seuratwrapper package)#####


read.loom.matrices_wt <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}




ReadVelocity_wt<-function (file, engine = "hdf5r", verbose = TRUE) 
{
  
  if (verbose) {
    sink(file = stderr(), type = "output")
    on.exit(expr = sink())
    ldat <- read.loom.matrices_wt(file = file, engine = engine)
  }
  else {
    invisible(x = capture.output(ldat <- read.loom.matrices_wt(file = file,  engine = engine)))
  }
  return(ldat)
}

load_STfun_v1<-function(
    matrix_dir="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/filtered_feature_bc_matrix/",
    path_position="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_positions_list.csv",
    path_json="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/scalefactors_json.json",
    image_path_high="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_hires_image.png",
    image_path_low="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_lowres_image.png",
    quickcheck=FALSE
){
  #define for loading Tirosh's ST
  #spaceranger 1.3.1
  tissue_positions_list <- read.csv(path_position, header=FALSE)
  colnames(tissue_positions_list)<-c(
    'barcode',
    'in_tissue',
    'array_row',
    'array_col',
    'pxl_col_in_fullres',
    'pxl_row_in_fullres'
  )
  rownames(tissue_positions_list)<-tissue_positions_list$barcode
  
  
  
  scale_coor<- fromJSON(file = path_json)
  
  #scale to low resolution version's image
  tissue_positions_list['pxl_col_low']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_lowres_scalef
  tissue_positions_list['pxl_row_low']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_lowres_scalef
  
  #scale to high resolution version's image
  tissue_positions_list['pxl_col_high']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_hires_scalef
  tissue_positions_list['pxl_row_high']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_hires_scalef
  
  # barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  # features.path <- paste0(matrix_dir, "features.tsv.gz")
  # matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  # mydata <- (readMM(file = matrix.path))
  # feature.names = read.delim(features.path, 
  #                            header = FALSE,
  #                            stringsAsFactors = FALSE)
  # barcode.names = read.delim(barcode.path, 
  #                            header = FALSE,
  #                            stringsAsFactors = FALSE)
  # colnames(mydata) = barcode.names$V1
  # rownames(mydata) = feature.names$V2
  
  mydata<-Read10X_h5(paste(Data_dir,'filtered_feature_bc_matrix.h5',sep=''))
  
  
  
  if(quickcheck){
    mydata_sub<-mydata[,rownames(tissue_positions_list)[which(tissue_positions_list$in_tissue==1 )]]
    
    plot_dat<-data.frame(
      X=tissue_positions_list[colnames(mydata_sub),'pxl_row_low'],
      Y=tissue_positions_list[colnames(mydata_sub),'pxl_col_low'],
      Geneexp=mydata_sub[which.max(rowMeans(mydata_sub)),]
    )
    
    
    gout<-spatial_plot_self2(plot_dat,image_path=image_path_low,
                             x='X',
                             y='Y',
                             expres_choice='Geneexp',
                             point_size=2,alphaa=c(1,1))
    
    print(    gout$with_dots+
                labs(x='',y='',color='Geneexp')+
                scale_color_continuous(type = "viridis")+
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.text = element_blank(),
                      axis.ticks = element_blank()))
    #print(spatial_plot_wt(pltdat=t(plot_dat), igene='Geneexp', coordinates = plot_dat[,c(1,2)]) )
    
  }
  
  
  
  return(list(
    mydata=mydata,
    tissue_positions_list=tissue_positions_list,
    scale_coor=scale_coor,
    image_path_high=image_path_high,
    image_path_low=image_path_low
  ))
  
}


#load ST data fun (new)#########
load_STfun_v2<-function(
  matrix_dir="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/filtered_feature_bc_matrix/",
  path_position="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_positions_list.csv",
  path_json="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/scalefactors_json.json",
  image_path_high="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_hires_image.png",
  image_path_low="D:/Spatial_scRNAseq/BTC_data/Data_GRCh38/NH18-1181_2/spatial/tissue_lowres_image.png",
  quickcheck=FALSE
){
  #spaceranger v1.3
  # tissue_positions_list <- read.csv(path_position, header=FALSE)
  # colnames(tissue_positions_list)<-c(
  #   'barcode',
  #   'in_tissue',
  #   'array_row',
  #   'array_col',
  #   'pxl_col_in_fullres',
  #   'pxl_row_in_fullres'
  # )
  
  #now spaceranger v2.0
  tissue_positions_list <- read.csv(path_position, header=TRUE)
  rownames(tissue_positions_list)<-tissue_positions_list$barcode
  
  

  
  scale_coor<- fromJSON(file = path_json)
  
  #scale to low resolution version's image
  tissue_positions_list['pxl_col_low']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_lowres_scalef
  tissue_positions_list['pxl_row_low']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_lowres_scalef
  
  #scale to high resolution version's image
  tissue_positions_list['pxl_col_high']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_hires_scalef
  tissue_positions_list['pxl_row_high']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_hires_scalef
  
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mydata <- (readMM(file = matrix.path))
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mydata) = barcode.names$V1
  rownames(mydata) = feature.names$V2
  
  
  
  
  if(quickcheck){
    mydata_sub<-mydata[,rownames(tissue_positions_list)[which(tissue_positions_list$in_tissue==1 )]]
    
    plot_dat<-data.frame(
      X=tissue_positions_list[colnames(mydata_sub),'pxl_col_low'],
      Y=tissue_positions_list[colnames(mydata_sub),'pxl_row_low'],
      Geneexp=mydata_sub[which.max(rowMeans(mydata_sub)),]
    )
    
    
    gout<-spatial_plot_self2(plot_dat,image_path=image_path_low,
                             x='X',
                             y='Y',
                             expres_choice='Geneexp',
                             point_size=2,alphaa=c(1,1))

    print(    gout$with_dots+
                labs(x='',y='',color='Geneexp')+
                scale_color_continuous(type = "viridis")+
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.text = element_blank(),
                      axis.ticks = element_blank()))
    #print(spatial_plot_wt(pltdat=t(plot_dat), igene='Geneexp', coordinates = plot_dat[,c(1,2)]) )
    
  }
  
  
  
  return(list(
    mydata=mydata,
    tissue_positions_list=tissue_positions_list,
    scale_coor=scale_coor,
    image_path_high=image_path_high,
    image_path_low=image_path_low
  ))
  
}

load_STfun<-function(
  path="D:/Spatial_scRNAseq/spatialLIBD/Data/151673/",
  filename_h5='151673_filtered_feature_bc_matrix.h5',
  filename_tissue_positions_list="tissue_positions_list.txt",
  quickcheck=FALSE
){
  
  filee<-paste(path,'/',filename_h5,sep='')
  tissue_positions_list <- read.csv(paste(path,'/',filename_tissue_positions_list,sep=''), header=FALSE)
  rownames(tissue_positions_list)<-tissue_positions_list$V1
  colnames(tissue_positions_list)<-c(
    'barcode',
    'in_tissue',
    'array_row',
    'array_col',
    'pxl_col_in_fullres',
    'pxl_row_in_fullres'
  )
  #tissue_positions_list<-tissue_positions_list[rownames(tissue_positions_list)[which(tissue_positions_list$in_tissue==1 )],]
  scale_coor<- fromJSON(file = paste(path,"scalefactors_json.json",sep=''))
  
  #scale to low resolution version's image
  tissue_positions_list['pxl_col_low']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_lowres_scalef
  tissue_positions_list['pxl_row_low']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_lowres_scalef
  
  #scale to high resolution version's image
  tissue_positions_list['pxl_col_high']<-tissue_positions_list$pxl_col_in_fullres*scale_coor$tissue_hires_scalef
  tissue_positions_list['pxl_row_high']<-tissue_positions_list$pxl_row_in_fullres*scale_coor$tissue_hires_scalef
  
  mydata<-Read10X_h5(filee)
  
  
  image_path_high<-paste(path,"tissue_hires_image.png",sep='')
  image_path_low<-paste(path,"tissue_lowres_image.png",sep='')
  
  
  if(quickcheck){
    
    
    
    mydata_sub<-mydata[,rownames(tissue_positions_list)[which(tissue_positions_list$in_tissue==1 )]]
    
    plot_dat<-data.frame(
      X=tissue_positions_list[colnames(mydata_sub),'pxl_row_low'],
      Y=tissue_positions_list[colnames(mydata_sub),'pxl_col_low'],
      Geneexp=mydata_sub[which.max(rowMeans(mydata_sub)),]
    )
    
    
    gout<-spatial_plot_self2(plot_dat,image_path=image_path_low,
                             x='X',
                             y='Y',
                             expres_choice='Geneexp',
                             point_size=2,alphaa=c(1,1))
    
print(    gout$with_dots+
            labs(x='',y='',color='Geneexp')+
            scale_color_continuous(type = "viridis")+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text = element_blank(),
                  axis.ticks = element_blank()))
#print(spatial_plot_wt(pltdat=t(plot_dat), igene='Geneexp', coordinates = plot_dat[,c(1,2)]) )
    
  }
  

  
  return(list(
    mydata=mydata,
    tissue_positions_list=tissue_positions_list,
    scale_coor=scale_coor,
    image_path_high=image_path_high,
    image_path_low=image_path_low
  ))
  
}


##Adapted from SpaGCN paper########

#find neighbour cells for a specific cluster
st_find_neighbours<-function(cell_labels,target_cluster,coordinates,neighbours=50,plot_check=FALSE){
  target_cells<-names(cell_labels)[which(cell_labels==target_cluster)]
  distt<-as.matrix(dist(coordinates))
  distt_sub<-distt[target_cells,-which(colnames(distt) %in% target_cells)]
  
  for(x in 1:dim(distt_sub)[1]){
    sdf<-order(distt_sub[x,],decreasing=FALSE)[1:neighbours]
    distt_sub[x,sdf]=1
    distt_sub[x,-sdf]=0
  }
  neighbours_cells<-as.vector(apply(distt_sub,1,function(x){colnames(distt_sub)[which(x==1)]}))
  
  neighbours_cells<-neighbours_cells[!duplicated(neighbours_cells)]
  
  if(plot_check){
    plot(coordinates[,c(1,2)])
    see<-names(cell_labels)[which(cell_labels==target_cluster)]
    points(coordinates[see,c(1,2)],col=3,pch=16)
    points(coordinates[neighbours_cells,c(1,2)],col=4,pch=16)
    legend('topright',legend=c('target','neighbours'),col=c(3,4),pch=16,bty='n')
  }
  return(list(
    target_cells=target_cells,
    neighbours_cells=neighbours_cells
  ))
  
}

#run wilcox for certain genes between target and neighbours
wilcox_dat<-function(genesused,target_dat,neighbour_dat,  prop_target.threshold,
                     prop_target_neighbour,
                     p.threshold,
                     verbose=TRUE,
                     using_method='BioQC',
                     NCores=4){
if(length(genesused)>1){
  wout<-WilcoxDETest_wt(
    data.use=cbind(target_dat,neighbour_dat)[genesused,],
    cells.1=colnames(target_dat),
    cells.2=colnames(neighbour_dat),
    verbose = verbose,
    using_method=using_method,
    NCores=NCores
  )
  pval_vec<-wout$p_val
  pval_vec2<-p.adjust(pval_vec,method='BH')
  names(pval_vec)<-genesused
  names(pval_vec2)<-genesused
  

  
  logFC<-rowMeans(log2(neighbour_dat[genesused,]+1))-rowMeans(log2(target_dat[genesused,]+1))
  logFC<-logFC[intersect(names(logFC),names(which(pval_vec2<p.threshold)))]
  pval_vec2<-pval_vec2[pval_vec2<p.threshold]
  target_markers<-sort(logFC[which(logFC<mean(logFC))],decreasing=FALSE)
  
} else{
  wout<-wilcox.test(target_dat[genesused,],neighbour_dat[genesused,])
  pval_vec<-wout$p.value
  pval_vec2<-pval_vec
  names(pval_vec)<-genesused
  names(pval_vec2)<-genesused
  
  logFC<-mean(log2(neighbour_dat[genesused,]+1))-mean(log2(target_dat[genesused,]+1))
  target_markers<-genesused
  
}
 

  
  return(list(
    target_markers=target_markers,
    logFC=logFC,
    pval_ori=pval_vec,
    pval_adj=pval_vec2
  ))
}


#wilcox test between target cluster and its neighbours
st_find_markers<-function(
  Data,
  target_cells,
  neighbours_cells,
  prop_target.threshold=0.8,
  prop_target_neighbour=1,
  p.threshold=0.5,
  using_method='BioQC',
  NCores=4
){
  
  Data<-Check_input(Data)
  
  target_dat<-Data[,target_cells]
  neighbour_dat<-Data[,neighbours_cells]
  #proportion of non zero count
  prop_target<-diff(Matrix::t(target_dat)@p)/dim(target_dat)[2]
  prop_neighbour<-diff(Matrix::t(neighbour_dat)@p)/dim(neighbour_dat)[2]
  g_1<-rownames(target_dat)[which(prop_target/prop_neighbour>=prop_target_neighbour)]
  g_2<-rownames(target_dat)[which(prop_target>=prop_target.threshold)]
  
  g_1_down<-rownames(target_dat)[which((1-prop_target)/(1-prop_neighbour)>=prop_target_neighbour)]
  g_2_down<-rownames(target_dat)[which((1-prop_target)>=prop_target.threshold)]
  
  
  
  genesused<-intersect(g_1,g_2)
  genesused_down<-intersect(g_1_down,g_2_down)
  if(length(genesused)==0){
    out_up<-NULL
  }else{
    out_up<-wilcox_dat(genesused=genesused,target_dat=target_dat,neighbour_dat=neighbour_dat,  prop_target.threshold=prop_target.threshold,prop_target_neighbour=prop_target_neighbour,p.threshold=p.threshold,  using_method=using_method, NCores=NCores)
  }
  
  if(length(genesused_down)==0){
    out_down<-NULL
  }else{
    out_down<-wilcox_dat(genesused=genesused_down,target_dat=target_dat,neighbour_dat=neighbour_dat,prop_target.threshold=prop_target.threshold,prop_target_neighbour=prop_target_neighbour,p.threshold=p.threshold, using_method=using_method, NCores=NCores)
  }

  return(list(
    out_up=out_up,
    out_down=out_down
  ))
}






st_find_markers_full<-function(
  Data,
  cell_labels,
  coordinates,
  neighbours=30,
  plot_check=FALSE,
  prop_target.threshold=0.8,
  prop_target_neighbour=1,
  p.threshold=0.5,
  using_method='BioQC',
  NCores=4
){
  Data<-Check_input(Data)
  
  
  unii<-unique(cell_labels)
  markers_list<-list()
  
  for(cluind in seq(1,length(unii))){
    print(paste('Find markers in cluster: ',cluind,sep=''))
    
    target_cluster=unii[cluind]
    
    stneiout<-st_find_neighbours(cell_labels=cell_labels,target_cluster=target_cluster,coordinates=coordinates,neighbours=neighbours,plot_check=plot_check)
    
    
    markers<-st_find_markers(
      Data=Data,
      target_cells=stneiout$target_cells,
      neighbours_cells=stneiout$neighbours_cells,
      prop_target.threshold=prop_target.threshold,
      prop_target_neighbour=prop_target_neighbour,
      p.threshold=p.threshold,
      using_method=using_method, 
      NCores=NCores
    )
    markers_list[[cluind]]<-markers
    
  }
  
  names(markers_list)<-unique(cell_labels)
  
  return(markers_list)
}





#load ST data fun (old)#########
load_sp<-function(path_tsv){
  
  dat_st<-(fread(path_tsv))
  vecc<-dat_st[,1]
  dat_st<-dat_st[,-1]
  dat_st<-t(dat_st)
  colnames(dat_st)<-as.character(vecc$V1)
  
  dat_coor<-do.call(rbind,strsplit(vecc$V1,'x'))
  dat_coor2<-cbind(as.numeric((dat_coor[,1])),as.numeric((dat_coor[,2])))
  rownames(dat_coor2)<-vecc$V1
  
  # vecc$V1[1:10]
  # dat_coor2[1:10,]

  return(list(
    Data=dat_st,
    Corrdinates=dat_coor2
  ))
}




getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


match_annotations<-function(x_target,x){
  #make sure:
  #all.equal(names(x_target),names(x))
  
  library(riverplot)
  message(
    cyan$bold(     '
#plot(rp, plot_area = 0.95, yscale=0.06)'))
  
  
  levels(x_target)<-paste('left:',levels(x_target),sep='')
  levels(x)<-paste('right:',levels(x),sep='')
  
  # ID1<-paste('target:',levels(x_target),sep='')
  # ID2<-paste('compare:',levels(x),sep='')
  ID1<-unique(x_target)
  ID2<-unique(x)
  
  nodes=data.frame(ID=c(as.character(ID1),as.character(ID2)),
                   x=c(rep(1,length(unique(x_target))),rep(2,length(unique(x)))),
                   y=c(seq(1,length(unique(x_target))),seq(1,length(unique(x)))))

  
  edges<-expand.grid(ID1,ID2)
  colnames(edges)<-c('N1','N2')
  edges['Value']<-rep(0,dim(edges)[1])
  
  
  
  for(i in 1:dim(edges)[1]){
    
    edges[i,3]<-length(intersect(
      names(x[x==edges[i,2]]),
      names(x_target[x_target==edges[i,1]])
    ))
  
  }
  
  
  
  palette = SpatialColors(max(length(ID1),length(ID2)))
  styles = lapply(nodes$y, function(n) {
    list(col = palette[n+1], lty = 0, textcol = "black")
  })
  names(styles) = nodes$ID
  
  rp <- list(nodes =nodes, edges = edges, styles = styles)
  
  class(rp) <- c(class(rp), "riverplot")
  #plot(rp, plot_area = 0.95, yscale=0.06)
  
  
  return(list(nodes=nodes,edges=edges,rp=rp))
  
}


mse_fun<-function(x,y){
  return(mean(((x-y)^2),na.rm =TRUE))
}





##spatialDE##########
spatialDE_wt<-function(Data,spatial_locs){
  library(Giotto)
  library(reticulate)
  use_condaenv("giotto_env")
  checkGiottoEnvironment()
  
  
  gobject<-createGiottoObject(
    raw_exprs=Data,
    spatial_locs = spatial_locs,
    norm_expr = NULL,
    norm_scaled_expr = NULL,
    custom_expr = NULL,
    cell_metadata = NULL,
    gene_metadata = NULL,
    spatial_network = NULL,
    spatial_network_name = NULL,
    spatial_grid = NULL,
    spatial_grid_name = NULL,
    spatial_enrichment = NULL,
    spatial_enrichment_name = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    images = NULL,
    offset_file = NULL,
    instructions = NULL,
    cores = NA
  )
  
  
  spaout<-spatialDE(
    gobject = gobject,
    expression_values = "raw",
    size = c(4, 2, 1),
    color = c("blue", "green", "red"),
    sig_alpha = 0.5,
    unsig_alpha = 0.5,
    python_path = NULL,
    show_plot = NA,
    return_plot = NA,
    save_plot = NA,
    save_param = list(),
    default_save_name = "SpatialDE"
  )
  
  

  spatialDE_summarytable<-spaout$results$results
  rownames(spatialDE_summarytable)<-spatialDE_summarytable$g

  spatialDE_summarytable<-spatialDE_summarytable[order(spatialDE_summarytable$BIC,decreasing=TRUE),]
  
  
  
  return(list(
    spaout=spaout,
    spatialDE_summarytable=spatialDE_summarytable
  ))
}


###BayesSpace########

stpreprocess_wt<-function(Data,spatial_locs,n.HVGs=2000,n.PCs=15,log.normalize=TRUE,assay.type='logcounts'){
  
  sce <- SingleCellExperiment(assays=list(counts=Data))
  if(!is.null(spatial_locs)){
    sce@colData$row<-spatial_locs[,1]
    sce@colData$col<-spatial_locs[,2]
  }

  if (log.normalize)
    sce <- logNormCounts(sce)
  
  dec <- modelGeneVar(sce, assay.type=assay.type)
  top <- getTopHVGs(dec, n=n.HVGs)
  sce <- runPCA(sce, subset_row=top, ncomponents=n.PCs, exprs_values=assay.type)
  rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
  return(sce)
}

BayesSpace_tune_wt<-function(
  Data,
  spatial_locs,
  num_clusters=7,
  n.PCs=15,
  n.HVGs=1000,
  init=NULL,
  log.normalize=TRUE,
  assay.type = "logcounts",
  rangee=seq(2,10)
){
  library(BayesSpace)
  
  sce <- SingleCellExperiment(assays=list(counts=Data))
  sce@colData$row<-spatial_locs[,1]
  sce@colData$col<-spatial_locs[,2]
  set.seed(102)
  sce <- spatialPreprocess(
    sce, 
    platform="Visium",
    n.PCs=n.PCs, 
    n.HVGs=n.HVGs, 
    log.normalize=log.normalize,
    assay.type = assay.type)
  
  sce <- qTune(sce, qs=rangee, platform="Visium", d=15)
  # qPlot(sce)

  
  message(cyan$bold('
qPlot(sce)'))
  return(sce)
}

#allow non gene count as input
spatialPreprocess_wt<-function (sce, platform = c("Visium", "ST"), n.PCs = 15, n.HVGs = 2000, skip.PCA = FALSE, log.normalize = TRUE, assay.type = "logcounts", BSPARAM = ExactParam()) 
{
  library(BiocSingular)
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- match.arg(platform)
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  if (!skip.PCA) {
    if (log.normalize) 
      sce <- logNormCounts(sce)
    
    if(!is.null(n.HVGs)){
      dec <- modelGeneVar(sce, assay.type = assay.type)
      top <- getTopHVGs(dec, n = n.HVGs)
      sce <- runPCA(sce, subset_row = top, ncomponents = n.PCs, 
                    exprs_values = assay.type, BSPARAM = BSPARAM)
      rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
    }else{
      sce <- runPCA(sce, ncomponents = n.PCs, 
                    exprs_values = assay.type, BSPARAM = BSPARAM)
    }
    
  }
  sce
}
BayesSpace_wt<-function(
  Data,
  spatial_locs,
  num_clusters=7,
  n.PCs=15,
  n.HVGs=1000,#set to NULL if non gene count
  init=NULL,
  log.normalize=TRUE,
  assay.type = "logcounts"
  ){
  library(BayesSpace)
  
  sce <- SingleCellExperiment(assays=list(counts=Data))
  sce@colData$row<-spatial_locs[,1]
  sce@colData$col<-spatial_locs[,2]
  set.seed(102)
  sce <- spatialPreprocess_wt(
    sce, 
    platform="Visium",
    n.PCs=n.PCs, 
    n.HVGs=n.HVGs, 
    log.normalize=log.normalize,
    assay.type = assay.type)
  
  # sce <- qTune(sce, qs=seq(2, 10), platform="Visium", d=15)
  # qPlot(sce)
  set.seed(12300)
  sce <- spatialCluster(
    sce, 
    q=num_clusters, 
    platform="Visium", 
    d=15,
    init = init,
    init.method="mclust", 
    model="normal", 
    gamma=3,
    nrep=1000, 
    burn.in=100,
    save.chain=FALSE)
  names(sce$spatial.cluster)<-colnames(Data)
  
  message(cyan$bold('
cluster_idents<-sce$spatial.cluster'))
  return(sce)
}



platform_neighbors_wt <- function(coordinates, platform,index0=FALSE) {
  if (platform == "Visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (platform == "ST") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  # colnames(stRNA_list$tissue_positions_list)
  # coordinates<-stRNA_list$tissue_positions_list[,c('array_col','array_row')]
  # colnames(coordinates)<-c('col','row')
  # coordinates<-coordinates[colnames(mydata_sub_new),]
  colnames(coordinates)<-c('col','row')
  spot.positions <- coordinates
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  if(index0==TRUE){
    neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  }

  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  
  df_j <- purrr::map(df_j, function(nbrs) purrr::discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(purrr::keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          nrow(coordinates), " spots.")
  names(df_j)<-rownames(coordinates)
  return(df_j)
  
  #debug
  # names(df_j)<-rownames(coordinates)
  # indd<-30
  # celll<-colnames(mydata_sub_new)[indd]
  # plot(coordinates[colnames(mydata_sub_new),],col='lightgrey')
  # points(coordinates[celll,],col=3,pch=16)
  # points(coordinates[df_j[celll][[1]]+1,],col=4)
}


##infCNV##########
infCNV_wt<-function(
    Data,
    annotations_file,
    gene_order_file,
    ref_group_names,
    cluster_by_groups=FALSE,
    cluster_references=FALSE,
    BayesMaxPNormal=0,
    denoise=TRUE,
    HMM=FALSE,
    cutoff=0.1,
    chr_exclude=c('chrX', 'chrY', 'chrM'),
    num_threads=8,
    no_plot=TRUE,
    no_prelim_plot=TRUE,
    plot_probabilities=FALSE,
    inspect_subclusters=FALSE,
    scale_data=FALSE,
    out_dir=tempfile()){
  
  
  infercnv_obj = infercnv::CreateInfercnvObject(
    raw_counts_matrix=Data,
    annotations_file=annotations_file,
    gene_order_file=gene_order_file,
    ref_group_names=ref_group_names,
    min_max_counts_per_cell=c(1,+Inf),
    chr_exclude = chr_exclude)
  
  
  infercnv_obj= infercnv::run(
    infercnv_obj,
    cutoff=cutoff,
    out_dir=out_dir,
    cluster_by_groups=cluster_by_groups,
    cluster_references=cluster_references,
    BayesMaxPNormal=BayesMaxPNormal,
    denoise=denoise,
    HMM=HMM,resume_mode=FALSE,
    num_threads=num_threads,
    no_plot=no_plot,
    no_prelim_plot=no_prelim_plot,
    plot_probabilities=plot_probabilities,
    inspect_subclusters=inspect_subclusters,
    scale_data = scale_data
  )
  
  return(infercnv_obj)
}

#cellchat#######
cellchat_wt<-function(object,meta,group.by='type',species='H'){
  cellchat <- createCellChat(object = object, meta = meta, group.by = group.by)
  if(species=='H'){
    CellChatDB <- CellChatDB.human
    PPI_used<-PPI.human
  }else{
    CellChatDB <- CellChatDB.mouse
    PPI_used<-PPI.mouse
  }
  
  
  #showDatabaseCategory(CellChatDB)
  #CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) 
  #future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.p=0.1)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI_used)
  #inference
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat,thresh=0.1)
  cellchat <- aggregateNet(cellchat,thresh=0.1)
  #groupSize <- as.numeric(table(cellchat@idents))
  return(cellchat)
}








###check input#########
as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}






#' @title Check input 
#'
#' @description  Check input 
#'
#' @param Data Input data.
#' @return Matrix (of type matrix in R).
#'
#' @examples
#' aa<-matrix(seq(1,6),nrow=2,ncol=3)
#' Check_input(aa)
#'
#' @export
Check_input <- function(Data){
  
  # Adapted from SCnorm
  if (methods::is(Data, "SummarizedExperiment") | methods::is(Data, "SingleCellExperiment")) {
    if (is.null(SummarizedExperiment::assayNames(Data))
        || SummarizedExperiment::assayNames(Data)[1] !=
        "Counts") {
      message("Renaming the first element in assays(Data) to 'Counts'")
      SummarizedExperiment::assayNames(Data)[1] <- "Counts"
      if (
        is.null(
          colnames(SummarizedExperiment::assays(Data)[["Counts"]]))) {
        stop("Must supply sample/cell names!")
      }
    }
    Data <- SummarizedExperiment::assays(Data)[["Counts"]]
  }
  
  
  #use dgCMatrix
  if(is(Data, 'sparseMatrix')){
    if(!is(Data, 'dgCMatrix')){
      Data <- as(Data, "dgCMatrix")
      #Data<-as_matrix(Data)
    } 
  } else{
    Data <- as(as.matrix(Data), "dgCMatrix")
  }
  
  
  return(Data)
}





#Some other tips###########

##change order of factor x axis in ggplot2########
#plot_z$groupp<-factor(plot_z$groupp,levels = c("1","2","3","0"))


#######这个方式漂亮！！！！！！！！！！！！！
#shifted_df <- tidyr::pivot_longer(data = data, cols = dplyr::all_of(states),  names_to = "gene_set", values_to = "gene_set_expr")




firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}



##gene signature from Neftel 2019 (scalop)#####
#counts_vec should have gene as names
SC_fun<-function(gene_set,counts_vec,nbins=30,rn=100){
  binss<-OneR::bin(data=log(counts_vec+1), nbins = nbins)
  names(binss)<-names(counts_vec)
  ctr_genes<-NULL

  for(i in 1:length(gene_set)){
    poolt<-names(binss)[which(binss==binss[gene_set[i]])]
    sii<-ifelse(length(poolt)>rn,rn,length(poolt))
    sg<-sample(poolt,sii)
    ctr_genes<-c(ctr_genes,sg)
  }
  SC<-mean(counts_vec[gene_set])-mean(counts_vec[ctr_genes])
  return(SC)
}

SC_data_fun<-function(inputdat,gene_set,nbins=30,rn=100,NCores=4){

  cluster <- makeCluster(NCores, type = "SOCK")
  registerDoSNOW(cluster)
  getDoParWorkers()
  
  iterations <- dim(inputdat)[2]
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  SC_out<-foreach(i=1:dim(inputdat)[2],.combine=c,.options.snow = opts,.export = c('SC_fun'))%dopar%{
    scall<-SC_fun(gene_set=gene_set,counts_vec=inputdat[,i],nbins=nbins,rn=rn)
    return(scall)
  }
  
  close(pb)
  stopCluster(cluster)

  
  return(SC_out)
}



library(SpatialPack)
modified.ttest_fun<-function(x,y,coords){
  mtout=modified.ttest(x=x, y=y, coords=coords, nclass = 20)
  return(c(mtout$p.value))
}
#sub_input_1 (H, 2 for mouse): rows scores, columns spots
corPval_fun<-function(sub_input_1,sub_input_2, coords, NCores=6,slide_ind=1){
  
  cluster <- makeCluster(NCores, type = "SOCK")
  registerDoSNOW(cluster)
  getDoParWorkers()
  
  iterations <- dim(sub_input_1)[1]
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  BB_parmat <- foreach(
    Geneind = seq_len(iterations),
    .combine = rbind,
    .options.snow = opts,
    .packages = c('SpatialPack')) %dopar% {
      
      tempdat<-NULL
      for(jj in 1:dim(sub_input_2)[1]){
        
        mtout=modified.ttest(x=sub_input_1[Geneind,], y=sub_input_2[jj,], coords=coords, nclass = 20)
        tempdat<-rbind(tempdat,c(mtout$corr,mtout$p.value))
      }
      return( tempdat)
    }
  
  close(pb)
  stopCluster(cluster)
  
  cornames<-NULL
  for(i in 1:dim(sub_input_1)[1]){
    for(j in 1:dim(sub_input_2)[1]){
      cornames<-rbind(cornames,c(rownames(sub_input_1)[i],rownames(sub_input_2)[j]))
    }
    
  }
  corout<-cbind(BB_parmat,cornames)
  colnames(corout)<-c('cor','pval','H','M')
  corout<-as.data.frame(corout)
  corout$pval<-as.numeric(corout$pval)
  corout$cor<-as.numeric(corout$cor)
  corout['adjpval']<-p.adjust(corout$pval,method='BH')
  corout['pairs']<-paste(corout$H,'-',corout$M,sep='')
  corout['slides']<-slide_ind

  return(corout)
}

#Seurat WNN: weighted nearest neighbour######
WNN_wt<-function(Hdat,Mdat,resolution=c(0.5,0.7,0.9),pcadims=20){
  S_core<-CreateSeuratObject(counts=Hdat, project = "SeuratProject", assay = "H")
  adt_assay <- CreateAssayObject(counts = Mdat)
  S_core[["M"]] <- adt_assay
  
  DefaultAssay(S_core) <- 'H'
  S_core <- NormalizeData(S_core) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'Hpca')
  
  DefaultAssay(S_core) <- 'M'
  # we will use all M features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the 
  VariableFeatures(S_core) <- rownames(S_core[["M"]])
  S_core <- NormalizeData(S_core) %>% FindVariableFeatures() %>%  ScaleData() %>% RunPCA(reduction.name = 'Mpca')
  S_core <- FindMultiModalNeighbors(
    S_core, reduction.list = list("Hpca", "Mpca"), 
    dims.list = list(1:pcadims, 1:pcadims), modality.weight.name = "RNA.weight"
  )
  
  S_core <- RunUMAP(S_core, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  S_core<- FindClusters(S_core, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
  
  return(S_core)
}


#fgsea####

fgsea_wt<-function(fgsea_sets,stats,minSize=10,nproc=4,...){
  fgseaRes<- fgsea(pathways=fgsea_sets, stats =stats,minSize=minSize,nproc=nproc,...)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(dplyr::desc(NES))
  fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES, -log2err) %>%
    arrange(padj) %>%
    head()
  fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
  fgseaResTidy<-as.data.frame(fgseaResTidy)
  rownames(fgseaResTidy)<-fgseaResTidy$pathway
  return(fgseaResTidy)
}


#run on output of wilcoxauc from presto package

fgsea_wilcoxauc<-function(wout,Sample,fgsea_sets,thres_padj=0.05,n_top=10,n_bot=10,nproc=4){
  Genes<-wout %>%
    dplyr::filter(group == Sample) %>%
    arrange(dplyr::desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- deframe(Genes)
  head(ranks)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks,minSize=10,nproc=nproc)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(dplyr::desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -log2err) %>% 
    arrange(padj) %>% 
    head()
  fgseaResTidy$Enrichment = ifelse(fgseaResTidy$NES > 0, "Up-regulated", "Down-regulated")
  Filtidy<-fgseaResTidy %>% filter(padj < thres_padj) 
  filtRes = rbind(Filtidy %>% filter(Enrichment == "Up-regulated") %>% head(n= n_top),
                  Filtidy %>% filter(Enrichment == "Down-regulated") %>% head(n= n_bot))
  Y<-fgseaResTidy
  names(Y)[names(Y)=="NES"]=paste(Sample)
  fgsea_out<-Y[,c("pathway",paste(Sample))]
  return(fgsea_out)
}





compute_ustat <- function(Xr, cols, n1n2, group.size) {
  grs <- sumGroups(Xr, cols)
  
  if (is(Xr, "dgCMatrix")) {
    gnz <- (group.size - nnzeroGroups(Xr, cols))
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - group.size *
      (group.size + 1 ) / 2
  } else {
    ustat <- grs - group.size * (group.size + 1 ) / 2
  }
  return(ustat)
}

compute_pval <- function(ustat, ties, N, n1n2) {
  z <- ustat - .5 * n1n2
  z <- z - sign(z) * .5
  .x1 <- N ^ 3 - N
  .x2 <- 1 / (12 * (N^2 - N))
  rhs <- lapply(ties, function(tvals) {
    (.x1 - sum(tvals ^ 3 - tvals)) * .x2
  }) %>% unlist
  usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
  z <- t(z / usigma)
  
  pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
  return(pvals)
}
tidy_results <- function(wide_res, features, groups) {
  res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame()
  colnames(res) <- names(wide_res)
  res$feature <- rep(features, times = length(groups))
  res$group <- rep(groups, each = length(features))
  res %>% dplyr::select(
    .data$feature,
    .data$group,
    .data$avgExpr,
    .data$logFC,
    .data$logFC_avg,
    .data$statistic,
    .data$auc,
    .data$pval,
    .data$padj,
    .data$pct_in,
    .data$pct_out
  )
}



wilcoxauc_wtlfc<-function(X, y, groups_use = NULL, verbose = TRUE, ...) {
  ## Check and possibly correct input values
  if (is(X, "dgeMatrix")) X <- as.matrix(X)
  if (is(X, "data.frame")) X <- as.matrix(X)
  if (is(X, "dgTMatrix")) X <- as(X, "dgCMatrix")
  if (is(X, "TsparseMatrix")) X <- as(X, "dgCMatrix")
  if (ncol(X) != length(y)) stop("number of columns of X does not
                                match length of y")
  if (!is.null(groups_use)) {
    idx_use <- which(y %in% intersect(groups_use, y))
    y <- y[idx_use]
    X <- X[, idx_use]
  }
  
  y <- factor(y)
  idx_use <- which(!is.na(y))
  if (length(idx_use) < length(y)) {
    y <- y[idx_use]
    X <- X[, idx_use]
    if (verbose)
      message("Removing NA values from labels")
  }
  
  group.size <- as.numeric(table(y))
  if (length(group.size[group.size > 0]) < 2) {
    stop("Must have at least 2 groups defined.")
  }
  
  if (is.null(row.names(X))) {
    row.names(X) <- paste0("Feature", seq_len(nrow(X)))
  }
  
  ## Compute primary statistics
  group.size <- as.numeric(table(y))
  n1n2 <- group.size * (ncol(X) - group.size)
  if (is(X, "dgCMatrix")) {
    rank_res <- rank_matrix(Matrix::t(X))
  } else {
    rank_res <- rank_matrix(X)
  }
  
  ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
  auc <- t(ustat / n1n2)
  pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2)
  fdr <- apply(pvals, 2, function(x) p.adjust(x, "BH"))
  
  ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
  group_sums <- sumGroups(X, y, 1)
  group_nnz <- nnzeroGroups(X, y, 1)
  group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
  group_pct_out <- -group_nnz %>%
    sweep(2, colSums(group_nnz) , "+") %>% 
    sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
  group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
  cs <- colSums(group_sums)
  gs <- as.numeric(table(y))
  lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
    group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
  }))
  
  lfc_avg <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
    group_means[, g] / ((cs - group_sums[g, ]) / (length(y) - gs[g]))
  })) %>% log()
  
  
  
  res_list <- list(auc = auc,
                   pval = pvals,
                   padj = fdr,
                   pct_in = 100 * group_pct,
                   pct_out = 100 * group_pct_out,
                   avgExpr = group_means,
                   statistic = t(ustat),
                   logFC = lfc,
                   logFC_avg=lfc_avg)
  return(tidy_results(res_list, row.names(X), levels(y)))
}






wilcoxauc_wt<-function(X,y,...){
  
  #wout<-wilcoxauc(X=X,y=y,...)
  wout<-wilcoxauc_wtlfc(X=X,y=y,...)
  if(is.factor(y)){
    Y=levels(y)
  }else{
    Y=unique(y)
  }

  
  delist<-lapply(seq(1,length(Y)),function(x){
    g=Y[x]
    degenes<-wout %>%
      dplyr::filter(group == g) %>%
      dplyr::arrange(dplyr::desc(auc)) %>%
      dplyr::select(feature, auc,pval,padj,logFC,logFC_avg)
    rownames(degenes)<-degenes$feature
    
    return(degenes)
    
  })
  names(delist)<-Y
  return(delist)
}



#Moran I test#######
library(pbmcapply)
#library(spdep) #install error
library(plotrix)

ST_sc_wt<-function(s_dat,datainput,feature,radius=5){
  tissue.pos <- data.frame(
    barcodes=colnames(s_dat),
    x=s_dat$array_row,
    y=s_dat$array_col)
  
  plot.new()
  segments <- pbmcapply::pbmclapply(1:nrow(tissue.pos), function(xx){
    
    segment_c <-plotrix::draw.circle(x=tissue.pos$x[xx], 
                                     y=tissue.pos$y[xx],
                                     radius=radius*2) 
    
    segment <- as.data.frame(do.call(cbind, segment_c))
    
    segment <- sp::Polygon(cbind(segment$x, segment$y))
    segment <- sp::Polygons(list(segment), tissue.pos$barcodes[xx])
    return(segment)
  }, mc.cores = 8)
  
  
  SpP = sp::SpatialPolygons(segments, 1:length(segments))
  
  
  SrDf = sp::SpatialPolygonsDataFrame(SpP, as.data.frame(t(datainput[feature,])))
  
  
  
  nb <- spdep::poly2nb(SrDf, queen=F, snap=25)
  lw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE)
  
  message(paste0(Sys.time(), "  Computing the Moran’s I statistic"))
  
  # stat <- purrr::map(.x=feature, .f=function(x){
  #   print(x)
  #   model <- spdep::moran.test(SrDf@data[,x],lw, zero.policy=T)
  #   model$estimate[[1]]
  # })
  # 
  
  stat <- pbmcapply::pbmclapply(1:length(feature), function(xx){
    x=feature[xx]
    print(x)
    model <- spdep::moran.test(SrDf@data[,x],lw, zero.policy=T)
    model$estimate[[1]]
  }, mc.cores = 8)
  
  names(stat) <- feature
  return(stat)
}

#AUCcell#####
library(AUCell)
library(DelayedArray)

AUCell_batch <- function(inp_data, genesets, num_batches = 100) {
  ## Scores a data matrix with AUCell in batches. Idea is to limit memory consumption when
  ## scoring with AUCell
  ## INPUTS:
  ##    inp_data = input data, either a dxn matrix of d features, n samples or a Seurat object
  ##                containing such a matrix
  ##    genesets = named list of character vectors, each consisting of a set of gene symbols
  ##    num_batches = number of batches to run AUCell for. More batches = fewer cells (observations)
  ##                  for each batch used for scoring
  ##    slot = slot to use if using a Seurat object
  ##    assay = assay to use if using a Seurat object
  ## RETURNS:
  ##  either an nxp matrix (samples x scores)
  if (is.matrix(inp_data) || is(inp_data, 'dgCMatrix')) {
    num_cells <- ncol(inp_data)
    batch_size <- ceiling(num_cells/num_batches)
    score_mat <- c()
    print('Running AUCell scoring')
    Sys.time()
    for (i in 1:num_batches) {
      print(paste('batch', i, Sys.time()))
      ind1 <- (i-1)*batch_size + 1
      ind2 <- i*batch_size
      if (ind2 > num_cells) {
        ind2 <- num_cells
      }
      gene_rankings <- AUCell::AUCell_buildRankings(inp_data[,ind1:ind2], plotStats = FALSE)
      score_mat_i <- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings)
      score_mat_i <- t(assay(score_mat_i, 'AUC'))
      score_mat <- rbind(score_mat, score_mat_i)
      gc(full = TRUE, verbose = TRUE)
    }
    print('Finished Scoring')
    print(Sys.time())
    return(as.data.frame(score_mat))
  } else if (class(inp_data) == 'seurat') {
    score_mat <- AUCell_batch(inp_data = GetAssayData(inp_data, slot = 'raw.data', assay = 'RNA'), 
                              genesets = genesets, 
                              num_batches = num_batches)
    colnames(score_mat) <- paste0(colnames(score_mat), '_AUC')
    inp_data <- AddMetaData(inp_data, as.data.frame(score_mat))
  }
}




AUCell_simple <- function(inp_data, genesets) {
  gene_rankings <- AUCell::AUCell_buildRankings(inp_data, plotStats = FALSE)
  score_mat<- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings)
  score_mat<- as.data.frame(t(assay(score_mat, 'AUC')))

  return(score_mat)
  gc(full = TRUE, verbose = TRUE)
  
}



AUCell_run_wt <- function(exprMatrix, 
                       geneSets, 
                       featureType="genes", 
                       keepZeroesAsNA=FALSE, 
                       normAUC=TRUE, 
                       aucMaxRank=ceiling(0.05*nrow(exprMatrix)), 
                       BPPARAM=MulticoreParam(6)) 
{
  collected <- DelayedArray::blockApply(DelayedArray(exprMatrix), 
                          FUN=.AUCell_run_internal, 
                          geneSets=geneSets,
                          featureType=featureType,
                          keepZeroesAsNA=keepZeroesAsNA,
                          normAUC=normAUC,
                          aucMaxRank=aucMaxRank,
                          BPPARAM=BPPARAM,
                          grid=colAutoGrid(exprMatrix)
  )
  cells_AUC=do.call(cbind, collected)
  score_mat<- as.data.frame(t(assay(cells_AUC, 'AUC')))
  
  return(list(score_mat=score_mat,
              cells_AUC=cells_AUC))
}

.AUCell_run_internal <- function(exprMat_block, geneSets, featureType, keepZeroesAsNA, normAUC, aucMaxRank) 
{
  ranked <- .AUCell_buildRankings(exprMat_block, splitByBlocks=TRUE, 
                                 featureType=featureType,
                                 keepZeroesAsNA=keepZeroesAsNA, 
                                 plotStats=FALSE, verbose=FALSE)
  AUCell_calcAUC(geneSets, ranked, normAUC=normAUC, aucMaxRank=aucMaxRank, verbose=FALSE)
}


.AUCell_buildRankings <- function(exprMat, featureType="genes",
                                  splitByBlocks=FALSE,
                                  keepZeroesAsNA=FALSE, 
                                  BPPARAM=NULL,
                                  plotStats=TRUE, verbose=TRUE,
                                  nCores=NULL, mctype=NULL)
{
  if((!splitByBlocks) && ("dgCMatrix" %in% class(exprMat)))
    stop("To use a dgCMatrix as input set 'splitByBlocks=TRUE'.")
  if(!is.null(nCores))
    warning("nCores is no longer used. It will be deprecated in the next AUCell version.")
  
  ## Still needed? 
  if(keepZeroesAsNA){
    zeroesCoords <- which(exprMat==0, arr.ind=TRUE)
  }
  
  ## Calculate expression stats?
  nGenesDetected <- numeric(0)
  if(plotStats)
  {
    msg <- tryCatch(plotGeneCount(exprMat, plotStats=plotStats, verbose=verbose),
                    error = function(e) {
                      return(e)
                    })
    if(methods::is(msg,"error")) {
      warning("There has been an error in plotGeneCount() [Message: ",
              msg$message, "]. Proceeding to calculate the rankings...", sep="")
    }else{
      if(is.numeric(nGenesDetected))
        nGenesDetected <- msg
    }
  }
  
  ## Rank genes
  rowNames <- rownames(exprMat)
  colNames <- colnames(exprMat)
  exprMat <- -exprMat # ro rank decreasingly
  
  # Would prefer not to use the block inside the function... but
  # the method for sparse matrices does not allow ties.method='random'
  if(splitByBlocks){
    exprMat <- do.call(cbind,
                       blockApply(DelayedArray(exprMat),
                                  FUN=DelayedMatrixStats::colRanks,
                                  ties.method="random", preserveShape=TRUE,
                                  BPPARAM=BPPARAM,
                                  grid=colAutoGrid(exprMat))) 
  }else{
    exprMat <- DelayedMatrixStats::colRanks(exprMat, ties.method="random", preserveShape=TRUE, decreasing=TRUE) 
  }
  
  rownames(exprMat) <- rowNames
  colnames(exprMat) <- colNames
  
  if(keepZeroesAsNA){
    exprMat[which(zeroesCoords==0, arr.ind=TRUE)] <- NA
  }
  
  ## Format & return
  names(dimnames(exprMat)) <- c(featureType, "cells")
  new("aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays=list(ranking=exprMat)),
      nGenesDetected=nGenesDetected)
}




#Assign cell labels based on AUC value and back ground genes######

Assign_WT<-function(
  Data,
  gene_sets,
  nbins=30,
  rn=5,
  numboot=20,
  thres=0.99,
  boot=TRUE, 
  BPPARAM=MulticoreParam(6)
){
  #Data=bayXIAO_TUMOR_TAM$Bay_out[geneused,]
  geneused<-rownames(Data)
  
  
  message('Run signature using original data')
  temp<-AUCell_run_wt(Data,gene_sets, 
                      BPPARAM=BPPARAM)
  AUCout<-temp$score_mat
  AUCoutS<-apply(AUCout,2,scale)
  rownames(AUCoutS)<-rownames(AUCout)
  
  
  if(boot){
    
  
  
  counts_vec<-rowMeans(Data)
  set.seed(12300)
  
  
  binss<-OneR::bin(data=log(counts_vec+1), nbins = nbins)
  names(binss)<-names(counts_vec)
  
  ctr_genes_list_list<-list()
  for(ib in 1:numboot){
    ctr_genes_list<-list()
    for(gind in 1:length(gene_sets)){
      gene_set=gene_sets[[gind]]
      
      ctr_genes<-NULL
      for(i in 1:length(gene_set)){
        poolt<-names(binss)[which(binss==binss[gene_set[i]])]
        poolgenes=setdiff(poolt,gene_set[i])
        sii<-ifelse(length(poolgenes)>rn,rn,length(poolgenes))
        sg<-sample(poolgenes,sii)
        ctr_genes<-c(ctr_genes,sg)
      }
      ctr_genes_list[[gind]]<-ctr_genes
    }
    names(ctr_genes_list)<-names(gene_sets)
    ctr_genes_list_list[[ib]]<-ctr_genes_list
  }
  
  
  
  
  coo<-list()
  for(ib in 1:numboot){
    ff<-ctr_genes_list_list[[ib]]
    names(ff)<-paste('boot',ib,'_',names(ff),sep='')
    coo<-c(coo,ff)
  }
  
  #count time
  
  message('Run signature using background genes')
  start_time <- Sys.time()
  bootall<-AUCell_run_wt(Data,coo, 
                         BPPARAM=BPPARAM)
  AUCoutboot<-bootall$score_mat
  AUCoutbootS<-apply(AUCoutboot,2,scale)
  end_time <- Sys.time()
  
  f2=end_time - start_time
  print(f2)
  
  

  
  
  AUCbootarray_list<-list()
  AUCbootSarray_list<-list()
  
  for(ib in 1:numboot){
    l=dim(AUCout)[2]
    subind<-seq((ib-1)*l+1,ib*l)
    sub<-AUCoutboot[,subind]
    subS<-AUCoutbootS[,subind]
    colnames(sub)<-colnames(AUCout)
    colnames(subS)<-colnames(AUCout)
    
    AUCbootarray_list[[ib]]<-sub
    AUCbootSarray_list[[ib]]<-subS
  }
  
  
  AUCbootarrayLIST<-list(
    AUCbootarray_list=AUCbootarray_list,
    AUCbootSarray_list=AUCbootSarray_list
  )
  
  
  
  message('Now, annotate cell types')
  #annotation cell types#####
  
  mat_thres<-matrix(0,nrow=dim(Data)[2],ncol=length(gene_sets))
  mat_flag<-matrix(0,nrow=dim(Data)[2],ncol=length(gene_sets))
  mat_flagS<-matrix(0,nrow=dim(Data)[2],ncol=length(gene_sets))
  rownames(mat_thres)<-colnames(Data)
  rownames(mat_flag)<-colnames(Data)
  rownames(mat_flagS)<-colnames(Data)
  
  colnames(mat_thres)<-names(gene_sets)
  colnames(mat_flag)<-names(gene_sets)
  colnames(mat_flagS)<-names(gene_sets)
  
  
  vec_anno_list<-list()
  for(indspot in 1:dim(Data)[2]){
    vec_anno<-NULL
    vec_score<-NULL
    nn<-NULL
    for(indsig in 1:dim(AUCoutS)[2]){
      vec_boot=unlist(lapply(AUCbootSarray_list,function(x){return(x[indspot,indsig])}))
      f1n <- fitdistr(vec_boot,"normal")
      thre<-qnorm(thres,mean=f1n$estimate[1],sd=f1n$estimate[2])
      anno<-ifelse(AUCoutS[indspot,indsig]>thre,colnames(AUCout)[indsig],NA)
      scoree<-ifelse(AUCoutS[indspot,indsig]>thre,AUCoutS[indspot,indsig],NA)
      nn<-c(nn,ifelse(AUCoutS[indspot,indsig]>thre,colnames(AUCoutS)[indsig],NA))
      vec_anno<-c(vec_anno,anno)
      vec_score<-c(vec_score,scoree)
      mat_thres[indspot,indsig]<-thre
      mat_flag[indspot,indsig]<-ifelse(AUCoutS[indspot,indsig]>thre,AUCout[indspot,indsig],NA)
      mat_flagS[indspot,indsig]<-ifelse(AUCoutS[indspot,indsig]>thre,AUCoutS[indspot,indsig],NA)
      
    }
    names(vec_anno)<-nn
    names(vec_score)<-nn
    
    vec_anno<-vec_anno[!is.na(vec_anno)]
    vec_score<-vec_score[!is.na(vec_score)]
    vec_anno_list[[indspot]]<-sort(vec_score,decreasing=TRUE)
  }
  names(vec_anno_list)<-colnames(Data)
  
  numanno<-unlist(lapply(vec_anno_list,length))
  table(numanno)
  ANNO<-rep('nosig',length(numanno))
  names(ANNO)<-names(numanno)
  ANNO[numanno==0]='nosig'
  ANNO[numanno==1]=names(unlist(unname(vec_anno_list[numanno==1])))
  
  k=names(which(numanno>=2))
  for(j in 1:length(k)){
    v=vec_anno_list[k[j]][[1]]
    #r=abs((v[1]-v[2])/v[2])
    r=v[1]
    ANNO[k[j]]<-ifelse((r>0 ),names(v)[1],'mixture')
    #ANNO[k[j]]<-'mixture'
  }
  print(table(ANNO))
  

  return(
    list(
      AUCout=AUCout,
      AUCoutS=AUCoutS,
      vec_anno_list=vec_anno_list,
      AUCbootarrayLIST=AUCbootarrayLIST,
      ANNOTATION=ANNO,
      mat_thres=mat_thres,
      mat_flag=mat_flag,
      mat_flagS=mat_flagS
    )
  )
  
  }else{
    #no boot
    
    return(
      list(
        AUCout=AUCout,
        AUCoutS=AUCoutS
      )
    )
  }
  
  
}


#AUC assign label based on spatial location###########
annofun<-function(Sample,cluster_AUCsmooth,Assign,pick){
  
  A=as.data.frame(Assign$mat_flagS[pick,])
  A['cluster']<-cluster_AUCsmooth[pick]
  A['spots']<-pick
  A2=reshape2::melt(A)
  atemp<-as.data.frame(A2)
  
  
  mdatall<-NULL
  for(Sample_INDEX in 1:4){
    s_dat<-get(paste("sdat_",Sample,'_',Sample_INDEX,sep ='' ))
    image_path<-s_dat@misc$IMAGE$image_path_low
    
    library(moranfast)
    pltdat<-data.frame(
      X=s_dat$array_row,
      Y=s_dat$array_col)
    
    
    
    l<-levels(cluster_AUCsmooth)
    g<-names(gene_sets)
    
    ip<-intersect(rownames(pltdat),atemp$spots)
    atemp2<-atemp[which(atemp$spots %in% ip),]
    
    l2<-unique(as.character(cluster_AUCsmooth[ip]))
    
    
    for(il in 1:length(l2)){
      pp<-NULL
      for(ig in 1:length(g)){
        sl<-l2[il]
        sg<-g[ig]
        p<-atemp2$value[which(atemp2$cluster==sl & atemp2$variable==sg)]
        names(p)<-atemp2$spots[which(atemp2$cluster==sl & atemp2$variable==sg)]
        pp<-cbind(pp,p)
      }
      colnames(pp)<-g
      f<-which(colSums(pp)==0)
      if(length(f)>0){
        pp<-pp[,-f]
      }
      
      
      
      
      mdat<-NULL
      ft<-NULL
      for(j in 1:dim(pp)[2]){
        
        inn<-pp[,j]
        inn<-inn[!is.na(inn)]
        ip2<-intersect(names(inn),(ip))
        if(length(ip2)>10){
          
          moran<-moranfast(inn[ip2], pltdat[ip2,1], pltdat[ip2,2])
          vec<-c(moran$observed,moran$p.value)
          names(vec)<-c('obs','pval')
          mdat<-rbind(mdat,vec)
          ft<-c(ft,colnames(pp)[j])
        }
      }
      
      
      if(!is.null(mdat)){
        
        
        rownames(mdat)<-ft
        mdat<-as.data.frame(mdat)
        mdat['cluster']<-l[il]
        mdat['section']<-Sample_INDEX
        mdat['type']<-rownames(mdat)
        mdatall<-rbind(mdatall,mdat)
      }
    }
    
  }
  
  
  head(mdatall)
  
  l<-levels(cluster_AUCsmooth)
  anno<-NULL
  
  for(il in 1:length(l)){
    
    sub<-mdatall[which(mdatall$cluster==l[il]),]
    mm<-NULL
    for(i in 1:4){
      sub2<-sub[which(sub$section==i),]
      
      kk=intersect(which.min(sub2$pval),which(sub2$obs>0))
      mm<-c(mm,sub2$type[kk])
    }
    anno<-c(anno,getmode(mm))
    
    
  }
  names(anno)<-l
  anno
  fineanno=mapvalues(cluster_AUCsmooth,from=names(anno),to=as.character(anno))
  return(fineanno)
}



annofun_cell2loc<-function(Sample,A_cell2loc,cluster_AUCsmooth,pick){
  
  A=as.data.frame(A_cell2loc[pick,])
  A['cluster']<-cluster_AUCsmooth[pick]
  A['spots']<-pick
  A2=reshape2::melt(A)
  atemp<-as.data.frame(A2)
  
  
  mdatall<-NULL
  for(Sample_INDEX in 1:4){
    s_dat<-get(paste("sdat_",Sample,'_',Sample_INDEX,sep ='' ))
    image_path<-s_dat@misc$IMAGE$image_path_low
    
    library(moranfast)
    pltdat<-data.frame(
      X=s_dat$array_row,
      Y=s_dat$array_col)
    
    
    
    l<-levels(cluster_AUCsmooth)
    g<-colnames(A_cell2loc)
    
    ip<-intersect(rownames(pltdat),atemp$spots)
    atemp2<-atemp[which(atemp$spots %in% ip),]
    
    l2<-unique(as.character(cluster_AUCsmooth[ip]))
    
    
    for(il in 1:length(l2)){
      pp<-NULL
      for(ig in 1:length(g)){
        sl<-l2[il]
        sg<-g[ig]
        p<-atemp2$value[which(atemp2$cluster==sl & atemp2$variable==sg)]
        names(p)<-atemp2$spots[which(atemp2$cluster==sl & atemp2$variable==sg)]
        pp<-cbind(pp,p)
      }
      colnames(pp)<-g
      f<-which(colSums(pp)==0)
      if(length(f)>0){
        pp<-pp[,-f]
      }
      
      
      
      
      mdat<-NULL
      ft<-NULL
      for(j in 1:dim(pp)[2]){
        
        inn<-pp[,j]
        inn<-inn[!is.na(inn)]
        ip2<-intersect(names(inn),(ip))
        if(length(ip2)>10){
          
          moran<-moranfast(inn[ip2], pltdat[ip2,1], pltdat[ip2,2])
          vec<-c(moran$observed,moran$p.value)
          names(vec)<-c('obs','pval')
          mdat<-rbind(mdat,vec)
          ft<-c(ft,colnames(pp)[j])
        }
      }
      
      
      if(!is.null(mdat)){
        
        
        rownames(mdat)<-ft
        mdat<-as.data.frame(mdat)
        mdat['cluster']<-l[il]
        mdat['section']<-Sample_INDEX
        mdat['type']<-rownames(mdat)
        mdatall<-rbind(mdatall,mdat)
      }
    }
    
  }
  
  
  
  head(mdatall)
  
  l<-levels(cluster_AUCsmooth)
  anno<-NULL
  
  for(il in 1:length(l)){
    
    sub<-mdatall[which(mdatall$cluster==l[il]),]
    mm<-NULL
    for(i in 1:4){
      sub2<-sub[which(sub$section==i),]
      
      kk=intersect(which.min(sub2$pval),which(sub2$obs>0))
      mm<-c(mm,sub2$type[kk])
    }
    anno<-c(anno,getmode(mm))
    
    
  }
  names(anno)<-l
  anno
  fineanno=mapvalues(cluster_AUCsmooth,from=names(anno),to=as.character(anno))
  table(fineanno)
  
  
  return(fineanno)
}


#A_cell2loc: cell by cell types
#cluster_AUCsmooth： clusters
library(moranfast)
annofun_gsig<-function(s_dat,A_cell2loc,clusters){
  
  pick<-Reduce(intersect,list(
    names(clusters),
    colnames(s_dat),
    rownames(A_cell2loc))
    )
  
  A=as.data.frame(A_cell2loc[pick,])
  A['cluster']<-clusters[pick]
  A['spots']<-pick
  A2=reshape2::melt(A)
  atemp<-as.data.frame(A2)
  
  
  mdatall<-NULL

  pltdat<-data.frame(
    X=s_dat$array_row,
    Y=s_dat$array_col)
    
    
    
    l<-levels(clusters)
    g<-colnames(A_cell2loc)
    
    atemp2<-atemp
    
    l2<-unique(as.character(clusters[pick]))
    
    
    for(il in 1:length(l2)){
      pp<-NULL
      for(ig in 1:length(g)){
        sl<-l2[il]
        sg<-g[ig]
        p<-atemp2$value[which(atemp2$cluster==sl & atemp2$variable==sg)]
        names(p)<-atemp2$spots[which(atemp2$cluster==sl & atemp2$variable==sg)]
        pp<-cbind(pp,p)
      }
      colnames(pp)<-g
      
      
      f<-which(colSums(pp)==0)
      if(length(f)>0){
        pp<-pp[,-f]
      }
      
      mdat<-NULL
      ft<-NULL
      for(j in 1:dim(pp)[2]){
        
        inn<-pp[,j]
        inn<-inn[!is.na(inn)]
        ip2<-intersect(names(inn),(pick))
        if(length(ip2)>10){
          
          moran<-moranfast(inn[ip2], pltdat[ip2,1], pltdat[ip2,2])
          vec<-c(moran$observed,moran$p.value)
          names(vec)<-c('obs','pval')
          mdat<-rbind(mdat,vec)
          ft<-c(ft,colnames(pp)[j])
        }
      }
      
      
      if(!is.null(mdat)){
        rownames(mdat)<-ft
        mdat<-as.data.frame(mdat)
        mdat['cluster']<-l[il]
        mdat['type']<-rownames(mdat)
        mdatall<-rbind(mdatall,mdat)
      }
    }

  l<-levels(clusters)
  anno<-NULL
  
  for(il in 1:length(l)){
    
    sub<-mdatall[which(mdatall$cluster==l[il]),]
    mm<-NULL
    sub2<-sub
    kk=intersect(which.min(sub2$pval),which(sub2$obs>0))
    kk2=kk[which.max(sub2[kk,'obs'])]

    anno<-c(anno,sub2$type[kk2])

  }
  names(anno)<-l

  fineanno=mapvalues(clusters,from=names(anno),to=as.character(anno))
  table(fineanno)
  
  
  return(fineanno)
}







#leiden clustering from inferCNV#########
cluster_seurat_preprocess_routine_wt <- function(expr_data, k_nn, resolution_parameter, objective_function='modularity',variable.features.n=3000,algorithm=1) {
  message('algorithm: 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm')
  
  seurat_obs = CreateSeuratObject(expr_data, "assay" = "infercnv", project = "infercnv", names.field = 1)
  seurat_obs=SCTransform(seurat_obs, verbose = TRUE,method = "glmGamPoi",assay="infercnv",variable.features.n =variable.features.n)
  # seurat_obs = FindVariableFeatures(seurat_obs) # , selection.method = "vst", nfeatures = 2000
  # seurat_obs = tryCatch(FindVariableFeatures(seurat_obs), 
  #                       warning=function(w) {
  #                         flog.info(paste0("Got a warning:\n\t", w$message, "\n\nFalling back to simple Leiden clustering for this chromosome.\n"))
  #                       })
  
  if ("Seurat" %in% is(seurat_obs)) {
    all.genes <- rownames(seurat_obs)
    seurat_obs<- RunPCA(object = seurat_obs, assay = "SCT", npcs = 20)
    seurat_obs<- FindNeighbors(object = seurat_obs, assay = "SCT", reduction = "pca", dims = 1:20)
    seurat_obs <- FindClusters(object = seurat_obs, resolution = resolution_parameter,algorithm=algorithm)
    
    
    # seurat_obs <- ScaleData(seurat_obs, features = all.genes)
    # seurat_obs = RunPCA(seurat_obs, npcs=10) # only settings dims to 10 since FindNeighbors only uses 1:10 by default, if needed, could add optional settings for npcs and dims
    # seurat_obs = FindNeighbors(seurat_obs, k.param=k_nn)
    
    #graph_obj = graph_from_adjacency_matrix(seurat_obs@graphs$SCT_snn, mode="min", weighted=TRUE)
    
    # parlist<-list()
    # for(re in 1:length(resolution_parameter)){
    #   partition_obj = cluster_leiden(graph_obj, resolution_parameter=resolution_parameter[re], objective_function=objective_function)
    #   partition = partition_obj$membership
    #   parlist[[re]]<-factor(partition)
    # }
    # names(parlist)<-resolution_parameter
    
    cldat<-seurat_obs@meta.data[,grep(colnames(seurat_obs@meta.data),pattern='SCT_snn_res.')]
    cldat<-as.list(cldat)
    #cldat<-lapply(cldat,function(x){names(x)=colnames(seurat_obs)})
    parlist=cldat

  }
  else {
    parlist = leiden_simple_snn_wt(expr_data, k_nn, resolution_parameter, objective_function)
  }
  
  # names(partition)<-colnames(expr_data)
  # partition<-factor(partition)
  
  return(list(partition=parlist,seurat_obs=seurat_obs))
}
library(RANN)

leiden_simple_snn_wt <- function(expr_data, k_nn, resolution_parameter, objective_function='modularity') {
  snn <- nn2(t(expr_data), k=k_nn)$nn.idx
  
  sparse_adjacency_matrix <- sparseMatrix(
    i = rep(seq_len(ncol(expr_data)), each=k_nn), 
    j = t(snn),
    x = rep(1, ncol(expr_data) * k_nn),
    dims = c(ncol(expr_data), ncol(expr_data)),
    dimnames = list(colnames(expr_data), colnames(expr_data))
  )
  
  graph_obj = graph_from_adjacency_matrix(sparse_adjacency_matrix, mode="undirected")
  partition_obj = cluster_leiden(graph_obj, resolution_parameter=resolution_parameter, objective_function=objective_function)
  partition = partition_obj$membership
  
  names(partition)<-colnames(expr_data)
  partition<-factor(partition)
  return(partition)
}



cluster_simple_snn_wt <- function(expr_data, k_nn, resolution_parameter, objective_function='modularity',method='Louvain') {
  
  #expr_data: column sample, row feature
  
  snn <- nn2(t(expr_data), k=k_nn)$nn.idx
  
  sparse_adjacency_matrix <- sparseMatrix(
    i = rep(seq_len(ncol(expr_data)), each=k_nn), 
    j = t(snn),
    x = rep(1, ncol(expr_data) * k_nn),
    dims = c(ncol(expr_data), ncol(expr_data)),
    dimnames = list(colnames(expr_data), colnames(expr_data))
  )
  
  graph_obj = graph_from_adjacency_matrix(sparse_adjacency_matrix, mode="undirected")
  if(length(resolution_parameter)==1){
    if(method=='Leiden'){
      partition_obj = cluster_leiden(graph_obj, resolution_parameter=resolution_parameter, objective_function=objective_function)
      partition = partition_obj$membership
    }else if(method=='Louvain'){
      partition_obj = cluster_louvain(graph_obj, resolution=resolution_parameter)
      partition = partition_obj$membership
    }
    
    
    
    names(partition)<-colnames(expr_data)
    partition<-factor(partition)
    return(partition)
  }else{
    partition_list<-list()
    for(ind in 1:length(resolution_parameter)){
      re<-resolution_parameter[ind]
      
      if(method=='Leiden'){
        partition_obj = cluster_leiden(graph_obj, resolution_parameter=re, objective_function=objective_function)
        partition = partition_obj$membership
      }else if(method=='Louvain'){
        partition_obj = cluster_louvain(graph_obj, resolution=re)
        partition = partition_obj$membership
      }
      names(partition)<-colnames(expr_data)
      partition<-factor(partition)
      
      partition_list[[ind]]<-partition
      
    }
    names(partition_list)<-resolution_parameter
    
    return(partition_list)
  }
 
}



cluster_seurat_wt <- function(expr_data, k_nn, resolution_parameter, modularity.fxn=1,algorithm=1,distance.matrix=FALSE,...) {
  #expr_data: row as samples
  nnout<-FindNeighbors(expr_data,k.param=k_nn,distance.matrix=distance.matrix,...)
  
  clu<-FindClusters(nnout$snn,resolution=resolution_parameter,modularity.fxn=modularity.fxn,algorithm = algorithm,...)

  return(clu)
}


#copykat (remove heatmap plotting)##########

copykat_wt <- function(rawmat=rawdata, id.type="S", cell.line="no", ngene.chr=5,LOW.DR=0.05, UP.DR=0.1, win.size=25, norm.cell.names="", KS.cut=0.1, sam.name="", distance="euclidean", output.seg="FALSE", plot.genes="TRUE", genome="hg20", n.cores=1){
  
  start_time <- Sys.time()
  set.seed(1234)
  sample.name <- paste(sam.name,"_copykat_", sep="")
  
  print("running copykat v1.1.0")
  
  print("step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))
  
  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
  
  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
  
  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }
  
  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  
  print("step 2: annotations gene coordinates ...")
  if(genome=="hg20"){
    anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  } else if(genome=="mm10"){
    anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
    dim(rawmat)
  }
  anno.mat <- anno.mat[order(as.numeric(anno.mat$abspos), decreasing = FALSE),]
  
  # print(paste(nrow(anno.mat)," genes annotated", sep=""))
  
  ### module 3 removing genes that are involved in cell cycling
  
  if(genome=="hg20"){
    HLAs <- anno.mat$hgnc_symbol[grep("^HLA-", anno.mat$hgnc_symbol)]
    toRev <- which(anno.mat$hgnc_symbol %in% c(as.vector(cyclegenes[[1]]), HLAs))
    if(length(toRev)>0){
      anno.mat <- anno.mat[-toRev, ]
    }
  }
  #  print(paste(nrow(anno.mat)," genes after rm cell cycle genes", sep=""))
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }
  
  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")
  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }
  
  # print(paste("filtered out ", length(ToRemov2), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)
  
  #print(paste("A total of ", ncol(norm.mat), " cells, ", nrow(norm.mat), " genes after preprocessing", sep=""))
  
  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }
  
  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
  
  colnames(norm.mat.smooth) <- colnames(norm.mat)
  
  print("step 4: measuring baselines ...")
  if (cell.line=="yes"){
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
    
  } else if(length(norm.cell.names)>1){
    
    #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", sep=""))
    
    if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")
    
    d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells
    
    km <- 6
    fit <- hclust(d, method="ward.D2")
    CL <- cutree(fit, km)
    
    while(!all(table(CL)>5)){
      km <- km -1
      CL <- cutree(fit, k=km)
      if(km==2){
        break
      }
    }
    
    WNS <- "run with known normal"
    preN <- norm.cell.names
    ##relative expression using pred.normal cells
    norm.mat.relat <- norm.mat.smooth-basel
    
  }else {
    basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if (WNS =="unclassified.prediction"){
      
      basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
      basel <-basa$basel
      WNS <- basa$WNS
      
      preN <- basa$preN
      
    }
    ##relative expression using pred.normal cells
    norm.mat.relat <- norm.mat.smooth-basel
    
  }
  
  ###use a smaller set of genes to perform segmentation
  DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
  ##relative expression using pred.normal cells
  norm.mat.relat <- norm.mat.relat[which(DR2>=UP.DR),]
  
  ###filter cells
  anno.mat2 <- anno.mat[which(DR2>=UP.DR), ]
  
  ToRemov3 <- NULL
  for(i in 8:ncol(anno.mat2)){
    cell <- cbind(anno.mat2$chromosome_name, anno.mat2[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat2)[i]
      ToRemov3 <- c(ToRemov3, rm)
    }
    i<- i+1
  }
  
  if(length(ToRemov3)==ncol(norm.mat.relat)) stop ("all cells are filtered")
  
  if(length(ToRemov3)>0){
    norm.mat.relat <-norm.mat.relat[, -which(colnames(norm.mat.relat) %in% ToRemov3)]
    #print(paste("filtered out ", length(ToRemov3), " cells with less than ",ngene.chr, " genes per chr", sep=""))
  }
  
  #print(paste("final segmentation: ", nrow(norm.mat.relat), " genes; ", ncol(norm.mat.relat), " cells", sep=""))
  
  CL <- CL[which(names(CL) %in% colnames(norm.mat.relat))]
  CL <- CL[order(match(names(CL), colnames(norm.mat.relat)))]
  
  print("step 5: segmentation...")
  results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = KS.cut, n.cores=n.cores)
  
  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 50%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*KS.cut, n.cores=n.cores)
  }
  
  if(length(results$breaks)<25){
    print("too few breakpoints detected; decreased KS.cut to 75%")
    results <- CNA.MCMC(clu=CL, fttmat=norm.mat.relat, bins=win.size, cut.cor = 0.5*0.5*KS.cut, n.cores=n.cores)
  }
  
  if(length(results$breaks)<25) stop ("too few segments; try to decrease KS.cut; or improve data")
  
  colnames(results$logCNA) <- colnames(norm.mat.relat)
  results.com <- apply(results$logCNA,2, function(x)(x <- x-mean(x)))
  RNA.copycat <- cbind(anno.mat2[, 1:7], results.com)
  
  #write.table(RNA.copycat, paste(sample.name, "CNA_raw_results_gene_by_cell.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
  
  if(genome=="hg20"){
    print("step 6: convert to genomic bins...") ###need multi-core
    Aj <- convert.all.bins.hg20(DNA.mat = DNA.hg20, RNA.mat=RNA.copycat, n.cores = n.cores)
    
    uber.mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
    
    print("step 7: adjust baseline ...")
    
    if(cell.line=="yes"){
      
      mat.adj <- data.matrix(Aj$RNA.adj[, 4:ncol(Aj$RNA.adj)])
      #write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
      
      if(distance=="euclidean"){
        hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
      }else {
        hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
      }
      
      
      saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))
      
      #plot heatmap
      print("step 8: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
      
      chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
      rbPal1 <- colorRampPalette(c('black','grey'))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR,CHR)
      
      
      if (ncol(mat.adj)< 3000){
        h <- 10
      } else {
        h <- 15
      }
      
      col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
      #library(parallelDist)
      
      if(distance=="euclidean"){
        # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
        # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
        #           ColSideColors=chr1,Colv=NA, Rowv=TRUE,
        #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
        #           keysize=1, density.info="none", trace="none",
        #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
        #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        # dev.off()
        ### add a step to plot out gene by cell matrix
        if(plot.genes=="TRUE"){
          
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
          rbPal1g <- colorRampPalette(c('black','grey'))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg,CHRg)
          
          # pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
          # heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
          #           ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
          #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          #           keysize=1, density.info="none", trace="none",
          #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
          # dev.off()
        }
        #end of ploting gene by cell matrix
        
      } else {
        # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
        # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
        #           ColSideColors=chr1,Colv=NA, Rowv=TRUE,
        #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
        #           keysize=1, density.info="none", trace="none",
        #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
        #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        # dev.off()
        ### add a step to plot out gene by cell matrix
        if(plot.genes=="TRUE"){
          
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
          rbPal1g <- colorRampPalette(c('black','grey'))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg,CHRg)
          
          # pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
          # heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
          #           ColSideColors=chr1g,Colv=NA, Rowv=TRUE,
          #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          #           keysize=1, density.info="none", trace="none",
          #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
          # dev.off()
        }
        #end of ploting gene by cell matrix
      }
      
      end_time<- Sys.time()
      print(end_time -start_time)
      
      reslts <- list(cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
      names(reslts) <- c("CNAmat","hclustering")
      return(reslts)
    } else {
      ########## cell line mode ends here ####################
      
      #removed baseline adjustment
      if(distance=="euclidean"){
        hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
      }else {
        hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
      }
      hc.umap <- cutree(hcc,2)
      names(hc.umap) <- colnames(results.com)
      
      cl.ID <- NULL
      for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
      }
      
      com.pred <- names(hc.umap)
      com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
      com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
      names(com.pred) <- names(hc.umap)
      
      ################removed baseline adjustment
      results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
      results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
      results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)
      
      cf.h <- apply(results.com.rat.norm, 1, sd)
      base <- apply(results.com.rat.norm, 1, mean)
      
      adjN <- function(j){
        a <- results.com.rat[, j]
        a[abs(a-base) <= 0.25*cf.h] <- mean(a)
        a
      }
      
      
      mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
      adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
      colnames(adj.results) <- colnames(results.com.rat)
      
      #rang <- 0.5*(max(adj.results)-min(adj.results))
      #mat.adj <- adj.results/rang
      mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))
      
      print("step 8: final prediction ...")
      
      if(distance=="euclidean"){
        hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
      }else {
        hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
      }
      
      hc.umap <- cutree(hcc,2)
      names(hc.umap) <- colnames(results.com)
      
      saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))
      
      cl.ID <- NULL
      for(i in 1:max(hc.umap)){
        cli <- names(hc.umap)[which(hc.umap==i)]
        pid <- length(intersect(cli, preN))/length(cli)
        cl.ID <- c(cl.ID, pid)
        i<- i+1
      }
      
      com.preN <- names(hc.umap)
      com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
      com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
      names(com.preN) <- names(hc.umap)
      
      if(WNS=="unclassified.prediction"){
        com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
        com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
      }
      
      print("step 9: saving results...")
      
      ##add back filtered cells as not defined in prediction results
      '%!in%' <- function(x,y)!('%in%'(x,y))
      
      ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
      if(length(ndef)>0){
        res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
        colnames(res) <- c("cell.names", "copykat.pred")
      } else {
        res <- data.frame(cbind(names(com.preN), com.preN))
        colnames(res) <- c("cell.names", "copykat.pred")
      }
      ##end
      #write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
      
      ####save copycat CNA
      #write.table(cbind(Aj$RNA.adj[, 1:3], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
      
      ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
      print("step 10: ploting heatmap ...")
      my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
      
      chr <- as.numeric(Aj$DNA.adj$chrom) %% 2+1
      rbPal1 <- colorRampPalette(c('black','grey'))
      CHR <- rbPal1(2)[as.numeric(chr)]
      chr1 <- cbind(CHR,CHR)
      
      rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
      compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
      
      cells <- rbind(compreN_pred,compreN_pred)
      
      if (ncol(mat.adj)< 3000){
        h <- 10
      } else {
        h <- 15
      }
      
      col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
      
      if(distance=="euclidean"){
        # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
        # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
        #           ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
        #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
        #           keysize=1, density.info="none", trace="none",
        #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
        #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        # 
        # legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
        # dev.off()
        
        ### add a step to plot out gene by cell matrix
        if(plot.genes=="TRUE"){
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
          rbPal1g <- colorRampPalette(c('black','grey'))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg,CHRg)
          
          # pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
          # heatmap.3(t(results.com),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
          #           ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
          #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          #           keysize=1, density.info="none", trace="none",
          #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
          # dev.off()
        }
        #end of ploting gene by cell matrix
        
        
        
      } else {
        # jpeg(paste(sample.name,"heatmap.jpeg",sep=""), height=h*250, width=4000, res=100)
        # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
        #           ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
        #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
        #           keysize=1, density.info="none", trace="none",
        #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
        #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
        # 
        # legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=1)
        # 
        # dev.off()
        ### add a step to plot out gene by cell matrix
        if(plot.genes=="TRUE"){
          dim(results.com)
          rownames(results.com) <- anno.mat2$hgnc_symbol
          chrg <- as.numeric(anno.mat2$chrom) %% 2+1
          rbPal1g <- colorRampPalette(c('black','grey'))
          CHRg <- rbPal1(2)[as.numeric(chrg)]
          chr1g <- cbind(CHRg,CHRg)
          
          # pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
          # heatmap.3(t(results.com),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
          #           ColSideColors=chr1g,RowSideColors=cells, Colv=NA, Rowv=TRUE,
          #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
          #           keysize=1, density.info="none", trace="none",
          #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
          #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
          # dev.off()
        }
        #end of ploting gene by cell matrix
      }
      
      if(output.seg=="TRUE"){
        print("generating seg files for IGV viewer")
        
        thisRatio <- cbind(Aj$RNA.adj[, 1:3], mat.adj)
        Short <- NULL
        chr <- rle(thisRatio$chrom)[[2]]
        
        for(c in 4:ncol(thisRatio))
        {
          for (x in 1:length(chr)){
            thisRatio.sub <- thisRatio[which(thisRatio$chrom==chr[x]), ]
            seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]
            
            rle.length.sub <- rle(thisRatio.sub[,c])[[1]]
            
            num.mark.sub <- seq(1,length(rle.length.sub),1)
            loc.start.sub <-seq(1,length(rle.length.sub),1)
            loc.end.sub <- seq(1,length(rle.length.sub),1)
            
            len <-0
            j <-1
            
            for (j in 1: length(rle.length.sub)){
              num.mark.sub[j] <- rle.length.sub[j]
              loc.start.sub[j] <- thisRatio.sub$chrompos[len+1]
              len <- num.mark.sub[j]+len
              loc.end.sub[j] <- thisRatio.sub$chrompos[len]
              j <- j+1
            }
            
            ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
            chrom <- rep(chr[x], times=length(rle.length.sub))
            Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
            Short <- rbind(Short, Short.sub)
            x <- x+1
          }
          c<- c+1
        }
        
        colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
        head(Short)
        write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")
        
      }
      end_time<- Sys.time()
      print(end_time -start_time)
      
      reslts <- list(res, cbind(Aj$RNA.adj[, 1:3], mat.adj), hcc)
      names(reslts) <- c("prediction", "CNAmat","hclustering")
      return(reslts)
    }
    
  }
  
  if(genome=="mm10") {
    uber.mat.adj <- data.matrix(results.com)
    dim(uber.mat.adj)
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(uber.mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(uber.mat.adj, method = distance)), method = "ward.D")
    }
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }
    
    com.pred <- names(hc.umap)
    com.pred[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.pred[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.pred) <- names(hc.umap)
    
    ################removed baseline adjustment
    results.com.rat <- uber.mat.adj-apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
    
    results.com.rat <- apply(results.com.rat,2,function(x)(x <- x-mean(x)))
    results.com.rat.norm <- results.com.rat[,which(com.pred=="diploid")]; dim(results.com.rat.norm)
    
    cf.h <- apply(results.com.rat.norm, 1, sd)
    base <- apply(results.com.rat.norm, 1, mean)
    
    adjN <- function(j){
      a <- results.com.rat[, j]
      a[abs(a-base) <= 0.25*cf.h] <- mean(a)
      a
    }
    
    
    mc.adjN <-  parallel::mclapply(1:ncol(results.com.rat),adjN, mc.cores = n.cores)
    adj.results <- matrix(unlist(mc.adjN), ncol = ncol(results.com.rat), byrow = FALSE)
    colnames(adj.results) <- colnames(results.com.rat)
    
    #rang <- 0.5*(max(adj.results)-min(adj.results))
    #mat.adj <- adj.results/rang
    mat.adj <- t(t(adj.results)-apply(adj.results,2,mean))
    
    print("step 8: final prediction ...")
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(mat.adj),threads =n.cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(mat.adj, method = distance)), method = "ward.D")
    }
    
    hc.umap <- cutree(hcc,2)
    names(hc.umap) <- colnames(results.com)
    
    saveRDS(hcc, file = paste(sample.name,"clustering_results.rds",sep=""))
    
    cl.ID <- NULL
    for(i in 1:max(hc.umap)){
      cli <- names(hc.umap)[which(hc.umap==i)]
      pid <- length(intersect(cli, preN))/length(cli)
      cl.ID <- c(cl.ID, pid)
      i<- i+1
    }
    
    com.preN <- names(hc.umap)
    com.preN[which(hc.umap == which(cl.ID==max(cl.ID)))] <- "diploid"
    com.preN[which(hc.umap == which(cl.ID==min(cl.ID)))] <- "aneuploid"
    names(com.preN) <- names(hc.umap)
    
    if(WNS=="unclassified.prediction"){
      com.preN[which(com.preN == "diploid")] <- "c1:diploid:low.conf"
      com.preN[which(com.preN == "aneuploid")] <- "c2:aneuploid:low.conf"
    }
    
    print("step 9: saving results...")
    
    ##add back filtered cells as not defined in prediction results
    '%!in%' <- function(x,y)!('%in%'(x,y))
    ndef <- colnames(rawmat)[which(colnames(rawmat) %!in% names(com.preN))]
    if(length(ndef)>0){
      res <- data.frame(cbind(c(names(com.preN),ndef), c(com.preN, rep("not.defined",length(ndef)))))
      colnames(res) <- c("cell.names", "copykat.pred")
    } else {
      res <- data.frame(cbind(names(com.preN), com.preN))
      colnames(res) <- c("cell.names", "copykat.pred")
    }
    ##end
    write.table(res, paste(sample.name, "prediction.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
    
    ####save copycat CNA
    #write.table(cbind(anno.mat2[, 1:7], mat.adj), paste(sample.name, "CNA_results.txt", sep=""), sep="\t", row.names = FALSE, quote = F)
    
    ####%%%%%%%%%%%%%%%%%next heatmaps, subpopulations and tSNE overlay
    print("step 10: ploting heatmap ...")
    my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
    
    rownames(mat.adj) <- anno.mat2$mgi_symbol
    chrg <- as.numeric(anno.mat2$chromosome_name) %% 2+1
    rle(as.numeric(anno.mat2$chromosome_name))
    rbPal1g <- colorRampPalette(c('black','grey'))
    CHRg <- rbPal1g(2)[as.numeric(chrg)]
    chr1g <- cbind(CHRg,CHRg)
    
    
    rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
    compreN_pred <- rbPal5(2)[as.numeric(factor(com.preN))]
    
    cells <- rbind(compreN_pred,compreN_pred)
    
    if (ncol(mat.adj)< 3000){
      h <- 10
    } else {
      h <- 15
    }
    
    col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
    
    if(distance=="euclidean"){
      
      # pdf(paste(sample.name,"with_genes_heatmap.pdf",sep=""), height=h*2.5, width=40)
      # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =n.cores, method = distance), hclustfun = function(x) hclust(x, method="ward.D"),
      #           ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
      #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
      #           keysize=1, density.info="none", trace="none",
      #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
      #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
      # dev.off()
      # 
      
    } else {
      
      # pdf(paste(sample.name,"with_genes_heatmap1.pdf",sep=""), height=h*2.5, width=40)
      # heatmap.3(t(mat.adj),dendrogram="r", distfun = function(x) as.dist(1-cor(t(x), method = distance)), hclustfun = function(x) hclust(x, method="ward.D"),
      #           ColSideColors=chr1g,RowSideColors=cells,Colv=NA, Rowv=TRUE,
      #           notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
      #           keysize=1, density.info="none", trace="none",
      #           cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
      #           symm=F,symkey=F,symbreaks=T,cex=1, main=paste(WNS1,"; ",WNS, sep=""), cex.main=4, margins=c(10,10))
      # dev.off()
      
      #end of ploting gene by cell matrix
    }
    
    if(output.seg=="TRUE"){
      print("generating seg files for IGV viewer")
      
      thisRatio <- cbind(anno.mat2[, c(2,3,1)], mat.adj)
      Short <- NULL
      chr <- rle(thisRatio$chromosome_name)[[2]]
      
      for(c in 4:ncol(thisRatio))
      {
        for (x in 1:length(chr)){
          thisRatio.sub <- thisRatio[which(thisRatio$chromosome_name==chr[x]), ]
          seg.mean.sub <- rle(thisRatio.sub[,c])[[2]]
          
          rle.length.sub <- rle(thisRatio.sub[,c])[[1]]
          
          num.mark.sub <- seq(1,length(rle.length.sub),1)
          loc.start.sub <-seq(1,length(rle.length.sub),1)
          loc.end.sub <- seq(1,length(rle.length.sub),1)
          
          len <-0
          j <-1
          
          for (j in 1: length(rle.length.sub)){
            num.mark.sub[j] <- rle.length.sub[j]
            loc.start.sub[j] <- thisRatio.sub$start_position[len+1]
            len <- num.mark.sub[j]+len
            loc.end.sub[j] <- thisRatio.sub$start_position[len]
            j <- j+1
          }
          
          ID <- rep(colnames(thisRatio[c]), times=length(rle.length.sub))
          chrom <- rep(chr[x], times=length(rle.length.sub))
          Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
          Short <- rbind(Short, Short.sub)
          x <- x+1
        }
        c<- c+1
      }
      
      colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","seg.mean")
      
      #write.table(Short, paste(sample.name, "CNA_results.seg", sep=""), row.names = FALSE, quote=FALSE, sep="\t")
      
    }
    end_time<- Sys.time()
    print(end_time -start_time)
    
    reslts <- list(res, cbind(anno.mat2[, 1:7], mat.adj), hcc)
    names(reslts) <- c("prediction", "CNAmat","hclustering")
    return(reslts)
    
  }
}

#enricher GO terms########
ENRICHER_WT<-function(allgenes,topgenes,org=org.Mm.eg.db,DROP_GO_PATTERN=DROP_GO_PATTERN){
  anno <- AnnotationDbi::select(x=org,
                                keys=allgenes,
                                columns=c("ENSEMBL","SYMBOL", "GENENAME","GO"), keytype="SYMBOL")
  gmt_in <-anno[which(anno$ONTOLOGY=='BP'),c(4,1)]
  colnames(gmt_in)<-c('term','gene')
  gmt_in$term=goIdToTerm(gmt_in$term)
  gmt_in2=gmt_in[-grep(gmt_in$term,pattern=DROP_GO_PATTERN),]
  
  enriched <- clusterProfiler::enricher(gene = topgenes,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        universe = allgenes,
                                        TERM2GENE = gmt_in2)
  qq=enriched@result
  return(qq)
}
ENRICHER_WT_SELF<-function(allgenes,topgenes,org=org.Mm.eg.db,DROP_GO_PATTERN=DROP_GO_PATTERN,list_marker){
  
  vecc<-unlist(list_marker,use.names = FALSE)
  vecname<-rep(names(list_marker),unlist(lapply(list_marker,length)))
  
  
  gmt_in2=data.frame(
    term=vecname,
    gene=vecc
  )
  
  enriched <- clusterProfiler::enricher(gene = topgenes,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        universe = allgenes,
                                        TERM2GENE = gmt_in2)
  
  if(!is.null(enriched)){
    qq=enriched@result
  }else{
    qq=NA
  }

  return(qq)
}



#from pRoloc PACKAGE, fail install in lovelace#######
goIdToTerm <- function(x, names = TRUE, keepNA = TRUE) {
  stopifnot(requireNamespace("GO.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ans <- rep(NA_character_, length(x))
  names(ans) <- x
  ids <- AnnotationDbi::GOID(GO.db::GOTERM)
  i <- match(x, ids)
  k <- which(!is.na(i))
  res <- AnnotationDbi::Term(GO.db::GOTERM[i[k]])
  ans[k] <- res
  if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
  if (!names) names(ans) <- NULL
  return(ans)
}

##' @rdname goIdToTerm
goTermToId <- function(x, names = TRUE, keepNA = TRUE) {
  stopifnot(requireNamespace("GO.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ans <- rep(NA_character_, length(x))
  names(ans) <- x
  terms <- AnnotationDbi::Term(GO.db::GOTERM)
  i <- match(x, terms)
  k <- which(!is.na(i))
  res <- AnnotationDbi::GOID(GO.db::GOTERM[i[k]])
  ans[k] <- res
  if (!keepNA) ans[is.na(ans)] <- names(ans[is.na(ans)])
  if (!names) names(ans) <- NULL
  return(ans)
}
#SCENIC wt (mayby use pyscenic)#############
# library(SCENIC)
# library(SCopeLoomR)
# 
# 
# initializeScenic_wt<-function (org = NULL, dbDir = "databases", dbs = NULL, datasetTitle = "", 
#                                nCores = 4, dbIndexCol = "features") 
# {
#   inputDataset <- list(org = org, datasetTitle = datasetTitle, 
#                        cellInfo = c("int/cellInfo.Rds", NA), colVars = c("int/colVars.Rds", 
#                                                                          NA), int_01 = c("int/cellColorNgenes.Rds", NA))
#   if (!org %in% c("mgi", "hgnc", "dmel")) 
#     stop("'org' should be one of: mgi, hgnc, dmel.")
#   defaultDBs <- FALSE
#   if (is.null(dbs)) {
#     data(defaultDbNames)
#     dbs <- defaultDbNames[[org]]
#   }
#   dbsFound <- unlist(unname(lapply(dbs, function(x) setNames(file.exists(file.path(dbDir, 
#                                                                                    x)), unname(file.path(dbDir, x))))))
#   if (any(!dbsFound)) {
#     stop("The following RcisTarget databases were not found: ", 
#          paste(paste0("\n- ", names(dbsFound[which(!dbsFound)])), 
#                collapse = " "), "\nMake sure the arguments 'dbDir' and 'dbs' are correct.")
#     dbs <- NULL
#   }
#   else {
#     message("Motif databases selected: ", paste(paste0("\n  ", 
#                                                        dbs, collapse = " ")))
#   }
#   loadAttempt <- sapply(dbs, function(x) dbLoadingAttempt(file.path(dbDir, 
#                                                                     x), indexCol = dbIndexCol))
#   if (any(!loadAttempt)) 
#     warning("It was not possible to load the following databses; check whether they are downloaded correctly: \n", 
#             paste(dbs[which(!loadAttempt)], collapse = "\n"))
#   db_mcVersion <- dbVersion(dbs)
#   scenicSettings = list(dbs = dbs, dbDir = dbDir, db_mcVersion = db_mcVersion, 
#                         db_annotFiles = NULL, verbose = TRUE, nCores = nCores, 
#                         seed = 123, devType = "pdf", modules = list(weightThreshold = 0.001), 
#                         regulons = list(), aucell = list(smallestPopPercent = 0.25), 
#                         defaultTsne = list(dims = 50, perpl = 30, aucType = "AUC"), 
#                         tSNE_filePrefix = "int/tSNE")
#   scenicFiles <- list(output = c(s2_motifEnrichment = "output/Step2_MotifEnrichment.tsv", 
#                                  s2_motifEnrichmentHtml = "output/Step2_MotifEnrichment_preview.html", 
#                                  s2_regulonTargetsInfo = "output/Step2_regulonTargetsInfo.tsv", 
#                                  s3_AUCheatmap = "output/Step3_RegulonActivity_heatmap", 
#                                  s3_AUCtSNE_colAct = "output/Step3_RegulonActivity_tSNE_colByActivity", 
#                                  s3_AUCtSNE_colProps = "output/Step3_RegulonActivity_tSNE_colByCellProps", 
#                                  s4_boxplotBinaryActivity = "output/Step4_BoxplotActiveCellsRegulon", 
#                                  s4_binaryActivityHeatmap = "output/Step4_BinaryRegulonActivity_Heatmap_", 
#                                  s4_binarytSNE_colAct = "output/Step4_BinaryRegulonActivity_tSNE_colByActivity", 
#                                  s4_binarytSNE_colProps = "output/Step4_BinaryRegulonActivity_tSNE_colByCellProps", 
#                                  loomFile = "output/scenic.loom"), int = list(genesKept = c("int/1.1_genesKept.Rds", 
#                                                                                             TRUE), corrMat = c("int/1.2_corrMat.Rds", TRUE), genie3wm = c("int/1.3_GENIE3_weightMatrix.Rds", 
#                                                                                                                                                           FALSE), genie3ll = c("int/1.4_GENIE3_linkList.Rds", TRUE), 
#                                                                               genie3weighPlot = c("int/1.5_weightPlot", TRUE), tfModules_asDF = c("int/1.6_tfModules_asDF.Rds", 
#                                                                                                                                                   TRUE), tfModules_forEnrichment = c("int/2.1_tfModules_forMotifEnrichmet.Rds", 
#                                                                                                                                                                                      FALSE), motifs_AUC = c("int/2.2_motifs_AUC.Rds", 
#                                                                                                                                                                                                             FALSE), motifEnrichment_full = c("int/2.3_motifEnrichment.Rds", 
#                                                                                                                                                                                                                                              FALSE), motifEnrichment_selfMotifs_wGenes = c("int/2.4_motifEnrichment_selfMotifs_wGenes.Rds", 
#                                                                                                                                                                                                                                                                                            FALSE), regulonTargetsInfo = c("int/2.5_regulonTargetsInfo.Rds", 
#                                                                                                                                                                                                                                                                                                                           NA), regulons = c("int/2.6_regulons_asGeneSet.Rds", 
#                                                                                                                                                                                                                                                                                                                                             NA), regulons_incidMat = c("int/2.6_regulons_asIncidMat.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                        NA), aucell_regulons = c("int/3.1_regulons_forAUCell.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                 NA), aucell_genesStatsPlot = c("int/3.2_aucellGenesStats", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                NA), aucell_rankings = c("int/3.3_aucellRankings.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                         NA), aucell_regulonAUC = c("int/3.4_regulonAUC.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    NA), aucell_thresholds = c("int/3.5_AUCellThresholds.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               NA), aucell_thresholdsTxt = c("int/3.5_AUCellThresholds_Info.tsv", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             NA), aucell_binary_full = c("int/4.1_binaryRegulonActivity.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         NA), aucell_binary_nonDupl = c("int/4.2_binaryRegulonActivity_nonDupl.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NA), aucell_regulonSelection = c("int/4.3_regulonSelections.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         NA), aucell_binaryRegulonOrder = c("int/4.4_binaryRegulonOrder.Rds", 
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            NA)))
#   scenicFiles$output <- cbind(fileName = scenicFiles$output)
#   scenicFiles$int <- do.call(rbind, scenicFiles$int)
#   colnames(scenicFiles$int) <- c("fileName", "keep")
#   scenicFiles$int <- scenicFiles$int[, "fileName", drop = FALSE]
#   dir.create("int", showWarnings = FALSE)
#   dir.create("output", showWarnings = FALSE)
#   object <- new("ScenicOptions", inputDatasetInfo = inputDataset, 
#                 status = list(current = 0, values = c(`0` = "SCENIC initialized", 
#                                                       `1` = "Co-expression modules", `2` = "Regulons", 
#                                                       `3` = "Cells scored", `4` = "SCENIC run completed")), 
#                 settings = scenicSettings, fileNames = scenicFiles)
#   motifAnnot <- getDbAnnotations(object)
#   featuresWithAnnot <- checkAnnots_wt(object, motifAnnot,dbIndexCol=dbIndexCol)
#   if (any(featuresWithAnnot == 0)) 
#     message("Missing annotations for: \n", paste("\t", names(which(featuresWithAnnot == 
#                                                                      0))))
#   return(object)
# }
# 
# 
# checkAnnots_wt<-function (object, motifAnnot,dbIndexCol) 
# {
#   allFeaturesInAnnot <- unlist(motifAnnot[, 1])
#   featuresWithAnnot <- lapply(getDatabases(object), function(dbFile) {
#     rnktype = dbIndexCol
#     nRnks <- getRanking(RcisTarget::importRankings(dbFile, 
#                                                    columns = rnktype))
#     nRnks <- dplyr::pull(nRnks, rnktype)
#     length(intersect(allFeaturesInAnnot, nRnks))/length(unique(c(allFeaturesInAnnot, 
#                                                                  nRnks)))
#   })
#   return(featuresWithAnnot)
# }
# 
# regulon_specificity_scores<-function(mat_AUC,cluster){
#   mat_AUC_short <- sapply(split(names(cluster), cluster),
#                           function(cells) rowMeans(mat_AUC[,cells]))
#   return(mat_AUC_short)
# }


#fun z score######
fun_z<-function(data){
  tempp<-as.data.frame(apply(data,2,scale))
  rownames(tempp)<-rownames(data)
  return(tempp)
}