#load packages########
library(org.Mm.eg.db)
library(org.Hs.eg.db)


library(genefilter)
library(topGO)

library(babelgene)

DROP_GO_PATTERN='hair|sex|heart|liver|gonad|lung|sexual|skin|cardiac|molting|behavior|bone|animal|urogenital|renal|kidney|digestive|digestion|ear|light|salt|pregnancy|ovulation|eye|sperm|odontogenesis|penile'

Map_H_to_M<-function(gh){
  
  outs=orthologs(genes = gh, species = "Mus musculus",human=TRUE)
  return(data.frame(H=outs$human_symbol,M=outs$symbol))
}

Map_M_to_H<-function(gm){
  outs=orthologs(genes = gm, species = "Mus musculus",human=FALSE)
  return(data.frame(H=outs$human_symbol,M=outs$symbol))
}




#Prepare fgsea database#######

# m_df_H<- msigdbr(species = "Homo sapiens", category = cat) 
# Annot.pathway2_H=as.data.frame(levels(as.factor(m_df_H$gs_name)))
# names(Annot.pathway2_H)="pathway"
# fgsea_sets_H<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
# 
# m_df_M<- msigdbr(species = "Mus musculus", category = cat) 
# Annot.pathway2_M=as.data.frame(levels(as.factor(m_df_M$gs_name)))
# names(Annot.pathway2_M)="pathway"
# fgsea_sets_M<- m_df_M %>% split(x = .$gene_symbol, f = .$gs_name)



#The text, condition treated vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated). source:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#note-on-factor-levels
###extract significant genes#########
extract_results<-function(DESeq2Res,padj_threshold=0.1){
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$padj), ]
  DESeq2Res <- DESeq2Res[ !is.na(DESeq2Res$pvalue), ]
  DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]
  library(fdrtool)
  FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, statistic= "normal", plot = FALSE)
  FDR.DESeq2Res$param[1, "sd"]
  DESeq2Res[,"padj"]  <- p.adjust(FDR.DESeq2Res$pval, method = "BH")
  ## Gene ontology enrichment analysis######
  sigGenes <- rownames(subset(DESeq2Res, padj < padj_threshold))

  return(list(
    DESeq2Res=DESeq2Res,
    sigGenes=sigGenes
  ))
  
}

###extract background genes#######

extract_bc<-function(AnnotationDb=org.Mm.eg.db,input_genes,sigGenes,overallBaseMean,keytype="SYMBOL"){
  #columns(org.Mm.eg.db)
  anno <- AnnotationDbi::select(x=AnnotationDb, 
                                keys=input_genes, 
                                columns=c("ENSEMBL","SYMBOL", "GENENAME","GO"), keytype=keytype)
  #head( keys(org.Mm.eg.db, keytype="SYMBOL") )
  anSig<-anno[anno[,keytype] %in% sigGenes,]
  
  #overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
  sig_idx <- match(anSig[,keytype], rownames(overallBaseMean))
  
  #find background genes for every significant gene
  backG <- c()
  
  for(i in sig_idx){
    ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
    backG <- c(backG, ind)
  }
  
  backG <- unique(backG)
  backG <- rownames(overallBaseMean)[backG]
  #remove significant genes from background genes
  backG <- setdiff(backG,  anSig[,keytype])
  
  return(list(
    backG=backG,
    anSig=anSig,
    overallBaseMean=overallBaseMean
  ))
}


####run topGO#########

extract_GO<-function(overallBaseMean,anSig,backG,keytype='SYMBOL',topNodes=200,nodeSize=5,mapping="org.Mm.eg.db"){
  onts = c( "MF", "BP", "CC" )
  
  geneIDs = rownames(overallBaseMean)
  inUniverse = geneIDs %in% c(anSig[,keytype],  backG) 
  inSelection =  geneIDs %in% anSig[,keytype] 
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- geneIDs[inUniverse]
  
  
  tab = as.list(onts)
  names(tab) = onts
  for(i in 1:3){
    
    ## prepare data
    tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=nodeSize, annot=annFUN.org, mapping=mapping, ID = keytype)
    
    ## run tests
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
    ## look at results
    tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                          Fisher.classic = resultTopGO.classic,
                          orderBy = "Fisher.classic" , topNodes = topNodes)
    
  }
  names(tab)<-onts
  return(tab)
}

###wrapper for any two groups#######
wrapper_deseq<-function(AnnotationDb=org.Mm.eg.db,DESeq2Res,keytype="SYMBOL",padj_threshold=0.1,topNodes=200,nodeSize=5,mapping="org.Mm.eg.db"){
  
  exout<-extract_results(DESeq2Res=DESeq2Res,padj_threshold=padj_threshold)
  message(blue('Finished: extract results'))
  
  DESeq2Res=exout$DESeq2Res
  sigGenes<-exout$sigGenes
  ## Gene ontology enrichment analysis######
  overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
  
  
  bcout<-extract_bc(AnnotationDb=AnnotationDb,input_genes=rownames(DESeq2Res),sigGenes=sigGenes,overallBaseMean=overallBaseMean,keytype=keytype)
  
  message(blue('Finished: extract background genes'))
  
  anSig=bcout$anSig
  backG=bcout$backG
  
  ##Running topGO#######
  tab<-extract_GO(overallBaseMean=overallBaseMean,anSig=anSig,backG=backG,keytype=keytype,topNodes=topNodes,nodeSize=nodeSize,mapping=mapping)
  
  return(
    list(tab=tab,
         sigGenes=sigGenes,
         overallBaseMean=overallBaseMean,
         bcout=bcout)
  )
  
}

wrapper_deseq_more<-function(paris_used,data_meta_merged,counts_dat_used,condition="type",AnnotationDb=org.Mm.eg.db,keytype="SYMBOL",padj_threshold=0.1,topNodes=200,nodeSize=5,mapping="org.Mm.eg.db"){
  subsamples<-which(data_meta_merged[,condition] %in% paris_used)
  
  submeta<-data_meta_merged[subsamples,]

  
  DE_1 <- DESeqDataSetFromMatrix(countData=counts_dat_used[,subsamples], colData=submeta, design=as.formula(paste("~",condition,sep='')))
  
  
  levels(DE_1$type)<-paris_used
  DE_1=DESeq(DE_1)
  
  #resultsNames(DE_1)
  RE_1<-results(DE_1, c(condition, paris_used)) 
  tab_1<-wrapper_deseq(AnnotationDb=AnnotationDb,DESeq2Res=RE_1,keytype=keytype,padj_threshold=padj_threshold,topNodes=topNodes,nodeSize=nodeSize,mapping=mapping)

  
  return(
    list(
      DE=DE_1,
      RE=RE_1,
      sigGenes=tab_1[[2]],
      overallBaseMean=tab_1[[3]],
      bcout=tab_1[[4]],
      tab=tab_1[[1]])
  )
}

wrapper_deseq_up_down<-function(AnnotationDb=org.Mm.eg.db,DESeq2Res,keytype="SYMBOL",padj_threshold=0.1,topNodes=200,nodeSize=5,mapping="org.Mm.eg.db",lf_up_thres=1,lf_down_thres=-1,bc_genes=NULL){
  
  exout<-extract_results(DESeq2Res=DESeq2Res,padj_threshold=padj_threshold)
  message(blue('Finished: extract results'))
  
  DESeq2Res=exout$DESeq2Res
  sigGenes<-exout$sigGenes
  
  DESeq2Res_sub<-as.data.frame(DESeq2Res)[sigGenes,]
  sigGenes_up<-rownames(DESeq2Res_sub)[which(DESeq2Res_sub$log2FoldChange>(lf_up_thres))]
  sigGenes_down<-rownames(DESeq2Res_sub)[which(DESeq2Res_sub$log2FoldChange<(lf_down_thres))]

  
  
  
  ## Gene ontology enrichment analysis######
  overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])
  
  
  bcout_up<-extract_bc(AnnotationDb=AnnotationDb,input_genes=rownames(DESeq2Res),sigGenes=sigGenes_up,overallBaseMean=overallBaseMean,keytype=keytype)
  bcout_down<-extract_bc(AnnotationDb=AnnotationDb,input_genes=rownames(DESeq2Res),sigGenes=sigGenes_down,overallBaseMean=overallBaseMean,keytype=keytype)
  
  message(blue('Finished: extract background genes'))
  
  anSig_up=bcout_up$anSig
  anSig_down=bcout_down$anSig

  if(is.null(bc_genes)){
    backG_up=bcout_up$backG
    backG_down=bcout_down$backG
  }else{
    backG_up=bc_genes
    backG_down=bc_genes
  }
  
  
  ##Running topGO#######
  tab_up<-extract_GO(overallBaseMean=overallBaseMean,anSig=anSig_up,backG=backG_up,keytype=keytype,topNodes=topNodes,nodeSize=nodeSize,mapping=mapping)
  tab_down<-extract_GO(overallBaseMean=overallBaseMean,anSig=anSig_down,backG=backG_down,keytype=keytype,topNodes=topNodes,nodeSize=nodeSize,mapping=mapping)
  
  return(
    list(
      tab_up=tab_up,
      tab_down=tab_down,
      overallBaseMean=overallBaseMean,
      bcout_up=bcout_up,
      bcout_down=bcout_down)
  )
  
}



##extract go terms (CEMiTool)#########
#prepare annotation file
# anno <- AnnotationDbi::select(x=org.Mm.eg.db, 
#                               keys=rownames(DESeq2Table), 
#                               columns=c("ENSEMBL","SYMBOL", "GENENAME","GO"), keytype="SYMBOL")
# gmt_in <-anno[which(anno$ONTOLOGY=='BP'),c(4,1)]
# colnames(gmt_in)<-c('term','gene')
# gmt_in$term=goIdToTerm(gmt_in$term)

extractGO_CEMiTool_simple<-function(topgenes,allgenes,gmt_list){
  enriched <- clusterProfiler::enricher(gene = topgenes,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        universe = allgenes,
                                        TERM2GENE = gmt_list)
  
  
  if(!is.null(enriched)){
    results <- enriched@result
    results['Annotated']= as.numeric(do.call(rbind,str_split(results$BgRatio,pattern='/'))[,1])
    results['Significant']= as.numeric(do.call(rbind,str_split(results$GeneRatio,pattern='/'))[,1])
    
    
    colnames(results)[1]<-'Term'
    results$Term<-factor(results$Term,levels=rev(results$Term))
  } else{
    results=NULL
  }
  
  
  return(results)
}


extractGO_CEMiTool<-function(
  DESeq2Res,
  bc_genes=NULL,#use all genes from DESeq2Res as background if is NULL
  padj_threshold=0.1,
  lf_thres=1,
  gmt_list
  ){
  exout<-extract_results(DESeq2Res=DESeq2Res,padj_threshold=padj_threshold)
  message(blue('Finished: extract results'))
  
  DESeq2Res=exout$DESeq2Res
  sigGenes<-exout$sigGenes
  
  DESeq2Res_sub<-as.data.frame(DESeq2Res)[sigGenes,]
  
  if(lf_thres>=0){
    topgenes<-rownames(DESeq2Res_sub)[which(DESeq2Res_sub$log2FoldChange>(lf_thres))]
  }else{
    topgenes<-rownames(DESeq2Res_sub)[which(DESeq2Res_sub$log2FoldChange<(lf_thres))]
  }

  #use all genes as background genes
  if(!is.null(bc_genes)){
    allgenes=bc_genes
  }else{
    allgenes=rownames(DESeq2Res)
  }
  
  
  enriched <- clusterProfiler::enricher(gene = topgenes,
                                        pvalueCutoff = 1,
                                        qvalueCutoff = 1,
                                        universe = allgenes,
                                        TERM2GENE = gmt_list)
  
  
  if(!is.null(enriched)){
    results <- enriched@result
    results['Annotated']= as.numeric(do.call(rbind,str_split(results$BgRatio,pattern='/'))[,1])
    results['Significant']= as.numeric(do.call(rbind,str_split(results$GeneRatio,pattern='/'))[,1])
    
    
    colnames(results)[1]<-'Term'
    results$Term<-factor(results$Term,levels=rev(results$Term))
  } else{
    results=NULL
  }

  
  return(results)
}









#single data#######
##input data, all genes, significant genes
## return GO terms


GOout<-function(
  DATA,
  GENES,
  sigGenes,
  typee='mouse',
  AnnotationDb=org.Mm.eg.db,
  mapping="org.Mm.eg.db",
  topNodes=100
){
  if(typee=='mouse'){
    sigGenes=firstup(sigGenes)
  } else if(typee=='human'){
    sigGenes=toupper(sigGenes)
  }
  
  
  anno <- AnnotationDbi::select(x=AnnotationDb, 
                                keys=(GENES), 
                                columns=c("ENTREZID","ENSEMBL","SYMBOL", "GENENAME","GO"), keytype="SYMBOL")
  anno<-na.omit(anno)

  gee<-intersect(anno$SYMBOL,sigGenes)
  
  
  overallBasemean<-matrix(rowMeans(DATA))
  rownames(overallBasemean)<-(rownames(DATA))
  bcout<-extract_bc(AnnotationDb=AnnotationDb,input_genes=GENES,sigGenes=sigGenes,overallBaseMean=overallBasemean,keytype="SYMBOL")
  
  GOout<-extract_GO(overallBaseMean=bcout$overallBaseMean,anSig=bcout$anSig,backG=bcout$backG,keytype='SYMBOL',topNodes=topNodes,nodeSize=5,mapping=mapping)
  
  return(GOout)
}






###plots##########

dotplot_wt_dat<-function(GO_input,topNodes){
  GO_input$Fisher.classic<-as.numeric(GO_input$Fisher.classic)
  GO_input$Fisher.classic[is.na(GO_input$Fisher.classic)]=  min(GO_input$Fisher.classic[!is.na(GO_input$Fisher.classic)])/2
  
  
  if(sum(duplicated(GO_input$Term))>0){
    GO_input<-GO_input[!duplicated(GO_input$Term),]
  }
  GO_input$Term<-factor(GO_input$Term,levels=rev(GO_input$Term))
  return(GO_input[1:topNodes,])
}


dotplot_wt<-function(GO_input,topNodes=15){
  #GO_input: output from GOout
  GO_input<-dotplot_wt_dat(GO_input=GO_input,topNodes=topNodes)
  
  gout<-ggplot(data=GO_input[1:topNodes,],aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(Fisher.classic)))+
    scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
    geom_point()+
    labs(x='Significant/Annotated',y='',color='-log(pval)')+
    theme_bw()
  return(gout)
  
}


dotplot_enricher_wt<-function(GO_input,topNodes=15){
  #GO_input: output from GOout

  gout<-ggplot(data=GO_input[1:topNodes,],aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(qvalue)))+
    scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
    geom_point()+
    labs(x='Significant/Annotated',y='',color='-log(qvalue)')+
    theme_bw()
  return(gout)
  
}




dotplot_wt_up_down<-function(GO_input_up,GO_input_down,topNodes=15){
  GO_input_up<-dotplot_wt_dat(GO_input=GO_input_up,topNodes=topNodes)
  GO_input_down<-dotplot_wt_dat(GO_input=GO_input_down,topNodes=topNodes)

  GO_input<-rbind(
    cbind(GO_input_up,l2fold_direction=rep('up',dim(GO_input_up)[1])),
    cbind(GO_input_down,l2fold_direction=rep('down',dim(GO_input_down)[1]))
  )
  GO_input$l2fold_direction<-factor(GO_input$l2fold_direction,levels=c('up','down'))
  
  gout<-ggplot(data=GO_input,aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(Fisher.classic)))+
    facet_grid(rows = vars(l2fold_direction),scales = "free_y")+
    scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
    geom_point()+
    labs(x='Significant/Annotated',y='',color='-log(pval)')+
    theme_bw()
  
  return(gout)
  
}


dotplot_enricher_wt_up_down<-function(GO_input_up,GO_input_down,topNodes_up=15,topNodes_down=5){
  if(!is.null(GO_input_up) & !is.null(GO_input_down)){
    GO_input<-rbind(
      cbind(GO_input_up[1:topNodes_up,],l2fold_direction=rep('up',dim(GO_input_up[1:topNodes_up,])[1])),
      cbind(GO_input_down[1:topNodes_down,],l2fold_direction=rep('down',dim(GO_input_down[1:topNodes_down,])[1]))
    )
    GO_input$l2fold_direction<-factor(GO_input$l2fold_direction,levels=c('up','down'))
    
    gout<-ggplot(data=GO_input,aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(qvalue)))+
      facet_grid(rows = vars(l2fold_direction),scales = "free_y")+
      scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
      geom_point()+
      labs(x='Significant/Annotated',y='',color='-log(qvalue)')+
      theme_bw()
    
  } else if(!is.null(GO_input_up) & is.null(GO_input_down)){
    GO_input=GO_input_up
    gout<-ggplot(data=GO_input[1:topNodes_up,],aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(qvalue)))+
      scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
      geom_point()+
      labs(x='Significant/Annotated',y='',color='-log(qvalue), up')+
      theme_bw()
  } else if(!is.null(GO_input_down) & is.null(GO_input_up)){
    GO_input=GO_input_down
    gout<-ggplot(data=GO_input[1:topNodes_down,],aes(x=Significant/Annotated,y=Term,size=Significant,color=-log(qvalue)))+
      scale_color_gradient(low = "yellow", high = "red", na.value = NA)+
      geom_point()+
      labs(x='Significant/Annotated',y='',color='-log(qvalue), down')+
      theme_bw()
  }


  

  
  return(gout)
  
}


