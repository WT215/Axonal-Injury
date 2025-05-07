load('~/RData/PSEDUOBULK_MYELINHIGH.RData')
load('~/RData/MYELIN_PDX_LATE.RData')
load('~/RData/MYELIN_PDX_EARLY.RData')

#all genes
allgenes<-Reduce(intersect,list(
  results_NSG$gene,results_PDX$gene,results_PDX_E$gene
))

#Select DE genes

#NSG
de_nsg<-results_NSG$gene[which(abs(results_NSG$log2FoldChange)>=0.5 & results_NSG$padj<0.05)]

#late
de_pdx<-results_PDX$gene[which(abs(results_PDX$log2FoldChange)>=0.5 & results_PDX$padj<0.05)]

#early
de_pdx_e<-results_PDX_E$gene[which(abs(results_PDX_E$log2FoldChange)>=0.5 & results_PDX_E$padj<0.05)]


degenes<-Reduce(union,list(de_nsg,de_pdx,de_pdx_e))
length(degenes)

logFCdat<-cbind(
  results_NSG[degenes,'log2FoldChange'],
  results_PDX_E[degenes,'log2FoldChange'],
  results_PDX[degenes,'log2FoldChange']
)
rownames(logFCdat)<-degenes
colnames(logFCdat)<-c('NSG','Early','Late')
logFCdat[is.na(logFCdat)]<-0


#cluster genes based on logFC######
centers=5
set.seed(12300)
kout<-kmeans(logFCdat,centers = centers)
unnii<-sort(unique(kout$cluster))



hdat<-NULL
for(j in 1:length(unnii)){
  gg<-names(kout$cluster)[which(kout$cluster==unnii[j])]
  #hdat<-rbind(hdat,apply(logFCdat[gg,],2,median))
  hdat<-rbind(hdat,colMeans(logFCdat[gg,],na.rm=TRUE))
}
rownames(hdat)<-unnii
colnames(hdat)<-c('NSG','Early','Late')

cellsize=10
phee<-Heatmap(MinMax(hdat,-2,2),
              cluster_columns = FALSE,
              col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              heatmap_legend_param = list( title = paste('log2FC',sep=''),direction = "horizontal"),
              width = ncol(hdat)*unit(cellsize, "mm"), 
              height = nrow(hdat)*unit(cellsize, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", hdat[i, j]), x, y, gp = gpar(fontsize = 10))
              })

phee<-draw(phee,
           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(draw(phee))


#enrichment analysis######
enrich_list<-list()
for(i in 1:centers){
  gg<-names(kout$cluster)[which(kout$cluster==i)]
  
  length(gg)
  enrich<-ENRICHER_WT(allgenes=allgenes,topgenes=gg,org=org.Mm.eg.db,DROP_GO_PATTERN = DROP_GO_PATTERN)
  #enrich<-enrich[enrich$p.adjust<0.1,]
  enrich_list[[i]]<-enrich
}
names(enrich_list)<-paste('cluster',1:centers)

datsup_fig2b<-NULL
for(i in 1:length(enrich_list)){
  dat<-enrich_list[[i]]
  dat['cluster']<-names(enrich_list)[i]
  datsup_fig2b<-rbind(datsup_fig2b,dat)
}

data<-NULL
for(i in 1:length(enrich_list)){
  tempp<-enrich_list[[i]]
  tempp['cluster']<-names(enrich_list)[i]
  data<-rbind(data,tempp)
}



enrich<-ENRICHER_WT(allgenes=allgenes,topgenes=gg,org=org.Mm.eg.db,DROP_GO_PATTERN = DROP_GO_PATTERN)

vec_terms<-c('angiogenesis','immune response','cell adhesion','wound healing','extracellular matrix organization','cell migration')




list_genes<-list()
for(i in 1:length(vec_terms)){
  list_genes[[i]]<-str_split(enrich[vec_terms[i],'geneID'],pattern='/')[[1]]
}
names(list_genes)<-vec_terms

heat_gene<-NULL
for(i in 1:length(list_genes)){
  gg<-list_genes[[i]]
  
  intersect(de_nsg,gg)
  heat_gene<-rbind(heat_gene,colMeans(logFCdat[intersect(gg,rownames(logFCdat)),]))
}
rownames(heat_gene)<-names(list_genes)
summary(as.vector(heat_gene))
phee_sub<-Heatmap(MinMax(heat_gene,0.1,0.7),
                  #cluster_rows = FALSE,
                  col=colorRamp2(c(0.1,0.7), c("white", "red")),
                  heatmap_legend_param = list( title = paste('cluster 2: log2FC',sep=''),direction = "horizontal"),
                  cluster_columns = FALSE,
                  width = ncol(heat_gene)*unit(cellsize, "mm"), 
                  height = nrow(heat_gene)*unit(cellsize, "mm"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", heat_gene[i, j]), x, y, gp = gpar(fontsize = 10))
                  },
                  show_row_dend = FALSE)

phee_sub<-draw(phee_sub,
               heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat_sub<-grid.grabExpr(draw(phee_sub))


gall<-ggarrange(plotlist = list(gout_heat,gout_heat_sub),nrow=1,ncol=2)
gall