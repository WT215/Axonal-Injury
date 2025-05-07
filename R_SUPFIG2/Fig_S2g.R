#load library#########
source("~/utils/ANA_SOURCE_SHORT.R", echo=FALSE)


#Figure S2g#####
#each sample: GO terms by c(S, pval, padj) matrix
load('/home/clustor2/ma/w/wt215/PROJECT_ST/AUC_GOBP/LIST_MK_M.RData')


keepgo<-Reduce(union,lapply(LIST_MK,function(x){rownames(x)[which(x$pval<0.1)]}))

#matrix of S values
Sheat<-Reduce(cbind,lapply(LIST_MK,function(x){return(x[keepgo,'S'])}))
colnames(Sheat)<-names(LIST_MK)
rownames(Sheat)<-keepgo

Heat<-t(Sheat)

set.seed(12300)
cluu<-kmeans(Sheat,centers=8)
unii<-sort(unique(cluu$cluster))

CLUSTER_GO<-list()
for(i in 1:8){
  CLUSTER_GO[[i]]<-names(cluu$cluster)[which(cluu$cluster==i)]
}
names(CLUSTER_GO)<-paste('cluster',seq(1,8))



Heat<-NULL
nn<-NULL
num_split<-NULL
list_cluster<-list()
for(i in 1:length(unii)){
  gg<-names(cluu$cluster)[which(cluu$cluster==unii[i])]
  Heat<-cbind(Heat,t(Sheat[gg,]))
  nn<-c(nn,gg)
  num_split<-c(num_split,rep(unii[i],length(gg)))
  
  rr<-rowMeans(Sheat[gg,])
  
  
  tempp<-as.data.frame(Sheat[names(sort(rr,decreasing=TRUE)),])
  tempp['GO']<-gg
  list_cluster[[i]]<-tempp
  
}
names(list_cluster)<-paste('cluster',unii)


qq<-NULL
for(i in 1:length(list_cluster)){
  tempp<-list_cluster[[i]]
  tempp['GO']<-rownames(tempp)
  tempp['cluster']<-names(list_cluster)[i]
  qq<-rbind(qq,tempp)
}
TABLES5=list(`Supplementary Table 5`=qq)


COL<-SpatialColors(length(unii))


topanno<-HeatmapAnnotation(
  clusters=anno_block(gp = gpar(fill = COL ),
                      labels =unii, 
                      labels_gp = gpar(col = "black", fontsize = 8)),
  which='column',
  simple_anno_size = unit(5, "cm"))

phee<-Heatmap(Heat,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_column_names = FALSE,
              show_column_dend = FALSE,
              column_split = num_split,
              top_annotation = topanno,
              heatmap_legend_param = list( title = paste('Mann−Kendall test',sep=''),direction = "horizontal"),
              #column_km = 8,
              column_title=NULL)


phee<-draw(phee,
           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(draw(phee))
