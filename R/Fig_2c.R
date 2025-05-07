#load source data####
source("~/ANA_SOURCE_SHORT.R", echo=FALSE)
load('~/Fig2d_time_L24E2.RData')

#select important GO terms####
terms<-Reduce(union,lapply(fglist,function(x){return(x$pathway[which(x$padj<0.1)])}))
datNES<-Reduce(cbind,lapply(fglist,function(x){return(x[terms,'NES'])}))
rownames(datNES)<-terms
colnames(datNES)<-names(fglist)
datNES[is.na(datNES)]<-0


#cluster based on NES#######
centers=3
set.seed(12300)
kout<-kmeans(datNES,centers=centers)
clu<-kout$cluster
names(clu)<-rownames(datNES)
table(clu)

listgenes<-list()
for(i in 1:centers){
  listgenes[[i]]<-names(clu)[which(clu==i)]
  
}
names(listgenes)<-paste('cluster',1:centers)
unlist(lapply(listgenes,length))

meandat<-NULL
for(i in 1:length(listgenes)){
  
  meandat<-rbind(meandat,colMeans(datNES[listgenes[[i]],]))
}
rownames(meandat)<-names(listgenes)
cellsize=10
phee<-Heatmap(MinMax(meandat,-2,2),
              col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              cluster_columns = FALSE,
              heatmap_legend_param = list( title = paste('NES',sep=''),direction = "horizontal"),
              width = ncol(meandat)*unit(cellsize, "mm"), 
              height = nrow(meandat)*unit(cellsize, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", meandat[i, j]), x, y, gp = gpar(fontsize = 10))
              })

phee<-ComplexHeatmap::draw(phee,
                           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(ComplexHeatmap::draw(phee))




#Fig2d heatmap########

inputmarkers<-Markers_M_GO_1
ternss=intersect(inputmarkers,rownames(datNES))
names(ternss)<-names(inputmarkers[match(ternss,inputmarkers)])
heat<-datNES[ternss,]
rownames(heat)<-names(ternss)

cellsize=10
phee<-Heatmap(MinMax(heat,-2,2),
              col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              cluster_columns = FALSE,
              cluster_rows=FALSE,
              heatmap_legend_param = list( title = paste('NES',sep=''),direction = "horizontal"),
              width = ncol(heat)*unit(cellsize, "mm"), 
              height = nrow(heat)*unit(cellsize, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", heat[i, j]), x, y, gp = gpar(fontsize = 10))
              })


phee<-ComplexHeatmap::draw(phee,
                           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(ComplexHeatmap::draw(phee))

