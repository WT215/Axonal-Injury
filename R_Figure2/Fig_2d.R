load('~/RData/PSEDUOBULK_MYELINHIGH.RData')
load('~/RData/MYELIN_PDX_LATE.RData')
load('~/RData/MYELIN_PDX_EARLY.RData')
RAW<-cbind(RAW_NSG_HIGH,RAW_PDX_HIGH_E,RAW_PDX_HIGH)
stats=log2(rowMeans(RAW))
RAW_SUB<-RAW[names(stats)[which(stats>0)],]
colData_SUB<-data.frame(
  sample=colnames(RAW_SUB),
  stage=rep(c('NSG','Early','Late'),c(2,dim(RAW_PDX_HIGH_E)[2],dim(RAW_PDX_HIGH)[2]))
)
rownames(colData_SUB)<-colData_SUB$sample





load('/home/clustor2/ma/w/wt215/PROJECT_ST/FIGURES/Fig2d_time_L24E2.RData')
normcounts<-counts(DE_STAGE,norm=TRUE)


#select cell type marker genes from Ximerakis study
inputmarkers<-Marker_ximer[c('ASC','DC','MG','MAC','EC','Monocytes','OPC','OLG','Neurons')] #remove NEUT 03/02/2025


mdat<-NULL
for(i in 1:length(inputmarkers)){
  subdat<-normcounts[intersect(inputmarkers[[i]],rownames(normcounts)),]
  mdat<-rbind(mdat,apply(subdat,2,median))
}
rownames(mdat)<-names(inputmarkers)


pdat<-reshape2::melt(mdat)
pdat['stage']<-colData_SUB[pdat$Var2,'stage']


#run AUCell
AUC<-AUCell_run_wt(exprMatrix = normcounts,geneSets=inputmarkers)$score_mat
AUCZ<-as.data.frame(apply(AUC,2,scale))
rownames(AUCZ)<-rownames(AUC)



qdat<-AUCZ
qdat['stage']<-colData_SUB[rownames(qdat),'stage']
qdat2<-reshape2::melt(qdat)
qdat2$stage<-factor(qdat2$stage,levels=c('NSG','Early','Late'))
qdat2$variable<-factor(qdat2$variable,levels=rev(c('Neurons','ASC','OLG','OPC','Monocytes','DC','MAC','MG','EC')))


qdat3<-as.data.frame(
  qdat2 %>% dplyr::group_by(stage,variable)
  %>% dplyr::summarise(median=median(value),
                       sd=sd(value))
)

head(qdat3)
qdat3$stage<-factor(qdat3$stage,levels=c('NSG','Early','Late'))


hh<-reshape2::dcast(qdat3[,c('stage','variable','median')],variable ~stage)
rownames(hh)<-hh[,1]
hh<-hh[rev(levels(qdat2$variable)),-1]

rownames(hh)<-c("Neurons",'Astro','Oligo','OPC','Monocyte','DC','Macrophage','Microglia','EC')


cellsize=10
phee<-Heatmap(MinMax(hh,-0.5,0.5),
              col=colorRamp2(c( -0.5,0, 0.5), c('blue',"white", "red")),
              cluster_columns = FALSE,
              cluster_rows=FALSE,
              heatmap_legend_param = list( title = paste('median of gene signature',sep=''),direction = "horizontal"),
              width = ncol(hh)*unit(cellsize, "mm"), 
              height = nrow(hh)*unit(cellsize, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", hh[i, j]), x, y, gp = gpar(fontsize = 10))
              })


phee<-ComplexHeatmap::draw(phee,
                           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(ComplexHeatmap::draw(phee))
