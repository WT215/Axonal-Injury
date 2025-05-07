#load deconvolution resultst

LIST_DECON<-readRDS('/home/clustor2/ma/w/wt215/Axonal-Injury/RData/LIST_DECON_GITHUB.rds')


#load myelin labels#####
load('/home/clustor2/ma/w/wt215/PROJECT_ST/R/LABEL_MYELIN.RData')
spots_myelinhigh<-names(LABEL_MYELIN)[which(LABEL_MYELIN=='Myelin high')]
spots_myelinlow<-names(LABEL_MYELIN)[which(LABEL_MYELIN=='Myelin low')]
species='M'

#cell type change as function of time#####
TYPES<-c('ASC','DC','MG','MAC','EC','Monocytes','OPC','OLG','Neurons')

SAMPLES<-c('NSG',PICK_SAMPLES_INV,PICK_SAMPLES)
dat<-NULL
for(i in 1:length(SAMPLES)){
  Sample=SAMPLES[i]
  decon<-LIST_DECON[[Sample]][,TYPES]
  decon<-decon[grep(rownames(decon),pattern='S1.|S2.'),]
  spots<-intersect(rownames(decon),spots_myelinhigh)
  deconsub<-decon[spots,]
  dat<-cbind(dat,colMeans(deconsub))
}
colnames(dat)<-SAMPLES

datmean<-cbind(
  dat[,'NSG'],
  rowMeans(dat[,PICK_SAMPLES_INV]),
  rowMeans(dat[,PICK_SAMPLES])
)
colnames(datmean)<-c('NSG','Early','Late')

datmeanz<-t(apply(datmean,1,scale))
colnames(datmeanz)<-colnames(datmean)


#choose samples#######
Samples=PICK_SAMPLES
LIST_DECON_SUB<-LIST_DECON[Samples]

NUMSPOTS_DAT<-NULL
LIST_MEANDAT<-list()
for(i in 1:length(LIST_DECON_SUB)){
  tempp<-LIST_DECON_SUB[[i]]
  tempp<-tempp[grep(rownames(tempp),pattern='S1.|S2.'),]
  tempp_norm<-tempp/rowSums(tempp)
  
  Sample<-names(LIST_DECON_SUB)[i]
  list_spots_ratio<-readRDS(paste('~/utils/LIST_SPOT_RATIO/',Sample,'.rds',sep=''))
  
  tempmat<-NULL
  numspots<-NULL
  
  for(indbin in 1:length(list_spots_ratio)){
    binspots<-list_spots_ratio[[indbin]]
    spots<-intersect(rownames(tempp_norm),binspots)
    LL<-LABELS_3ROIS[spots]
    spots_tumor<-names(LL)[which(LL %in% c('bulk','mar','inv'))]
    numspots<-c(numspots,length(spots_tumor))
    
    
    if(length(spots)>1){
      tempmat<-rbind(tempmat,colMeans(tempp_norm[spots,]))
    }else if(length(spots)==1){
      tempmat<-rbind(tempmat,tempp_norm[spots,])
    }else{
      tempmat<-rbind(tempmat,rep(NA,dim(tempp_norm)[2]))
    }
    
  }
  rownames(tempmat)<-names(list_spots_ratio)
  LIST_MEANDAT[[i]]<-t(tempmat)
  
  NUMSPOTS_DAT<-rbind(NUMSPOTS_DAT,numspots)
  
}
names(LIST_MEANDAT)<-names(LIST_DECON_SUB)

colnames(NUMSPOTS_DAT)<-colnames(LIST_MEANDAT$`GCGR-L15`)
rownames(NUMSPOTS_DAT)<-names(LIST_MEANDAT)


NUMSPOTS_DAT_NORM<-NUMSPOTS_DAT/rowSums(NUMSPOTS_DAT)
dat_num_late<-data.frame(
  mean=colMeans(NUMSPOTS_DAT_NORM),
  sd=apply(NUMSPOTS_DAT_NORM,2,sd)
)

COL_S<-brewer.pal(n=dim(NUMSPOTS_DAT)[1],name='Set1')
names(COL_S)<-rownames(NUMSPOTS_DAT)

ha_top = HeatmapAnnotation(`Prop of tumor spots` = anno_lines(t(NUMSPOTS_DAT_NORM), pch = 16, gp = gpar(col = COL_S),pt_gp = gpar(col = COL_S),add_points = TRUE),height = unit(3, "cm"),show_legend = TRUE)
lgd_top = Legend(labels = names(COL_S), title = "Prop of tumor spots", legend_gp = gpar(fill = COL_S),title_position = "lefttop")




mean_matrix <- Reduce("+", lapply(LIST_MEANDAT, function(x) replace(x, is.na(x), 0))) / 
  Reduce("+", lapply(LIST_MEANDAT, function(x) !is.na(x)))

heatz<-t(apply(mean_matrix,1,scale))
colnames(heatz)<-unlist(lapply(str_split(names(list_spots_ratio),pattern=','),function(x){
  str_sub(x[2],0,-2)
}))

cellsize=5
phee<-ha_top %v%
  Heatmap(MinMax(heatz,-2,2),
          col=magma(10),
          #col=colorRamp2(c(-2,0,2), c( 'darkgreen',"white", "darkred")),
          cluster_columns = FALSE,
          heatmap_legend_param = list( title = paste('Proportion of cell types',sep=''),direction = "horizontal"),
          width = ncol(heatz)*unit(cellsize, "mm"), 
          height = nrow(heatz)*unit(cellsize, "mm"))



phee<-draw(phee,
           heatmap_legend_side = "top",annotation_legend_side = "top",
           annotation_legend_list = list(lgd_top))




gout_heat<-grid.grabExpr(draw(phee))






