#load library#########
source("~/ANA_SOURCE_SHORT.R", echo=FALSE)
#load AUC GOBP M#####

INPUT_SAMPLES=PICK_SAMPLES


LIST_GOBP_M=readRDS('/home/clustor2/ma/w/wt215/PROJECT_ST/AUC_GOBP/LIST_GOBP_M.rds')
LIST_GOBP_M_SUB<-LIST_GOBP_M[INPUT_SAMPLES]
LIST_GOBP_M_SUB<-lapply(LIST_GOBP_M_SUB,function(x){return(x[,Markers_M_GO_1])})




#AUC heatmap#######
#load list of heatmap data: list_heat, mean AUCell scores
load('~/Fig_2e.RData')


#mean except NA values
mean_matrix <- Reduce("+", lapply(list_heat, function(x) replace(x, is.na(x), 0))) / 
  Reduce("+", lapply(list_heat, function(x) !is.na(x)))

rownames(mean_matrix)<-rownames(list_heat[[1]])

heat<-t(apply(mean_matrix,1,scale))

rownames(heat)<-rownames(list_heat[[1]])
colnames(heat)<-unlist(lapply(str_split(colnames(list_heat[[1]]),pattern=','),function(x){
  str_sub(x[2],0,-2)
}))





phee=Heatmap(MinMax(heat,-2,2),
             col=magma(10),
             cluster_columns = FALSE,
             cluster_rows=FALSE,
             row_split = c(rep(1,6),rep(2,2)),
             row_title = NULL,
             heatmap_legend_param = list( title = paste('Average expression',sep=''),direction = "horizontal"))


phee<-draw(phee,
           heatmap_legend_side = "top",annotation_legend_side = "top")


gout_heat<-grid.grabExpr(draw(phee))


#PROP ANATOMICAL#######
REGIONS=c('CC','Cortex','septal nucleus','Hypothalamus','STR')
list_prop_anatomical<-list()
for(i in 1:length(INPUT_SAMPLES)){
  Sample=INPUT_SAMPLES[i]
  list_spots_ratio<-readRDS(paste('~/LIST_SPOT_RATIO/',Sample,'.rds',sep=''))
  
  
  LL<-Harmony_labels_prob_m[grep(names(Harmony_labels_prob_m),pattern=paste(Sample,"[.]",sep=''))]
  LL<-LL[grep(names(LL),pattern='S1.')]
  LL<-LL[which(LL %in% REGIONS)]
  
  vecc<-matrix(NA,nrow=length(REGIONS),ncol=length(list_spots_ratio))
  colnames(vecc)<-names(list_spots_ratio)
  rownames(vecc)<-REGIONS
  names(list_spots_ratio)
  for(indbin in 1:length(list_spots_ratio)){
    spots=intersect(list_spots_ratio[[indbin]],names(LL))
    
    ta<-table(LL[spots])
    ta<-ta/sum(ta)
    vecc[names(ta),indbin]<-ta
  }
  
  list_prop_anatomical[[i]]<-vecc
  
}
names(list_prop_anatomical)<-INPUT_SAMPLES


colSums(list_prop_anatomical$`GCGR-L15`,na.rm=TRUE)

mean_matrix_an <-Reduce("+", lapply(list_prop_anatomical, function(x) replace(x, is.na(x), 0)))/length(list_prop_anatomical)

rownames(mean_matrix_an)<-rownames(list_prop_anatomical[[1]])

plot(mean_matrix_an['CC',])
#plot(list_prop_anatomical[[3]]['CC',])






COL_ANOTAMICAL<-c('red','grey','lightgrey','black','skyblue')
names(COL_ANOTAMICAL)<-REGIONS

ha_anotomical = HeatmapAnnotation(`Prop of spots V2` = anno_lines(t(mean_matrix_an), 
                                                                  pch = 16, 
                                                                  gp = gpar(col = COL_ANOTAMICAL),
                                                                  pt_gp = gpar(col = COL_ANOTAMICAL),
                                                                  add_points = TRUE),
                                  height = unit(3, "cm"),
                                  show_legend = TRUE)

lgd_anotomical = Legend(labels = names(COL_ANOTAMICAL), title = "Prop of spots V2", legend_gp = gpar(fill = COL_ANOTAMICAL),
                        title_position = "lefttop")
#draw(ha_anotomical)







#PROP MYELIN#####
REGIONS_WM=c('Myelin low','Myelin high')
list_prop_wm<-list()
for(i in 1:length(INPUT_SAMPLES)){
  Sample=INPUT_SAMPLES[i]
  list_spots_ratio<-readRDS(paste('~/LIST_SPOT_RATIO/',Sample,'.rds',sep=''))
  
  
  LL<-LABEL_MYELIN[grep(names(LABEL_MYELIN),pattern=paste(Sample,"[.]",sep=''))]
  LL<-LL[grep(names(LL),pattern='S1.')]
  vecc<-matrix(NA,nrow=length(REGIONS_WM),ncol=length(list_spots_ratio))
  colnames(vecc)<-names(list_spots_ratio)
  rownames(vecc)<-REGIONS_WM
  names(list_spots_ratio)
  for(indbin in 1:length(list_spots_ratio)){
    spots=intersect(list_spots_ratio[[indbin]],names(LL))
    
    ta<-table(LL[spots])
    ta<-ta/sum(ta)
    vecc[names(ta),indbin]<-ta
  }
  
  list_prop_wm[[i]]<-vecc
  
}
names(list_prop_wm)<-INPUT_SAMPLES


mean_matrix_wm <-Reduce("+", lapply(list_prop_wm, function(x) replace(x, is.na(x), 0)))/length(list_prop_wm)

rownames(mean_matrix_wm)<-rownames(list_prop_wm[[1]])


COL_WM<-c('lightgrey','darkgreen')
names(COL_WM)<-REGIONS_WM

ha_wm = HeatmapAnnotation(`Prop of spots V1` = anno_lines(t(mean_matrix_wm), 
                                                          pch = 16, 
                                                          gp = gpar(col = COL_WM),
                                                          pt_gp = gpar(col = COL_WM),
                                                          add_points = TRUE),
                          height = unit(3, "cm"),
                          show_legend = TRUE)

lgd_wm = Legend(labels = names(COL_WM), title = "Prop of spots V1", legend_gp = gpar(fill = COL_WM),
                title_position = "lefttop")




#FINAL HEATMAP: FIG2e#########

cellsize=5


phee=ha_wm %v%
  ha_anotomical %v%
  Heatmap(MinMax(heat,-2,2),
          col=magma(10),
          #col=colorRamp2(c(-2,0,2), c( 'darkgreen',"white", "darkred")),
          cluster_columns = FALSE,
          cluster_rows=FALSE,
          row_split = c(rep(1,6),rep(2,2)),
          row_title = NULL,
          heatmap_legend_param = list( title = paste('Average expression',sep=''),direction = "horizontal"),
          width = ncol(heat)*unit(cellsize, "mm"), 
          height = nrow(heat)*unit(cellsize, "mm")
  )






phee<-draw(phee,
           heatmap_legend_side = "top",
           annotation_legend_side = "top",
           annotation_legend_list = list(lgd_wm,lgd_anotomical))


gout_heat<-grid.grabExpr(draw(phee))
