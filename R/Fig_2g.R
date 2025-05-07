source("~/ANA_SOURCE_SHORT.R", echo=FALSE)
#HEATMAP Fig2g##########
#quantified nuclei area from H&E images
load('~/RData/Fig2g_RaviST/LIST_RAVI_NUCLEI.RData')
#binning of Ravi spots based on tumor density inferred from nuclei area
load('~/RData/Fig2g_RaviST/LIST_RAVI_SPOTS.RData')
#AUCell data
load('~/RData/Fig2g_RaviST/LIST_RAVI_AUC.RData')

#number of ROI spots in each bin:
load('~/RData/Fig2g_RaviST/LIST_TABLE.RData')

#summarize AUCell score in each tumor bin
load('~/LIST_AUC.RData')

#Define region of interests for analysis
VEC_ROI<-c('Infiltrative','Cellular','Necrotic_Edge','Vascular_Hyper','Necrosis')

##heat proportion########
meanmat =Reduce("+", lapply(LIST_TABLE, function(x) replace(x, is.na(x), 0)))
meanmat <-meanmat /rowSums(meanmat )



COL_RAVI<-c('black','#E41A1C','lightgrey','#377EB8',"#FCCDE5")
names(COL_RAVI)<-VEC_ROI

ha_ravi = HeatmapAnnotation(`Proportion of spots in bins` = anno_lines((meanmat), pch = 16, gp = gpar(col = COL_RAVI),pt_gp = gpar(col = COL_RAVI),add_points = TRUE),height = unit(4.5, "cm"),show_legend = TRUE)
plot(ha_ravi)

lgd_anotomical = Legend(labels = names(COL_RAVI), title = "Proportion of spots in bins", legend_gp = gpar(fill = COL_RAVI),
                        title_position = "lefttop")





##Fig2g based on Ravi ST data########
LIST_AUC2=lapply(LIST_AUC,function(x){return(apply(x,2,scale))})


meanauc <- Reduce("+", lapply(LIST_AUC2, function(x) replace(x, is.na(x), 0)))/length(LIST_AUC)

rownames(meanauc)<-names(LIST_RAVI_SPOTS)

heat<-t(meanauc)

cellsize=5
colnames(heat)<-str_sub(unlist(lapply(str_split(colnames(heat),pattern=','),function(x){return(x[2])})),0,-2)


phee=ha_ravi %v%
  Heatmap(MinMax(heat,-0.5,0.5),
          #Heatmap(heat,
          col=magma(10),
          #col=colorRamp2(c(-0.5,0,0.5), c( 'darkgreen',"white", "darkred")),
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          heatmap_legend_param = list( title = paste('Average expression',sep=''),direction = "horizontal"),
          width = ncol(heat)*unit(cellsize, "mm"), 
          height = nrow(heat)*unit(cellsize, "mm"))

phee<-draw(phee,
           heatmap_legend_side = "top",
           annotation_legend_side = "top",
           annotation_legend_list = list(lgd_anotomical))


gout_heat<-grid.grabExpr(draw(phee))
