
#load data######

#summarized  AUCell score of GO terms in each one of twenty tumor bins
load('~/Axonal-Injury/RData/FigS2h-o/LIST_MEDIAN.RData')

#summarized  AUCell score of GO terms in each one of twenty tumor bins (control data, pseudo spot, see Method)
load('~/Axonal-Injury/RData/FigS2h-o/LIST_MEDIAN_CONTROL.RData')


#VEC GO in paper####
VEC_GO<-c(
  `Axon injury`=	"GOBP_NEURON_PROJECTION_REGENERATION",
  `Immune response`=	'GOBP_ACTIVATION_OF_IMMUNE_RESPONSE',
  `Cell death`	='GOBP_CELL_KILLING',
  `Wound healing`=	'GOBP_REGENERATION',
  `Inflammation`	='GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
  `Angiogenesis`	='GOBP_BLOOD_VESSEL_MORPHOGENESIS',
  `Neuronal activity`	='GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC' ,
  `Cell signalling`	='GOBP_REGULATION_OF_SIGNALING_RECEPTOR_ACTIVITY'
)


#FigS2: line charts##########
SAMPLES<-PICK_SAMPLES
#SAMPLES=c(PICK_SAMPLES,PICK_SAMPLES_INV)


COL_ST<-brewer.pal(length(SAMPLES),name = 'Set1')
names(COL_ST)<-SAMPLES
COL_ST[which(COL_ST=='#FFFF33')]='grey'


#just need the information of the 20 bins names: "(0,0.05]"   "(0.05,0.1]",...
list_spots_ratio<- readRDS(paste('~/utils/LIST_SPOT_RATIO/','GCGR-L15','.rds',sep=''))


#for each GO term, make the line chart
for(i in 1:length(VEC_GO)){
  B=VEC_GO[i]
  
  cdat<-do.call(rbind,lapply(LIST_MEDIAN_CONTROL[SAMPLES],function(x){
    if(length(intersect(B,colnames(x)))==1)
      return(x[,B])
  }))
  
  if(is.null(dim(cdat))){
    cdat2=NULL
  }else if(dim(cdat)[1]>=2){
    cdat2<-cbind(apply(cdat,2,mean,na.rm=TRUE),
                 apply(cdat,2,sd,na.rm=TRUE))
    colnames(cdat2)<-c('c_mean','c_sd')
    cdat2<-as.data.frame(cdat2)
    cdat2['density']<-rownames(cdat2)
    cdat2['c_mean2']<-cdat2$c_mean
    cdat2$c_mean2<-cdat2$c_mean2/cdat2$c_mean[!is.na(cdat2$c_mean)][1]
  }
  
  

  dat<-do.call(rbind,lapply(LIST_MEDIAN,function(x){return(x[,B])})) 
  
  dat<-dat/apply(dat,1,function(x){x[!is.na(x) & x>0][1]})
  
  pdat<-reshape2::melt(dat)
  head(pdat)
  
  pdat$Var2<-factor(pdat$Var2,levels=names(list_spots_ratio))
  
  vecy<-NULL
  for(jj in 1:length(SAMPLES)){
    pdatsub<-pdat[which(pdat$Var1==SAMPLES[jj]),]
    vecy<-c(vecy,pdatsub$value[max(which(!is.na(pdatsub$value)))])
  }
  
  labeldat<-data.frame(
    Var1=SAMPLES,
    Var2=factor('(0.95,1]'),
    value=vecy
  )
  
  
  
  p<-ggplot()+
    geom_point(data=pdat[which(pdat$Var1 %in% SAMPLES),],aes(x=Var2,y=value,group=Var1,color=Var1))+
    geom_line(data=pdat[which(pdat$Var1 %in% SAMPLES),],aes(x=Var2,y=value,group=Var1,color=Var1))+
    labs(x='',y='gene signature',title=names(VEC_GO)[which(VEC_GO==B)])+
    scale_color_manual(values=COL_ST,limits=names(COL_ST))+
    geom_text_repel(data=labeldat,aes(x=Var2,y=value,label=Var1,color=Var1))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = 'none')
  
  if(!is.null(cdat2)){
    p=p+
      geom_point(data=cdat2,aes(x=density,y=c_mean2))+
      geom_line(data=cdat2,aes(x=density,y=c_mean2,group=1))+
      geom_errorbar(data=cdat2,aes(x=density,ymin=c_mean2-c_sd, ymax=c_mean2+c_sd), width=.2,
                    position=position_dodge(0.05))
  }
  
  
  p

  
}




