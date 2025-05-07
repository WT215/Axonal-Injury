#load python environment####
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
library(rjson)
library(rhdf5)
library(readr)
library(reticulate)
#use_virtualenv('virtual_R42', required = TRUE)

#load data#####
patht="/home/clustor2/ma/w/wt215/RFILES/"


source(paste(patht,"utils_bayNorm.R",sep=''), echo=FALSE)
source(paste(patht,"utils_BTC.R",sep=''), echo=FALSE)
source(paste(patht,"deseq2_sup.R",sep=''), echo=FALSE)
source(paste(patht,"utils_plot.R",sep=''), echo=FALSE)

#load markers#######
source("/home/clustor2/ma/w/wt215/PROJECT_ST/R/LOAD_MARKERS.R", echo=TRUE)

load('/home/clustor2/ma/w/wt215/PROJECT_ZAN/GEO_SC_OLD/RData/Markers_Zan_OLD.RData')


#LOAD ST related #######
DROP_SAMPLES=c('GCGR-L3','GCGR-L22','G144')
PICK_SAMPLES=c("GCGR-L15", "GL23" ,    "GCGR-E43" ,"GCGR-L21", "GCGR-L24", "GCGR-L5"  ,"GCGR-L1" ,'GCGR-E20')
PICK_SAMPLES_INV<-c('GCGR-L22','GBM2','GCGR-L24E','GL23E','GCGR-E43E','GCGR-L5E')
names(PICK_SAMPLES_INV)<-c("NH18-742",'GBM2','GCGR-L24E','GL23E','GCGR-E43E','GCGR-L5E')




##load M tissue annotation based on NSG##########
#load('/home/clustor2/ma/w/wt215/BTC_project/BTC_R/NSG/RData/ANNOTATION_NSG.RData')
load('/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/Rdata_pre/ANNOTATION_NSG.RData')
ANNOTATION_NSG_SHORT<-as.character(ANNOTATION_NSG)
ANNOTATION_NSG_SHORT[which(ANNOTATION_NSG_SHORT %in% c('L1','L2-3','L4','L5','L6'))]='Cortex'
names(ANNOTATION_NSG_SHORT)<-names(ANNOTATION_NSG)



load("/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/Rdata_pre/Harmony_labels_prob_m_GCGR.RData")
ANNO_L24E_E43E_L5E_GL23E<-readRDS(paste("/home/clustor2/ma/w/wt215/PROJECT_ST/HAR_NSG/HarmonyOut/ANNO_L24E_E43E_L5E_GL23E",'.rds',sep=''))


Harmony_labels_prob_m=c(Harmony_labels_prob_m_GCGR,ANNO_L24E_E43E_L5E_GL23E) #ANNOTATION_NSG_SHORT has already been included in *_GCGR
#Harmony_labels_prob_m<-c(Harmony_labels_prob_m_GCGR,ANNOTATION_NSG_SHORT)
#grep(names(Harmony_labels_prob_m),pattern='742')





##IMAGE PATH from visium########
path_image<-'/home/clustor2/ma/w/wt215/PROJECT_ST/10X_OUTS_SUB'
files<-list.files(path_image)
PATH_IMAGES_LOWRES<-paste(path_image,'/',files,'/spatial/tissue_lowres_image.png',sep='')
names(PATH_IMAGES_LOWRES)<-unlist(lapply(str_split(files,pattern='_'),function(x){return(x[1])}))




##load ROI (LABEL_ROI.R)#########
load('/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/Rdata_pre/LABELS_3ROIS_GCGR.RData')
load('/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/Rdata_pre/LL_L3_L3T_G144.RData')#I did not do ROI for these 3 samples
load('/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/Rdata_pre/LL_L24E_E43E_L5E_GL23E.RData')

LABELS_3ROIS=c(LABELS_3ROIS_GCGR,LL_L3_L3T_G144,LL_L24E_E43E_L5E_GL23E)

#write.csv(LABELS_3ROIS,file='/home/clustor2/ma/w/wt215/PROJECT_ST/LABELS_3ROIS.csv')

#a<-unique(unlist(lapply(str_split(names(LABELS_3ROIS),pattern='[.]'),function(x){return(x[2])})))

species_vec<-c('H','M')


##LOAD MYELIN LABEL (LABEL_MYELIN.R)#####
load('/home/clustor2/ma/w/wt215/PROJECT_ST/R/LABEL_MYELIN.RData')


##define samples for loading#########
unique_samples<-c(
  # "NH16-2725"="GCGR-L3",
  # "G144"='G144',
  "NH17-2957"="GCGR-L15",
  
  "NH10-1123"="GL23",
  "G343"="GCGR-E43",
  
  "NH17-556"="GCGR-L21",
  "NH18-1181"='GCGR-L24',
  
  "NH17-427"="GCGR-L5",
  "NH16-1240"="GCGR-L1",
  "NH18-742"="GCGR-L22",
  
  #2023.12.06
  'GBM2'='GBM2',
  'GCGR-E20'='GCGR-E20',
  # 'GCGR-L19'='GCGR-L19',
  
  #"NH16-2725T"='GCGR-L3T',
  # "NH10-1123T"="GL23T",
  # 
  # 'GCGR-L5T'='GCGR-L5T',
  # 'GCGR-L21T'='GCGR-L21T',
  # 'GCGR-L24T1'='GCGR-L24T1',
  # 'GCGR-L24T2'='GCGR-L24T2',
  # 'GCGR-E43T'='GCGR-E43T',
  # 'GCGR-L1T'='GCGR-L1T',
  
  "NSG"='NSG'
)

ALLSAMPLES<-c(
  "NH16-2725"="GCGR-L3",
  "G144"='G144',
  "NH17-2957"="GCGR-L15",
  
  "NH10-1123"="GL23",
  "G343"="GCGR-E43",
  
  "NH17-556"="GCGR-L21",
  "NH18-1181"='GCGR-L24',
  
  "NH17-427"="GCGR-L5",
  "NH16-1240"="GCGR-L1",
  "NH18-742"="GCGR-L22",
  
  #2023.12.06
  'GBM2'='GBM2',
  'GCGR-E20'='GCGR-E20',
  # 'GCGR-L19'='GCGR-L19',#not use L19, no tumor
  
  "NH16-2725T"='GCGR-L3T',
  "NH10-1123T"="GL23T",
  
  'GCGR-L5T'='GCGR-L5T',
  'GCGR-L21T'='GCGR-L21T',
  'GCGR-L24T1'='GCGR-L24T1',
  'GCGR-L24T2'='GCGR-L24T2',
  'GCGR-E43T'='GCGR-E43T',
  'GCGR-L1T'='GCGR-L1T',
  
  #2024.11.22 4 early stage sampels, only S1 and S2
  'GCGR-L24E'='GCGR-L24E',
  'GL23E'='GL23E',
  'GCGR-E43E'='GCGR-E43E',
  'GCGR-L5E'='GCGR-L5E',
  
  
  "NSG"='NSG'
)

unique_samplesU<-c(
  "NH16-2725"='GCGR-L3',
  "NH10-1123"="GL23",
  "NH17-427"='GCGR-L5',
  "NH17-556"='GCGR-L21',
  "NH18-1181"='GCGR-L24',
  "NH18-1181"='GCGR-L24',
  
  'G343'='GCGR-E43',
  'NH16-1240'='GCGR-L1'
)

unique_samplesT<-c(
  "NH16-2725T"='GCGR-L3T',
  "NH10-1123T"="GL23T",
  'GCGR-L5T'='GCGR-L5T',
  'GCGR-L21T'='GCGR-L21T',
  'GCGR-L24T1'='GCGR-L24T1',
  'GCGR-L24T2'='GCGR-L24T2',
  
  'GCGR-E43T'='GCGR-E43T',
  'GCGR-L1T'='GCGR-L1T'
)


#just keep L24T1 for analysis
META_TREATED<-data.frame(
  U=unname(unique_samplesU)[-6],
  T=unname(unique_samplesT)[-6]
)

unique_samplesE<-c(
  'GCGR-L24E'='GCGR-L24E',
  'GL23E'='GL23E',
  'GCGR-E43E'='GCGR-E43E',
  'GCGR-L5E'='GCGR-L5E'
)

unique_samplesL<-c(
  'GCGR-L24'='GCGR-L24',
  'GL23'='GL23',
  'GCGR-E43'='GCGR-E43',
  'GCGR-L5'='GCGR-L5'
)

META_STAGE<-data.frame(
  Early=unname(unique_samplesE),
  Late=unname(unique_samplesL)
)

META_TIME<-rbind(
  c('GCGR-E43E','GCGR-E43','GCGR-E43T'),
  c('GCGR-L24E','GCGR-L24','GCGR-L24T1'),
  c('GL23E','GL23','GL23T'),
  c('GCGR-L5E','GCGR-L5','GCGR-L5T')
)
SUV_ST<-c(
  'GCGR-L24E'=75,
  'GL23E'=74,
  'GCGR-E43E'=35,
  'GCGR-L5E'=60
)


#gene locations######
# library(readr)
# genes_loc_H <- as.data.frame(read_csv("PROJECT_ST/GENE_LOCATION/genes_loc_H.csv"))
# genes_loc_M <- as.data.frame(read_csv("PROJECT_ST/GENE_LOCATION/genes_loc_M.csv"))
# genes_loc<-rbind(genes_loc_H,genes_loc_M)
# save(genes_loc,file=paste(patht,"genes_loc.RData",sep=''))


load(paste(patht,"genes_loc.RData",sep=''))
genes_loc_H<-genes_loc[grep(genes_loc$gene_id,pattern='GRCh38'),]

genes_loc_H$gene_name<-genes_loc_H$gene_name
genes_loc_H<-genes_loc_H[!duplicated(genes_loc_H$gene_name),]
rownames(genes_loc_H)<-genes_loc_H$gene_name
dim(genes_loc_H)

list_chr_H<-list()
for(i in 1:22){
  list_chr_H[[i]]<-rownames(genes_loc_H)[genes_loc_H$seqnames==paste('chr',i,sep='')]
}
names(list_chr_H)<-paste('chr',1:22,sep='')



genes_loc_M<-genes_loc[grep(genes_loc$gene_id,pattern='mm10'),]

# genes_loc_M$gene_name<-str_sub(genes_loc_M$gene_name,8)
# genes_loc_M<-genes_loc_M[!duplicated(genes_loc_M$gene_name),]
rownames(genes_loc_M)<-genes_loc_M$gene_name
# dim(genes_loc_M)

list_chr_M<-list()
for(i in seq(1,19)){
  print(i)
  list_chr_M[[i]]<-rownames(genes_loc_M)[genes_loc_M$seqnames==paste('chr',i,sep='')]
}
length(list_chr_M)

names(list_chr_M)<-paste('chr',seq(1,19),sep='')


#save(genes_loc_H,genes_loc_M,file='/home/clustor22/ma/w/wt215/BTC_project/genes_loc.RData')
#write_csv(genes_loc_H,file='/home/clustor22/ma/w/wt215/BTC_project/genes_loc_H.csv')
#write_csv(genes_loc_M,file='/home/clustor22/ma/w/wt215/BTC_project/genes_loc_M.csv')
# 
# 


#SOME COLORS##########
COL_DECON_XIMER<-c(
  'MAC'="#F781BF", 
  'MG'="#4DAF4A", 
  'DC'="#A65628",
  'Monocytes'="#999999",
  
  'ASC'="#E41A1C", 
  'Neurons'="#377EB8", 
  
  'OPC'="#984EA3", 
  'OLG'="#FFFF33", 
  
  'EC'="#FF7F00", 
  'Pericytes'='lightgrey',
  
  'Choroid plexus epithelial'='black',
  'Ependymocytes'='green'
)


COL_LUCY<-c(
  'NPC-like'= '#145fc8',
  'OPC-like'= '#875bc2',
  'iOPC-like' ='#74369b',
  'NCC-like' ='#5f6669',
  'AC-like' ='#db3f33',
  'Reactive AC-like'= '#ec8833',
  'MES-like'='#f1cb57',
  'Cycling PC-like'= '#a4a3a3'
)

COL_GBM3<-c(
  'NPC-like'= '#145fc8',
  'OPC-like'= '#875bc2',
  'iOPC-like' ='#74369b',
  'MES-Hyp' ='black',
  #'AC-like' ='#db3f33',
  'Reactive AC-like'= '#ec8833',
  'MES-like'='#f1cb57',
  'Cycling PC-like'= '#a4a3a3'
)
