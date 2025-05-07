#genenames<-rownames(mydata_sub)
#geeout<-H_M_names(genenames=genenames)
#source('/home/clustor2/ma/w/wt215/BTC_project/BTC_R/deseq2_sup.R', echo=FALSE)
#source('D:/Miriam_project/deseq2_sup.R', echo=FALSE)
source('/home/clustor2/ma/w/wt215/RFILES/deseq2_sup.R', echo=FALSE)

# library(One2One)
# # Download and format the homolog data from MGI
# allHomologs = load.homologs()
# species1="human"
# species2="mouse"
# ortholog_data = analyse.orthology(species1,species2,allHomologs)
# #print(ortholog_data$orthologs_one2one[1:10,])
# 
# orthologs_one2one=ortholog_data$orthologs_one2one
# 
# 
# 
# Map_H_to_M<-function(gh){
#   M=orthologs_one2one$mouse.symbol[which(orthologs_one2one$human.symbol %in% gh )]
#   
#   H=orthologs_one2one$human.symbol[which(orthologs_one2one$human.symbol %in% gh )]
#   return(cbind(H,M))
# }
# 
# Map_M_to_H<-function(gm){
#   M=orthologs_one2one$mouse.symbol[which(orthologs_one2one$mouse.symbol %in% gm )]
#   H=orthologs_one2one$human.symbol[which(orthologs_one2one$mouse.symbol %in% gm )]
#   return(cbind(H,M))
# }

#convert data for MA plot######
MAplot_convert<-function(R,G){
  M=log2(R/G)
  A=1/2*log2(R*G)
  d=data.frame(M=M,A=A)
  return(d)
}





#RUN AUCell across cell lines#####
AUCfun_baywsas<-function(species="M",OUTDIR=NULL,SAMPLES,inputmarkers,typee='bay'){
  AUCOUT_LIST<-list()

  
  for(ind_sample in 1:length(SAMPLES)){
    print(ind_sample)
    if(typee=='bay'){
      bayout=readRDS(paste('/home/clustor2/ma/w/wt215/PROJECT_ST/bayNorm/bay_',species,'_',SAMPLES[ind_sample],'.rds',sep=''))
      inputdat<-bayout$Bay_out
    }else if(typee=='rawTC'){
      aa<-mget(ls(.GlobalEnv,pattern = paste("raw_",species,'_',names(SAMPLES)[ind_sample],"_",sep='')), envir = .GlobalEnv)
      #print(length(aa))
      inputdat<-do.call(cbind,aa)
      ll<-LABELS_3ROIS[colnames(inputdat)]

      if(species=='M' & unique_samples[ind_sample]!='NSG'){
        inputdat<-inputdat[,names(ll)[-which(ll %in% c('necrotic/lowq'))]]
      }else if(species=='H'){
        inputdat<-inputdat[,names(ll)[-which(ll %in% c('normal','necrotic/lowq'))]]
      }
      inputdat<-t(t(inputdat)/colSums(inputdat))

      
    }
    
    
    
    tempdat<-tryCatch(
      expr=AUCell_run_wt(inputdat,inputmarkers,  BPPARAM=MulticoreParam(20)),
      error=function(e){
        message('Caught an error!')
        print(e)
        return(list(score_mat=NA))
      }               )
    #AUCOUT<-rbind(AUCOUT,tempdat$score_mat)
    output=tempdat$score_mat
    
    if(!is.null(OUTDIR)){
      dir.create(file.path(OUTDIR), showWarnings = FALSE)
      saveRDS(output,file=paste(OUTDIR,SAMPLES[ind_sample],'.rds',sep=''))
    }
    AUCOUT_LIST[[ind_sample]]<-output
    
  }
  names(AUCOUT_LIST)<-SAMPLES
  return(AUCOUT_LIST)
}



#convert gene symbols (adapted from https://github.com/zhanghao-njmu/SCP)#####
library(biomaRt)
GeneConvert <- function(geneID, geneID_from_IDtype = "symbol", geneID_to_IDtype = "entrez_id",
                        species_from = "Homo_sapiens", species_to = NULL,
                        Ensembl_version = 103, try_times = 5, mirror = NULL) {
  if (require("httr", quietly = TRUE)) {
    httr::set_config(httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))
  }
  
  if (missing(geneID)) {
    stop("geneID must be provided.")
  }
  if (is.null(species_to)) {
    species_to <- species_from
  }
  if (is.null(Ensembl_version)) {
    Ensembl_version <- "current_release"
  }
  species_from_split <- unlist(strsplit(species_from, split = "_"))
  species_to_split <- unlist(strsplit(species_to, split = "_"))
  species_from_simp <- paste0(tolower(substring(species_from_split[1], 1, 1)), species_from_split[2])
  species_to_simp <- paste0(tolower(substring(species_to_split[1], 1, 1)), species_to_split[2])
  geneID_from_IDtype <- sapply(geneID_from_IDtype, tolower)
  geneID_to_IDtype <- sapply(unique(geneID_to_IDtype), tolower)
  
  if ("symbol" %in% geneID_from_IDtype) {
    geneID_from_IDtype <- geneID_from_IDtype[geneID_from_IDtype != "symbol"]
    geneID_from_IDtype <- c("ensembl_symbol", "entrez_symbol", "uniprot_symbol", "wiki_symbol", geneID_from_IDtype)
  }
  from_IDtype <- sapply(geneID_from_IDtype, function(x) {
    switch(x,
           "ensembl_symbol" = "external_gene_name",
           "ensembl_id" = "ensembl_gene_id",
           "entrez_symbol" = "entrezgene_accession",
           "entrez_id" = "entrezgene_id",
           "uniprot_symbol" = "uniprot_gn_symbol",
           "wiki_symbol" = "wikigene_name"
    )
  })
  names(from_IDtype) <- geneID_from_IDtype
  
  to_IDtype <- sapply(geneID_to_IDtype, function(x) {
    switch(x,
           "symbol" = "external_gene_name",
           "ensembl_symbol" = "external_gene_name",
           "entrez_symbol" = "external_gene_name",
           "ensembl_id" = "ensembl_gene_id",
           "entrez_id" = "entrezgene_id"
    )
  })
  
  if (species_from != species_to && all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    to_IDtype <- sapply(to_IDtype, function(x) {
      switch(x,
             "external_gene_name" = "associated_gene_name",
             "ensembl_gene_id" = "ensembl_gene"
      )
    })
    to_attr <- paste(species_to_simp, to_IDtype, sep = "_homolog_")
    names(to_attr) <- geneID_to_IDtype
  } else {
    to_attr <- to_IDtype
    names(to_attr) <- geneID_to_IDtype
  }
  
  message("Connect to the Ensembl archives...")
  archives <- NULL
  ntry <- 0
  while (is.null(archives)) {
    ntry <- ntry + 1
    archives <- tryCatch(expr = {
      listEnsemblArchives()
    }, error = function(e) {
      message(e)
      message("Get errors when connecting with EnsemblArchives...\nRetrying...")
      Sys.sleep(1)
      return(NULL)
    })
    if (is.null(archives) && ntry >= try_times) {
      stop("Stop connecting...")
    }
  }
  
  Ensembl_version <- as.character(Ensembl_version)
  if (Ensembl_version == "current_release") {
    url <- archives[which(archives$current_release == "*"), "url"]
    version <- as.character(archives[which(archives$current_release == "*"), "version"])
    message("Using the ", Ensembl_version, "(", version, ")", " version of biomart...")
  } else if (Ensembl_version %in% archives$version) {
    url <- archives[which(archives$version == Ensembl_version), "url"]
    version <- as.character(archives[which(archives$version == Ensembl_version), "version"])
    message("Using the ", version, " version of biomart...")
  } else {
    stop("Ensembl_version is invalid. Must be one of current_release,", paste0(archives$version, collapse = ","))
  }
  
  message("Connect to the biomart...")
  mart <- NULL
  ntry <- 0
  
  while (is.null(mart)) {
    ntry <- ntry + 1
    if (!is.null(mirror)) {
      mart <- tryCatch(
        expr = {
          useEnsembl("ensembl", mirror = mirror)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with mirror...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
    }
    if (is.null(mart)) {
      mart <- tryCatch(
        expr = {
          useMart("ensembl", host = url)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with ensembl mart...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
    }
    if (is.null(mart) && ntry >= try_times) {
      stop("Stop connecting...")
    }
  }
  Datasets <- listDatasets(mart)
  
  dataset <- paste0(species_from_simp, "_gene_ensembl")
  if (!dataset %in% Datasets$dataset) {
    warning(paste0("Can not find the dataset for the species: ", species_from, " (", dataset, ")"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = NULL))
  }
  message("Connect to the dataset ", dataset, " ...")
  mart1 <- NULL
  ntry <- 0
  while (is.null(mart1)) {
    ntry <- ntry + 1
    mart1 <- tryCatch(
      expr = {
        useDataset(dataset = dataset, mart = mart)
      },
      error = function(e) {
        message(e)
        message("Get errors when connecting with ensembl mart...\nRetrying...")
        Sys.sleep(1)
        return(NULL)
      }
    )
    if (is.null(mart1) && ntry >= try_times) {
      stop("Stop connecting...")
    }
  }
  
  if (species_from != species_to && any(!geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
    dataset2 <- paste0(species_to_simp, "_gene_ensembl")
    if (!dataset2 %in% Datasets$dataset) {
      warning(paste0("Can not find the dataset for the species: ", species_from, " (", dataset2, ")"), immediate. = TRUE)
      return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = NULL))
    }
    message("Connect to the dataset ", dataset2, " ...")
    mart2 <- NULL
    ntry <- 0
    while (is.null(mart2)) {
      ntry <- ntry + 1
      mart2 <- tryCatch(
        expr = {
          useDataset(dataset = dataset2, mart = mart)
        },
        error = function(e) {
          message(e)
          message("Get errors when connecting with ensembl mart...\nRetrying...")
          Sys.sleep(1)
          return(NULL)
        }
      )
      if (is.null(mart2) && ntry >= try_times) {
        stop("Stop connecting...")
      }
    }
  }
  
  Attributes <- listAttributes(mart1)
  from_IDtype <- from_IDtype[from_IDtype %in% Attributes[, "name"]]
  geneID_res_list <- list()
  total <- length(geneID)
  geneID_res <- NULL
  if (any(!to_attr %in% Attributes$name)) {
    to_attr_drop <- to_attr[!to_attr %in% Attributes$name]
    to_attr <- to_attr[to_attr %in% Attributes$name]
    warning(paste0("Can not find the attributes for the species ", species_from, ": ", paste(to_attr_drop, collapse = ", ")), immediate. = TRUE)
    if (length(to_attr) == 0) {
      warning("No attribute found for the species ", species_from, ". Please check the 'Attributes' in the result.", immediate. = TRUE)
      return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Ensembl_version = version, Datasets = Datasets, Attributes = Attributes))
    }
  }
  
  message("Match and convert the geneID...\n")
  if (species_from != species_to) {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- getBM(
          mart = mart1,
          attributes = c(from_attr, "ensembl_gene_id"),
          filters = from_attr,
          values = list(geneID)
        )
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        colnames(geneID_res1) <- c("from_geneID", "ensembl_gene_id_tmp")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1[, "from_IDtype"] <- from_name
        
        if (all(geneID_to_IDtype %in% c("symbol", "ensembl_id"))) {
          geneID_res2 <- getBM(
            mart = mart1,
            attributes = unique(c("ensembl_gene_id", to_attr)),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id_tmp"])
          )
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          geneID_res2[, "ensembl_gene_id_tmp"] <- geneID_res2[, "ensembl_gene_id"]
          geneID_res2 <- geneID_res2[, c("ensembl_gene_id_tmp", to_attr)]
          geneID_res2 <- unique(reshape2::melt(geneID_res2, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
          geneID_res2$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res2$to_IDtype]
          geneID_res_merge <- merge(x = geneID_res1, y = geneID_res2, by = "ensembl_gene_id_tmp")
        } else {
          homolog_ensembl_gene <- paste(species_to_simp, "ensembl_gene", sep = "_homolog_")
          geneID_res2 <- getBM(
            mart = mart1,
            attributes = c("ensembl_gene_id", homolog_ensembl_gene),
            filters = "ensembl_gene_id",
            values = list(geneID_res1[, "ensembl_gene_id_tmp"])
          )
          geneID_res2 <- geneID_res2[geneID_res2[, "ensembl_gene_id"] %in% geneID_res1[, "ensembl_gene_id_tmp"], ]
          if (nrow(geneID_res2) == 0) {
            next
          }
          colnames(geneID_res2) <- c("ensembl_gene_id_tmp", homolog_ensembl_gene)
          
          geneID_res3 <- getBM(
            mart = mart2,
            attributes = unique(c("ensembl_gene_id", to_attr)),
            filters = "ensembl_gene_id",
            values = list(geneID_res2[, homolog_ensembl_gene])
          )
          geneID_res3 <- geneID_res3[geneID_res3[, "ensembl_gene_id"] %in% geneID_res2[, homolog_ensembl_gene], ]
          if (nrow(geneID_res3) == 0) {
            next
          }
          geneID_res3[, "ensembl_gene_id_tmp2"] <- geneID_res3[, "ensembl_gene_id"]
          geneID_res3 <- geneID_res3[, c("ensembl_gene_id_tmp2", to_attr)]
          geneID_res3 <- unique(reshape2::melt(geneID_res3, id.vars = "ensembl_gene_id_tmp2", variable.name = "to_IDtype", value.name = "to_geneID"))
          colnames(geneID_res3)[1] <- homolog_ensembl_gene
          geneID_res3$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res3$to_IDtype]
          geneID_res_merge <- Reduce(merge, list(geneID_res1, geneID_res2, geneID_res3))
        }
        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c("from_IDtype", "from_geneID", "to_IDtype", "to_geneID")]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  } else {
    for (from_attr in from_IDtype) {
      if (length(geneID) > 0) {
        geneID_res1 <- getBM(
          mart = mart1,
          attributes = unique(c("ensembl_gene_id", from_attr, to_attr)),
          filters = from_attr,
          values = list(geneID)
        )
        geneID_res1 <- geneID_res1[geneID_res1[, from_attr] %in% geneID, ]
        if (nrow(geneID_res1) == 0) {
          next
        }
        geneID_res1[, "ensembl_gene_id_tmp"] <- geneID_res1[, "ensembl_gene_id"]
        geneID_res1_from <- unique(geneID_res1[, c("ensembl_gene_id_tmp", from_attr)])
        colnames(geneID_res1_from) <- c("ensembl_gene_id_tmp", "from_geneID")
        from_name <- geneID_from_IDtype[which(from_IDtype == from_attr)]
        geneID_res1_from[, "from_IDtype"] <- from_name
        
        geneID_res1_to <- unique(geneID_res1[, c("ensembl_gene_id_tmp", to_attr)])
        geneID_res1_to <- unique(reshape2::melt(geneID_res1_to, id.vars = "ensembl_gene_id_tmp", variable.name = "to_IDtype", value.name = "to_geneID"))
        geneID_res1_to$to_IDtype <- setNames(names(to_attr), nm = to_attr)[geneID_res1_to$to_IDtype]
        
        geneID_res_merge <- merge(x = geneID_res1_from, y = geneID_res1_to, by = "ensembl_gene_id_tmp")
        
        geneID_res_merge[geneID_res_merge == ""] <- NA
        geneID_res_list[[from_attr]] <- geneID_res_merge[, c("from_IDtype", "from_geneID", "to_IDtype", "to_geneID")]
        ismap <- geneID %in% geneID_res_list[[from_attr]][, "from_geneID"]
        message(paste(sum(ismap), "genes mapped with", from_name))
        geneID <- geneID[!ismap]
      }
    }
  }
  message(
    paste0(
      paste0(rep("=", 30), collapse = ""), "\n",
      total - length(geneID), " genes mapped\n",
      length(geneID), " genes unmapped"
    ), "\n",
    paste0(rep("=", 30), collapse = ""), "\n"
  )
  geneID_res <- bind_rows(geneID_res_list) %>% unique()
  if (is.null(geneID_res) || nrow(geneID_res) == 0) {
    warning(paste0("No gene mapped"), immediate. = TRUE)
    return(list(geneID_res = NULL, geneID_collapse = NULL, geneID_expand = NULL, Datasets = Datasets, Attributes = Attributes))
  }
  geneID_collapse <- geneID_res %>%
    group_by(.data[["from_geneID"]], .data[["to_IDtype"]]) %>%
    mutate(
      from_geneID = unique(.data[["from_geneID"]]),
      to_IDtype = unique(.data[["to_IDtype"]]),
      to_geneID = list(unique(.data[["to_geneID"]][!.data[["to_geneID"]] %in% c("", NA)]))
    )
  geneID_collapse <- unique(as.data.frame(geneID_collapse[, c("from_geneID", "to_IDtype", "to_geneID")]))
  geneID_collapse <- geneID_collapse[sapply(geneID_collapse$to_geneID, length) > 0, ]
  geneID_collapse <- reshape2::dcast(geneID_collapse, formula = from_geneID ~ to_IDtype, value.var = "to_geneID")
  rownames(geneID_collapse) <- geneID_collapse[, "from_geneID"]
  geneID_expand <- unnest(data = geneID_collapse, cols = colnames(geneID_collapse)[sapply(geneID_collapse, class) == "list"])
  
  return(list(geneID_res = geneID_res, geneID_collapse = geneID_collapse, geneID_expand = geneID_expand, Ensembl_version = version, Datasets = Datasets, Attributes = Attributes, geneID_unmapped = geneID))
}



#read10xMolInfo .h5####
read10xMolInfo_wt<-function (sample, barcode.length = NULL, keep.unmapped = FALSE, 
                             get.cell = TRUE, get.umi = TRUE, get.gem = TRUE, get.gene = TRUE, 
                             get.reads = TRUE, get.library = TRUE, extract.library.info = FALSE, 
                             version = c("auto", "2", "3")) 
{
  sample <- path.expand(sample)
  version <- match.arg(version)
  if (version == "auto") {
    available <- h5ls(sample, recursive = FALSE)
    version <- if ("barcode_idx" %in% available$name) 
      "3"
    else "2"
  }
  data <- list()
  if (get.cell) {
    if (version == "3") {
      all.barcodes <- as.vector(h5read(sample, "/barcodes"))
      all.barcodes <- sub("-[0-9]+", "", all.barcodes)
      data$cell <- all.barcodes[as.vector(h5read(sample, 
                                                 "/barcode_idx")) + 1L]
    }
    else {
      data$cell <- get_cell_barcodes(sample, "barcode", 
                                     barcode.length)
    }
  }
  if (get.umi) {
    data$umi <- as.vector(h5read(sample, "/umi"))
  }
  if (get.gem) {
    data$gem_group <- as.vector(h5read(sample, "/gem_group"))
  }
  if (get.gene || !keep.unmapped) {
    path <- if (version == "3") 
      "/feature_idx"
    else "/gene"
    data$gene <- as.vector(h5read(sample, path)) + 1L
  }
  if (get.reads) {
    path <- if (version == "3") 
      "/count"
    else "/reads"
    data$reads <- as.vector(h5read(sample, path))
  }
  if (version == "3" && get.library) {
    data$library <- as.vector(h5read(sample, "/library_idx")) + 
      1L
  }
  if (length(data) == 0) {
    fhandle <- H5Fopen(sample)
    on.exit(H5Fclose(fhandle))
    dhandle <- H5Dopen(fhandle, "/umi")
    on.exit(H5Dclose(dhandle), add = TRUE, after = FALSE)
    space <- H5Dget_space(dhandle)
    on.exit(H5Sclose(space), add = TRUE, after = FALSE)
    N <- H5Sget_simple_extent_dims(space)
    stopifnot(N$rank == 1L)
    data <- make_zero_col_DFrame(N$size)
  }
  else {
    data <- do.call(DataFrame, data)
  }
  path <- if (version == "3") 
    "/features/name"
  else "/gene_ids"
  gene.ids <- as.vector(h5read(sample, path))
  if (!keep.unmapped) {
    keep <- data$gene <= length(gene.ids)
    if (!get.gene) {
      data$gene <- NULL
    }
    data <- data[keep, , drop = FALSE]
  }
  output <- list(data = data, genes = gene.ids)
  if (version == "3") {
    output$feature.type <- as.vector(h5read(sample, "/features/feature_type"))
    if (extract.library.info) {
      lib.info <- h5read(sample, "/library_info")
      output$library.info <- jsonlite::fromJSON(lib.info, 
                                                simplifyVector = FALSE)
    }
  }
  output
}

#function#####


H_M_names<-function(genenames){
  
  i1<-grep('GRCh38',genenames)
  s1<-str_sub(genenames[i1],8)
  i2<-grep('mm10',genenames)

  s2<-str_sub(genenames[i2],8)
  
  
  genenames_re<-genenames
  names(genenames_re)<-rep('H',length(genenames))
  genenames_re[i1]<-s1
  genenames_re[i2]<-s2
  names(genenames_re)[i2]<-rep('M',length(i2))
  
  outlist<-list(
    genes_h=s1,
    genes_m=s2,
    genenames_re=genenames_re
  )
  return(outlist)
}



GO_dots<-function(
  DATA,
  wout_list,
  MARKERS=NULL,
  HVG,
  target=0,
  num_markers=200,
  speciess='M',
  onts='BP'
  ){

  if(is.null(MARKERS)){
    MARKERS<-rownames(wout_list[[which(as.numeric(names(wout_list))==target)]])[1:num_markers]
  }

  #Background genes
  bc<-H_M_names(HVG)
  selggg_bc<-bc$genenames_re[which(names(bc$genenames_re)==speciess)]
  
  #Significant genes of background genes
  #selggg<-H_M_names(MARKERS)
  #selggg<-selggg$genenames_re[which(names(selggg$genenames_re)==speciess)]
  selggg<-MARKERS


  if(speciess=='M'){
    g_out<-GOout(
      DATA=DATA,
      GENES=bc$genenames_re,
      sigGenes=selggg,
      typee='mouse',
      AnnotationDb=org.Mm.eg.db,
      mapping="org.Mm.eg.db",
      topNodes=50
    )
  }else{
    g_out<-GOout(
      DATA=DATA,
      GENES=bc$genenames_re,
      sigGenes=selggg,
      typee='human',
      AnnotationDb=org.Hs.eg.db,
      mapping="org.Hs.eg.db",
      topNodes=50
    )
  }


  GO_out=g_out[onts][[1]]
  
  
  GO_out<-GO_out[!duplicated(GO_out$Term),]
  
  dots_out<-dotplot_wt(GO_input=GO_out,topNodes=30)+
    theme(text = element_text(size=20))
  
  return(
    list(dots_out=dots_out,
         GO_out=GO_out)
  )
}

#onts = c( "MF", "BP", "CC" )
##RcisTarget#########
# data(motifAnnotations_mgi)
# data(motifAnnotations_hgnc)
# 
# 
# RcisTarget_wt<-function(
#   geneLists,
#   typee='human',
#   top_motif=3
#   ){
#   if(typee=='human'){
#     # Search space: 500 bp around TSS - HUMAN
#     motifRankings <- importRankings("D:/Spatial_scRNAseq/RcisTarget/hg19-500bp-upstream-7species.mc9nr.feather")
# 
#     motifAnnot=motifAnnotations_hgnc
#   }else{
#     motifRankings <- importRankings("D:/Spatial_scRNAseq/RcisTarget/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
# 
#     motifAnnot=motifAnnotations_mgi
#   }
#   
#   
#   motifEnrichmentTable <- cisTarget(geneLists, 
#                                     motifRankings,
#                                     motifAnnot=motifAnnotations_mgi)
#   
#   #motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable)
#   
#   #resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]
#   
#   library(DT)
#   # datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
#   #           escape = FALSE, # To show the logo
#   #           filter="top", options=list(pageLength=5))
#   
# 
#   
#   signifMotifNames <- motifEnrichmentTable$motif[1:top_motif]
#   
#   incidenceMatrix <- getSignificantGenes(geneLists[[1]], 
#                                          motifRankings,
#                                          signifRankingNames=signifMotifNames,
#                                          plotCurve=FALSE, maxRank=5000, 
#                                          genesFormat="incidMatrix",
#                                          method="aprox")$incidMatrix
#   
#   return(list(
#     motifEnrichmentTable=motifEnrichmentTable,
#     incidenceMatrix=incidenceMatrix
#   ))
# }
# 


#KEGG pathway gene sets#########
KEGG_wt<-function(ENTRY=c("mmu04630")){
  outobject<-keggGet(ENTRY)
  genes_list<-list()
  for(i in 1:length(outobject)){
    names <- outobject[[i]]$GENE
    #Delete the gene number by deleting every other line
    namesodd <-  names[seq(0,length(names),2)]
    namestrue <- gsub("\\;.*","",namesodd)
    genes_list[[i]]<-namestrue
  }
  names(genes_list)<-names(ENTRY)
  return(genes_list)
}
#qq=KEGG_wt(c("mmu04630","mmu04659"))


#neighbours (function)########
# inds=19
# NOVASEQ_META$Sample_ID
# 
# inputind=4
# inputind=which(NOVASEQ_META$Sample=='NH18-742' & NOVASEQ_META$Index2==2)

neighbor_list_fun<-function(inds=1){
  Sample_ID=NOVASEQ_META$Sample_ID[inds]
  Sample=NOVASEQ_META$Sample[inds]
  Sample_INDEX=NOVASEQ_META$Index[inds]
  Sample_INDEX2=NOVASEQ_META$Index2[inds]
  
  GCGR_ID=NOVASEQ_META$GCCR_ID[inds]
  print(Sample_ID)
  
  labels_spots_sample<-get(paste("labels_spots_",Sample,'_',Sample_INDEX2,sep ='' ))
  
  ratio_vec<-get(paste("ratio_h_m_",Sample,'_',Sample_INDEX2,sep ='' ))
  total_h<-get(paste("total_h_",Sample,'_',Sample_INDEX2,sep ='' ))
  total_m<-get(paste("total_m_",Sample,'_',Sample_INDEX2,sep ='' ))
  s_dat<-get(paste("sdat_",Sample,'_',Sample_INDEX2,sep ='' ))
  
  coor<-data.frame(
    s_dat$array_col,
    s_dat$array_row
  )
  names_core<-names(labels_spots_sample[which(labels_spots_sample=='core')])
  a<-coor[names_core,]
  b<-colMeans(coor[names_core,])
  
  
  #centerpoint<-which.max(ratio_vec[names_core])
  
  
  
  if(dim(a)[1]!=0){
    
    distttomean<-sqrt(rowSums((a-b)^2))
    
    centerpoint<-which.min(distttomean)
    
    index_neighbours<-platform_neighbors_wt(coordinates=coor, platform='Visium',index0=FALSE)
    
    #find neighbours, exclude previous
    neighbor_list<-list()
    neighbor_list[[1]]<-names(centerpoint)
    qf=1
    while(length(neighbor_list[[qf]])>0){
      pick<-names(index_neighbours)[unlist(index_neighbours[neighbor_list[[qf]]])]
      #neighbor_list[[qf+1]]<-pick
      neighbor_list[[qf+1]]<-setdiff(pick,unlist(neighbor_list[1:qf]))
      qf<-qf+1
    }
    neighbor_list<-neighbor_list[-length(neighbor_list)]
    length(neighbor_list)
    return(list(
      neighbor_list=neighbor_list,
      sdat=s_dat,
      centerpoint=centerpoint,
      Sample_ID=Sample_ID,
      GCGR_ID=GCGR_ID,
      Sample=Sample,
      total_h=total_h,
      total_m=total_m
    ))
  } else{return(NULL)}
  
  
}

neighbor_list_fun_v2<-function(Sample,ind_slide,ind_center=1){
  #ind_center=1: pick the center of tumor, with minimum distance to centroid of bulk or mar or inv
  Sample_ID=paste(Sample,'_',ind_slide,sep='')
  
  GCGR_ID=paste(unique_samples[Sample],'_',ind_slide,sep='')
  print(Sample_ID)
  

  
  ratio_vec<-get(paste("ratio_h_m_",Sample,'_',ind_slide,sep ='' ))
  total_h<-get(paste("total_h_",Sample,'_',ind_slide,sep ='' ))
  total_m<-get(paste("total_m_",Sample,'_',ind_slide,sep ='' ))
  s_dat<-get(paste("sdat_",Sample,'_',ind_slide,sep ='' ))
  labels_spots_sample<-LABELS_3ROIS[colnames(s_dat)]
  
  coor<-data.frame(
    s_dat$array_col,
    s_dat$array_row
  )
  
  coor<-coor[names(labels_spots_sample)[which(labels_spots_sample!='necrotic/lowq')],]
  
  
  names_core<-names(labels_spots_sample[which(labels_spots_sample=='bulk')])
  a<-coor[names_core,]
  b<-colMeans(coor[names_core,])
  
  if(dim(a)[1]==0){
    names_core<-names(labels_spots_sample[which(labels_spots_sample=='mar')])
    a<-coor[names_core,]
    b<-colMeans(coor[names_core,])
    
    if(dim(a)[1]==0){
      names_core<-names(labels_spots_sample[which(labels_spots_sample=='inv')])
      a<-coor[names_core,]
      b<-colMeans(coor[names_core,])
    }
    
  }
  
  
  if(dim(a)[1]!=0){
    
    
    distttomean<-sort(sqrt(rowSums((a-b)^2)),decreasing=FALSE)
    
    centerpoint<-distttomean[ind_center]
    
    index_neighbours<-platform_neighbors_wt(coordinates=coor, platform='Visium',index0=FALSE)
    
    #find neighbours, exclude previous
    neighbor_list<-list()
    neighbor_list[[1]]<-names(centerpoint)
    qf=1
    while(length(neighbor_list[[qf]])>0){
      pick<-names(index_neighbours)[unlist(index_neighbours[neighbor_list[[qf]]])]
      #neighbor_list[[qf+1]]<-pick
      neighbor_list[[qf+1]]<-setdiff(pick,unlist(neighbor_list[1:qf]))
      qf<-qf+1
    }
    neighbor_list<-neighbor_list[-length(neighbor_list)]
    length(neighbor_list)
    return(list(
      index_neighbours=index_neighbours,
      neighbor_list=neighbor_list,
      sdat=s_dat,
      centerpoint=centerpoint,
      Sample_ID=Sample_ID,
      GCGR_ID=GCGR_ID,
      Sample=Sample,
      total_h=total_h,
      total_m=total_m
    ))
  } else{return(NULL)}
}


neighbourfun<-function(index_neighbours,index){
  nei=unique(unlist(index_neighbours[index]))
  return(nei)
}

layersfun<-function(index_neighbours,index,num_layers=5){
  layers_list<-list()
  for(j in 1:num_layers){
    if(j==1){
      layers_list[[j]]<-neighbourfun(index_neighbours,index)
    }
    else{
      layers_list[[j]]<-setdiff(neighbourfun(index_neighbours,layers_list[[j-1]]),index)
    }
    
  }
  names(layers_list)<-paste('layer:',seq(1,num_layers))
  
  layers_setdiff_list<-layers_list
  if(num_layers>1){
    for(j in num_layers:2){
      layers_setdiff_list[[j]]<-setdiff(layers_setdiff_list[[j]],c(Reduce(union,layers_setdiff_list[(j-1):1]),i))
    }
  }
  
  return(list(
    layers_list=layers_list,
    layers_setdiff_list=layers_setdiff_list
  ))
  
}



LOAD_ST_10Xfun_hm<-function(Sample,IND_SECTION,
                         path_sout="/home/clustor2/ma/w/wt215/BTC_project/BTC_R/Results_H_M/SEURAT_OBJECTS"){
  Sample_ID=paste(Sample,'-',IND_SECTION,sep='')
  

  seuratpath=paste(path_sout,'/',Sample_ID,'_',reftype,'.rds',sep='')
  
  assign(paste("sdat_",Sample,'_',IND_SECTION,sep ='' ),readRDS(seuratpath), envir = .GlobalEnv)
  sdat<-get(paste("sdat_",Sample,'_',IND_SECTION,sep ='' ))
  
  
  
  assign(paste("GENES_H_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$RNA@meta.features$H)
  assign(paste("GENES_M_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$RNA@meta.features$M)
  
  GENES_H<-get(paste("GENES_H_",Sample,'_',IND_SECTION,sep ='' ))
  GENES_M<-get(paste("GENES_M_",Sample,'_',IND_SECTION,sep ='' ))
  
  assign(paste("raw_H_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$RNA@counts[GENES_H,])
  assign(paste("raw_M_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$RNA@counts[GENES_M,])
  
  raw_H<-get(paste("raw_H_",Sample,'_',IND_SECTION,sep ='' ))
  raw_M<-get(paste("raw_M_",Sample,'_',IND_SECTION,sep ='' ))
  rownames(raw_H)<-str_sub(rownames(raw_H),8)
  rownames(raw_M)<-str_sub(rownames(raw_M),8)
  
  assign(paste("raw_H_",Sample,'_',IND_SECTION,sep ='' ),raw_H, envir = .GlobalEnv)
  assign(paste("raw_M_",Sample,'_',IND_SECTION,sep ='' ),raw_M, envir = .GlobalEnv)
  
  assign(paste("raw_mt_H_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$mt[grep(rownames(sdat@assays$mt),pattern='GRCh38'),], envir = .GlobalEnv)
  assign(paste("raw_mt_M_",Sample,'_',IND_SECTION,sep ='' ),sdat@assays$mt[grep(rownames(sdat@assays$mt),pattern='mm10'),], envir = .GlobalEnv)
  

  
  ratio_h_m<-colSums(raw_H)/colSums(raw_M)
  assign(paste("ratio_h_m_",Sample,'_',IND_SECTION,sep ='' ),ratio_h_m, envir = .GlobalEnv)
  assign(paste("total_h_",Sample,'_',IND_SECTION,sep ='' ),colSums(raw_H), envir = .GlobalEnv)
  assign(paste("total_m_",Sample,'_',IND_SECTION,sep ='' ),colSums(raw_M), envir = .GlobalEnv)

}

LOAD_ST_10Xfun_m<-function(Sample,IND_SECTION,
                            path_sout="/home/clustor2/ma/w/wt215/BTC_project/BTC_R/Results_H_M/SEURAT_OBJECTS"){
  Sample_ID=paste(Sample,'-',IND_SECTION,sep='')
  seuratpath=paste(path_sout,'/',Sample_ID,'_',reftype,'.rds',sep='')
  
  assign(paste("sdat_",Sample,'_',IND_SECTION,sep ='' ),readRDS(seuratpath), envir = .GlobalEnv)
  sdat<-get(paste("sdat_",Sample,'_',IND_SECTION,sep ='' ))
  
}
#load cell2loc output #########
loadcell2loc<-function(species='H',path,kr=vec_inds){
  A_q05_list<-list()
  # if(species=='H'){
  #   kr=(length(unique_samples)-1)
  # }else{
  #   kr=length(unique_samples)
  # }
  # 
  
  for(i in 1:length(kr)){
    ind_sample=kr[i]
    Sample<-names(unique_samples)[ind_sample]
    
    A_q05 <- as.data.frame(read_delim(paste(path,"A_q05_",Sample,'.csv',sep=''), delim = "\t", escape_double = FALSE, trim_ws = TRUE))
    
    rownames(A_q05)<-A_q05[,1]
    A_q05<-A_q05[,-1]
    colnames(A_q05)<-unlist(lapply(str_split(colnames(A_q05),pattern='_'),function(x){
      paste(x[-c(1,2,3,4)],collapse='_')
    }))
    A_q05_list[[i]]<-A_q05
    
    
  }
  names(A_q05_list)<-unique_samples[kr]
  return(A_q05_list)
}

#remove species prefix in gene name####
getHMgene=function(vec){
  temh<-str_sub(vec[grep(vec,pattern='GRCh38')],8)
  temm<-str_sub(vec[grep(vec,pattern='mm10')],8)
  return(list(
    Hgene=temh,
    Mgene=temm,
    Hind=grep(vec,pattern='GRCh38'),
    Mind=grep(vec,pattern='mm10')
  ))
}




#extract top 2 statistics########
#https://stackoverflow.com/questions/58504930/fastest-way-to-select-i-th-highest-value-from-row-and-assign-to-new-column/58505175#58505175
nth <- function(x, nth_largest=c(0,1)) {
  
  n=length(x) - (nth_largest - 1L)
  tempp<-c(sort(x,partial=n)[n],colnames(AUCZ)[order(x)[n]])
  return(tempp)
}

#my diy, use order once, maybe faster?
nth_v2<-function(x, nth_largest=c(0,1),namess) {
  
  n=length(x) - (nth_largest - 1L)
  
  orrr<-order(x)
  
  tempp<-c(x[orrr][n],namess[orrr[n]])
  return(tempp)
}


#for each row, extract the top 2 highest value and names
Extract_Top2<-function(AUCZ){
  qq=t(apply(AUCZ, 1, nth_v2, nth_largest = c(1,2),namess=colnames(AUCZ)))
  qq<-as.data.frame(qq)
  for(i in 1:2){
    qq[,i]<-as.numeric(qq[,i])
  }
  qq['FC']<-qq[,1]/qq[,2]
  colnames(qq)<-c('TOP1','TOP2','TYPE1','TYPE2','FC')
  return(qq)
}


#co-localization: adjacency: Tirosh 2024: https://github.com/tiroshlab/Spatial_Glioma/blob/f38f9b60ddbeefd9739e1483cbc66b2650f931d3/Module3_spatial_associations.R#L447########

#score_vecs: row: spots, column cell type, each column z-score norm

Tirosh_colocal_adj<-function(scorer_vecs,num_boot=200){
  message('score_vecs: row: spots, column cell type, each column z-score norm')
  
  scored_filt <- apply(score_vecs, 2, function(x) {
    if(length(x[x > 0]) == 0) {prop_vec <- setNames(rep(0, nrow(score_vecs)), rownames(score_vecs))}
    
    get_positives <- x[x > 0.1] #all scores >=0.1 - consider cell state present in that spot
    prop_vec <- setNames(rep(0, nrow(score_vecs)), rownames(score_vecs))
    prop_vec <- setNames(ifelse(rownames(score_vecs) %in% names(get_positives), yes =  get_positives[match(names(prop_vec), names( get_positives))], no = 0), rownames(score_vecs))
    return(prop_vec)
  })
  
  #convert to df
  new_decon_df <- as.data.frame(t(scored_filt)) %>% rownames_to_column(var = "barcodes")
  new_decon_mtrx <- new_decon_df %>%
    tibble::column_to_rownames("barcodes") %>%
    as.matrix()
  
  
  decon=t(new_decon_mtrx)
  obs_coloc <- colocal_mat((decon))
  
  
  #they shuffle number in score mat in each column
  spot_list <- lapply(c(1:num_boot),function(x){ #num. of shuffled matrices
    all_row <- sapply(c(1:dim(decon)[2]),function(c){
      new_row <- sample(c(1:length(row.names(decon))), length(row.names(decon)), replace = FALSE)
      return(new_row) 
    })
    new_mat <- sapply(c(1:dim(decon)[2]),function(c){
      sh_row <- decon[,c][all_row[,c]]#generate many deconv matrices by shuffling cell types
      return(sh_row)
    })
    colnames(new_mat) <- colnames(decon)
    return(new_mat)
  }) 
  
  #run colocal_mat_fun -- perform colocalization on the shuffled decon matrices
  mats_list<- lapply(spot_list, function(new_decon){
    colocal_mats<-colocal_mat((new_decon))
    return(colocal_mats)
  })
  mean_mat <- sapply(mats_list,function(x){ #expected colocalization calculated from shuffled decon matrices
    return(x[,"ab_comb"])
  })
  pairs_shuf_mean <- apply(mean_mat,1,mean)
  names(pairs_shuf_mean) <- mats_list[[1]][,"pairs"]
  rownames(mean_mat) <- mats_list[[1]][,"pairs"]
  
  
  r1 = obs_coloc$ab_comb #observed colocal
  r2 = pairs_shuf_mean #expected colocal
  
  n = dim(decon)[1] #num total spots
  
  fisher = ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n-3))+(1/(n-3)))^0.5
  
  p.value = (2*(1-pnorm(abs(as.matrix(fisher)))))
  effect_size <- r1/r2 #observed/expected
  results<-data.frame(pvalue = p.value,
                      effect_size = effect_size)
  
  
  return(results)
  
  
}




# # co-localization functions
# Count interactions between pairs 
colocal_mat <- function(x) {
  stopifnot(
    is.matrix(x), is.numeric(x),
    all(dim(x) > 0), ncol(x) > 1)
  
  if (is.null(colnames(x))) {
    colnames(x) <- seq_len(ncol(x))
  }
  df <- calc_pairs(x)
  df <- compute_interactions(x, df)
  return(df)
}
calc_pairs <- function(x) {
  x <- x > 0
  ab <- combn(colnames(x), 2)
  #y: number of spots with both cell types present
  y <- apply(ab, 2, function(.) sum(matrixStats::rowAlls(x[, ., drop = FALSE])))
  df <- data.frame(t(ab), y)
  names(df) <- c("pair1", "pair2", "n")
  return(df)
}
compute_interactions <- function(x, df) {
  y <- colnames(x)
  df$a <- factor(df$pair1, levels = y)
  df$b <- factor(df$pair2, levels = rev(y))
  t <- colSums(x > 0)
  a <- match(df$pair1, y)
  b <- match(df$pair2, y)
  df$ta <- t[a]
  df$tb <- t[b]
  df$pa <- df$n / df$ta
  df$pb <- df$n / df$tb
  df$ab_mean <- (df$pa + df$pb) / 2
  df$ab_comb <- df$n / (df$ta + df$tb)
  df$ab_comb2 <- df$n / (df$ta * df$tb)
  df$pairs <- paste(df$a ,'_', df$b,sep='')
  return(df)
}

# a1<-x[,1]>0
# a2<-x[,2]>0
# 
# jaccard::jaccard(as.numeric(a1),as.numeric(a2))
# 
# jaccard2 <- function(a, b) {
#   intersection = length(intersect(a, b))
#   union = length(a) + length(b) - intersection
#   return (intersection/union)
# }
# jaccard2(as.numeric(a1),as.numeric(a2))
# sum(as.numeric(a1),as.numeric(a2))
# length(which(a1==a2))
# 
# sumxy=sum(a1  & a2)
# unionxy <- sum(a1) + sum(a2) - sumxy
# sumxy/unionxy
# 
# a1 <- c(1,5,8,10,22,14,15,16,2,7)
# a2 <- c(10,12,13,2,7,9,2,7,23,15)
# 
# a1<-c(1,0,1,1,0)
# a2<-c(1,0,1,1,0)
# 
# jaccard::jaccard(as.numeric(a1),as.numeric(a2))
# jaccard2(as.numeric(a1),as.numeric(a2))
# 
# 
# Binary_A <- c(0,1,0,0,0,1,0,0,1,1) 
# Binary_B <- c(0,0,1,0,0,0,0,0,1,1) 
# jaccard::jaccard(Binary_A,Binary_B)


#get raw count data for H and M 
getRaw_WT<-function(Sample,Section,path_STsout="/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/SeuratObject"){
  (Sample_ID=paste(Sample,'-',Section,sep=''))
  image_path=PATH_IMAGES_LOWRES[paste(Sample,'-',Section,sep='')]
  
  
  
  if(Sample=='NSG'){
    seuratpath=paste(path_STsout,'/',Sample_ID,'_m','.rds',sep='')
    assign(paste("sdat",sep ='' ),readRDS(seuratpath), envir = .GlobalEnv)
    counts<-sdat@assays$RNA$counts

    rawM<-counts
    plot_dat<-data.frame(
      X=sdat$pxl_col_low,
      Y=sdat$pxl_row_low
    )

    return(list(
      image_path=image_path,
      sdat=sdat,
      plot_dat=plot_dat,
      rawM=rawM
    ))
    
  }else{
    seuratpath=paste(path_STsout,'/',Sample_ID,'_hm','.rds',sep='')
    assign(paste("sdat",sep ='' ),readRDS(seuratpath), envir = .GlobalEnv)
    plot_dat<-data.frame(
      X=sdat$pxl_col_low,
      Y=sdat$pxl_row_low
    )
    counts<-sdat@assays$RNA$counts
    rawH<-counts[grep(rownames(counts),pattern='GRCh38-'),]
    rawM<-counts[grep(rownames(counts),pattern='mm10-'),]
    rownames( rawH)<-str_sub(rownames( rawH),8)
    rownames( rawM)<-str_sub(rownames( rawM),8)
    return(list(
      image_path=image_path,
      sdat=sdat,
      rawH=rawH,
      rawM=rawM,
      plot_dat=plot_dat
    ))
  }
}

getSP_WT<-function(Sample,Section,speciessup='hm'){
  (Sample_ID=paste(Sample,'-',Section,sep=''))
  image_path=PATH_IMAGES_LOWRES[paste(Sample,'-',Section,sep='')]
  seuratpath=paste(path_STsout,'/',Sample_ID,'_',speciessup,'.rds',sep='')
  assign(paste("sdat",sep ='' ),readRDS(seuratpath), envir = .GlobalEnv)
  plot_dat<-data.frame(
    X=sdat$pxl_col_low,
    Y=sdat$pxl_row_low
  )
  return(
    list(
      sdat=sdat,
      plot_dat=plot_dat,
      image_path=image_path
    )
  )
}
#psudo spot, create data and then bayNorm#########
CONTROL_DATA<-function(
    GCGR_ID=GCGR_ID,
    species='M',
    list_spots=list_spots,
    seeds=12300){
  
  kk<-str_split(names(list_spots),pattern=',')
  unique_intervals=cbind(
    as.numeric(unlist(lapply(kk,function(x){str_sub(x[1],2)}))),
    as.numeric(unlist(lapply(kk,function(x){str_sub(x[2],0,-2)})))
  )

  # ind_sample=2
  # species='M'
  #Sample<-names(unique_samples)[ind_sample]
  #bayout=readRDS(paste('/home/clustor2/ma/w/wt215/BTC_project/BTC_R/Results_H_M/bay_objects/bayNorm_wsas/bay_',species,'_',unique_samples[ind_sample],'.rds',sep=''))
  
  # aa<-mget(ls(.GlobalEnv,pattern = paste("raw_",species,'_',GCGR_ID,"_",sep='')),envir = .GlobalEnv)
  # rawh<-do.call(cbind,aa)
  
  
  #ff<-intersect(GCGR_ID,c("GCGR-L24E","GL23E","GCGR-E43E","GCGR-L5E"))
  
  listr<-list()
  ratio<-NULL
  # if(GCGR_ID=="GCGR-L24E"){
  #   Total=1
  # }else{
  #   
  # }
  Total=2 #just use the first two sections

  for(Section in 1:Total){
    rdat<-getRaw_WT(Sample=GCGR_ID,Section=Section,path_STsout="/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/SeuratObject")
    listr[[Section]]=rdat$rawM
    
    ratio<-c(ratio,(colSums(rdat$rawH)/colSums(rdat$rawM)))
  }
  gg<-Reduce(intersect,lapply(listr,rownames))
  
  rawh<-Reduce(cbind,lapply(listr,function(x){return(x[gg,])}))
  
  
  
  
  lroi<-LABELS_3ROIS[colnames(rawh)]
  
  subb<- names(lroi)[!(lroi %in% c('necrotic/lowq'))]
  rawh<-rawh[,subb]
  
  numzero<-tabulate(rawh@i + 1)
  #hist(log2(numzero),breaks=200)
  #keep genes which have non zero counts in at least 1% of spots
  keep_genes<-rownames(rawh)[which(numzero>(dim(rawh)[2]*0.01))]
  
  
  #log2TC<-log2(colSums(rawh))
  #hist(log2TC,breaks=200)
  #rawh<-rawh[,names(log2TC)[which(log2TC>=6)]]
  #ratio<-Reduce(c,mget(ls(.GlobalEnv,pattern = paste("ratio_h_m_",GCGR_ID,"_",sep='')),envir = .GlobalEnv))
  if(GCGR_ID=="GCGR-L24E"){
    aa<-names(ratio[which((ratio)>=(0.01))])
  }else{
    aa<-names(ratio[which((ratio)>=(0.5) & (ratio)<=(1.5))])
  }

  #aa<-aa[grep(aa,pattern='S1.|S2.|S3.')]
  #table(Harmony_labels_m[aa])
  #aa<-intersect(aa,colnames(bayout$Bay_out))
  aa<-intersect(aa,colnames(rawh))
  
  s<-rawh[keep_genes,aa]
  p=apply(s,2,function(x){length(which(x!=0))/length(x)})
  
  #generate 500 pseudo spots
  num_pseudospots<-500
  #pick a spot which has the maximum number of non-zero count genes
  ps<-names(p[which.max(p)])
  
  set.seed(seeds)
  sim_beta<-runif(num_pseudospots,0.05,0.95)
  
  
  #print(length(rawh[keep_genes,ps]))
  
  simvec<-rawh[keep_genes,ps]
  #simvec=rowSums(rawh[rownames(bayout$Bay_out),aa])
  #which(simvec==0)
  simvec<-simvec[simvec>0]
  print(length(simvec))
  
  pseudodat=matrix(rep(simvec,num_pseudospots),nrow=length(simvec),ncol=num_pseudospots)
  pseudodat<-DownSampling(pseudodat,sim_beta)
  rownames(pseudodat)<-names(simvec)
  colnames(pseudodat)<-paste('c',seq(1,num_pseudospots),sep='')
  names(sim_beta)<-colnames(pseudodat)
  qb<-Beta_default(pseudodat,mean_beta = 0.5)
  
  #bin the data
  if(species=='M'){
    #for M, downsample M, low subsample prob means low tumor density
    list_pseudo<-list()
    for(i in 1:dim(unique_intervals)[1]){
      lin<-names(which(sim_beta<=unique_intervals[i,2] & sim_beta > unique_intervals[i,1]))
      list_pseudo[[i]]<-lin
    }
  } else if (species=='H'){
    #for H, downsample H, then low subsample prob means high tumor density
    list_pseudo<-list()
    for(i in dim(unique_intervals)[1]:1){
      lin<-names(which(sim_beta<=unique_intervals[i,2] & sim_beta > unique_intervals[i,1]))
      list_pseudo[[i]]<-lin
    }
  }
  
  
  names(list_pseudo)<-names(list_spots)
  
  baysub<-bayNorm(pseudodat,BETA_vec = qb,mean_version = TRUE,NCores = 16,BB_SIZE = FALSE)
  #tc<-t(t(pseudodat)/sim_beta)
  return(
    list(
      pseudodat=pseudodat,
      sim_beta=sim_beta,
      list_pseudo=list_pseudo,
      baysub=baysub
    )
  )
}



#GBM labels?#########

GBMlabel_WT<-function(
    Sample='GCGR-L3',
    LIST_DECON=LIST_DECON,
    Section=1,
    VEC_LUCY=c('AC-like','MES-like','NPC-like','OPC-like','Reactive AC-like','iOPC-like'),
    thres_quantile=NA,
    thres_prob=NA){
  
  deconU<-LIST_DECON[[Sample]][,VEC_LUCY]
  # deconU<-readRDS(paste('/home/clustor2/ma/w/wt215/PROJECT_ST/AUC_H/AUC_LucyZan_H_7/',Sample,'.rds',sep=''))[,VEC_LUCY]
  deconU<-deconU[grep(rownames(deconU),pattern=paste('S',Section,'.',sep='')),]
  deconUZ<-as.data.frame(apply(deconU,2,scale))
  rownames(deconUZ)<-rownames(deconU)
  TOP2_U<-Extract_Top2(deconUZ)
  

  L_U<-TOP2_U$TYPE1
  names(L_U)<-rownames(TOP2_U)
  
  # if(!is.na(thres_quantile)){
  #   spots_unclear<-rownames(TOP2_U)[which(TOP2_U$TOP1<max(c(quantile(TOP2_U$TOP1,thres_quantile),0)))]
  #   L_U[spots_unclear]<-'Unclear'
  # }
  if(!is.na(thres_prob)){
    aa<-apply(deconUZ,1,function(x){length(which(x<0))/dim(deconUZ)[2]})
    spots_unclear=names(aa)[which(aa>=thres_prob)]
    L_U[spots_unclear]<-'Unclear'
  }

  

  
  return(
    list(
      SCORE=deconU,
      TOP2=TOP2_U,
      Labels=L_U)
  )
}

#micronenvironment label
MElabel_WT<-function(
    Sample='GCGR-L3',
    LIST_DECON=LIST_DECON,
    Section=1,
    VEC_TYPES=c("ASC","Neurons" ,"MG" ,"OPC" ,"EC","OLG" ,"DC" ,"MAC","Pericytes" ,"Choroid plexus epithelial","Ependymocytes","Monocytes" )){
  
  rr<-getRaw_WT(Sample=Sample,Section=Section,path_STsout="/home/clustor2/ma/w/wt215/PROJECT_ST/RECOVER/SeuratObject")
  #check if choroid plexus exist
  prop_ttr<-length(which(rr$rawM[c('Ttr'),]>0))/dim(rr$rawM)[2]
  if(prop_ttr<0.2){
    VEC_TYPES<-setdiff(VEC_TYPES,c("Choroid plexus epithelial"))
  }
  
  deconU<-LIST_DECON[[Sample]][,VEC_TYPES]
  # deconU<-readRDS(paste('/home/clustor2/ma/w/wt215/PROJECT_ST/AUC_H/AUC_LucyZan_H_7/',Sample,'.rds',sep=''))[,VEC_LUCY]
  deconU<-deconU[grep(rownames(deconU),pattern=paste('S',Section,'.',sep='')),]
  deconUZ<-as.data.frame(apply(deconU,2,scale))
  rownames(deconUZ)<-rownames(deconU)
  TOP2_U<-Extract_Top2(deconUZ)
  L_U<-TOP2_U$TYPE1
  names(L_U)<-rownames(TOP2_U)
  
  return(
    list(
      SCORE=deconU,
      TOP2=TOP2_U,
      Labels=L_U)
  )
}

#by checking Ttr non-zero prop, I found that below 0.2, there are no Choroid plexus
# vecc<-NULL
# for(ind in 1:length(ALLSAMPLES)){
#   print(ALLSAMPLES[ind])
#   rr<-getRaw_WT(Sample=ALLSAMPLES[ind],Section=1)
#   vecc<-c(vecc,length(which(rr$rawM[c('Ttr'),]>0))/dim(rr$rawM)[2])
# }
# names(vecc)<-ALLSAMPLES



