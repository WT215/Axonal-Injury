
library(ggplotify)
library(directlabels)

library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
#frequent used colors##########
# message(
#   cyan$bold(
#     '
# Run plot_grid(qq2) to plot \n
# DIY color:                  \n
#   Discrete: scale_color_manual(values=SpatialColors(7)) \n
#   or: colorRampPalette(brewer.pal(12, "Set3"))(7)
#   or: colorRampPalette(palette_scvelo)(7)
#   or: scale_colour_gradient2(low="mediumseagreen",mid="lightyellow2",high="deeppink",na.value = "lightgrey")
#   Continuous: scale_color_continuous(type = "viridis")
#   '))

library(PBSmapping)
palette_scvelo<-c("#FFD0A6", "#B0CEE4", "#DDCECA" ,"#F1AAAB","#BFE2BF","#EFF0CD","#F6D5EC","#DED1EB","#D8D8D8")


circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  #https://stackoverflow.com/questions/6862742/draw-a-circle-with-ggplot2
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


#RNA velocity related######
#copied from https://github.com/velocyto-team/velocyto.R/blob/master/R/momentum_routines.R
read.loom.matrices <- function(file, engine='hdf5r') {
  if (engine == 'h5'){
    cat('reading loom file via h5...\n')
    f <- h5::h5file(file,mode='r');
    cells <- f["col_attrs/CellID"][];
    genes <- f["row_attrs/Gene"][];
    dl <- c(spliced="/layers/spliced",unspliced="/layers/unspliced",ambiguous="/layers/ambiguous");
    if("/layers/spanning" %in% h5::list.datasets(f)) {
      dl <- c(dl,c(spanning="/layers/spanning"))
    }
    dlist <- lapply(dl,function(path) {
      m <- as(f[path][],'dgCMatrix'); rownames(m) <- genes; colnames(m) <- cells; return(m)
    })
    h5::h5close(f)
    return(dlist)
  } else if (engine == 'hdf5r') {
    cat('reading loom file via hdf5r...\n')
    f <- hdf5r::H5File$new(file, mode='r')
    cells <- f[["col_attrs/CellID"]][]
    genes <- f[["row_attrs/Gene"]][]
    dl <- c(spliced="layers/spliced",
            unspliced="layers/unspliced",
            ambiguous="layers/ambiguous")
    if("layers/spanning" %in% hdf5r::list.datasets(f)) {
      dl <- c(dl, c(spanning="layers/spanning"))
    }
    dlist <- lapply(dl, function(path) {
      m <- as(t(f[[path]][,]),'dgCMatrix')
      rownames(m) <- genes; colnames(m) <- cells;
      return(m)
    })
    f$close_all()
    return(dlist)
  }
  else {
    warning('Unknown engine. Use hdf5r or h5 to import loom file.')
    return(list())
  }
}


##RNA velocity plot########
extract_velo_arrows<-function(velo_out,grid.n=20){
  gv<-velo_out$gvel
  gs<-velo_out$geshifts
  es<-velo_out$eshifts
  nd<-velo_out$vel
  cc<-velo_out$cc
  
  
  garrows=velo_out$garrows
  ars<-as.data.frame(velo_out$arrows)
  
  rx <- range(c(range(ars$x0),range(ars$x1)))
  ry <- range(c(range(ars$y0),range(ars$y1)))
  gx <- seq(rx[1],rx[2],length.out=grid.n)
  gy <- seq(ry[1],ry[2],length.out=grid.n)
  
  max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
  alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
  
  
  w <- expand.grid( gx,gy )
  plot_grid<-data.frame(
    X_grid=w[,1],
    Y_grid=w[,2]
  )
  plot_arrows<-as.data.frame(garrows)
  return(list(
    plot_arrows=plot_arrows,
    plot_grid=plot_grid
  ))
}


######FeaturePlot_wt###########
FeaturePlot_wt<-function(
  Data,
  Seurat_DimReduc,
  some_genes_used,
  cols =c('lightgrey', 'blue'),
  textsize=14,
  point.shape=20,
  point.size=1
  ){
  
  x.seurat <- CreateSeuratObject(counts =Data,assay = 'bayNorm')
  x.seurat@reductions$tsne<-Seurat_DimReduc
  reduction='tsne'

  cells<-colnames(x.seurat)
  

  pt.size = NULL
  min.cutoff = NA
  max.cutoff = NA
  order = FALSE
  
  dims<-c(1,2)
  dims <- paste0(Key(object = x.seurat[[reduction]]), dims)
  features=some_genes_used
  Features=some_genes_used
  

  
  data <- FetchData(
    object = x.seurat,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = 'data'
  )
  
  
  #Seurat's default
  # min.cutoff <- mapply(
  #   FUN = function(cutoff, feature) {
  #     return(ifelse(
  #       test = is.na(x = cutoff),
  #       yes = min(data[, feature]),
  #       no = cutoff
  #     ))
  #   },
  #   cutoff = min.cutoff,
  #   feature = features
  # )
  # max.cutoff <- mapply(
  #   FUN = function(cutoff, feature) {
  #     return(ifelse(
  #       test = is.na(x = cutoff),
  #       yes = max(data[, feature]),
  #       no = cutoff
  #     ))
  #   },
  #   cutoff = max.cutoff,
  #   feature = features
  # )
  
  
  #modified by Wenhao
  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = quantile(data[, feature],0.01),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = quantile(data[, feature],0.99),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )
  
  brewer.gran <- ifelse(
    test = length(x = cols) == 1,
    yes = brewer.pal.info[cols, ]$maxcolors,
    no = length(x = cols)
  )
  
  
  # Apply cutoffs
  data[, 4:ncol(x = data)] <- sapply(
    X = 4:ncol(x = data),
    FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index - 3], data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index - 3], data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      if (brewer.gran == 2) {
        return(data.feature)
      }
      data.cut <- if (all(data.feature == 0)) {
        0
      }
      else {
        as.numeric(x = as.factor(x = cut(
          x = as.numeric(x = data.feature),
          breaks = brewer.gran
        )))
      }
      return(data.cut)
    }
  )
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells

  split.by = NULL
  shape.by = NULL
  blend = FALSE
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells, drop = TRUE],
      object[[split.by, drop = TRUE]][cells, drop = TRUE]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  # Make list of plots
  plots <- vector(
    mode = "list",
    length = ifelse(
      test = blend,
      yes = 4,
      no = length(x = features) * length(x = levels(x = data$split))
    )
  )

  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  
  
  gout_list<-list()
  for(ii in 1:length(Features)){
    feature=Features[ii]
    ident <- levels(x = data$split)[1]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    
    data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
    
    cols.use <- NULL
    
    cols.grad<-cols
    
    
    data.single$orderr=rank(data.single[,feature],ties.method="first")
    kfe<-'orderr'
    
    
    plot<-ggplot(data = data.single,
                 mapping = aes_string(
                   x = dims[1],
                   y = dims[2],
                   color = paste0("`", feature, "`")
                 )) +
      geom_point(
        shape=point.shape,
        size=point.size
      ) +
      guides(color ='none',size='none') +
      scale_color_gradientn(
        colors = cols.grad,
        guide = "colorbar"
      )+
      labs(color = NULL,title=paste('Gene:',feature,sep=''))+
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            text = element_text(size=textsize))
    
    gout_list[[ii]]<-plot
  }
  

  
  return(gout_list)
  
}


#spatial plotting########
##from BayesSpace########
adjust_hex_centers <- function(spot_positions,r=1/2) {
  ## R = circumradius, distance from center to vertex
  ## r = inradius, distance from center to edge midpoint
  R <- (2 / sqrt(3)) * r
  
  ## Start at (1-indexed origin)
  spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos) + 1
  spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos) + 1
  
  ## Shift centers up so rows are adjacent
  spot_positions$y.pos <- spot_positions$y.pos * R * (3/2)
  
  ## Spot columns are offset by row
  ## (i.e. odd rows have odd numbered columns, even rows have even)
  ## Shift centers to the left so columns are adjacent (but hexes stay offset)
  spot_positions$x.pos <- (spot_positions$x.pos + 1) / 2
  
  spot_positions
}

adjust_hex_centers_px <- function(plot_dat,x='X',y='Y',spot_diameter_fullres,tissue_lowres_scalef,rotate='FALSE') {
  r=spot_diameter_fullres/2*tissue_lowres_scalef*3/2
  R <- (2 / sqrt(3)) * r
  
  spot_positions <- plot_dat[,c(x,y)]
  colnames(spot_positions)<-c("x.pos", "y.pos") 
  # spot_positions$x.pos <- spot_positions$x.pos - min(spot_positions$x.pos)+1
  # spot_positions$y.pos <- spot_positions$y.pos - min(spot_positions$y.pos)+1
  

  #rotate the hexagon?
  if(rotate){
    R2<-sqrt(R*R-r*r)
    vertex_offsets <- data.frame(x.offset=c(-R2,-R,-R2,R2,R,R2),
                                 y.offset=c(r,0,-r,-r,0,r))
  }else{
    ## vertices of each hex (with respect to center coordinates)
    ## start at top center, loop clockwise
    vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                                 y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  }

  
  
  
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  spot_vertices$spot<-rownames(plot_dat)
  
  #spot_vertices['expres_choice']=plot_dat[spot_vertices$spot,expres_choice]
  
  return(list(
    spot_positions=spot_positions,
    spot_vertices=spot_vertices
  ))
}


spatial_plot_self2_px<-function(plot_dat,
                                image_path=NULL,
                                images_cl=NULL,
                                x='X',
                                y='Y',
                                
                                #hexagonal params
                                rotate=FALSE,
                                hexagonal=FALSE, #TRUE for points
                                spot_diameter_fullres='',
                                tissue_lowres_scalef='',
                                
                                
                                expres_choice='clusters',
                                point_size=2,
                                alpha_image=1,
                                alpha_point=1,
                                image_w_h=c(1,1),
                                piecontents=NULL,
                                pie_scale=0.5,
                                shape=21,
                                color_pie=NA,
                                show.legend=TRUE,
                                ...){
  
  
  if(!is.null(image_path)){
    images_cl <- list()
    images_cl[[1]]=read.bitmap(image_path)
    #Transparency, learnt from https://stackoverflow.com/questions/11357926/r-add-alpha-value-to-png-image
    if(dim(images_cl[[1]])[3]==3){
      images_cl[[1]]<-abind(images_cl[[1]],matrix(1,nrow=dim(images_cl[[1]])[1],ncol=dim(images_cl[[1]])[2]),along=3)
    }
    
    
    images_cl[[1]]<- matrix(
      rgb(images_cl[[1]][,,1],images_cl[[1]][,,2],images_cl[[1]][,,3], images_cl[[1]][,,4] *alpha_image),
      nrow=dim(images_cl[[1]])[1],
      ncol=dim(images_cl[[1]])[2],
      byrow=FALSE)
    
  }
  
  
  grobs <- list()
  grobs[[1]]=rasterGrob(images_cl[[1]], width=unit(image_w_h[1],"npc"), height=unit(image_w_h[2],"npc"))
  height <- list()
  width <- list()
  height[[1]] <-  data.frame(height = nrow(images_cl[[1]]))
  width[[1]] <- data.frame(width = ncol(images_cl[[1]]))
  height <- bind_rows(height)
  width <- bind_rows(width)
  
  
  images_tibble <- tibble(sample=factor(c('figure1')), grob=grobs)
  images_tibble$height <- height$height
  images_tibble$width <- width$width
  
  qq_histo<-ggplot()+
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)
  
  
  if(hexagonal){
    adout<-adjust_hex_centers_px(plot_dat,x=x,y=y,
                                 spot_diameter_fullres=spot_diameter_fullres,
                                 tissue_lowres_scalef=tissue_lowres_scalef,rotate=rotate)
    spot_vertices=adout$spot_vertices
    spot_vertices['expres_choice']=plot_dat[spot_vertices$spot,expres_choice]
    spot_vertices['alpha_point']=plot_dat[spot_vertices$spot,alpha_point]
    
    ## Flip to match image orientation
    spot_vertices$x.vertex <- -spot_vertices$x.vertex
    
    qq_od<-ggplot()+
      geom_polygon(data=spot_vertices, aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~expres_choice, alpha=~(alpha_point)),show.legend=show.legend) +
      labs(fill="") +
      xlim(max(spot_vertices$x.vertex),min(spot_vertices$x.vertex))+
      ylim(max(spot_vertices$y.vertex),min(spot_vertices$y.vertex))+
      theme_void()
    
    
    qq_oo<-qq_histo+
      geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
      geom_polygon(data=spot_vertices, 
                   aes_(x=~(-x.vertex), y=~(y.vertex), group=~spot, fill=~(expres_choice), alpha=~(alpha_point)),show.legend=show.legend) +
      labs(fill="")+
      coord_cartesian(expand = FALSE) +
      xlim(0,max(width[[1]]))+
      ylim(max(height[[1]]),0)
    
    
    
    coord_settings <- list(coord_cartesian(expand = FALSE),
                           xlim(0, max(width[[1]])),
                           ylim(max(height[[1]]), 0))
    
    return(list(
      coord_settings=coord_settings,
      spot_vertices=spot_vertices,
      listwh=list(width=width,height=height),
      only_histo=qq_histo,
      only_dots=qq_od,
      with_dots=qq_oo))
    
  }else{
    qq_oo<- qq_histo+
      geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
      geom_point(plot_dat,mapping=aes(x=(!! rlang::sym(x)),y=(!! rlang::sym(y)),color=!! rlang::sym(expres_choice)),size = point_size,shape = shape, show.legend = show.legend,alpha=alpha_point)+
      coord_cartesian(expand = FALSE) +
      xlim(0,max(width[[1]]))+
      ylim(max(height[[1]]),0)
    
    qq_od<-ggplot()+
      geom_point(plot_dat,mapping=aes(x=(!! rlang::sym(x)),y=(!! rlang::sym(y)),color=!! rlang::sym(expres_choice)),size = point_size, stroke = 0.25,alpha=alpha_point,shape = shape)+
      coord_cartesian(expand = FALSE) +
      xlim(0,max(width[[1]]))+
      ylim(max(height[[1]]),0)+
      theme_bw()
    
    if(!is.null(piecontents)){
      qq_pie<-ggplot()+
        geom_scatterpie(data=plot_dat,aes(x=X,y=Y),
                        cols=piecontents, color=color_pie, alpha=.8,pie_scale=pie_scale,...)+
        coord_cartesian(expand = FALSE) +
        xlim(0,max(width[[1]]))+
        ylim(max(height[[1]]),0)+
        theme_bw()
      
      
      return(list(
        only_histo=qq_histo,
        only_histo_pie=qq_pie,
        only_dots=qq_od,
        with_dots=qq_oo))
    }
    
    
    return(list(
      only_histo=qq_histo,
      only_dots=qq_od,
      with_dots=qq_oo))
  }
}

spatial_plot_BayesSpace<-function(plot_dat,x='X',y='Y',expres_choice='clusters',r = 1/2){
  R = (2 / sqrt(3)) * r
  
  spot_positions <- plot_dat[,c(x,y)]
  colnames(spot_positions)<-c("x.pos", "y.pos")   
  
  spot_positions <- adjust_hex_centers(spot_positions,r=r)
  
  ## vertices of each hex (with respect to center coordinates)
  ## start at top center, loop clockwise
  vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
                               y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  
  spot_vertices <- merge(spot_positions, vertex_offsets)
  spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  spot_vertices$spot<-rownames(plot_dat)
  
  spot_vertices['expres_choice']=plot_dat[spot_vertices$spot,expres_choice]
  
  
  ## Flip to match image orientation
  spot_vertices$x.vertex <- -spot_vertices$x.vertex
  
  gout<-ggplot(data=spot_vertices, 
               aes_(x=~x.vertex, y=~y.vertex, group=~spot, fill=~(expres_choice))) +
    geom_polygon() +
    labs(fill="") +
    coord_equal() +
    theme_void()
  
  
  return(list(
    spot_vertices=spot_vertices,
    gout=gout
  ))
}



## From SPARK####
spatial_plot_wt <- function(pltdat, igene, coordinates = NULL, main = F, titlesize = 2,  pointsize = 3, xpand = 0, ypand = 1, title = NULL) {
  
  pltdat<-as.data.frame(t(pltdat))


  
  if (is.null(coordinates)) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[,1]), split = "x"))), ncol = 2)
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    if(dim(coordinates)[2]==2){
      coordinates<-as.data.frame(coordinates)
      colnames(coordinates)<-c("x", "y")
    }
    pd <- cbind(coordinates,pltdat)
  }
  
  igene_used<-which(colnames(pd)==igene)
  
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = pd[, igene_used])) + 
    geom_point(size = pointsize) + 
    scale_color_gradientn(colours = pal(5)) + 
    scale_x_discrete(expand = c(xpand, ypand)) + 
    scale_y_discrete(expand = c(xpand, ypand)) + 
    coord_equal() + 
    theme_bw()
  if (main) {
    if (is.null(title)) {
      title = colnames(pd)[igene_used]
    }
    out = gpt + 
      labs(title = title, x = NULL, y = NULL) + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
  }
  return(out)
}# end func



geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



##self made plot with figures##########

images_tibble_fun<-function(image_path="",alpha_image=1,image_w_h=c(1,1)){
  images_cl <- list()
  images_cl[[1]]=read.bitmap(image_path)
  #Transparency, learnt from https://stackoverflow.com/questions/11357926/r-add-alpha-value-to-png-image
  if(dim(images_cl[[1]])[3]==3){
    images_cl[[1]]<-abind(images_cl[[1]],matrix(1,nrow=dim(images_cl[[1]])[1],ncol=dim(images_cl[[1]])[2]),along=3)
  }
  
  
  images_cl[[1]]<- matrix(
    rgb(images_cl[[1]][,,1],images_cl[[1]][,,2],images_cl[[1]][,,3], images_cl[[1]][,,4] *alpha_image),
    nrow=dim(images_cl[[1]])[1],
    ncol=dim(images_cl[[1]])[2],
    byrow=FALSE)
  
  
  grobs <- list()
  grobs[[1]]=rasterGrob(images_cl[[1]], width=unit(image_w_h[1],"npc"), height=unit(image_w_h[2],"npc"))
  height <- list()
  width <- list()
  height[[1]] <-  data.frame(height = nrow(images_cl[[1]]))
  width[[1]] <- data.frame(width = ncol(images_cl[[1]]))
  height <- bind_rows(height)
  width <- bind_rows(width)
  
  
  images_tibble <- tibble(sample=factor(c('figure1')), grob=grobs)
  images_tibble$height <- height$height
  images_tibble$width <- width$width
  return(images_tibble)
}


spatial_plot_self2<-function(plot_dat,
                             image_path=NULL,
                             images_cl=NULL,
                             x='X',
                             y='Y',
                             expres_choice='clusters',
                             point_size=2,
                             alpha_image=1,
                             alpha_point=1,
                             image_w_h=c(1,1),
                             piecontents=NULL,
                             pie_scale=0.5,
                             shape=21,
                             color_pie=NA,
                             show.legend=FALSE,
                             ...){
  if(!is.null(image_path)){
    images_cl <- list()
    images_cl[[1]]=read.bitmap(image_path)
    #Transparency, learnt from https://stackoverflow.com/questions/11357926/r-add-alpha-value-to-png-image
    if(dim(images_cl[[1]])[3]==3){
      images_cl[[1]]<-abind(images_cl[[1]],matrix(1,nrow=dim(images_cl[[1]])[1],ncol=dim(images_cl[[1]])[2]),along=3)
    }
    
    
    images_cl[[1]]<- matrix(
      rgb(images_cl[[1]][,,1],images_cl[[1]][,,2],images_cl[[1]][,,3], images_cl[[1]][,,4] *alpha_image),
      nrow=dim(images_cl[[1]])[1],
      ncol=dim(images_cl[[1]])[2],
      byrow=FALSE)
    
  }
  
  
  grobs <- list()
  grobs[[1]]=rasterGrob(images_cl[[1]], width=unit(image_w_h[1],"npc"), height=unit(image_w_h[2],"npc"))
  height <- list()
  width <- list()
  height[[1]] <-  data.frame(height = nrow(images_cl[[1]]))
  width[[1]] <- data.frame(width = ncol(images_cl[[1]]))
  height <- bind_rows(height)
  width <- bind_rows(width)
  
  
  images_tibble <- tibble(sample=factor(c('figure1')), grob=grobs)
  images_tibble$height <- height$height
  images_tibble$width <- width$width
  
  qq_histo<-ggplot()+
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)
  
  
  qq_oo<- qq_histo+
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(plot_dat,mapping=aes(x=(!! rlang::sym(x)),y=(!! rlang::sym(y)),color=!! rlang::sym(expres_choice)),size = point_size,shape = shape, show.legend = show.legend,alpha=alpha_point)+
    coord_cartesian(expand = FALSE) +
    xlim(0,max(width[[1]]))+
    ylim(max(height[[1]]),0)
  
  qq_od<-ggplot()+
    geom_point(plot_dat,mapping=aes(x=(!! rlang::sym(x)),y=(!! rlang::sym(y)),color=!! rlang::sym(expres_choice)),size = point_size, stroke = 0.25,alpha=alpha_point,shape = shape)+
    coord_cartesian(expand = FALSE) +
    xlim(0,max(width[[1]]))+
    ylim(max(height[[1]]),0)+
    theme_bw()
  
  if(!is.null(piecontents)){
    qq_pie<-ggplot()+
      geom_scatterpie(data=plot_dat,aes(x=X,y=Y),
                      cols=piecontents, color=color_pie, alpha=.8,pie_scale=pie_scale,...)+
      coord_cartesian(expand = FALSE) +
      xlim(0,max(width[[1]]))+
      ylim(max(height[[1]]),0)+
      theme_bw()
    
    
    return(list(
      only_histo=qq_histo,
      only_histo_pie=qq_pie,
      only_dots=qq_od,
      with_dots=qq_oo))
  }
  
  
  return(list(
    only_histo=qq_histo,
    only_dots=qq_od,
    with_dots=qq_oo))
  
}

readimage_wt<-function(image_path=NULL,alpha_image=1){
  images_cl <- list()
  images_cl[[1]]=read.bitmap(image_path)
  #Transparency, learnt from https://stackoverflow.com/questions/11357926/r-add-alpha-value-to-png-image
  if(dim(images_cl[[1]])[3]==3){
    mat_image<-abind(images_cl[[1]],matrix(1,nrow=dim(images_cl[[1]])[1],ncol=dim(images_cl[[1]])[2]),along=3)
  }

  return(mat_image)
}

printimage_wt<-function(mat_image,image_path=NULL){
  if(!is.null(image_path)){
    image_path=NULL
    
    mat_image=read.bitmap(image_path)
    if(dim(mat_image)[3]==3){
      mat_image<-abind(mat_image,matrix(1,nrow=dim(mat_image)[1],ncol=dim(mat_image)[2]),along=3)
    }
  }
  
  images_cl <- list()
  images_cl[[1]]<- matrix(
    rgb(mat_image[,,1],mat_image[,,2],mat_image[,,3], mat_image[,,4] *1),
    nrow=dim(mat_image)[1],
    ncol=dim(mat_image)[2],
    byrow=FALSE)
  
  image_w_h=c(1,1)
  grobs <- list()
  
  grobs[[1]]=rasterGrob(mat_image, width=unit(image_w_h[1],"npc"), height=unit(image_w_h[2],"npc"))
  height <- list()
  width <- list()
  height[[1]] <-  data.frame(height = nrow(mat_image))
  width[[1]] <- data.frame(width = ncol(mat_image))
  height <- bind_rows(height)
  width <- bind_rows(width)
  
  images_tibble <- tibble(sample=factor(c('figure1')), grob=grobs)
  images_tibble$height <- height$height
  images_tibble$width <- width$width
  
  
  image_ori<-ggplot()+
    geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)
  return(image_ori)
  
}


#other utils from Seurat############
SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

RandomName <- function(length = 5L, ...) {
  CheckDots(..., fxns = 'sample')
  return(paste(sample(x = letters, size = length, ...), collapse = ''))
}
CheckDots <- function(..., fxns = NULL) {
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(
    X = fxns,
    FUN = function(x) {
      x <- tryCatch(
        expr = if (isS3stdGeneric(f = x)) {
          as.character(x = methods(generic.function = x))
        } else {
          x
        },
        error = function(...) {
          return(x)
        }
      )
      x <- if (is.character(x = x)) {
        sapply(X = x, FUN = argsAnywhere, simplify = FALSE, USE.NAMES = TRUE)
      } else if (length(x = x) <= 1) {
        list(x)
      }
      return(sapply(
        X = x,
        FUN = function(f) {
          return(names(x = formals(fun = f)))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  ))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  } else if (any(fxn.null)) {
    warning(
      "The following functions passed could not be found: ",
      paste(names(x = which(x = fxn.null)), collapse = ', '),
      call. = FALSE,
      immediate. = TRUE
    )
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = 'logical', length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = '...', x = fxn.args[[i]], fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(
        X = Filter(f = nchar, x = dfxns),
        FUN = function(x) {
          if (isS3method(method = x)) {
            x <- unlist(x = strsplit(x = x, split = '\\.'))
            x <- x[length(x = x) - 1L]
          }
          return(x)
        },
        FUN.VALUE = character(length = 1L)
      )
      message(
        "The following functions and any applicable methods accept the dots: ",
        paste(unique(x = fx), collapse = ', ')
      )
      if (any(nchar(x = dfxns) < 1)) {
        message(
          "In addition, there is/are ",
          length(x = Filter(f = Negate(f = nchar), x = dfxns)),
          " other function(s) that accept(s) the dots"
        )
      }
    } else {
      message("There is/are ", length(x = dfxns), 'function(s) that accept(s) the dots')
    }
  } else {
    unused <- Filter(
      f = function(x) {
        return(!x %in% unlist(x = fxn.args))
      },
      x = args.names
    )
    if (length(x = unused) > 0) {
      msg <- paste0(
        "The following arguments are not used: ",
        paste(unused, collapse = ', ')
      )
      switch(
        EXPR = getOption(x = "Seurat.checkdots"),
        "warn" = warning(msg, call. = FALSE, immediate. = TRUE),
        "stop" = stop(msg),
        "silent" = NULL,
        stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent")
      )
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message(
          "Suggested parameter: ",
          paste(unused.hints, "instead of", names(x = unused.hints), collapse = '; '),
          "\n"
        )
      }
    }
  }
}

####bindSC##########
UMAP_wt_plot<-function (meta = NULL, color = NULL, xlim = NULL, alpha = 0.1,size=1,ylim = NULL, showLabel = TRUE, title = NULL, mylabel = NULL) 
{
  library(ggrepel)
  dt <- data.frame(UMAP1 = meta[, "UMAP1"], UMAP2 = meta[, "UMAP2"], label = meta[, color])
  dt$label <- as.factor(dt$label)
  x_min <- xlim[1]
  x_max <- xlim[2]
  y_min <- ylim[1]
  y_max <- ylim[2]
  label_pos <- aggregate(. ~ label, dt, median)
  p1 <- ggplot(dt, aes(x = UMAP1, y = UMAP2, color = label)) + 
    geom_point(alpha = alpha, size = size) + scale_colour_manual(values = mylabel) + 
    xlim(x_min, x_max) + ylim(y_min, y_max) + 
    theme_classic() + 
    theme(legend.position = "none") + 
    geom_text_repel(data = label_pos, 
                    aes(label = label), 
                    color = "black", 
                    fontface = "bold", 
                    alpha = 0.75, 
                    box.padding = 0.5, 
                    point.padding = 0.1) + 
    NoLegend() + 
    theme(axis.text = element_blank(), 
          axis.title = element_blank(), 
          axis.ticks = element_blank())
  return(p1)
}



#complex heatmap: interval breaks from pheatmap#########

scale_vec_colours <- function(x, col = rainbow(10), breaks = NA, na_col){
  res <- col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]
  #res[is.na(res)] <- na_col
  res[which(x>=max(breaks))]<-col[length(col)]
  res[which(x<=min(breaks))]<-col[1]
  
  return(res)
}

scale_colours = function(mat, col = rainbow(10), breaks = NA, na_col){
  mat = as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks, na_col = na_col), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

colorRamp2_interval<-function(data,col,breaks,na_col="#DDDDDD"){
  qq<-scale_colours(data, col = col, breaks =breaks, na_col="#DDDDDD")

  col_fun = colorRamp2(as.vector(as.matrix(data)), as.vector(qq))
  return(col_fun)
  
}



#complex heatmap : marker genes: sep rows sep cols######
#some_genes:  0        
#  "AQP1"

#Labels_order: #5E4FA2        
#  "0"
#cluster_labels: S1.GCTAGTAGAGCTTGTA-1
#0

HEAT_WT<-function(heat_dat,cluster_labels,Labels_order,some_genes){
  
  topanno<-HeatmapAnnotation(
    clusters=anno_block(gp = gpar(fill =names(Labels_order) ),
                        labels = Labels_order, 
                        labels_gp = gpar(col = "black", fontsize = 3)), which='column',simple_anno_size = unit(1, "cm"))
  
  
  
  heat_labels<-factor(cluster_labels,levels=Labels_order)
  heat_dat_used<-heat_dat[,names(heat_labels)]
  heat_dat_used<-apply(heat_dat_used,1,function(x){(x-mean(x))/sd(x)})
  heat_dat_used<-t(heat_dat_used)
  
  Location=do.call(rbind,str_split(colnames(heat_dat_used),pattern='\\.'))[,1]
  Location<-paste(Location,'-',cluster_labels[colnames(heat_dat_used)],sep='')
  names(Location)<-colnames(heat_dat_used)
  
  numm_1<-cluster_labels[colnames(heat_dat_used)]
  numm_1<-droplevels(numm_1)
  
  numm_2<-factor(names(some_genes),levels=levels(numm_1))
  numm_3=Location[colnames(heat_dat_used)]
  
  
  
  Labels_order[match(names(some_genes),Labels_order)]
  
  
  
  
  leftanno<-HeatmapAnnotation(
    Genes=names(some_genes),
    which='row',
    col=list(
      Genes=setNames(names(Labels_order), Labels_order)),
    annotation_name_side ='top',
    show_legend =FALSE,
    show_annotation_name=FALSE
  )
  
  
  
  
  topanno<-HeatmapAnnotation(
    clusters=anno_block(gp = gpar(fill =setNames(names(Labels_order), Labels_order) ),
                        labels = Labels_order, 
                        labels_gp = gpar(col = "black", fontsize = 12)),
    which='column',
    simple_anno_size = unit(1, "cm"))
  
  numm_1<-factor(numm_1,levels=Labels_order)
  numm_2<-factor(numm_2,levels=Labels_order)
  
  
  phee<-Heatmap(matrix=MinMax(data = heat_dat_used, min = (-2.5), max = 2.5),
                
                col=  PurpleAndYellow(),
                cluster_columns=FALSE,
                cluster_rows=FALSE,
                column_names_rot=45,
                name='z-score of norm count',
                heatmap_legend_param = list(direction = "horizontal"),
                show_column_names = FALSE,
                top_annotation = topanno,
                left_annotation=leftanno,
                row_split=numm_2,
                column_split=numm_1,
                row_gap=unit(0.15, "cm"),
                column_gap=unit(0.15, "cm"),
                column_title=NULL,
                row_title=NULL
  )
  #mm font 3
  phee<-draw(phee,
             heatmap_legend_side = "top",annotation_legend_side = "right")
  
  
  gout_heat<-grid.grabExpr(draw(phee))
  return(gout_heat)
  
  
}

##extract_cluster_ComplexHeat####
extract_cluster_ComplexHeat<-function(ht,Data,direction='row'){
  #It is more suggested to do as `ht = draw(ht); row_order(ht)`. 
  if(direction=='row'){
    suppressWarnings(roworder<-row_order(ht))
  }else if(direction=='column'){
    suppressWarnings(roworder<-column_order(ht))
    Data=t(Data)
  }
  qq=names(roworder)
  for (i in 1:length(roworder)){

    if (i == 1) {
      clu <- t(t(row.names(Data)[roworder[[i]]]))
      out <- cbind(clu, paste("cluster", qq[i], sep=""))
      colnames(out) <- c("GeneID", "Cluster")
    } else {
      clu <- t(t(row.names(Data)[roworder[[i]]]))
      clu <- cbind(clu, paste("cluster", qq[i], sep=""))
      out <- rbind(out, clu)
    }
  }
  
  
  clu<-unique(out[,2])
  list_genes<-lapply(clu,function(x){out[which(out[,2]==x),1]})
  names(list_genes)<-clu
  
  return(list_genes)
}





#4 quadrant plot from Neftel 2019 and Dr. Dieter Henrik Heiland#########
hlpr_gene_set_name <- function(string){
  
  stringr::str_remove(string = string, pattern = "^.+?_")
  
}
#Labels_order: names: label name, values: color code
quadrant4_plot_fun<-function(data,color_to='clusters',states,Labels_order,pt_size=2,pt_alpha=1,font_corner=8,pos_corner=0.8,pt_shape=19,maxx=NULL,label_prop=FALSE){
  shifted_df <-
    tidyr::pivot_longer(
      data = data,
      cols = dplyr::all_of(states),
      names_to = "gene_set",
      values_to = "gene_set_expr"
    )
  
  max_localisation <-
    dplyr::group_by(shifted_df, barcodes) %>%
    dplyr::filter(gene_set_expr == max(gene_set_expr)) %>%
    dplyr::ungroup() %>%
    # rename the remaining gene sets to 'max_gene_set'
    dplyr::select(barcodes, max_gene_set = gene_set, max_expr = gene_set_expr) %>%
    # assign the vertical localistion of the state plot depending on where the maximum occured
    dplyr::mutate(max_loc = dplyr::if_else(max_gene_set %in% states[1:2], true = "top", false = "bottom"))
  # if(sum(duplicated(max_localisation$barcodes)>0)){
  #   max_localisation <-max_localisation [!duplicated(max_localisation$barcodes),]
  # }
  
  
  # calculate the x-position
  with_x_positions <-
    dplyr::left_join(x = data, y = max_localisation, by = "barcodes") %>%
    dplyr::mutate(
      pos_x = dplyr::case_when(
        max_loc == "top" & !!sym(states[1]) > !!sym(states[2]) ~ (log2(abs((!!sym(states[1]) - !!sym(states[2])) + 1)) * -1),
        max_loc == "top" & !!sym(states[2]) > !!sym(states[1]) ~ log2(abs((!!sym(states[2]) - !!sym(states[1])) + 1)),
        max_loc == "bottom" & !!sym(states[3]) > !!sym(states[4]) ~ (log2(abs((!!sym(states[3]) - !!sym(states[4])) + 1)) * -1),
        max_loc == "bottom" & !!sym(states[4]) > !!sym(states[3]) ~ log2(abs((!!sym(states[4]) - !!sym(states[3])) + 1)))
    )
  
  
  # calculate the y-position
  plot_df<-with_x_positions
  plot_df['pos_y']<-0
  
  for(i in 1:dim(with_x_positions)[1]){
    if(with_x_positions[i,'max_loc']=='bottom'){
      plot_df[i,'pos_y']<-(log2(abs(max(with_x_positions[i,states[c(3,4)]]) - max(with_x_positions[i,states[c(1,2)]]) + 1)) * -1)
    }else if(with_x_positions[i,'max_loc']=='top'){
      plot_df[i,'pos_y']<-(log2(abs(max(with_x_positions[i,states[c(1,2)]]) - max(with_x_positions[i,states[c(3,4)]]) + 1)))
    }
  }
  aa<-union(which(is.na(plot_df$pos_x)),which(is.na(plot_df$pos_y)))
  if(length(aa)>0){
    plot_df<-plot_df[-aa,]
  }

  # plot_df <-
  #   dplyr::group_by(with_x_positions, barcodes) %>%
  #   dplyr::mutate(
  #     pos_y = dplyr::case_when(
  #       max_loc == "bottom" ~ (log2(abs(max(c(!!sym(states[3]), !!sym(states[4]))) - max(!!sym(states[1]), !!sym(states[2])) + 1)) * -1),
  #       max_loc == "top" ~ log2(abs(max(c(!!sym(states[1]), !!sym(states[2]))) - max(!!sym(states[3]), !!sym(states[4])) + 1))
  #     )
  #   ) %>%
  #   dplyr::filter(!base::is.na(pos_x) & !is.na(pos_y))
  
  

  
  
  states <- hlpr_gene_set_name(states)
  color_to_lab <- hlpr_gene_set_name(color_to)
  
  
  xlab <- base::bquote(paste("log2(Score "[.(states[3])]*" - Score "[.(states[4])]*")"))
  ylab <- base::bquote(paste("log2(Score "[.(states[2])]*" - Score "[.(states[1])]*")"))
  
  if(is.null(maxx)){
    max <- base::max(base::abs(plot_df$pos_x), base::abs(plot_df$pos_y),na.rm = TRUE)
  }else{
    max=maxx
  }

  
  
  #quadrant 1->4
  quadrantlabel=data.frame(
    X=c(max*pos_corner,-max*pos_corner,-max*pos_corner,max*pos_corner),
    Y=c(max*pos_corner,max*pos_corner,-max*pos_corner,-max*pos_corner),
    label=states[c(2,1,3,4)]
    #label=states
  )
  if(label_prop){
    prop1<-round(length(which(plot_df$pos_x>0 & plot_df$pos_y>0))/dim(plot_df)[1],2)
    prop2<-round(length(which(plot_df$pos_x<0 & plot_df$pos_y>0))/dim(plot_df)[1],2)
    prop3<-round(length(which(plot_df$pos_x<0 & plot_df$pos_y<0))/dim(plot_df)[1],2)
    prop4<-round(length(which(plot_df$pos_x>0 & plot_df$pos_y<0))/dim(plot_df)[1],2)
    
    indd<-which(quadrantlabel$X>0 & quadrantlabel$Y>0)
    quadrantlabel$label[indd]<-paste(quadrantlabel$label[indd],' (',prop1,')',sep='')
    indd<-which(quadrantlabel$X<0 & quadrantlabel$Y>0)
    quadrantlabel$label[indd]<-paste(quadrantlabel$label[indd],' (',prop2,')',sep='')
    indd<-which(quadrantlabel$X<0 & quadrantlabel$Y<0)
    quadrantlabel$label[indd]<-paste(quadrantlabel$label[indd],' (',prop3,')',sep='')
    indd<-which(quadrantlabel$X>0 & quadrantlabel$Y<0)
    quadrantlabel$label[indd]<-paste(quadrantlabel$label[indd],' (',prop4,')',sep='')
  }
  

  
  
  
  
  if(!is.null(color_to)){
    if(length(Labels_order)>1){
      gout=ggplot2::ggplot(data = plot_df) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
        ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
        ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_to), size = pt_size, alpha = pt_alpha, data = plot_df,shape=pt_shape) +
        #scale_color_continuous(type = "viridis")+
        scale_color_manual(values=(Labels_order),limits=names(Labels_order),na.value = 'lightgrey')+
        ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(x = xlab, y = ylab, color = color_to_lab)+
        #ggplot2::labs(x = '', y = '', color = color_to_lab)+
        geom_text(data=quadrantlabel,aes(x=X,y=Y,label=label),size=font_corner)
    }else{
      gout=ggplot2::ggplot(data = plot_df) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
        ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
        ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y", color = color_to), size = pt_size, alpha = pt_alpha, data = plot_df,shape=pt_shape) +
        scale_color_viridis(option='magma')+
        ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
        ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(x = xlab, y = ylab, color = color_to_lab)+
        #ggplot2::labs(x = '', y = '', color = color_to_lab)+
        geom_text(data=quadrantlabel,aes(x=X,y=Y,label=label),size=font_corner)
    }

    
  }else{
    gout=ggplot2::ggplot(data = plot_df) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_hline(yintercept = 0,  linetype = "dashed", color = "lightgrey") +
      ggplot2::geom_point(mapping = ggplot2::aes_string(x = "pos_x", y = "pos_y"), size = pt_size, alpha = pt_alpha, data = plot_df,shape=pt_shape) +
      ggplot2::scale_x_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
      ggplot2::scale_y_continuous(limits = c(-max*1.1, max*1.1), expand = c(0,0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = xlab, y = ylab, color = color_to_lab)+
      #ggplot2::labs(x = '', y = '', color = color_to_lab)+
      geom_text(data=quadrantlabel,aes(x=X,y=Y,label=label),size=font_corner)
    
  }

  
  return(list(gout=gout,plot_df=plot_df))
}

#SCP blend color #########
library(ggnewscale)
library(scales)

Blend2Color <- function(C1, C2, mode = "blend") {
  c1 <- C1[[1]]
  c1a <- C1[[2]]
  c2 <- C2[[1]]
  c2a <- C2[[2]]
  A <- 1 - (1 - c1a) * (1 - c2a)
  if (A < 1.0e-6) {
    return(list(c(0, 0, 0), 1))
  }
  if (mode == "mix") {
    Result <- (c1 + c2) / 2
    Result[Result > 1] <- 1
  }
  if (mode == "blend") {
    Result <- (c1 * c1a + c2 * c2a * (1 - c1a)) / A
    A <- 1
  }
  if (mode == "screen") {
    Result <- 1 - (1 - c1) * (1 - c2)
  }
  if (mode == "multiply") {
    Result <- c1 * c2
  }
  
  return(list(Result, A))
}
BlendRGBList <- function(Clist, mode = "blend", RGB_BackGround = c(1, 1, 1)) {
  N <- length(Clist)
  ClistUse <- Clist
  while (N != 1) {
    temp <- ClistUse
    ClistUse <- list()
    for (C in temp[1:(length(temp) - 1)]) {
      c1 <- C[[1]]
      a1 <- C[[2]]
      c2 <- temp[[length(temp)]][[1]]
      a2 <- temp[[length(temp)]][[2]]
      ClistUse <- append(ClistUse, list(Blend2Color(C1 = list(c1, a1 * (1 - 1 / N)), C2 = list(c2, a2 * 1 / N), mode = mode)))
    }
    N <- length(ClistUse)
  }
  Result <- list(ClistUse[[1]][[1]], ClistUse[[1]][[2]])
  Result <- RGBA2RGB(Result, BackGround = RGB_BackGround)
  return(Result)
}

blendcolors <- function(colors, mode = "blend") {
  colors <- colors[!is.na(colors)]
  if (length(colors) == 0) {
    return(NA)
  }
  if (length(colors) == 1) {
    return(colors)
  }
  rgb <- as.list(as.data.frame(col2rgb(colors) / 255))
  Clist <- lapply(rgb, function(x) {
    list(x, 1)
  })
  blend_color <- BlendRGBList(Clist, mode = mode)
  blend_color <- rgb(blend_color[1], blend_color[2], blend_color[3])
  return(blend_color)
}
RGBA2RGB <- function(RGBA, BackGround = c(1, 1, 1)) {
  A <- RGBA[[length(RGBA)]]
  RGB <- RGBA[[-length(RGBA)]] * A + BackGround * (1 - A)
  return(RGB)
}
adjcolors <- function(colors, alpha) {
  color_df <- as.data.frame(col2rgb(colors) / 255)
  colors_out <- sapply(color_df, function(color) {
    color_rgb <- RGBA2RGB(list(color, alpha))
    return(rgb(color_rgb[1], color_rgb[2], color_rgb[3]))
  })
  return(colors_out)
}
bg_color='lightgrey80'


palette_scp <- function(x, n = 100, palette = "Paired", palcolor = NULL, type = "auto",
                        matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80",palette_list) {
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }
  if (!palette %in% names(palette_list)) {
    stop("The palette name is invalid! You can check the available palette names with 'show_palettes()'. Or pass palette colors via the 'palcolor' parameter.")
  }
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) {
    palcolor <- palette_list[[palette]]
  }
  pal_n <- length(palcolor)
  
  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }
  
  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x))
    }
    n_x <- nlevels(x)
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- colorRampPalette(palcolor)(n_x)
    } else {
      color <- ifelse(rep(n_x, n_x) <= pal_n,
                      palcolor[1:n_x],
                      colorRampPalette(palcolor)(n_x)
      )
    }
    names(color) <- levels(x)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x) && all(!is.na(x))) {
      stop("'x' must be type of numeric when use continuous color palettes.")
    }
    if (all(is.na(x))) {
      values <- as.factor(rep(0, n))
    } else if (length(unique(na.omit(as.numeric(x)))) == 1) {
      values <- as.factor(rep(unique(na.omit(as.numeric(x))), n))
    } else {
      values <- cut(x, breaks = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = n + 1), include.lowest = TRUE)
    }
    
    n_x <- nlevels(values)
    color <- ifelse(rep(n_x, n_x) <= pal_n,
                    palcolor[1:n_x],
                    colorRampPalette(palcolor)(n_x)
    )
    names(color) <- levels(values)
    if (any(is.na(x))) {
      color <- c(color, setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      if (all(is.na(x))) {
        color <- NA_color
      } else if (length(unique(na.omit(x))) == 1) {
        color <- color[as.character(unique(na.omit(x)))]
        color[is.na(color)] <- NA_color
      } else {
        color <- color[as.character(values)]
        color[is.na(color)] <- NA_color
      }
    }
  }
  
  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}


SCP_blend_fun<-function(pltdat,
                        features,
                        color_blend_mode = 'multiply',
                        label_point_color='black',
                        label_point_size=1,
                        label_repulsion=20,
                        label_segment_color = "black",
                        label.fg='white',
                        label.bg='black',
                        label.bg.r=0.1,
                        label.size=4,
                        palette="Set1",
                        label = TRUE,
                        palcolor=NULL){
  
  dat=pltdat


  palette_list <- list(brewer.pal(length(features), palette))
  names(palette_list )<-palette
  colors <- palette_scp(features, type = "discrete", palette = palette, palcolor = palcolor,palette_list=palette_list)
  
  
  #color_blend_mode = c("blend", "mix", "screen", "multiply")
  #color_blend_mode = 'multiply'
  
  colors_list <- list()
  value_list <- list()
  pal_list <- list()
  temp_geom <- list()
  legend_list <- list()
  for (i in seq_along(colors)) {
    lower_col <- adjcolors(colors[i], 0.1)
    colors_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = TRUE, matched = TRUE, palcolor = c(lower_col, colors[i]),palette = palette,palette_list=palette_list)
    pal_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = FALSE, matched = FALSE, palcolor = c(lower_col, colors[i]),palette = palette,palette_list=palette_list)
    
    
    value_list[[i]] <- seq(min(dat[, names(colors)[i]], na.rm = TRUE), max(dat[, names(colors)[i]], na.rm = TRUE), length.out = 100)
    
    
    temp_geom[[i]] <- list(
      suppressWarnings(geom_point(data = dat, mapping = aes(x = .data[["X"]], y = .data[["Y"]], color = .data[[names(colors)[i]]]))),
      scale_color_gradientn(
        colours = pal_list[[i]],
        values = scales::rescale(value_list[[i]]), na.value = bg_color,
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
      ),
      new_scale_color()
    )
    legend_list[[i]] <- get_legend(ggplot(dat, aes(x = X, y = Y)) +
                                     temp_geom[[i]])
    
    
  }
  # plot(legend_list[[4]])
  # length(legend_list)
  
  for (j in seq_len(nrow(dat))) {
    dat[j, "color_blend"] <- blendcolors(sapply(colors_list, function(x) x[j]), mode = color_blend_mode)
  }
  
  
  dat["color_value"] <- colSums(col2rgb(dat[, "color_blend"]))
  dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_value"] <- NA
  dat <- dat[order(dat[, "color_value"], decreasing = TRUE, na.last = FALSE), ]
  dat[rowSums(is.na(dat[, names(colors)])) == length(colors), "color_blend"] <- bg_color
  
  
  dat[, "features"] <- paste(features, collapse = "|")
  
  p <- ggplot(dat)
  
  
  p<-p + suppressWarnings(geom_point(
    mapping = aes(x = .data[["X"]], y = .data[["Y"]], color = .data[["color_blend"]]))) +
    scale_color_identity() +
    theme(legend.position ='none')+
    new_scale_color()+
    theme_classic()
  
  
  if(label){
    
    label_df <- reshape2::melt(p$data, measure.vars = features)
    
    label_df <- label_df %>%
      dplyr::group_by(variable) %>%
      dplyr::filter(value >= quantile(value, 0.985, na.rm = TRUE) & value <= quantile(value, 0.99, na.rm = TRUE)) %>%
      dplyr::reframe(x = median(.data[["X"]]), y = median(.data[["Y"]])) %>%
      as.data.frame()
    colnames(label_df)[1] <- "label"
    label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
    label_df[, "rank"] <- seq_len(nrow(label_df))
    
    p <- p + geom_point(
      data = label_df, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
      color = label_point_color, size = label_point_size
    ) + geom_text_repel(
      data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]], color = .data[["label"]]),
      fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
      point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
      color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, show.legend = FALSE
    )
  }
  
  
  
  legend_nrow <- min(ceiling(sqrt(length(legend_list))), 3)
  total <- length(legend_list)
  leg_list <- list()
  n <- 1
  for (i in 1:total) {
    if (i == 1 || is.null(leg)) {
      leg <- legend_list[[i]]
    } else {
      leg <- cbind(leg, legend_list[[i]])
    }
    if (i %% legend_nrow == 0) {
      leg_list[[n]] <- leg
      leg <- NULL
      n <- n + 1
    }
    if (i %% legend_nrow != 0 && i == total) {
      ncol_insert <- dim(leg_list[[n - 1]])[2] - dim(leg)[2]
      for (col_insert in 1:ncol_insert) {
        leg <- gtable_add_cols(leg, sum(leg_list[[n - 1]]$widths) / ncol_insert, -1)
      }
      leg_list[[n]] <- leg
    }
  }
  legend <- do.call(rbind, leg_list)
  #plot(legend)
  
  legend2=legend
  
  grob <- ggplotGrob(p)
  grob <- gtable_add_cols(grob, sum(legend$widths), -1)
  grob <- gtable_add_grob(grob, legend, t = grob$layout[grepl(pattern = "panel", x = grob$layout$name), "t"], l = dim(grob)[2])
  
  gout<-plot_grid(grob)
  return(gout)
}

#cell2loc:https://github.com/BayraktarLab/cell2location/blob/9577c1645a22e74af1246453772ed9942eaa11ce/cell2location/plt/plot_spatial.py#L179##########
create_colormap<-function(vec_RGB=c(240,228,66),white_spacing=20){
  spacing = as.integer(white_spacing * 2.55)
  N = 255
  M = 3
  alphas=c(rep(0,spacing * M),seq(0, 1.0, length.out=(N - spacing) * M))
  vals=matrix(1,nrow=N*M,ncol=4)
  for(i in 1:3){
    vals[,i]=vec_RGB[i]/255
  }
  vals[, 4] = alphas
  return(vals)
}






clipR <- function(x, min_val, max_val) {
  pmin(pmax(x, min_val), max_val)
}





matcol<-rbind(
  c(240, 228, 66),
  #c(213, 94, 0),#red
  c(214, 39, 40),#red
  
  c(86, 180, 233),
  c(0, 158, 115),
  c(90, 20, 165),
  c(255, 133, 26),
  c(200, 200, 200),
  c(50, 50, 50)
)
rownames(matcol)<-c('yellow','red','blue','green','purple','orange','grey','white')

coolist<-list()
for(i in 1:dim(matcol)[1]){
  coolist[[i]]<-create_colormap(vec_RGB=matcol[i,])
}
names(coolist)<-rownames(matcol)

coolist_rgb<-list()
for(i in 1:length(coolist)){
  xx=coolist[[i]]
  coolist_rgb[[i]]<-apply(xx,1,function(x){
    rgb(x[1],x[2],x[3],x[4])
  })
  
}
names(coolist_rgb)<-rownames(matcol)



coolist_fun<-list()
for(i in 1:dim(matcol)[1]){
  coolist_fun[[i]]<-colorRampPalette(rgb(matcol[i,1],matcol[i,2], matcol[i,3],maxColorValue=255))
}
names(coolist_fun)<-rownames(matcol)


SCP_blend_fun_WT<-function(pltdat,
                           features,
                           features_colors=NULL,
                           image_path=NULL,
                           title='',
                           alpha_image=1,
                           image_w_h=c(1,1),
                           label_point_color='black',
                           label_point_size=1,
                           label_repulsion=20,
                           label_segment_color = "black",
                           label.fg='black',
                           label.bg='white',
                           label.bg.r=0.1,
                           label.size=4,
                           palette="Set1",
                           label = TRUE,
                           palcolor=NULL,
                           offset=0.01,
                           max_color_quantile=0.98){
  # r=1/2
  # R = (2 / sqrt(3)) * r
  # spot_positions<-coordinates
  # colnames(spot_positions)<-c("x.pos", "y.pos")   
  # spot_positions <- adjust_hex_centers(spot_positions)
  # 
  # ## vertices of each hex (with respect to center coordinates)
  # ## start at top center, loop clockwise
  # vertex_offsets <- data.frame(x.offset=c(0, r, r, 0, -r, -r),
  #                              y.offset=c(-R, -R/2, R/2, R, R/2, -R/2))
  # 
  # spot_vertices <- merge(spot_positions, vertex_offsets)
  # spot_vertices$x.vertex <- spot_vertices$x.pos + spot_vertices$x.offset
  # spot_vertices$y.vertex <- spot_vertices$y.pos + spot_vertices$y.offset
  # spot_vertices<-spot_vertices*s_dat@misc$IMAGE$scale_coor$tissue_lowres_scalef
  # 
  # 
  # spot_vertices$spot<-rownames(spot_positions)
  # 
  # 
  # spot_vertices['color_blend']<-dat[spot_vertices$spot,'color_blend']
  # 
  # ggplot() +
  #   geom_polygon(data=spot_vertices, 
  #                aes_(x=~x.vertex, y=~y.vertex, group=~spot,fill=~color_blend),color="#d8dcd6") 
  
  
  
  
  
  if(!is.null(image_path)){
    
    images_cl <- list()
    images_cl[[1]]=read.bitmap(image_path)
    #Transparency, learnt from https://stackoverflow.com/questions/11357926/r-add-alpha-value-to-png-image
    if(dim(images_cl[[1]])[3]==3){
      images_cl[[1]]<-abind(images_cl[[1]],matrix(1,nrow=dim(images_cl[[1]])[1],ncol=dim(images_cl[[1]])[2]),along=3)
    }
    
    
    images_cl[[1]]<- matrix(
      rgb(images_cl[[1]][,,1],images_cl[[1]][,,2],images_cl[[1]][,,3], images_cl[[1]][,,4] *alpha_image),
      nrow=dim(images_cl[[1]])[1],
      ncol=dim(images_cl[[1]])[2],
      byrow=FALSE)
    
    
    grobs <- list()
    grobs[[1]]=rasterGrob(images_cl[[1]], width=unit(image_w_h[1],"npc"), height=unit(image_w_h[2],"npc"))
    height <- list()
    width <- list()
    height[[1]] <-  data.frame(height = nrow(images_cl[[1]]))
    width[[1]] <- data.frame(width = ncol(images_cl[[1]]))
    height <- bind_rows(height)
    width <- bind_rows(width)
    
    
    images_tibble <- tibble(sample=factor(c('figure1')), grob=grobs)
    images_tibble$height <- height$height
    images_tibble$width <- width$width
    
    qq_histo<-ggplot()+
      geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)
    
  }
  
  
  #color: c('yellow','red','blue','green','purple','grey','white')
  dat=pltdat

  
  counts<-dat[,features]
  
  weightedcolors<-weighted_colors_fun(counts=counts+offset,features=features,features_colors=features_colors,max_color_quantile=max_color_quantile)
  dat[,'color_blend']=weightedcolors
  
  if(length(features_colors)==length(features)){
    coolused<-coolist[features_colors]
  }else{
    coolused<-coolist[seq(1,length(features))]
  }
  
  
  
  
  
  colors_list <- list()
  value_list <- list()
  pal_list <- list()
  temp_geom <- list()
  legend_list <- list()
  for (i in seq_along(features)) {
    coo=coolused[[i]]
    vec_color=NULL
    for(ff in 1:dim(coo)[1]){
      vec_color<-c(vec_color,rgb(red=coo[ff,1],green=coo[ff,2],blue=coo[ff,3],alpha=coo[ff,4], maxColorValue = 1))
    }
    
    #pal_list[[i]] <- palette_scp(dat[, names(colors)[i]], type = "continuous", NA_color = NA, NA_keep = FALSE, matched = FALSE, palcolor = c(lower_col, colors[i]),palette = palette,palette_list=palette_list)
    pal_list[[i]] <-vec_color
    
    value_list[[i]] <- seq(min(dat[, features[i]], na.rm = TRUE), max(dat[, features[i]], na.rm = TRUE), length.out = dim(coo)[1])
    
    
    temp_geom[[i]] <- list(
      suppressWarnings(geom_point(data = dat, mapping = aes(x = .data[["X"]], y = .data[["Y"]], color = .data[[features[i]]]))),
      scale_color_gradientn(
        # colours = pal_list[[i]],
        # values = scales::rescale(value_list[[i]]), na.value = bg_color,
        colours = vec_color,
        values = scales::rescale(value_list[[i]]), na.value = bg_color,
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.hjust = 0)
      ),
      new_scale_color()
    )
    legend_list[[i]] <- get_legend(ggplot(dat, aes(x = X, y = Y)) +
                                     temp_geom[[i]])
    
    
  }
  
  
  
  # plot(legend_list[[1]])
  # length(legend_list)
  

  
  
  
  
  dat["color_value"] <- colSums(col2rgb(dat[, "color_blend"]))
  dat[rowSums(is.na(dat[, features])) == length(features), "color_value"] <- NA
  dat <- dat[order(dat[, "color_value"], decreasing = TRUE, na.last = FALSE), ]
  dat[rowSums(is.na(dat[, features])) == length(features), "color_blend"] <- bg_color
  
  
  dat[, "features"] <- paste(features, collapse = "|")
  
  if(!is.null(image_path)){
    p <- qq_histo
  }else{
    p <- ggplot(dat)
  }

  
  
  #p<-p + suppressWarnings(geom_point(
  p<-p+suppressWarnings(geom_point(data=dat,
    mapping = aes(x = .data[["X"]], y = .data[["Y"]], color = .data[["color_blend"]]),shape = 19, stroke = 0.5)) +
    scale_color_identity() +
    theme(legend.position ='none')+
    new_scale_color()+
    theme_classic()
  
  if(!is.null(image_path)){
    p <- p+    
      coord_cartesian(expand = FALSE) +
      labs(x='',y='',title=title)+
      xlim(0,max(width[[1]]))+
      ylim(max(height[[1]]),0)+
      theme_bw()+
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank())
  }
  

  
  
  
  
  if(label){
    
    #label_df <- reshape2::melt(p$data, measure.vars = features)
    label_df <- reshape2::melt(dat, measure.vars = features)
    
    label_df <- label_df %>%
      dplyr::group_by(variable) %>%
      dplyr::filter(value >= quantile(value, 0.985, na.rm = TRUE) & value <= quantile(value, 0.99, na.rm = TRUE)) %>%
      dplyr::reframe(x = median(.data[["X"]]), y = median(.data[["Y"]])) %>%
      as.data.frame()
    colnames(label_df)[1] <- "label"
    label_df <- label_df[!is.na(label_df[, "label"]), , drop = FALSE]
    label_df[, "rank"] <- seq_len(nrow(label_df))
    
    p <- p + geom_point(
      data = label_df, mapping = aes(x = .data[["x"]], y = .data[["y"]]),
      color = label_point_color, size = label_point_size
    ) + geom_text_repel(
      data = label_df, aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label"]], color = .data[["label"]]),
      fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
      point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
      color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE, show.legend = FALSE
    )
  }
  
  
  
  legend_nrow <- min(ceiling(sqrt(length(legend_list))), 3)
  total <- length(legend_list)
  leg_list <- list()
  n <- 1
  for (i in 1:total) {
    if (i == 1 || is.null(leg)) {
      leg <- legend_list[[i]]
    } else {
      leg <- cbind(leg, legend_list[[i]])
    }
    if (i %% legend_nrow == 0) {
      leg_list[[n]] <- leg
      leg <- NULL
      n <- n + 1
    }
    if (i %% legend_nrow != 0 && i == total) {
      ncol_insert <- dim(leg_list[[n - 1]])[2] - dim(leg)[2]
      for (col_insert in 1:ncol_insert) {
        leg <- gtable_add_cols(leg, sum(leg_list[[n - 1]]$widths) / ncol_insert, -1)
      }
      leg_list[[n]] <- leg
    }
  }
  legend <- do.call(rbind, leg_list)
  #plot(legend)
  
  legend2=legend
  

  
  
  grob <- ggplotGrob(p)
  grob <- gtable_add_cols(grob, sum(legend$widths), -1)
  grob <- gtable_add_grob(grob, legend, t = grob$layout[grepl(pattern = "panel", x = grob$layout$name), "t"], l = dim(grob)[2])
  
  gout<-plot_grid(grob)
  
  #gout
  
  return(list(gout=gout,p=p,legend=legend))
}





weighted_colors_fun<-function(counts,features,features_colors,max_color_quantile=0.98){
  
  weights_mat<-matrix(0,nrow=dim(counts)[1],ncol=length(features))
  colors_mat<-array(0,dim=c(dim(counts)[1],length(features),4))
  
  if(length(features_colors)==length(features)){
    coolist_rgb_used<-coolist_rgb[features_colors]
  }else{
    coolist_rgb_used<-coolist_rgb[seq(1,length(features))]
  }
  
  
  for(i in 1:length(features)){
    x=counts[, features[i]]
    min_color_intensity = min(x)
    
    max_color_intensity = min(c(quantile(x, max_color_quantile),Inf))
    
    
    min_value=min_color_intensity
    max_value=max_color_intensity
    
    aa=(clipR(x, min_value, max_value) - min_value) / (max_value - min_value)
    
    weights=clipR(x / (max_color_intensity + 1e-10), 0, 1)
    weights[x<min_color_intensity]=0
    
    
    
    xout<-colorRamp(c(coolist_rgb_used[[i]][1],coolist_rgb_used[[i]][length(coolist_rgb_used[[i]])]),alpha = TRUE)(aa) 
    
    colors_mat[,i,]=xout/255
    weights_mat[,i]=weights
  }
  
  
  
  
  colors_ryb =  array(0,dim=c(dim(weights_mat),3))
  dim(colors_ryb)
  
  
  for(i in 1:dim(colors_mat)[1]){
    colors_ryb[i,,]<-RGB2RYB(colors_mat[i,,1:3])
  }
  
  
  #kernel_weights<-array(0,dim=c(dim(dat)[1],length(features),1))
  kernel_weights=weights_mat^2
  
  tempp<-colors_ryb
  for(i in 1:length(features)){
    tempp[,i,]<-tempp[,i,]*kernel_weights[,i]
  }
  
  weighted_colors_ryb<-apply(tempp,c(1,3),sum,na.rm = TRUE)/rowSums(kernel_weights)
  
  
  
  #weighted_colors_ryb[is.na(weighted_colors_ryb)]<-1/3
  
  weighted_colors = array(0,dim=c(dim(weights_mat)[1],4))
  
  dim(weighted_colors)
  weighted_colors[,1:3]<-RYB2RGB(weighted_colors_ryb)
  
  weighted_colors[,4]<-apply(colors_mat[,,4],1,max)
  weighted_colors_rgb<-apply(weighted_colors,1,function(x){
    rgb(x[1],x[2],x[3],x[4])
  })
  
  
  return(weighted_colors_rgb)
}

#spot.theme#########
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 12, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 12)),
  #theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15))
  #scale_x_discrete(position = "top")
  )


#heatmap word cloud########
#https://jokergoo.github.io/2020/05/31/word-cloud-as-heatmap-annotation/
word_cloud_grob = function(text, fontsize, 
                           line_space = unit(4, "pt"), word_space = unit(4, "pt"), max_width = unit(80, "mm"), 
                           col = function(fs) circlize::rand_color(length(fs), luminosity = "dark"),
                           test = FALSE) { # width in mm
  
  if(length(text) != length(fontsize)) {
    stop("`text` and `fontsize` should the same length.")
  }
  
  od = order(fontsize, decreasing = TRUE)
  text = text[od]
  fontsize = fontsize[od]
  
  if(Sys.info()["sysname"] == "Darwin" && dev.interactive()) {
    ComplexHeatmap:::dev.null()
    on.exit(ComplexHeatmap:::dev.off2())
  }
  
  n = length(text)
  text_gb_lt = lapply(seq_len(n), function(i) textGrob(text[i], gp = gpar(fontsize = fontsize[i])))
  text_width = vapply(text_gb_lt, function(gb) convertWidth(grobWidth(gb), "mm", valueOnly = TRUE), 0)
  text_height = vapply(text_gb_lt, function(gb) convertHeight(grobHeight(gb), "mm", valueOnly = TRUE), 0)
  
  if(is.unit(line_space)) line_space = convertHeight(line_space, "mm", valueOnly = TRUE)
  if(is.unit(word_space)) word_space = convertWidth(word_space, "mm", valueOnly = TRUE)
  
  x = numeric(n)
  y = numeric(n)
  current_line_height = 0
  current_line_width = 0
  
  # the first text
  current_line_height = text_height[1]
  current_line_width = text_width[1]
  x[1] = 0
  y[1] = 0
  
  w = text_width[1]
  h = text_height[1]
  
  if(is.unit(max_width)) {
    max_width = convertWidth(max_width, "mm", valueOnly = TRUE)
  } 
  
  for(i in seq_len(n)[-1]) {
    # the next text can be put on the same line
    if(current_line_width + text_width[i] <= max_width) {
      x[i] = current_line_width + word_space
      y[i] = y[i-1] # same as previous one
      current_line_width = x[i] + text_width[i]
      w = max(w, current_line_width)
      h = max(h, y[i] + text_height[i])
    } else { # the next text need to be put on the next line
      x[i] = 0
      y[i] = current_line_height + line_space
      current_line_width = text_width[i]
      current_line_height = y[i] + text_height[i]
      w = max(w, current_line_width)
      h = max(h, current_line_height)
    }
  }
  
  if(is.character(col) || is.numeric(col)) {
    if(length(col) == 1) col = rep(col, n)
    col_fun = function(fontsize) return(col)
  } else if(is.function(col)) {
    col_fun = col
  } else {
    stop("`col` can only be a function or a character vector.")
  }
  
  if(test) {
    gl = gList(
      rectGrob(),
      textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col_fun(fontsize)), 
               default.units = "mm", just = c(0, 0)),
      rectGrob(x = x, y = y, width = text_width, height = text_height, default.units = "mm", just = c(0, 0))
      
    )
  } else {
    gl = gList(
      textGrob(text, x = x, y = y, gp = gpar(fontsize = fontsize, col = col_fun(fontsize)), 
               default.units = "mm", just = c(0, 0))
    )
  }
  
  gb = gTree(children = gl, cl = "word_cloud", vp = viewport(width = unit(w, "mm"), height = unit(h, "mm")))
  return(gb)
}

widthDetails.word_cloud = function(x) {
  x$vp$width
}

heightDetails.word_cloud = function(x) {
  x$vp$height
}

scale_fontsize = function(x, rg = c(1, 30), fs = c(4, 16)) {
  k = (fs[2] - fs[1])/(rg[2] - rg[1]) 
  b = fs[2] - k*rg[2]
  y = k*x + b
  y[y < fs[1]] = fs[1]
  y[y > fs[2]] = fs[2]
  round(y)
}

panel_fun = function(index, nm) {
  # background
  grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
  # border
  grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"), 
             default.units = "npc")
  gb = gbl[[nm]]
  # a viewport within the margins
  pushViewport(viewport(x = margin/2, y = margin/2, 
                        width = grobWidth(gb), height = grobHeight(gb),
                        just = c("left", "bottom")))
  grid.draw(gb)
  popViewport()
}




#FONT#########
library(extrafont)
# font_import()
# loadfonts(device = "all")

#Liana########
liana_get_freq <- function(liana_res){
  liana_res %>%
    dplyr::group_by(source, target) %>%
    dplyr::summarise(freq = n(), .groups = 'keep') %>%
    pivot_wider(id_cols = source,
                names_from = target,
                values_from = freq,
                values_fill = 0) %>%
    dplyr::arrange(source) %>%
    dplyr::ungroup() %>%
    as.data.frame() %>%
    column_to_rownames('source') %>%
    as.matrix()
}




#reproduce inferCNV heatmap###########
#adopted from https://jokergoo.github.io/2020/10/29/make-genome-scale-heatmap/
average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}


HEATPRE_CNV_WT<-function(genes_loc,Heat,w = 1e5){
  #genes_loc:seqnames, start and end
  #Heat: cells by genes
  
  chr_gr=GRanges(seqnames=genes_loc$seqnames,ranges = IRanges(genes_loc[, 'start'] , genes_loc[, 'end']))
  chr_window = makeWindows(chr_gr, w = w)

  genes<-colnames(Heat)
  gr1 = GRanges(seqnames = genes_loc[genes,'seqnames'], ranges = IRanges(genes_loc[genes,'start'], genes_loc[genes,'end']))
  
  num_mat = average_in_window(chr_window, gr1, t(Heat))
  
  chr = as.vector(seqnames(chr_window))
  chr_level = paste0("chr", 1:22)
  chr = factor(chr, levels = chr_level)
  num_mat=t(num_mat)
  colnames(num_mat)<-chr
  rownames(num_mat)<-rownames(Heat)
  # heat[is.na(heat)]<-1, make NA white color?
  return(
    list(
      num_mat=(num_mat),
      chr=chr
    )
  )
  
  
  # heat=heatpre$num_mat
  # chr=heatpre$chr
  # heat2<-heat
  # #heat2<-heat2[,-which(is.na(colSums(heat2)))]
  # #heat[is.na(heat)]=0
  # 
  # phee<-Heatmap(
  #   MinMax(heat2,-0.05,0.05), name = "mat",
  #   col = colorRamp2(c(-0.05, 0, 0.05), c("darkblue", "white", "darkred")),
  #   column_split = chr,
  #   heatmap_legend_param = list( title = paste(' CNV',sep=''),direction = "horizontal"),
  #   cluster_columns  = FALSE,
  #   show_row_dend = FALSE,
  #   #row_split = subcluster,
  #   cluster_rows = FALSE,
  #   show_column_dend = FALSE,
  #   show_column_names = FALSE,
  #   show_row_names = FALSE,
  #   column_title_gp = gpar(fontsize = 10), border = TRUE,
  #   column_gap = unit(0, "points"),
  #   column_title_rot = 90)
  # phee<-draw(phee,
  #            heatmap_legend_side = "top",annotation_legend_side = "top")
  # gout<-grid.grabExpr(draw(phee))
  
}

HEATPRE_CNV_WT_Mouse<-function(genes_loc,Heat,w = 1e5){
  #genes_loc:seqnames, start and end
  #Heat: cells by genes
  
  chr_gr=GRanges(seqnames=genes_loc$seqnames,ranges = IRanges(genes_loc[, 'start'] , genes_loc[, 'end']))
  chr_window = makeWindows(chr_gr, w = w)
  
  genes<-colnames(Heat)
  gr1 = GRanges(seqnames = genes_loc[genes,'seqnames'], ranges = IRanges(genes_loc[genes,'start'], genes_loc[genes,'end']))
  
  num_mat = average_in_window(chr_window, gr1, t(Heat))
  
  chr = as.vector(seqnames(chr_window))
  chr_level = paste0("chr", 1:19)
  chr = factor(chr, levels = chr_level)
  num_mat=t(num_mat)
  colnames(num_mat)<-chr
  rownames(num_mat)<-rownames(Heat)
  # heat[is.na(heat)]<-1, make NA white color?
  return(
    list(
      num_mat=(num_mat),
      chr=chr
    )
  )
}



#Then do the following
# max(num_mat,na.rm = TRUE)
# min(num_mat,na.rm=TRUE)
# 
# heat<-num_mat
# heat[is.na(heat)]<-1
# 
# Heatmap(
#   MinMax(heat,0.9,1.1), name = "mat",
#   col = colorRamp2(c(0.9, 1, 1.1), c("darkblue", "white", "darkred")),
#   column_split = chr,
#   cluster_columns  = FALSE,
#   show_column_dend = FALSE,
#   cluster_column_slices = FALSE,
#   column_title = "numeric matrix",
#   row_title_rot = 0, row_title_gp = gpar(fontsize = 10), border = TRUE,
#   row_gap = unit(0, "points"))