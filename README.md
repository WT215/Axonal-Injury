# Axonal-Injury
This resource provides the R code to reproduce key results described in Clements, Tang, Baronik, Simpson Ragdale et al. "Axonal injury is a targetable driver of glioblastoma progression".

## Getting started
1. Load all R code/RData from ./utils;
2. ./RData contains RData files needed for reproducing corresponding figures from the paper;
3. ./R_Figure2 contains code for reproducing main figure 2;
4. ./R_SUPFIG2 contains code for reproducing supplementary figure 2.
5. Please modify the path in the R code where there is needed.


## Session information

```
R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/London
tzcode source: system (glibc)

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
  [1] DESeq2_1.46.0               biomaRt_2.62.1             
  [3] babelgene_22.9              topGO_2.58.0               
  [5] SparseM_1.84-2              GO.db_3.20.0               
  [7] graph_1.84.1                genefilter_1.88.0          
  [9] org.Hs.eg.db_3.20.0         org.Mm.eg.db_3.20.0        
 [11] AnnotationDbi_1.68.0        RANN_2.6.2                 
 [13] moranfast_1.0               DelayedArray_0.32.0        
 [15] SparseArray_1.6.2           S4Arrays_1.6.0             
 [17] AUCell_1.28.0               plotrix_3.8-4              
 [19] pbmcapply_1.5.1             SpatialPack_0.4-1          
 [21] fastmatrix_0.5-7721         doSNOW_1.0.20              
 [23] snow_0.4-4                  iterators_1.0.14           
 [25] randomForest_4.7-1.2        extrafont_0.19             
 [27] scales_1.3.0                ggnewscale_0.5.1           
 [29] PBSmapping_2.74.1           EnrichedHeatmap_1.36.0     
 [31] directlabels_2024.1.21      ggplotify_0.1.2            
 [33] mclust_6.1.1                doRNG_1.8.6.1              
 [35] rngtools_1.5.2              foreach_1.5.2              
 [37] vcd_1.4-13                  scuttle_1.16.0             
 [39] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
 [41] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
 [43] IRanges_2.40.1              S4Vectors_0.44.0           
 [45] MatrixGenerics_1.18.1       matrixStats_1.5.0          
 [47] BioQC_1.34.0                Biobase_2.66.0             
 [49] BiocGenerics_0.52.0         future.apply_1.11.3        
 [51] pbapply_1.7-2               future_1.34.0              
 [53] abind_1.4-8                 circlize_0.4.16            
 [55] ComplexHeatmap_2.22.0       ggVennDiagram_1.5.2        
 [57] RColorBrewer_1.1-3          viridis_0.6.5              
 [59] viridisLite_0.4.2           crayon_1.5.3               
 [61] cowplot_1.1.3               plyr_1.8.9                 
 [63] rdist_0.0.5                 readbitmap_0.1.5           
 [65] igraph_2.1.4                irlba_2.3.5.1              
 [67] Rtsne_0.17                  uwot_0.2.3                 
 [69] Matrix_1.7-2                bayNorm_1.24.0             
 [71] fitdistrplus_1.2-2          survival_3.8-3             
 [73] MASS_7.3-65                 lubridate_1.9.4            
 [75] forcats_1.0.0               stringr_1.5.1              
 [77] dplyr_1.1.4                 purrr_1.0.4                
 [79] tidyr_1.3.1                 tibble_3.2.1               
 [81] tidyverse_2.0.0             readxl_1.4.4               
 [83] fgsea_1.32.2                presto_1.0.0               
 [85] data.table_1.17.0           Rcpp_1.0.14                
 [87] gridGraphics_0.5-1          gridExtra_2.3              
 [89] ggpubr_0.6.0                gtable_0.3.6               
 [91] ggrepel_0.9.6               ggplot2_3.5.1              
 [93] SeuratData_0.2.2.9002       Seurat_5.2.1               
 [95] SeuratObject_5.0.2          sp_2.2-0                   
 [97] trend_1.1.6                 reshape2_1.4.4             
 [99] BiocParallel_1.40.0         magrittr_2.0.3             
[101] OmnipathR_3.14.0            reticulate_1.41.0          
[103] readr_2.1.5                 rhdf5_2.50.2               
[105] rjson_0.2.23               

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2         vroom_1.6.5              
  [3] GSEABase_1.68.0           progress_1.2.3           
  [5] tiff_0.1-12               goftest_1.2-3            
  [7] Biostrings_2.74.1         vctrs_0.6.5              
  [9] spatstat.random_3.3-2     digest_0.6.37            
 [11] png_0.1-8                 shape_1.4.6.1            
 [13] deldir_2.0-4              parallelly_1.42.0        
 [15] magick_2.8.5              httpuv_1.6.15            
 [17] withr_3.0.2               xfun_0.51                
 [19] memoise_2.0.1             systemfonts_1.2.1        
 [21] ragg_1.3.3                zoo_1.8-13               
 [23] GlobalOptions_0.1.2       R.oo_1.27.0              
 [25] Formula_1.2-5             prettyunits_1.2.0        
 [27] KEGGREST_1.46.0           promises_1.3.2           
 [29] httr_1.4.7                rstatix_0.7.2            
 [31] globals_0.16.3            rhdf5filters_1.18.0      
 [33] rstudioapi_0.17.1         UCSC.utils_1.2.0         
 [35] miniUI_0.1.1.1            generics_0.1.3           
 [37] curl_6.2.1                zlibbioc_1.52.0          
 [39] polyclip_1.10-7           GenomeInfoDbData_1.2.13  
 [41] quadprog_1.5-8            xtable_1.8-4             
 [43] doParallel_1.0.17         evaluate_1.0.3           
 [45] BiocFileCache_2.14.0      hms_1.1.3                
 [47] filelock_1.0.3            colorspace_2.1-1         
 [49] ROCR_1.0-11               spatstat.data_3.1-4      
 [51] lmtest_0.9-40             later_1.4.1              
 [53] lattice_0.22-6            spatstat.geom_3.3-5      
 [55] scattermore_1.2           XML_3.99-0.18            
 [57] RcppAnnoy_0.0.22          pillar_1.10.1            
 [59] nlme_3.1-167              compiler_4.4.3           
 [61] beachmat_2.22.0           RSpectra_0.16-2          
 [63] stringi_1.8.4             tensor_1.5               
 [65] extraDistr_1.10.0         locfit_1.5-9.12          
 [67] bit_4.5.0.1               fastmatch_1.1-6          
 [69] textshaping_1.0.0         codetools_0.2-20         
 [71] GetoptLong_1.0.5          plotly_4.10.4            
 [73] mime_0.12                 splines_4.4.3            
 [75] fastDummies_1.7.5         dbplyr_2.5.0             
 [77] sparseMatrixStats_1.18.0  cellranger_1.1.0         
 [79] Rttf2pt1_1.3.12           knitr_1.49               
 [81] blob_1.2.4                clue_0.3-66              
 [83] fs_1.6.5                  listenv_0.9.1            
 [85] checkmate_2.3.2           DelayedMatrixStats_1.28.1
 [87] logger_0.4.0              ggsignif_0.6.4           
 [89] statmod_1.5.0             tzdb_0.4.0               
 [91] pkgconfig_2.0.3           tools_4.4.3              
 [93] cachem_1.1.0              RSQLite_2.3.9            
 [95] rvest_1.0.4               DBI_1.2.3                
 [97] fastmap_1.2.0             rmarkdown_2.29           
 [99] ica_1.0-3                 broom_1.0.7              
[101] patchwork_1.3.0           dotCall64_1.2            
[103] carData_3.0-5             farver_2.1.2             
[105] yaml_2.3.10               bmp_0.3                  
[107] cli_3.6.4                 lifecycle_1.0.4          
[109] backports_1.5.0           annotate_1.84.0          
[111] timechange_0.3.0          ggridges_0.5.6           
[113] progressr_0.15.1          limma_3.62.2             
[115] jsonlite_1.9.1            edgeR_4.4.2              
[117] RcppHNSW_0.6.0            bit64_4.6.0-1            
[119] yulab.utils_0.2.0         spatstat.utils_3.1-2     
[121] zip_2.3.2                 spatstat.univar_3.1-2    
[123] R.utils_2.13.0            lazyeval_0.2.2           
[125] shiny_1.10.0              htmltools_0.5.8.1        
[127] sctransform_0.4.1         rappdirs_0.3.3           
[129] glue_1.8.0                spam_2.11-1              
[131] httr2_1.1.0               XVector_0.46.0           
[133] jpeg_0.1-10               extrafontdb_1.0          
[135] R6_2.6.1                  labeling_0.4.3           
[137] cluster_2.1.8             Rhdf5lib_1.28.0          
[139] tidyselect_1.2.1          xml2_1.3.7               
[141] car_3.1-3                 munsell_0.5.1            
[143] KernSmooth_2.23-26        htmlwidgets_1.6.4        
[145] rlang_1.1.5               spatstat.sparse_3.1-0    
[147] spatstat.explore_3.3-4    Cairo_1.6-2   
```
