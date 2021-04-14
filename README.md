# Heterogeneity in the response of different subtypes of Drosophila melanogaster midgut cells to viral infections

### Description
This page contains all R scripts used for the analysis of viral infections on public single-cell RNA sequencing (scRNA-seq) data from Drosophila melanogaster enteroendocrine midgut cells.

### Overview
Raw data was downloaded from the sequence read archives (SRA) portal from NCBI (https://www.ncbi.nlm.nih.gov/sra; BioProject accession numbers: PRJNA547484 and PRJNA493298). Gene count matrices were generated with cellranger v3.1.0 using a concatenated reference genome containing the D. melanogaster, Thika virus (TV), Drosophila melanogaster Nora virus (DMelNV) genomic sequences and Drosophila C virus (DCV). Downstream analyses were performed with the R scripts present in the current project. 

#### Contents
##### source_me.R
Load a few objects and contains all functions used in the analysis.

##### pre_processing.R
Filtering, cluster generation and cell type annotation.

##### infection_dynamics.R
Determination of infected cells and ifection dynamics analysis.

##### DE_glmGamPoi_analysis.R
Differential expression analysis with glmGamPoi.

##### graph_glmGamPoi_DE.R
Network and reactome analysis.

##### gene_matrices/
Contains count matrices, as well as a list of DEGs from bulk RNA-seq data from DMelNV-infected flies.

##### network/
Contains the D. melanogaster interactome files.

##### name_id_symbol_list
File to convert gene IDs to symbol.

#### References
##### cellranger
10X Genomics, CA, USA

https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome

##### Seurat
Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W.M., Hao, Y., Stoeckius, M., Smibert, P. and Satija, R. (2019). Comprehensive integration of single-cell data. Cell 177, 1888-1902.

##### glmGamPoi
Ahlmann-Eltze, C., Huber, W. (2020) glmGamPoi: Fitting Gamma-Poisson Generalized Linear Models on Single Cell Count Data. Bioinformatics Dec 9, btaa1009. 

##### ReactomePA
Yu, G., He, QY. (2016) ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Mol. BioSyst. 12, 477-479. 

##### scRNA-seq data
Guo, X., Yin, C., Yang, F., Zhang, Y., Huang, H., Wang, J., Deng, B., Cai T., Rao, Y., and Xi, R. (2019). The cellular diversity and transcription factor code of Drosophila enteroendocrine cells. Cell Rep. 29, 4172-4185.

Hung, R.J., Hu, Y., Kirchner, R., Liu, Y., Xu, C., Comjean, A., Tattikota, S.G., Li, F., Song, W., Sui, S.H. and Perrimon, N. (2020). A cell atlas of the adult Drosophila midgut. Proc. Natl. Acad. Sci. USA. 21, 1514-1523.

##### bulk RNA-seq data
Lopez, W., Page, A.M., Carlson, D.J., Ericson, B.L., Cserhati, M.F., Guda, C. and Carlson, K.A. (2018). Analysis of immune-related genes during Nora virus infection of Drosophila melanogaster using next generation sequencing. AIMS Microbiol. 4, 123.

This data is available under the https://creativecommons.org/licenses/by/4.0/ license. Minor changes to the original file format were made.

##### D. melanogaster interactome
Ding, X.B., Jin, J., Tao, Y.T., Guo, W.P., Ruan, L., Yang, Q.L., Chen, P.C., Yao, H., Zhang, H.B. and Chen, X. (2020). Predicted Drosophila Interactome Resource and web tool for functional interpretation of differentially expressed genes. Database 2020, baaa005.

http://drosophila.biomedtzc.cn/#/download

```
> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] writexl_1.3.1               ggrepel_0.9.1               DropletUtils_1.10.3         SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
 [6] Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.2         IRanges_2.24.1              S4Vectors_0.28.1           
[11] BiocGenerics_0.36.0         MatrixGenerics_1.2.0        matrixStats_0.57.0          ReactomePA_1.34.0           glmGamPoi_1.2.0            
[16] broom_0.7.3                 car_3.0-10                  carData_3.0-4               data.table_1.13.6           forcats_0.5.0              
[21] stringr_1.4.0               purrr_0.3.4                 readr_1.4.0                 tidyr_1.1.2                 tidyverse_1.3.0            
[26] igraph_1.2.6                ggforce_0.3.2               agricolae_1.3-3             gridExtra_2.3               ggpubr_0.4.0               
[31] reshape2_1.4.4              limma_3.46.0                sctransform_0.3.2           ggplot2_3.3.3               viridis_0.5.1              
[36] viridisLite_0.3.0           tibble_3.0.5                dplyr_1.0.3                 Seurat_3.2.3               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.1            scattermore_0.7           R.methodsS3_1.8.1         bit64_4.0.5               irlba_2.3.3               DelayedArray_0.16.0      
  [7] R.utils_2.10.1            rpart_4.1-15              RCurl_1.98-1.2            generics_0.1.0            cowplot_1.1.1             RSQLite_2.2.2            
 [13] shadowtext_0.0.7          RANN_2.6.1                combinat_0.0-8            future_1.21.0             bit_4.0.4                 enrichplot_1.10.2        
 [19] spatstat.data_1.7-0       xml2_1.3.2                lubridate_1.7.9.2         httpuv_1.5.5              assertthat_0.2.1          hms_1.0.0                
 [25] promises_1.1.1            dbplyr_2.1.0              readxl_1.3.1              DBI_1.1.1                 htmlwidgets_1.5.3         ellipsis_0.3.1           
 [31] backports_1.2.1           deldir_0.2-9              sparseMatrixStats_1.2.0   vctrs_0.3.6               ROCR_1.0-11               abind_1.4-5              
 [37] cachem_1.0.2              withr_2.4.1               checkmate_2.0.0           goftest_1.2-2             cluster_2.1.0             DOSE_3.16.0              
 [43] lazyeval_0.2.2            crayon_1.4.0              edgeR_3.32.1              pkgconfig_2.0.3           tweenr_1.0.1              nlme_3.1-149             
 [49] rlang_0.4.10              globals_0.14.0            questionr_0.7.4           lifecycle_0.2.0           miniUI_0.1.1.1            modelr_0.1.8             
 [55] rsvd_1.0.3                cellranger_1.1.0          polyclip_1.10-0           lmtest_0.9-38             graph_1.68.0              Matrix_1.2-18            
 [61] Rhdf5lib_1.12.0           zoo_1.8-8                 reprex_1.0.0              ggridges_0.5.3            png_0.1-7                 bitops_1.0-6             
 [67] R.oo_1.24.0               KernSmooth_2.23-17        rhdf5filters_1.2.0        blob_1.2.1                DelayedMatrixStats_1.12.2 qvalue_2.22.0            
 [73] parallelly_1.23.0         rstatix_0.6.0             ggsignif_0.6.0            klaR_0.6-15               reactome.db_1.74.0        beachmat_2.6.4           
 [79] scales_1.1.1              memoise_2.0.0             graphite_1.36.0           magrittr_2.0.1            plyr_1.8.6                ica_1.0-2                
 [85] zlibbioc_1.36.0           compiler_4.0.3            scatterpie_0.1.5          dqrng_0.2.1               RColorBrewer_1.1-2        fitdistrplus_1.1-3       
 [91] cli_2.3.0                 XVector_0.30.0            listenv_0.8.0             patchwork_1.1.1           pbapply_1.4-3             MASS_7.3-53              
 [97] mgcv_1.8-33               tidyselect_1.1.0          stringi_1.5.3             highr_0.8                 GOSemSim_2.16.1           locfit_1.5-9.4           
[103] grid_4.0.3                fastmatch_1.1-0           tools_4.0.3               future.apply_1.7.0        rio_0.5.16                rstudioapi_0.13          
[109] foreign_0.8-80            farver_2.0.3              Rtsne_0.15                ggraph_2.0.4              digest_0.6.27             rvcheck_0.1.8            
[115] BiocManager_1.30.10       shiny_1.6.0               Rcpp_1.0.6                scuttle_1.0.4             later_1.1.0.1             RcppAnnoy_0.0.18         
[121] httr_1.4.2                AnnotationDbi_1.52.0      colorspace_2.0-0          rvest_0.3.6               fs_1.5.0                  tensor_1.5               
[127] reticulate_1.18           splines_4.0.3             uwot_0.1.10               spatstat.utils_1.20-2     graphlayouts_0.7.1        plotly_4.9.3             
[133] xtable_1.8-4              jsonlite_1.7.2            AlgDesign_1.2.0           spatstat_1.64-1           tidygraph_1.2.0           R6_2.5.0                 
[139] pillar_1.4.7              htmltools_0.5.1           mime_0.9                  glue_1.4.2                fastmap_1.0.1             BiocParallel_1.24.1      
[145] codetools_0.2-16          fgsea_1.16.0              lattice_0.20-41           curl_4.3                  leiden_0.3.7              zip_2.1.1                
[151] GO.db_3.12.1              openxlsx_4.2.3            survival_3.2-7            munsell_0.5.0             DO.db_2.9                 rhdf5_2.34.0             
[157] GenomeInfoDbData_1.2.4    HDF5Array_1.18.0          labelled_2.7.0            haven_2.3.1               gtable_0.3.0
```
