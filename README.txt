Mammola S. et al. (2026). The biodiversity we ignore: a global meta-analysis on taxonomic bias. Ambio

# Data sources:

Two datasets are provided:

-Database_meta_analysis_bias.xlsx: is the master database for analyses.

-geo.tsv: is used to generate the map in Figure 1.

# Code: 

-meta_analysis.R: R code to reproduce the analysis, annotated throughout. The session parameters used to run the code are as follows:

R version 4.3.0 (2023-04-21)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Rome
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] gtable_0.3.6        TMB_1.9.5           bayestestR_0.13.2   ggplot2_3.5.2      
 [5] insight_0.20.3      lattice_0.21-8      mathjaxr_1.6-0      numDeriv_2016.8-1.1
 [9] vctrs_0.6.5         tools_4.3.0         generics_0.1.3      datawizard_0.10.0  
[13] tibble_3.2.1        proxy_0.4-27        fansi_1.0.6         pacman_0.5.1       
[17] pkgconfig_2.0.3     Matrix_1.6-0        KernSmooth_2.23-20  DHARMa_0.4.6       
[21] RColorBrewer_1.1-3  ggridges_0.5.6      readxl_1.4.3        lifecycle_1.0.4    
[25] compiler_4.3.0      farver_2.1.1        terra_1.7-39        codetools_0.2-19   
[29] class_7.3-21        pillar_1.9.0        nloptr_2.0.3        tidyr_1.3.0        
[33] MASS_7.3-58.4       classInt_0.4-9      metadat_1.2-0       boot_1.3-28.1      
[37] nlme_3.1-162        tidyselect_1.2.0    ggh4x_0.3.1         performance_0.12.2 
[41] mvtnorm_1.2-4       sf_1.0-14           dplyr_1.1.2         purrr_1.0.1        
[45] forcats_1.0.0       splines_4.3.0       cowplot_1.1.1       rnaturalearth_1.0.1
[49] grid_4.3.0          cli_3.6.2           metafor_4.4-0       magrittr_2.0.3     
[53] utf8_1.2.4          e1071_1.7-13        scales_1.4.0        estimability_1.4.1 
[57] httr_1.4.6          emmeans_1.8.7       lme4_1.1-34         cellranger_1.1.0   
[61] coda_0.19-4         glmmTMB_1.1.7       parameters_0.21.6   gghalves_0.1.4     
[65] rlang_1.1.3         Rcpp_1.0.12         xtable_1.8-4        glue_1.7.0         
[69] DBI_1.1.3           rstudioapi_0.15.0   minqa_1.2.5         jsonlite_1.8.8     
[73] R6_2.5.1            units_0.8-2 
