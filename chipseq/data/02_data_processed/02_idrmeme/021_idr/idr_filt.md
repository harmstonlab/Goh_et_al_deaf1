idr_filt
================
Qian Hui TAN
2/14/23

- <a href="#setup" id="toc-setup"><span
  class="toc-section-number">1</span> Setup</a>
  - <a href="#load-libraries" id="toc-load-libraries"><span
    class="toc-section-number">1.1</span> Load libraries</a>
  - <a href="#load-data" id="toc-load-data"><span
    class="toc-section-number">1.2</span> Load data</a>
- <a href="#data-wrangling" id="toc-data-wrangling"><span
  class="toc-section-number">2</span> Data wrangling</a>
- <a href="#exporting" id="toc-exporting"><span
  class="toc-section-number">3</span> Exporting</a>

In this notebook, we filter the IDR peaks to only contain high-quality
peaks (IDR 0 - 0.05).

# Setup

## Load libraries

``` r
suppressPackageStartupMessages({
    # Load libraries
    library(ChIPseeker)
    library(GenomicRanges)
})
```

## Load data

``` r
deaf1_bed <- readPeakFile("deaf1_idr_all.bed")
```

``` r
length(deaf1_bed)
```

    [1] 1550

# Data wrangling

`deaf1_bed` has a few columns that are named V4 to V14. We rename them
so they make sense:

``` r
head(deaf1_bed)
```

    GRanges object with 6 ranges and 11 metadata columns:
          seqnames              ranges strand |          V4        V5          V6
             <Rle>           <IRanges>  <Rle> | <character> <integer> <character>
      [1]     chr1 166378462-166382216      * |           .       845           .
      [2]     chr7 126427901-126429745      * |           .       826           .
      [3]     chr5 143527422-143528800      * |           .       836           .
      [4]    chr15   95790481-95792623      * |           .       808           .
      [5]     chr8     8688634-8691843      * |           .       854           .
      [6]     chr7     4792394-4793284      * |           .       826           .
                 V7        V8        V9       V10       V11       V12       V13
          <numeric> <numeric> <integer> <integer> <numeric> <integer> <integer>
      [1]      1.99      2.04 166378501 166382203  186.8078 166378461 166382216
      [2]      1.90      1.99 126428279 126429503  131.4282 126427900 126429745
      [3]      1.94      2.02 143527539 143528761  119.3424 143527421 143528800
      [4]      1.82      1.95  95790480  95792531  109.1668  95790507  95792623
      [5]      2.03      2.06   8689367   8691612  109.1657   8688633   8691843
      [6]      1.89      1.99   4792415   4793271   96.2429   4792393   4793284
                V14
          <numeric>
      [1]   245.758
      [2]   164.870
      [3]   149.103
      [4]   109.693
      [5]   144.287
      [6]   102.300
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

We rename the metadata columns (mcols) with the correct names:

``` r
colnames(mcols(deaf1_bed)) <- c(
                  "idr_name", "scaled_idr", "idr_strand",
                   "local_idr", "global_idr",
                   "rep1_chromstart", "rep1_chromend", "rep1_signal",
                   "rep2_chromstart", "rep2_chromend", "rep2_signal"
                  )
```

``` r
head(deaf1_bed)
```

    GRanges object with 6 ranges and 11 metadata columns:
          seqnames              ranges strand |    idr_name scaled_idr  idr_strand
             <Rle>           <IRanges>  <Rle> | <character>  <integer> <character>
      [1]     chr1 166378462-166382216      * |           .        845           .
      [2]     chr7 126427901-126429745      * |           .        826           .
      [3]     chr5 143527422-143528800      * |           .        836           .
      [4]    chr15   95790481-95792623      * |           .        808           .
      [5]     chr8     8688634-8691843      * |           .        854           .
      [6]     chr7     4792394-4793284      * |           .        826           .
          local_idr global_idr rep1_chromstart rep1_chromend rep1_signal
          <numeric>  <numeric>       <integer>     <integer>   <numeric>
      [1]      1.99       2.04       166378501     166382203    186.8078
      [2]      1.90       1.99       126428279     126429503    131.4282
      [3]      1.94       2.02       143527539     143528761    119.3424
      [4]      1.82       1.95        95790480      95792531    109.1668
      [5]      2.03       2.06         8689367       8691612    109.1657
      [6]      1.89       1.99         4792415       4793271     96.2429
          rep2_chromstart rep2_chromend rep2_signal
                <integer>     <integer>   <numeric>
      [1]       166378461     166382216     245.758
      [2]       126427900     126429745     164.870
      [3]       143527421     143528800     149.103
      [4]        95790507      95792623     109.693
      [5]         8688633       8691843     144.287
      [6]         4792393       4793284     102.300
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

While this IDR file has IDR values, it has not been filtered. We know
this because we should only have 608/1550 peaks if we check the IDR
output, but here we have 1550 peaks.

To be sure, we do a quick summary of the scaled IDR values,
min(int(log2(-125IDR), 1000). Peaks with an IDR of 0 have a score of
1000, idr 0.05 have a score of int(-125log2(0.05)) = 540, and idr 1.0
has a score of 0.

``` r
summary(deaf1_bed$scaled_idr)
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      121.0   191.0   413.0   431.1   667.0   887.0 

Min score is 121 instead of 540 - that means they have not been filtered
out. We filter them out:

``` r
# Filter and keep only highly reproducible peaks:
# IDR 0 - 0.05 (scaled IDR 540 - 1000)
deaf1_filt <- deaf1_bed[deaf1_bed$scaled_idr >= 540, ]
length(deaf1_filt)
```

    [1] 608

``` r
head(deaf1_filt)
```

    GRanges object with 6 ranges and 11 metadata columns:
          seqnames              ranges strand |    idr_name scaled_idr  idr_strand
             <Rle>           <IRanges>  <Rle> | <character>  <integer> <character>
      [1]     chr1 166378462-166382216      * |           .        845           .
      [2]     chr7 126427901-126429745      * |           .        826           .
      [3]     chr5 143527422-143528800      * |           .        836           .
      [4]    chr15   95790481-95792623      * |           .        808           .
      [5]     chr8     8688634-8691843      * |           .        854           .
      [6]     chr7     4792394-4793284      * |           .        826           .
          local_idr global_idr rep1_chromstart rep1_chromend rep1_signal
          <numeric>  <numeric>       <integer>     <integer>   <numeric>
      [1]      1.99       2.04       166378501     166382203    186.8078
      [2]      1.90       1.99       126428279     126429503    131.4282
      [3]      1.94       2.02       143527539     143528761    119.3424
      [4]      1.82       1.95        95790480      95792531    109.1668
      [5]      2.03       2.06         8689367       8691612    109.1657
      [6]      1.89       1.99         4792415       4793271     96.2429
          rep2_chromstart rep2_chromend rep2_signal
                <integer>     <integer>   <numeric>
      [1]       166378461     166382216     245.758
      [2]       126427900     126429745     164.870
      [3]       143527421     143528800     149.103
      [4]        95790507      95792623     109.693
      [5]         8688633       8691843     144.287
      [6]         4792393       4793284     102.300
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

That’s strange. Why do we have 608 peaks here, but 607 from the cutoff?

Let’s take a look at the peaks that barely pass the cutoff:

``` r
deaf1_bed[deaf1_bed$scaled_idr >= 540 & deaf1_bed$scaled_idr < 542,]
```

    GRanges object with 3 ranges and 11 metadata columns:
          seqnames              ranges strand |    idr_name scaled_idr  idr_strand
             <Rle>           <IRanges>  <Rle> | <character>  <integer> <character>
      [1]     chr3   88579545-88579672      * |           .        540           .
      [2]     chr5 147860479-147860661      * |           .        541           .
      [3]     chr6   40471242-40471447      * |           .        540           .
          local_idr global_idr rep1_chromstart rep1_chromend rep1_signal
          <numeric>  <numeric>       <integer>     <integer>   <numeric>
      [1]      0.78        1.3        88579544      88579672     52.0412
      [2]      0.78        1.3       147860540     147860661     32.8254
      [3]      0.78        1.3        40471268      40471404     15.6725
          rep2_chromstart rep2_chromend rep2_signal
                <integer>     <integer>   <numeric>
      [1]        88579555      88579671     29.2992
      [2]       147860478     147860661     23.9655
      [3]        40471241      40471447     20.7535
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

We have 2 peaks with scaled_idr 540, and another 1 with scaled_idr =
541. Since they are `integers`, I’m guessing there was a rounding error.

In any case, we include both the peaks with scaled_idr = 540.

# Exporting

``` r
write.table(deaf1_filt, "deaf1_idr_filt.bed", quote = F, 
            col.names = FALSE, row.names = FALSE, sep = "\t")

saveRDS(deaf1_filt, "deaf1_idr_filt.RDS")
```

``` r
sessionInfo()
```

    R version 4.2.2 (2022-10-31)
    Platform: aarch64-apple-darwin20 (64-bit)
    Running under: macOS Ventura 13.1

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     

    other attached packages:
    [1] GenomicRanges_1.50.2 GenomeInfoDb_1.34.7  IRanges_2.32.0      
    [4] S4Vectors_0.36.1     BiocGenerics_0.44.0  ChIPseeker_1.34.1   

    loaded via a namespace (and not attached):
      [1] shadowtext_0.1.2                       
      [2] fastmatch_1.1-3                        
      [3] BiocFileCache_2.6.0                    
      [4] plyr_1.8.8                             
      [5] igraph_1.3.5                           
      [6] lazyeval_0.2.2                         
      [7] splines_4.2.2                          
      [8] BiocParallel_1.32.5                    
      [9] ggplot2_3.4.0                          
     [10] digest_0.6.31                          
     [11] yulab.utils_0.0.6                      
     [12] htmltools_0.5.4                        
     [13] GOSemSim_2.24.0                        
     [14] viridis_0.6.2                          
     [15] GO.db_3.16.0                           
     [16] fansi_1.0.4                            
     [17] magrittr_2.0.3                         
     [18] memoise_2.0.1                          
     [19] Biostrings_2.66.0                      
     [20] graphlayouts_0.8.4                     
     [21] matrixStats_0.63.0                     
     [22] enrichplot_1.18.3                      
     [23] prettyunits_1.1.1                      
     [24] colorspace_2.1-0                       
     [25] blob_1.2.3                             
     [26] rappdirs_0.3.3                         
     [27] ggrepel_0.9.2                          
     [28] xfun_0.37                              
     [29] dplyr_1.1.0                            
     [30] crayon_1.5.2                           
     [31] RCurl_1.98-1.10                        
     [32] jsonlite_1.8.4                         
     [33] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
     [34] scatterpie_0.1.8                       
     [35] ape_5.6-2                              
     [36] glue_1.6.2                             
     [37] polyclip_1.10-4                        
     [38] gtable_0.3.1                           
     [39] zlibbioc_1.44.0                        
     [40] XVector_0.38.0                         
     [41] DelayedArray_0.24.0                    
     [42] scales_1.2.1                           
     [43] DOSE_3.24.2                            
     [44] DBI_1.1.3                              
     [45] Rcpp_1.0.10                            
     [46] plotrix_3.8-2                          
     [47] viridisLite_0.4.1                      
     [48] progress_1.2.2                         
     [49] gridGraphics_0.5-1                     
     [50] tidytree_0.4.2                         
     [51] bit_4.0.5                              
     [52] httr_1.4.4                             
     [53] fgsea_1.24.0                           
     [54] gplots_3.1.3                           
     [55] RColorBrewer_1.1-3                     
     [56] ellipsis_0.3.2                         
     [57] pkgconfig_2.0.3                        
     [58] XML_3.99-0.13                          
     [59] farver_2.1.1                           
     [60] dbplyr_2.3.0                           
     [61] utf8_1.2.3                             
     [62] ggplotify_0.1.0                        
     [63] tidyselect_1.2.0                       
     [64] rlang_1.0.6                            
     [65] reshape2_1.4.4                         
     [66] AnnotationDbi_1.60.0                   
     [67] munsell_0.5.0                          
     [68] tools_4.2.2                            
     [69] cachem_1.0.6                           
     [70] cli_3.6.0                              
     [71] generics_0.1.3                         
     [72] RSQLite_2.2.20                         
     [73] evaluate_0.20                          
     [74] stringr_1.5.0                          
     [75] fastmap_1.1.0                          
     [76] yaml_2.3.7                             
     [77] ggtree_3.6.2                           
     [78] knitr_1.42                             
     [79] bit64_4.0.5                            
     [80] tidygraph_1.2.3                        
     [81] caTools_1.18.2                         
     [82] purrr_1.0.1                            
     [83] KEGGREST_1.38.0                        
     [84] ggraph_2.1.0                           
     [85] nlme_3.1-162                           
     [86] aplot_0.1.9                            
     [87] xml2_1.3.3                             
     [88] biomaRt_2.54.0                         
     [89] compiler_4.2.2                         
     [90] rstudioapi_0.14                        
     [91] filelock_1.0.2                         
     [92] curl_5.0.0                             
     [93] png_0.1-8                              
     [94] treeio_1.22.0                          
     [95] tibble_3.1.8                           
     [96] tweenr_2.0.2                           
     [97] stringi_1.7.12                         
     [98] GenomicFeatures_1.50.4                 
     [99] lattice_0.20-45                        
    [100] Matrix_1.5-3                           
    [101] vctrs_0.5.2                            
    [102] pillar_1.8.1                           
    [103] lifecycle_1.0.3                        
    [104] data.table_1.14.6                      
    [105] cowplot_1.1.1                          
    [106] bitops_1.0-7                           
    [107] patchwork_1.1.2                        
    [108] rtracklayer_1.58.0                     
    [109] qvalue_2.30.0                          
    [110] R6_2.5.1                               
    [111] BiocIO_1.8.0                           
    [112] KernSmooth_2.23-20                     
    [113] gridExtra_2.3                          
    [114] codetools_0.2-19                       
    [115] boot_1.3-28.1                          
    [116] MASS_7.3-58.2                          
    [117] gtools_3.9.4                           
    [118] assertthat_0.2.1                       
    [119] SummarizedExperiment_1.28.0            
    [120] rjson_0.2.21                           
    [121] withr_2.5.0                            
    [122] GenomicAlignments_1.34.0               
    [123] Rsamtools_2.14.0                       
    [124] GenomeInfoDbData_1.2.9                 
    [125] parallel_4.2.2                         
    [126] hms_1.1.2                              
    [127] grid_4.2.2                             
    [128] ggfun_0.0.9                            
    [129] tidyr_1.3.0                            
    [130] HDO.db_0.99.1                          
    [131] rmarkdown_2.20                         
    [132] MatrixGenerics_1.10.0                  
    [133] ggforce_0.4.1                          
    [134] Biobase_2.58.0                         
    [135] restfulr_0.0.15                        
