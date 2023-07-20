022_meme_prep
================
Qian Hui TAN
2/14/23

- <a href="#meme_prep" id="toc-meme_prep"><span
  class="toc-section-number">1</span> 07_MEME_Prep</a>
- <a href="#setup" id="toc-setup"><span
  class="toc-section-number">2</span> Setup</a>
  - <a href="#checking" id="toc-checking"><span
    class="toc-section-number">2.1</span> Checking</a>
- <a href="#resizing" id="toc-resizing"><span
  class="toc-section-number">3</span> Resizing</a>

# 07_MEME_Prep

In this notebook, we extract sequences from IDR peaks, and save them
into a fasta ile. This fasta file will be used to run MEME motif finding
enrichments on the server.

``` r
suppressPackageStartupMessages({
 library(GenomicRanges)
 library(BSgenome.Mmusculus.UCSC.mm10)
 library(Biostrings)
 library(GenomicFeatures)
  
})
```

# Setup

``` r
deaf1_bed <- read.table("../021_idr/deaf1_idr_filt.bed")
genome <- BSgenome.Mmusculus.UCSC.mm10
```

## Checking

``` r
# check that min IDR is 540 
summary(deaf1_bed$scaled_idr)
```

    Length  Class   Mode 
         0   NULL   NULL 

``` r
head(deaf1_bed)
```

         V1        V2        V3   V4 V5 V6  V7 V8   V9  V10       V11       V12
    1  chr1 166378462 166382216 3755  *  . 845  . 1.99 2.04 166378501 166382203
    2  chr7 126427901 126429745 1845  *  . 826  . 1.90 1.99 126428279 126429503
    3  chr5 143527422 143528800 1379  *  . 836  . 1.94 2.02 143527539 143528761
    4 chr15  95790481  95792623 2143  *  . 808  . 1.82 1.95  95790480  95792531
    5  chr8   8688634   8691843 3210  *  . 854  . 2.03 2.06   8689367   8691612
    6  chr7   4792394   4793284  891  *  . 826  . 1.89 1.99   4792415   4793271
            V13       V14       V15      V16
    1 186.80782 166378461 166382216 245.7580
    2 131.42822 126427900 126429745 164.8703
    3 119.34243 143527421 143528800 149.1025
    4 109.16676  95790507  95792623 109.6928
    5 109.16571   8688633   8691843 144.2870
    6  96.24286   4792393   4793284 102.2998

# Resizing

``` r
# Get a range within 500 bps of each chip seq peak
deaf1_peaks = GRanges(deaf1_bed$V1,
               IRanges(deaf1_bed$V2, deaf1_bed$V3))

deaf1_peaks = keepStandardChromosomes(deaf1_peaks, 
                                      species = "Mus musculus", 
                                      pruning.mode = "coarse")

deaf1_peaks = resize(deaf1_peaks, width = 500, fix = "center")

head(deaf1_peaks)
```

    GRanges object with 6 ranges and 0 metadata columns:
          seqnames              ranges strand
             <Rle>           <IRanges>  <Rle>
      [1]     chr1 166380089-166380588      *
      [2]     chr7 126428573-126429072      *
      [3]     chr5 143527861-143528360      *
      [4]    chr15   95791302-95791801      *
      [5]     chr8     8689989-8690488      *
      [6]     chr7     4792589-4793088      *
      -------
      seqinfo: 20 sequences from an unspecified genome; no seqlengths

``` r
# Retrieve the genome
deaf1_peak_seqs = getSeq(genome, deaf1_peaks)
names(deaf1_peak_seqs) = paste("seq", 1:length(deaf1_peak_seqs), sep = "_")
head(deaf1_peak_seqs)
```

    DNAStringSet object of length 6:
        width seq                                               names               
    [1]   500 CCACACCGCCGGGCTGTCAGCAG...GGTAGAGTAAAACCAGGCGACAT seq_1
    [2]   500 TGGAGAAAGAGATTTGCAGGAAG...GGATCCGGAAACAGCCGAGCGGA seq_2
    [3]   500 GGCTCGGGAAGCGCTGGGCGAGG...TTCCTGAGCCCCCTTGTTATGTG seq_3
    [4]   500 CGGGGCGACGCCGGGCGGAGGAC...GGAGCAGCGTTAAAAAGTTACAC seq_4
    [5]   500 CGCTCGAACTCCGCCTTCTTCTC...CCTAGCAGGCCCGGCCCGGAGTC seq_5
    [6]   500 CTAGACAAGGGCGCAAAACCGTT...TGGGCTTTCTCTCCTTGTACCCA seq_6

``` r
writeXStringSet(deaf1_peak_seqs, "deaf1_peak_seqs.fa")
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
     [1] GenomicFeatures_1.50.4             AnnotationDbi_1.60.0              
     [3] Biobase_2.58.0                     BSgenome.Mmusculus.UCSC.mm10_1.4.3
     [5] BSgenome_1.66.2                    rtracklayer_1.58.0                
     [7] Biostrings_2.66.0                  XVector_0.38.0                    
     [9] GenomicRanges_1.50.2               GenomeInfoDb_1.34.7               
    [11] IRanges_2.32.0                     S4Vectors_0.36.1                  
    [13] BiocGenerics_0.44.0               

    loaded via a namespace (and not attached):
     [1] MatrixGenerics_1.10.0       httr_1.4.4                 
     [3] bit64_4.0.5                 jsonlite_1.8.4             
     [5] assertthat_0.2.1            BiocFileCache_2.6.0        
     [7] blob_1.2.3                  GenomeInfoDbData_1.2.9     
     [9] Rsamtools_2.14.0            yaml_2.3.7                 
    [11] progress_1.2.2              pillar_1.8.1               
    [13] RSQLite_2.2.20              lattice_0.20-45            
    [15] glue_1.6.2                  digest_0.6.31              
    [17] htmltools_0.5.4             Matrix_1.5-3               
    [19] XML_3.99-0.13               pkgconfig_2.0.3            
    [21] biomaRt_2.54.0              zlibbioc_1.44.0            
    [23] BiocParallel_1.32.5         tibble_3.1.8               
    [25] KEGGREST_1.38.0             generics_0.1.3             
    [27] ellipsis_0.3.2              cachem_1.0.6               
    [29] SummarizedExperiment_1.28.0 cli_3.6.0                  
    [31] magrittr_2.0.3              crayon_1.5.2               
    [33] memoise_2.0.1               evaluate_0.20              
    [35] fansi_1.0.4                 xml2_1.3.3                 
    [37] tools_4.2.2                 prettyunits_1.1.1          
    [39] hms_1.1.2                   BiocIO_1.8.0               
    [41] lifecycle_1.0.3             matrixStats_0.63.0         
    [43] stringr_1.5.0               DelayedArray_0.24.0        
    [45] compiler_4.2.2              rlang_1.0.6                
    [47] grid_4.2.2                  RCurl_1.98-1.10            
    [49] rstudioapi_0.14             rjson_0.2.21               
    [51] rappdirs_0.3.3              bitops_1.0-7               
    [53] rmarkdown_2.20              restfulr_0.0.15            
    [55] codetools_0.2-19            DBI_1.1.3                  
    [57] curl_5.0.0                  R6_2.5.1                   
    [59] GenomicAlignments_1.34.0    knitr_1.42                 
    [61] dplyr_1.1.0                 fastmap_1.1.0              
    [63] bit_4.0.5                   utf8_1.2.3                 
    [65] filelock_1.0.2              stringi_1.7.12             
    [67] parallel_4.2.2              Rcpp_1.0.10                
    [69] vctrs_0.5.2                 png_0.1-8                  
    [71] dbplyr_2.3.0                tidyselect_1.2.0           
    [73] xfun_0.37                  
