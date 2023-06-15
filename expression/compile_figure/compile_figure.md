# compile_figure

- [<span class="toc-section-number">1</span> Intro](#intro)
  - [<span class="toc-section-number">1.1</span> Setup](#setup)
  - [<span class="toc-section-number">1.2</span> Read in
    files](#read-in-files)
- [<span class="toc-section-number">2</span> Individual
  plots](#individual-plots)
- [<span class="toc-section-number">3</span> Combined
  plot](#combined-plot)

# Intro

In this notebook, we compile the supplemental figure for the paper.

## Setup

``` r
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
```

## Read in files

``` r
df_pdac <- read_csv("../gse133979_pdac/output/01_QC/pdac_deaf1_norm.csv")
```

    New names:
    Rows: 34 Columns: 6
    ── Column specification
    ──────────────────────────────────────────────────────── Delimiter: "," chr
    (4): ...1, sample, condition, sex dbl (2): norm_counts, sample_no
    ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    Specify the column types or set `show_col_types = FALSE` to quiet this message.
    • `` -> `...1`

``` r
df_tscmko <- read_csv("../sarcoatlas_tscmko/output/01_QC/sarcoatlas_tscmko_deaf1_normcounts.csv")
```

    New names:
    Rows: 20 Columns: 3
    ── Column specification
    ──────────────────────────────────────────────────────── Delimiter: "," chr
    (2): ...1, sample_id dbl (1): counts
    ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    Specify the column types or set `show_col_types = FALSE` to quiet this message.
    • `` -> `...1`

``` r
df_agingcrrm <- read_csv("../sarcoatlas_agingcrrm/output/01_QC/sarcoatlas_agingcrrm_deaf1_normcounts.csv")
```

    New names:
    Rows: 46 Columns: 6
    ── Column specification
    ──────────────────────────────────────────────────────── Delimiter: "," chr
    (5): ...1, sample_id, age, muscle, rep dbl (1): counts
    ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    Specify the column types or set `show_col_types = FALSE` to quiet this message.
    • `` -> `...1`

``` r
df_time <- read_csv("../sarcoatlas_timecourse_mouse/output/01_QC/sarcoatlas_timeseries_deaf1_normcounts.csv")
```

    New names:
    Rows: 50 Columns: 5
    ── Column specification
    ──────────────────────────────────────────────────────── Delimiter: "," chr
    (4): ...1, age, muscle, rep dbl (1): norm_counts
    ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    Specify the column types or set `show_col_types = FALSE` to quiet this message.
    • `` -> `...1`

``` r
# global plot theme

th <- theme_classic() +
  theme(axis.title.y = element_text(face = "italic"),
        legend.position = "none")
```

# Individual plots

<div class="panel-tabset">

## PDAC

``` r
head(df_pdac)
```

| …1           | norm_counts | sample       | condition | sex | sample_no |
|:-------------|------------:|:-------------|:----------|:----|----------:|
| cachexia_f_1 |    55.68739 | cachexia_f_1 | cachexia  | f   |         1 |
| cachexia_f_2 |    68.72244 | cachexia_f_2 | cachexia  | f   |         2 |
| cachexia_f_3 |    58.80755 | cachexia_f_3 | cachexia  | f   |         3 |
| cachexia_f_4 |    58.57067 | cachexia_f_4 | cachexia  | f   |         4 |
| cachexia_f_5 |    78.69038 | cachexia_f_5 | cachexia  | f   |         5 |
| cachexia_f_6 |    77.57164 | cachexia_f_6 | cachexia  | f   |         6 |

Change ctrl to Control, cachexia to Cachexia:

``` r
df_pdac$condition <- ifelse(df_pdac$condition == "ctrl", 
                            "Ctrl", 
                            "Cachexia")

df_pdac$condition <- factor(df_pdac$condition, 
                            levels = c("Ctrl", "Cachexia"))
```

``` r
 #| label: plot-a-pdac

gse_lab <- "GSE133979"

p1 <- ggplot(df_pdac, aes(x = condition, y = norm_counts,
                          color = condition)) +
  geom_boxplot() + 
  geom_point() +
  scale_y_continuous(limits = c(40, 160)) +
  labs(
    title = gse_lab,
    y = "DEAF1 Expression") +
  scale_color_brewer(palette = "Set1") +
  th +
  theme(axis.title.x = element_blank())

p1
```

![](compile_figure_files/figure-commonmark/unnamed-chunk-6-1.png)

## TSCmKO

``` r
df_tscmko <- df_tscmko %>% 
  separate(sample_id, into = c("age", "condition", "rep"), 
           remove = FALSE)

df_tscmko$age <- ifelse(df_tscmko$age == "X3M", "3", "6")

df_tscmko$condition <- ifelse(df_tscmko$condition == "CTRL", 
                              "Control", "TSCmKO")

head(df_tscmko)
```

| …1              |   counts | sample_id       | age | condition | rep  |
|:----------------|---------:|:----------------|:----|:----------|:-----|
| X3M_CTRL_rep1   | 232.6514 | X3M_CTRL_rep1   | 3   | Control   | rep1 |
| X3M_CTRL_rep2   | 247.3036 | X3M_CTRL_rep2   | 3   | Control   | rep2 |
| X3M_CTRL_rep3   | 245.2994 | X3M_CTRL_rep3   | 3   | Control   | rep3 |
| X3M_CTRL_rep4   | 286.0115 | X3M_CTRL_rep4   | 3   | Control   | rep4 |
| X3M_CTRL_rep5   | 264.3815 | X3M_CTRL_rep5   | 3   | Control   | rep5 |
| X3M_TSCmKO_rep1 | 197.4112 | X3M_TSCmKO_rep1 | 3   | TSCmKO    | rep1 |

``` r
gse_lab <- "GSE139213"

p2 <- ggplot(df_tscmko, aes(x = age, y = counts, color = age)) +
  geom_boxplot() +
  geom_point() +
  scale_y_continuous(limits = c(100, 350)) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = gse_lab, 
    x = "Age (months)", 
    y = "DEAF1 Expression") +
  facet_wrap(~ condition) +
  theme_classic() +
  theme(legend.position = "none") + 
  th

p2
```

![](compile_figure_files/figure-commonmark/plot-b-tscmko-1.png)

## Aging-crrm

``` r
df_agingcrrm$age <- ifelse(df_agingcrrm$age == "X10m", "Young",
                           "Old")

df_agingcrrm$age <- factor(df_agingcrrm$age,
                           levels = c("Young", "Old"))

df_agingcrrm$muscle <- gsub("Sol", "SOL", df_agingcrrm$muscle)

df_agingcrrm$muscle <- factor(df_agingcrrm$muscle, 
                              levels = c("GAS", "TA", "SOL", "TRI"))

head(df_agingcrrm)
```

| …1            |   counts | sample_id     | age   | muscle | rep  |
|:--------------|---------:|:--------------|:------|:-------|:-----|
| X10m_Sol_rep1 | 117.8690 | X10m_Sol_rep1 | Young | SOL    | rep1 |
| X10m_Sol_rep2 | 105.9695 | X10m_Sol_rep2 | Young | SOL    | rep2 |
| X10m_Sol_rep3 | 111.7793 | X10m_Sol_rep3 | Young | SOL    | rep3 |
| X10m_Sol_rep4 | 119.0885 | X10m_Sol_rep4 | Young | SOL    | rep4 |
| X10m_Sol_rep5 | 107.6350 | X10m_Sol_rep5 | Young | SOL    | rep5 |
| X10m_Sol_rep6 | 106.0339 | X10m_Sol_rep6 | Young | SOL    | rep6 |

``` r
gse_lab <- "GSE139204"

p3 <- ggplot(df_agingcrrm, aes(x = age, y = counts,
                         color = age)) +
  geom_boxplot() +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(100, 250)) +
  labs(
    title = gse_lab,
    y = "DEAF1 Expression") +
  facet_wrap(~muscle, nrow = 1) +
  th +
  theme(axis.title.x = element_blank())

p3
```

![](compile_figure_files/figure-commonmark/plot-c-agingcrrm-1.png)

## Timecourse - mouse

``` r
unique(df_time$age)
```

    [1] "X8m"  "X18m" "X22m" "X24m" "X26m" "X28m"

``` r
df_time$age <- as.factor(gsub("X|m", "", df_time$age))

# arrange in correct order
df_time$age <- factor(df_time$age, 
                      levels = c("8", "18", "22", "24", 
                                 "26", "28"))

df_time
```

| …1            | norm_counts | age | muscle | rep  |
|:--------------|------------:|:----|:-------|:-----|
| X8m_GAS_rep1  |   161.15864 | 8   | GAS    | rep1 |
| X8m_GAS_rep2  |   172.90319 | 8   | GAS    | rep2 |
| X8m_GAS_rep3  |   159.38220 | 8   | GAS    | rep3 |
| X8m_GAS_rep4  |   170.57883 | 8   | GAS    | rep4 |
| X8m_GAS_rep5  |   166.86661 | 8   | GAS    | rep5 |
| X8m_GAS_rep6  |   157.23516 | 8   | GAS    | rep6 |
| X8m_GAS_rep7  |   167.42633 | 8   | GAS    | rep7 |
| X8m_GAS_rep8  |   157.40237 | 8   | GAS    | rep8 |
| X18m_GAS_rep1 |   164.86595 | 18  | GAS    | rep1 |
| X18m_GAS_rep2 |   205.43891 | 18  | GAS    | rep2 |
| X18m_GAS_rep3 |   126.70343 | 18  | GAS    | rep3 |
| X18m_GAS_rep4 |   156.24575 | 18  | GAS    | rep4 |
| X18m_GAS_rep5 |   154.00373 | 18  | GAS    | rep5 |
| X18m_GAS_rep6 |   118.77189 | 18  | GAS    | rep6 |
| X18m_GAS_rep7 |   136.22613 | 18  | GAS    | rep7 |
| X18m_GAS_rep8 |   143.21406 | 18  | GAS    | rep8 |
| X22m_GAS_rep1 |   169.14703 | 22  | GAS    | rep1 |
| X22m_GAS_rep2 |   149.85283 | 22  | GAS    | rep2 |
| X22m_GAS_rep3 |   175.60666 | 22  | GAS    | rep3 |
| X22m_GAS_rep4 |   170.10363 | 22  | GAS    | rep4 |
| X22m_GAS_rep5 |   151.26732 | 22  | GAS    | rep5 |
| X22m_GAS_rep6 |   176.35855 | 22  | GAS    | rep6 |
| X22m_GAS_rep7 |   169.59500 | 22  | GAS    | rep7 |
| X22m_GAS_rep8 |   141.80138 | 22  | GAS    | rep8 |
| X24m_GAS_rep1 |   132.20351 | 24  | GAS    | rep1 |
| X24m_GAS_rep2 |   151.10455 | 24  | GAS    | rep2 |
| X24m_GAS_rep3 |   149.25532 | 24  | GAS    | rep3 |
| X24m_GAS_rep4 |   119.58540 | 24  | GAS    | rep4 |
| X24m_GAS_rep5 |   171.95977 | 24  | GAS    | rep5 |
| X24m_GAS_rep6 |   171.24356 | 24  | GAS    | rep6 |
| X24m_GAS_rep7 |   131.63596 | 24  | GAS    | rep7 |
| X24m_GAS_rep8 |   166.97187 | 24  | GAS    | rep8 |
| X26m_GAS_rep1 |   108.51790 | 26  | GAS    | rep1 |
| X26m_GAS_rep2 |   164.95612 | 26  | GAS    | rep2 |
| X26m_GAS_rep3 |   128.44601 | 26  | GAS    | rep3 |
| X26m_GAS_rep4 |   151.23339 | 26  | GAS    | rep4 |
| X26m_GAS_rep5 |   134.03827 | 26  | GAS    | rep5 |
| X26m_GAS_rep6 |   126.09438 | 26  | GAS    | rep6 |
| X26m_GAS_rep7 |   126.12345 | 26  | GAS    | rep7 |
| X26m_GAS_rep8 |   122.17300 | 26  | GAS    | rep8 |
| X26m_GAS_rep9 |   131.52481 | 26  | GAS    | rep9 |
| X28m_GAS_rep1 |   154.74911 | 28  | GAS    | rep1 |
| X28m_GAS_rep2 |    83.93168 | 28  | GAS    | rep2 |
| X28m_GAS_rep3 |   137.40020 | 28  | GAS    | rep3 |
| X28m_GAS_rep4 |   149.69437 | 28  | GAS    | rep4 |
| X28m_GAS_rep5 |   141.86211 | 28  | GAS    | rep5 |
| X28m_GAS_rep6 |   154.92932 | 28  | GAS    | rep6 |
| X28m_GAS_rep7 |   103.52808 | 28  | GAS    | rep7 |
| X28m_GAS_rep8 |   124.84789 | 28  | GAS    | rep8 |
| X28m_GAS_rep9 |   158.25580 | 28  | GAS    | rep9 |

``` r
gse_lab <- "GSE145480"

p4 <- ggplot(df_time, aes(x = age, y = norm_counts)) +
  geom_boxplot() + 
  geom_point() +
  scale_y_continuous(limits = c(0, 400)) +
  labs(
    title = gse_lab,
    x = "Age (months)", 
    y = "DEAF1 Expression"
  ) +
  th

p4
```

![](compile_figure_files/figure-commonmark/plot-d-time-1.png)

</div>

# Combined plot

``` r
layout <- "
AACCCC
BBDDDD
"

combined_plot <- p1 + p2 + p3 + p4 +
  plot_layout(design = layout)

combined_plot
```

![](compile_figure_files/figure-commonmark/combined-plot-1.png)

``` r
ggsave(filename = "combined_plot.pdf", 
       plot = combined_plot, device = "pdf", 
       width = 20, height = 25,
       units = "cm")
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
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] patchwork_1.1.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.1    
     [5] purrr_1.0.1     readr_2.1.3     tidyr_1.3.0     tibble_3.2.1   
     [9] ggplot2_3.4.2   tidyverse_1.3.2

    loaded via a namespace (and not attached):
     [1] lubridate_1.9.1     assertthat_0.2.1    digest_0.6.31      
     [4] utf8_1.2.3          R6_2.5.1            cellranger_1.1.0   
     [7] backports_1.4.1     reprex_2.0.2        evaluate_0.20      
    [10] httr_1.4.5          pillar_1.9.0        rlang_1.1.0        
    [13] googlesheets4_1.0.1 readxl_1.4.1        rstudioapi_0.14    
    [16] rmarkdown_2.20      textshaping_0.3.6   labeling_0.4.2     
    [19] googledrive_2.0.0   bit_4.0.5           munsell_0.5.0      
    [22] broom_1.0.3         compiler_4.2.2      modelr_0.1.10      
    [25] xfun_0.37           systemfonts_1.0.4   pkgconfig_2.0.3    
    [28] htmltools_0.5.4     tidyselect_1.2.0    fansi_1.0.4        
    [31] crayon_1.5.2        tzdb_0.3.0          dbplyr_2.3.0       
    [34] withr_2.5.0         grid_4.2.2          jsonlite_1.8.4     
    [37] gtable_0.3.3        lifecycle_1.0.3     DBI_1.1.3          
    [40] magrittr_2.0.3      scales_1.2.1        cli_3.6.1          
    [43] stringi_1.7.12      vroom_1.6.1         farver_2.1.1       
    [46] fs_1.6.2            xml2_1.3.3          ragg_1.2.5         
    [49] ellipsis_0.3.2      generics_0.1.3      vctrs_0.6.1        
    [52] RColorBrewer_1.1-3  tools_4.2.2         bit64_4.0.5        
    [55] glue_1.6.2          hms_1.1.2           parallel_4.2.2     
    [58] fastmap_1.1.1       yaml_2.3.7          timechange_0.2.0   
    [61] colorspace_2.1-0    gargle_1.3.0        rvest_1.0.3        
    [64] knitr_1.42          haven_2.5.1        
