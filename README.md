# deaf1_paper
All scripts and code necessary to recreate analyses and figures for Hongwen's Deaf1 paper

# ChIP-seq

Code for ChIP-seq figure. Contains Deaf1-chipseq data. 

# Expression

Code for Supplementary Figure 1. Contains RNA-seq data from GEO. 

## Data

Counts for the three sarcoatlas (mouse) datasets were downloaded from GEO. For `gse133979_pdac` (human), raw data was downloaded from SRA, realigned and then analyzed. 

Title                       | GEO accession | Annotation file used 
--------------------------- | ------------- | ----------------------------------
sarcoatlas_agingcrrm        | [GSE139209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139209) | Mus_musculus.GRCm38.94.chr.gtf.gz
sarcoatlas_timecourse_mouse | [GSE145480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145480) | Mus_musculus.GRCm38.94.chr.gtf.gz
sarcoatlas_tscmko           | [GSE139213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139213) | Mus_musculus.GRCm38.94.chr.gtf.gz
gse133979_pdac              | [GSE133979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133979) | Homo_sapiens.GRCh38.109.chr.gtf.gz

## Analysis

**! Ensure that the annotation files are in the [expression/annotation](expression/annotation/) folder before running the analysis.**

If running from scratch, these scripts can be found as `.qmd` files in the [expression](expression) folder. If the repository is cloned as is, and the .Rproj file is used, there should have no issues running the scripts.

## Checksums

Annotation file                    | MD5 checksum
---------------------------------- | --------------------------------
Homo_sapiens.GRCh38.109.chr.gtf.gz | 4fbfbb5c5fadf35f50f8f7134d7a2412
Mus_musculus.GRCm38.94.chr.gtf.gz  | 9aa004b6c98fc1aec6973af98e22b822

## Version info

