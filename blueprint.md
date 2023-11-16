{{ template:title }}
{{ template:description}}
{{ template:badges }}
{{ bullets }}
{{ template:toc }}


## ChIP-seq

Code for ChIP-seq figure. Contains analysis of Deaf-1 ChIP-seq data for C2C12 mouse myoblasts.

### Important files
- **`.bed` file with highly reproducible IDR peaks**: [chipseq/data/02_data_processed/02_idrmeme/021_idr/deaf1_idr_filt.bed](chipseq/data/02_data_processed/02_idrmeme/021_idr/deaf1_idr_filt.bed)
- **`MEME` motif enrichment results**: [chipseq/data/02_data_processed/02_idrmeme/022_meme/](chipseq/data/02_data_processed/02_idrmeme/022_meme/)

### Code files
- [chipseq/data/02_data_processed/analysis.md](chipseq/data/02_data_processed/analysis.md): Code for upstream alignment and filtering (`bwa`, `samtools`, `sambamba`)
- [chipseq/data/02_data_processed/02_idrmeme/021_idr](chipseq/data/02_data_processed/02_idrmeme/021_idr): Code and results for `IDR` 
- [chipseq/data/02_data_processed/02_idrmeme/022_meme](chipseq/data/02_data_processed/02_idrmeme/022_meme): Code and results for `MEME` 
- [chipseq/analysis](chipseq/analysis/) folder: Code for downstream analysis in R. (QC, GO enrichment etc)

### Brief description of pipeline
1. Quick QC with `fastqc`
2. Align fastqc files to genome with `bwa`. 
3. Flag duplicates and unmapped reads with `samblaster`
4. Filter to remove duplicates and unmapped reads with `sambamba`
5. Check quality of filtered reads with `samtools stats` and `sambamba`
6. QC using ChIPQC in R
6. Call peaks with `macs2`
6. Download files for visualization in IGV 
7. Keep only high quality, reproducible peaks with `idr`. 
8. Download files for downstream analysis in R. 

### Peak stats

Using `macs2` q = 0.1: 

- IP1: 1860 peaks 
- IP2: 3126 peaks
- 1550 total peaks upon merging with IDR 
- 608 peaks that are 0 <= IDR <= 0.05
- 595 unique genes bound by 608 IDR peaks
- 12 genes have >1 binding site; max number of binding sites in a gene is 3. 

### SHA checksums
- `Mus_musculus.GRCm38.102.chr.gtf`: e211ecb4ee8b735630a57c32a18715d0
- `mm10.fa`: 7e329b2bf419a9f5a7dc42c625c884ac

### Version info 

- FastQC: **v0.11.8**
- BWA: **0.7.17-r1198-dirty**
- samblaster: **0.1.24**
- sambamba: **0.7.0**
- samtools: **1.9**
- macs2: **2.1.2**
- ChIPQC: **1.32.2**

## Expression (RNA-seq)

Code for Supplementary Figure 1. Contains RNA-seq data from GEO. 

### Data

Counts for the three sarcoatlas (mouse) datasets were downloaded from GEO. For `gse133979_pdac` (human), raw data was downloaded from SRA, realigned and then analyzed. 

Title                       | GEO accession | Annotation file used 
--------------------------- | ------------- | ----------------------------------
sarcoatlas_agingcrrm        | [GSE139209](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139209) | Mus_musculus.GRCm38.94.chr.gtf.gz
sarcoatlas_timecourse_mouse | [GSE145480](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145480) | Mus_musculus.GRCm38.94.chr.gtf.gz
sarcoatlas_tscmko           | [GSE139213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139213) | Mus_musculus.GRCm38.94.chr.gtf.gz
gse133979_pdac              | [GSE133979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133979) | Homo_sapiens.GRCh38.109.chr.gtf.gz

### Checksums

Annotation file                    | MD5 checksum
---------------------------------- | --------------------------------
Homo_sapiens.GRCh38.109.chr.gtf.gz | 4fbfbb5c5fadf35f50f8f7134d7a2412
Mus_musculus.GRCm38.94.chr.gtf.gz  | 9aa004b6c98fc1aec6973af98e22b822

### Version info
- FastQC: **v0.11.8**
- STAR: **2.7.1a**
- RSEM: **v1.3.1**
- MultiQC: **version 1.14 (931327d)**

## Analysis

**! Ensure that the annotation files are in the [expression/annotation](expression/annotation/) folder before running the analysis.**

If running from scratch, these scripts can be found as `.qmd` files in the [expression](expression) folder. If the repository is cloned as is, and the .Rproj file is used, there should have no issues running the scripts.

## Acknowledgements

This table of contents was built with https://github.com/andreasbm/readme




