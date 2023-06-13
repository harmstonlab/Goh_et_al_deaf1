# Deaf-1 ChIP-seq analysis

Analysis of Deaf-1 ChIP-seq data for C2C12 mouse myoblasts.

## Quick navigation 
- [notebooks](notebooks) folder: Code for downstream analysis in R. (QC, GO enrichment etc)
- [data/02_data_processed/01_align_callpeak](data/02_data_processed/01_align_callpeak): Code for alignment and filtering from the server (`bwa`, `samtools`, `sambamba`).
- [data/02_data_processed/02_idrmeme](data/02_data_processed/02_idrmeme): Code and results for `IDR` and `MEME`. 

## Important files
- **`.bed` file with highly reproducible IDR peaks**: [data/02_data_processed/02_idrmeme/021_idr/deaf1_idr_filt.bed](data/02_data_processed/02_idrmeme/021_idr/deaf1_idr_filt.bed)
- **`MEME` motif enrichment results**: [data/02_data_processed/02_idrmeme/022_meme/](data/02_data_processed/02_idrmeme/022_meme/)

## Sample info

Original name | New name | Description |
|---|---| --- |
Clean_IP-C2C12-1_1.fq.gz | IP1.fq.gz | Sample 1 
Clean_IP-C2C12-2_1.fq.gz | IP2.fq.gz | Sample 2 
Clean_Input-C2C12-1_1.fq.gz | IP1_input.fq.gz | Sample 1 input (non-IP) 
Clean_Input-C2C12-2_1.fq.gz | IP2_input.fq.gz | Sample 2 input (non-IP) 

## Brief description of pipeline
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


## Peak stats:

Using `macs2` q = 0.1: 

- IP1: 1860 peaks 
- IP2: 3126 peaks
- 1550 total peaks upon merging with IDR 
- 608 peaks that are 0 <= IDR <= 0.05
- 595 unique genes bound by 608 IDR peaks
- 12 genes have >1 binding site; max number of binding sites in a gene is 3. 

## SHA checksums
- `Mus_musculus.GRCm38.102.chr.gtf`: e211ecb4ee8b735630a57c32a18715d0
- `mm10.fa`: 7e329b2bf419a9f5a7dc42c625c884ac

## Version info 

- FastQC: **v0.11.8**
- BWA: **0.7.17-r1198-dirty**
- samblaster: **0.1.24**
- sambamba **0.7.0**
- samtools **1.9**
- macs2 **2.1.2**
- ChIPQC **1.32.2**






