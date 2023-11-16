# 03_sambamba_samplestats

This folder contains the flag statistics of aligned reads, after alignment with `BWA` and filtering with `sambamba`.

## Quick navigation

File name | Description
|---| --------|
*.srt | Unfiltered reads, sorted by genomic coordinates.
*.filtered.srt | Filtered using `sambamba` for read quality >0, not duplicate, not unmapped.
*_flagstats.txt | Quick flag tabulation by `sambamba`. Contains short table summary of flags (duplicated, mapped, etc)
*_stats.txt | Detailed flag info from `samtools`. 
folders | Plots for flag stats, generated from `plot-bamtools`