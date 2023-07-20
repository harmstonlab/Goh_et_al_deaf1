# 022_MEME

Run the MEME motif finding enrichment on the server to identify Deaf1-binding motif. This serves as a QC for ChIP-seq data.  

## Method
Peaks were extracted from `deaf1_idr_filt.bed` using the `meme_prep.Rmd` file, and saved as the `deaf1_peak_seqs.fa` fasta file. We vary the `minsites` threshold from 60 - 450 (10% to 75% of all sequences)

## Result
- 1 motif was found in all 608/608 IDR peaks; this is consistent throughout minsites 150 and above. 



## Code: 

Flags explaination: 
- nmotifs: Pick top 10 motifs
- maxw: Max motif width of 20 base pairs (expected of transcription factors)
- revcomp: allow matches on + and - strands
- minsites: Minimum number of sites that contain this motif

Least stringent: present in 10% of all sequences
```sh
meme deaf1_peak_seqs.fa -dna -nmotifs 10 -maxw 20 -revcomp -minsites 60 -o meme_60
```

Moderately stringent: present in 25% of all sequences

```sh
meme deaf1_peak_seqs.fa -dna -nmotifs 10 -maxw 20 -revcomp -minsites 150 -o meme_150
```


Moderately stringent: present in 50% of all sequences

```sh
meme deaf1_peak_seqs.fa -dna -nmotifs 10 -maxw 20 -revcomp -minsites 300 -o meme_300
```


More stringent: present in 75% of all sequences

```sh
meme deaf1_peak_seqs.fa -dna -nmotifs 10 -maxw 20 -revcomp -minsites 450 -o meme_450
```



