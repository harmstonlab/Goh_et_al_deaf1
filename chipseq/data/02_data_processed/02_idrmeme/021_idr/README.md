# IDR

Calculates the IDR score (how reproducible a peak is) between two replicates. 

Note: The original code only adds an IDR score to the .bed file - it contains ALL peaks, and is not automatically filtered according to the IDR threshold. 

## Files of interest: 

- `deaf1_idr_filt.bed`: Contains filtered IDR peaks passing threshold (IDR 0 - 0.05). 608 peaks total. 
- `deaf1_idr_all.bed`: Contains all peaks with IDR scores, non filtered. 1550 peaks total. 

- `idr_filt.md` and `.Rmd`: Contains code used to filter the IDR peaks. 


## Code: 

```sh
# Sort peak by IP1, IP2
sort -k8,8nr IP1_peaks.narrowPeak > IP1_sorted_peaks.narrowPeak
sort -k8,8nr IP2_peaks.narrowPeak > IP2_sorted_peaks.narrowPeak
```

```sh
# Run IDR
idr --samples IP1_sorted_peaks.narrowPeak IP2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file deaf1_idr_all.bed \
--output-file-type bed \
--plot \
--log-output-file deaf1_idr.log

```

