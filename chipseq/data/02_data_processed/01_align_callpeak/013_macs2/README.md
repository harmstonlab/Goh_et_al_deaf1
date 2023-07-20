
# Peak calling


## MACS2 peak calling

We call peaks with `macs2`. A brief explaination of flags used: 

- `-t`: treatment (the IP file)
- `-c`: control (non-IP file)
- `-f`: format of input (BAM file)
- `-bdg` output bedgraph file for control and treatment
- `--name`: prefix for files
- `q`: set q value of 0.1
- `2>`: redirects process messages to a log file, which we can watch using `tail -f IP1-macs2.log` 


We first call peaks for our filtered bam files (i.e. remove duplicates and unmapped reads)

```sh
#IP1
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_srt_filtered.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_input_srt_filtered.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP1 \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP1-macs2.log
```

```sh
#IP2
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_srt_filtered.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_input_srt_filtered.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP2 \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP2-macs2.log
```

For comparison, we also call peaks for the unfiltered bam files

```sh
#IP1, unfiltered
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_srt.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_input_srt.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP1_unfiltered \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP1_unfiltered_macs2.log
```

```sh
#IP2, unfiltered
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_srt.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_input_srt.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP2_unfiltered \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP2_unfiltered_macs2.log
```

Merging them all into a script 
```sh

echo "Processing IP2" \
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_srt_filtered.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_input_srt_filtered.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP2 \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP2-macs2.log \

echo "Processing IP1_input" \
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_srt.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP1_input_srt.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP1_unfiltered \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP1_unfiltered_macs2.log \

echo "Processing IP2 input" \
macs2 callpeak -t  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_srt.bam \
	-c  ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/IP2_input_srt.bam \
 	-f BAM -g 1.87e+9 \
	-q 0.1 \
    --bdg \
	--name IP2_unfiltered \
	--outdir ~/deaf1_chipseq/data/02_data_processed/04_macs2 2> ~/deaf1_chipseq/data/02_data_processed/04_macs2/IP2_unfiltered_macs2.log

```


