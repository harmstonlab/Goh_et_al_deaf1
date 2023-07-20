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

## 1. Fastqc ##

```sh
fastqc -t 4 *.fq.gz
```
All samples have 50 bps. Read quality is good (>28, green) for all bases; no trimming needed.


## 2. Align reads with bwa ##

Our reads are 50bp long, so we use bwa aln. 

```sh
# Alignment 
for i in IP1 IP1_input \
         IP2 IP2_input \

do
    bwa aln -t 20 \
            /home/shared/genomes/mm10/BWAIndex/mm10.fa \
            ~/deaf1_chipseq/data/01_data_raw/${i}.fq.gz > \
            ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.sai

done

```

```sh
# Output aligned files to a sam file 
for i in IP1 IP1_input \
         IP2 IP2_input \

do
    bwa samse /home/shared/genomes/mm10/BWAIndex/mm10.fa \
             ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.sai \
             ~/deaf1_chipseq/data/01_data_raw/${i}.fq.gz > ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.sam 

done
```


## 3.1 Sort files by read names

`samblaster` requires files to be sorted by read_id (aka the QNAME). We run a quick `samtools` sort with the `-n` parameter: 

```sh
for i in IP1_input IP2 IP2_input \

do
# Sort by read names 
samtools sort -o ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}_byqname.sam \
              -n -O sam -@20 \
              ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.sam
done
```

## 3.2 Mark duplicates with samblaster

Next, we mark duplicates with `samblaster`. The `-M` parameter is for compatibility (there's a flag that has changed across the versions)

```sh
for i in IP1_input IP2 IP2_input \

do 
    samblaster -M -i ${i}_byqname.sam -o ${i}_marked.sam
done

```

## 3.3 Output as bam
Output as `bam`. The `-S` parameter specifies that we're inputting a `.sam` file, the `-b` parameter specifies a `.bam` file for output. 

```sh
for i in IP1_input IP2 IP2_input \

do
    samtools view --threads 20 -Sb ${i}_marked.sam > ${i}_marked.bam
done
```

## 4 Sort BAM files 
We're almost there!

Lastly, we sort the bam files according to chromosomal coordinates, which will allow us to index them for downstream analysis. We do this in two ways: 
## 4.1 Sort all reads, without filtering

No filtering for low quality, no removal for duplicates. 

```sh
for i in IP1 IP1_input \
         IP2 IP2_input \

do
        # Sort alignments by leftmost coordinates (default). Use 16 threads.
        samtools sort -@ 16 ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam \
         -o ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt.bam \
         # Index output files. Use 16 threads
        samtools index -@ 16 ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt.bam
done
```

Remove all the sam files for now - they are way too large and take up too much space - if we need them we can always convert the bam to sam. 
``
rm *.sam
rm *.sai
``


## 4.2 Filter and sort 

I tried using samtools, but somehow was unable to remove the duplicates - see the bottom segment for failed attempts and the respective code - so I switched to `sambamba` `instead. 

A big issue with `sambamba sort` is that it writes a chunk of code to our `/tmp` directory, which doesn't have much space to begin with, and is already quite full. I dealt with it by going into `/tmp` and removing the previous `sambamba` runs that were taking up space, and it worked. 

If it doesn't work, a potential workaround: try adding `--tmpdir=$WORKDIR/tmp ` to `sambamba sort` if the `tmp` directory is out of space.

Note: if it fails to run, check whether last executed line of code has a `\`. If yes, program waits for extra input (since it indicates code continues on to another line).

```sh
for i in IP1 IP1_input IP2 IP2_input \

do
    echo "1/5 Sorting genomic coordinates for ${i}"
    sambamba sort -t 16 --show-progress ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam 

    echo "2/5 Filtering ${i}" 
    sambamba view -h -t 30 -f bam --show-progress \
    -F "mapping_quality > 0 and not unmapped and not duplicate" \
    ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.sorted.bam > ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt_filtered.bam 

    echo "3/5 Tabulating flagstats for ${i}" 
    sambamba flagstat -t 30 --show-progress \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt_filtered.bam > ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_filtered_flagstats.txt
    # Unfiltered:
    sambamba flagstat -t 30 --show-progress \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt.bam > ~/deaf1_chipseq/data/02_data_processed/samplestats/$#{i}_srt_flagstats.txt 

    echo "4/5 Calculating samtools stats for ${i}" 
    # Filtered:
    samtools stats -@30 \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt_filtered.bam > \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_filtered_stats.txt 
    # Unfiltered stats:
    samtools stats -@30 \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_srt.bam > \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_stats.txt 

    echo "5/5 Plotting stats for ${i}" 
    # Filtered:
    plot-bamstats -p ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_filtered_stats/ \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_filtered_stats.txt 
    # Unfiltered stats
    plot-bamstats -p ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_stats/ \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_srt_stats.txt 
done 

```
Create indexes for filtered: 

```sh
for i in IP1 IP1_input \
         IP2 IP2_input \

do
        # Sort alignments by leftmost coordinates (default). Use 30 threads.
        samtools sort -@ 30 ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/${i}_srt_filtered.bam \
         -o ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/${i}_srt_filtered.bam \
         # Index output files. Use 30 threads
        samtools index -@ 30 ~/deaf1_chipseq/data/02_data_processed/02_bwa_srt/${i}_srt_filtered.bam 
done
```





# Failed attempts at removing duplicates using samtools

All code below failed for various reasons. Not sure why - which is why I used sambamba instead. 

**Attempt 0: The default on the pipeline.** 


Something is wrong with the samblaster step - this doesn't give any output at all. 
`samblaster` marks duplicates in read-id grouped paired end sam files. 

```sh

for i in IP1 IP1_input \
         IP2 IP2_input \

do
            # Sort by read names (QNAMES) instead of chromosomal coordinates
            # Use 20 threads for sorting
            samtools sort -@20 ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.sam -n  -O sam | \   
            # Mark duplicates with samblaster
            samblaster -M -i stdin -o stdout | \ 
            #output as bam
            samtools view --threads 20 -Sb stdin > ~/deaf1_chipseq/data/02_data_processed/03_bwa/${i}.bam 
done
```

**Attempt 1: Filter all 3 flags together** 
-doesn't work, somehow - not sure why.

```sh
for i in IP1 IP1_input \
         IP2 IP2_input \

do
    # Exclude mapping quality below 1 (i.e 0, which are the ones that failed to align.) 
    # Exclude Flag 0x4 (unmapped), 0x100 (multimapper), 0x400 (duplicates). 
    #Output as bam. 
    samtools view -@ 16 -q 1 -bF 0x4,0x100,0x400 \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam > \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam
    # Sort and index
    samtools sort -@ 16 ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam \
                -o ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam
    samtools index -@ 16 ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam
done
```


**Attempt 2: filter separately then merge**
-this doesn't work because it's a MERGE combine, not a merge intersect - meaning that 
things we've filtered out are now added back in which is ridiculous but ok. 

Note: `bwa-aln` has scores ranging from 0-37. The default pipeline threshold of 40 does not work. 
We exclude unmapped reads, multimappers and duplicates. 

```sh
for i in IP1 IP1_input \
         IP2 IP2_input \

do
    # Exclude mapping quality below 1 (i.e 0, which are the ones that failed to align.) 
    # Exclude Flag 0x4 (unmapped), 0x100 (multimapper), 0x400 (duplicates). 
    # Filter separately and then merge bc it doesn't accept multiple flags at once. 
    # Output as bam. 
    samtools view -@ 16 -q 1 -bF 0x4 \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam > \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_tmp1.bam

    samtools view -@ 16 -q 1 -bF 0x100 \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam > \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_tmp2.bam

    samtools view -@ 16 -q 1 -bF 0x400 \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam > \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_tmp3.bam

    samtools merge -u \
     -o ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam \
     ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_tmp[123].bam

    # Sort and index
    samtools sort -@ 16 ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam \
                -o ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam \
    samtools index -@ 16 ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam \

    # Check out flagstats
    samtools flagstat -@16 \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam > \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_filtered_srt_flagstats.txt
done
```
This does not work - somehow the duplicate removal fails. 



**Attempt 3: We do this in two sequential steps:**

1. Remove duplicates (0x400)
2. Remove unmapped reads (0x4)


```sh
for i in IP1 IP1_input \
         IP2 IP2_input \

do
    # Exclude mapping quality below 1 (i.e 0, which are the ones that failed to align.) 
    # Exclude 0x400 (duplicates) first, then exclude 0x4 (unmapped reads)
    # Output as bam. 

    # 1. Remove duplicates
    samtools view -@ 16 -q 1 -bF 0x400 \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_marked.bam > \
                ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam \
    # 2. Remove unmapped reads
    #samtools view -@ 16 -q 1 -bF 0x4 \
    #            ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam > \
    #            ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam \

    # 3. Sort and index
    samtools sort -@ 16 ~/deaf1_chipseq/data/02_data_processed/03_bwa_cp/${i}_filtered.bam |
    samtools index -@ 16 ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam \

    # 4. Check out flagstats
    samtools flagstat -@16 \
    ~/deaf1_chipseq/data/02_data_processed/bwa_srt/${i}_filtered_srt.bam > \
    ~/deaf1_chipseq/data/02_data_processed/samplestats/${i}_filtered_srt_flagstats.txt
done
```
This still fails. Somehow samtools fails to remove the 0x400 flag even though we explicitly specify it. 


