#!/bin/sh
# Running STAR
# Input isn't in a gz format so I removed the readFiles zcat line

for i in *.fastq

do
  STAR  --genomeDir /home/shared/genomes/hg38/hg38_109/ \
        --sjdbGTFfile /home/shared/genomes/hg38/hg38_109/Homo_sapiens.GRCh38.109.chr.gtf \
        --readFilesIn /home/qianhui/deaf1/data/01_data_raw/${i}\
        --runThreadN 20 \
        --twopassMode Basic \
        --outWigType bedGraph \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --runDirPerm All_RWX \
        --outFileNamePrefix /home/qianhui/deaf1/data/02_aligned/01_STAR/${i}
done

# Indexing bam files

for i in *.sortedByCoord.out.bam

do
  echo "Indexing: "$i        
  samtools index $i $i".bai"
done
