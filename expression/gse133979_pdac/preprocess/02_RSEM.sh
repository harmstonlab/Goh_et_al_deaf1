for i in *.

do
  rsem-calculate-expression --alignments \
                            --no-bam-output \
                            --strandedness none \
                            --num-threads 20 \
                            /home/qianhui/deaf1/data/02_aligned/01_STAR/${i}.Aligned.toTranscriptome.out.bam \
                            /home/shared/genomes/hg38/hg38_109/hg38_109 \
                            $i 
done



while read i

do
  rsem-calculate-expression --alignments \
                            --no-bam-output \
                            --strandedness none \
                            --temporary-folder temp \
                            --num-threads 20 \
                            /home/qianhui/deaf1/data/02_aligned/01_STAR/${i}.Aligned.toTranscriptome.out.bam \
                            /home/shared/genomes/hg38/hg38_109/hg38_109 \
                            $i 
done < /home/qianhui/deaf1/data/01_data_raw/srr_acc_list.txt