
        set terminal png size 600,400 truecolor
        set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_srt_filtered_stats/indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "IP1_srt_filtered_stats.txt" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	122932
2	26864
3	3978
end
1	126116
2	29612
3	4192
end
1	0.974753
2	0.907200
3	0.948950
end
