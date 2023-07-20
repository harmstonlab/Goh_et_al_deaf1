
        set terminal png size 600,400 truecolor
        set output "IP1_srt_stats/indel-dist.png"
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
        set title "IP1_srt_stats.txt" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	142845
2	28884
3	4325
end
1	184444
2	33792
3	4465
end
1	0.774463
2	0.854759
3	0.968645
end
