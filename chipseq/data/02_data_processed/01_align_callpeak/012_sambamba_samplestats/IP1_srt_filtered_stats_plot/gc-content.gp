
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_srt_filtered_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP1_srt_filtered_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",42.71) at 42.71,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000021
2	0.000081
4	0.000244
6	0.000621
8	0.001476
10	0.003178
12	0.007061
14	0.014695
16	0.029186
18	0.053928
20	0.093139
22	0.149035
24	0.224851
26	0.316216
28	0.422181
30	0.535640
32	0.650075
34	0.757489
36	0.852055
38	0.925547
40	0.975417
42	1.000000
44	0.998280
46	0.968118
48	0.910728
50	0.824898
52	0.712487
54	0.584893
56	0.454536
58	0.333936
60	0.231560
62	0.152921
64	0.095983
66	0.058927
68	0.036051
70	0.022516
72	0.015027
74	0.010801
76	0.008031
78	0.005916
80	0.004386
82	0.003136
84	0.002140
86	0.001306
88	0.000789
90	0.000379
92	0.000173
94	0.000065
96	0.000025
98	0.000007
end
