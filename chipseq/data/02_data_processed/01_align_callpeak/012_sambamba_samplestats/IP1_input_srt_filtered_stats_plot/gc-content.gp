
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_input_srt_filtered_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP1_input_srt_filtered_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",42.71) at 42.71,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000038
2	0.000141
4	0.000409
6	0.000975
8	0.002167
10	0.004576
12	0.009547
14	0.018871
16	0.035624
18	0.064010
20	0.106271
22	0.166120
24	0.244031
26	0.337525
28	0.443542
30	0.556250
32	0.667460
34	0.771991
36	0.862697
38	0.934466
40	0.980063
42	1.000000
44	0.993335
46	0.956837
48	0.893834
50	0.802011
52	0.685193
54	0.554929
56	0.424221
58	0.305266
60	0.206180
62	0.131541
64	0.079161
66	0.045897
68	0.025784
70	0.014639
72	0.008774
74	0.005656
76	0.003850
78	0.002701
80	0.001915
82	0.001331
84	0.000843
86	0.000508
88	0.000297
90	0.000141
92	0.000055
94	0.000024
96	0.000010
98	0.000001
end
