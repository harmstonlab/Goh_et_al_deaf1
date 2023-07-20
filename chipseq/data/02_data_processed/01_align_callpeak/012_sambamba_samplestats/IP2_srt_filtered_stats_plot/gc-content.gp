
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_srt_filtered_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP2_srt_filtered_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",42.71) at 42.71,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000018
2	0.000084
4	0.000268
6	0.000618
8	0.001457
10	0.003243
12	0.007151
14	0.015102
16	0.030162
18	0.055518
20	0.096060
22	0.153695
24	0.230282
26	0.323081
28	0.429551
30	0.543298
32	0.655885
34	0.763024
36	0.856660
38	0.928709
40	0.977548
42	1.000000
44	0.995390
46	0.963632
48	0.906086
50	0.820405
52	0.709686
54	0.583113
56	0.454262
58	0.334331
60	0.233326
62	0.155092
64	0.098218
66	0.061104
68	0.038040
70	0.024336
72	0.016435
74	0.011969
76	0.009054
78	0.006772
80	0.005011
82	0.003700
84	0.002447
86	0.001543
88	0.000882
90	0.000481
92	0.000210
94	0.000080
96	0.000027
98	0.000007
end
