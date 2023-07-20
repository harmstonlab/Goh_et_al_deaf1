
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_input_srt_filtered_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP2_input_srt_filtered_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",42.71) at 42.71,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000037
2	0.000151
4	0.000463
6	0.001125
8	0.002465
10	0.005055
12	0.010374
14	0.020656
16	0.039047
18	0.069223
20	0.114721
22	0.178168
24	0.259580
26	0.355320
28	0.462447
30	0.573288
32	0.682100
34	0.783636
36	0.869677
38	0.937342
40	0.979678
42	1.000000
44	0.993402
46	0.960621
48	0.899315
50	0.810018
52	0.694029
54	0.565629
56	0.433696
58	0.313987
60	0.214425
62	0.137380
64	0.083304
66	0.048502
68	0.027625
70	0.015869
72	0.009715
74	0.006347
76	0.004306
78	0.003116
80	0.002272
82	0.001528
84	0.001027
86	0.000612
88	0.000346
90	0.000176
92	0.000070
94	0.000024
96	0.000011
98	0.000002
end
