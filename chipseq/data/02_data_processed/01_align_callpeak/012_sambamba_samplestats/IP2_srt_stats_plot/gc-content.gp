
            set terminal png size 600,400 truecolor
            set output "IP2_srt_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP2_srt_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",40.70) at 40.70,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000036
2	0.000121
4	0.000292
6	0.000591
8	0.001307
10	0.002788
12	0.006068
14	0.012674
16	0.025242
18	0.046388
20	0.080229
22	0.128467
24	0.194056
26	0.277584
28	0.379027
30	0.501286
32	0.635914
34	0.789238
36	0.926720
38	0.987144
40	1.000000
42	0.967708
44	0.941383
46	0.902113
48	0.837611
50	0.759075
52	0.654917
54	0.537961
56	0.415321
58	0.305656
60	0.216660
62	0.148047
64	0.096165
66	0.062562
68	0.041051
70	0.027862
72	0.020399
74	0.016257
76	0.012117
78	0.008900
80	0.006642
82	0.004967
84	0.003222
86	0.002114
88	0.001251
90	0.000764
92	0.000407
94	0.000193
96	0.000090
98	0.000036
end
