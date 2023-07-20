
            set terminal png size 600,400 truecolor
            set output "IP1_srt_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP1_srt_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",40.70) at 40.70,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000039
2	0.000124
4	0.000280
6	0.000597
8	0.001335
10	0.002766
12	0.005998
14	0.012337
16	0.024310
18	0.044784
20	0.077282
22	0.123982
24	0.188432
26	0.270220
28	0.371091
30	0.493231
32	0.630875
34	0.788436
36	0.930817
38	0.988732
40	1.000000
42	0.966167
44	0.940068
46	0.900751
48	0.835942
50	0.756886
52	0.651716
54	0.534405
56	0.411288
58	0.301950
60	0.212123
62	0.142989
64	0.091482
66	0.058196
68	0.036921
70	0.024252
72	0.017166
74	0.013229
76	0.009759
78	0.007175
80	0.005330
82	0.003909
84	0.002627
86	0.001662
88	0.001066
90	0.000623
92	0.000341
94	0.000145
96	0.000067
98	0.000027
end
