
            set terminal png size 600,400 truecolor
            set output "IP2_input_srt_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP2_input_srt_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",40.70) at 40.70,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000067
2	0.000213
4	0.000497
6	0.001100
8	0.002278
10	0.004558
12	0.009194
14	0.018062
16	0.033907
18	0.059891
20	0.099018
22	0.153970
24	0.226000
26	0.314879
28	0.420252
30	0.541690
32	0.673199
34	0.815098
36	0.935493
38	0.989938
40	1.000000
42	0.979590
44	0.956033
46	0.916099
48	0.847862
50	0.763850
52	0.651961
54	0.529776
56	0.402617
58	0.290979
60	0.199929
62	0.129432
64	0.079175
66	0.047122
68	0.027188
70	0.016051
72	0.010045
74	0.006808
76	0.004619
78	0.003375
80	0.002527
82	0.001782
84	0.001241
86	0.000792
88	0.000512
90	0.000307
92	0.000175
94	0.000074
96	0.000035
98	0.000011
end
