
            set terminal png size 600,400 truecolor
            set output "IP1_input_srt_stats/gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "IP1_input_srt_stats.txt" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",40.70) at 40.70,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        0	0.000065
2	0.000196
4	0.000450
6	0.000945
8	0.002000
10	0.004114
12	0.008430
14	0.016493
16	0.030991
18	0.055330
20	0.091999
22	0.144001
24	0.213116
26	0.299716
28	0.402994
30	0.523743
32	0.654330
34	0.794733
36	0.917058
38	0.980086
40	1.000000
42	0.985140
44	0.961274
46	0.917363
48	0.846476
50	0.758763
52	0.646604
54	0.522369
56	0.395102
58	0.283730
60	0.192855
62	0.124098
64	0.075430
66	0.044479
68	0.025379
70	0.014832
72	0.009085
74	0.006136
76	0.004203
78	0.003003
80	0.002197
82	0.001618
84	0.001109
86	0.000719
88	0.000467
90	0.000306
92	0.000160
94	0.000082
96	0.000036
98	0.000012
end
