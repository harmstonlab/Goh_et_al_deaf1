
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_input_srt_filtered_stats/coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "IP1_input_srt_filtered_stats.txt" noenhanced
            set xrange [:6]
            plot '-' with lines notitle
        1	775074760
2	364740844
3	128695723
4	37397746
5	9484983
6	2204500
7	504816
8	124166
9	38491
10	17885
11	11057
12	8485
13	7029
14	6003
15	5361
16	4754
17	4179
18	3581
19	3172
20	2663
21	2496
22	2143
23	1799
24	1506
25	1207
26	1042
27	833
28	735
29	598
30	562
31	487
32	404
33	323
34	297
35	236
36	200
37	180
38	150
39	162
40	166
41	161
42	153
43	174
44	126
45	121
46	119
47	127
48	83
49	111
50	99
51	117
52	83
53	109
54	93
55	84
56	86
57	95
58	72
59	77
60	77
61	58
62	50
63	78
64	47
65	37
66	34
67	50
68	31
69	51
70	24
71	31
72	28
73	34
74	23
75	48
76	29
77	26
78	21
79	26
80	33
81	39
82	37
83	18
84	10
85	3
86	6
87	2
88	8
89	4
90	3
91	5
92	8
93	1
94	2
end
