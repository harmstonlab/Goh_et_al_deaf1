
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_input_srt_filtered_stats/coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "IP2_input_srt_filtered_stats.txt" noenhanced
            set xrange [:5]
            plot '-' with lines notitle
        1	781090901
2	366680526
3	127497165
4	36130065
5	8879583
6	1975996
7	432074
8	103469
9	33067
10	16410
11	10904
12	8554
13	7421
14	6037
15	5374
16	4515
17	3768
18	3623
19	3150
20	2691
21	2316
22	1890
23	1728
24	1467
25	1220
26	1063
27	916
28	695
29	562
30	489
31	391
32	342
33	303
34	239
35	271
36	237
37	235
38	224
39	192
40	165
41	192
42	164
43	144
44	136
45	124
46	109
47	132
48	108
49	124
50	105
51	107
52	111
53	108
54	66
55	69
56	76
57	88
58	78
59	75
60	66
61	63
62	74
63	64
64	50
65	46
66	45
67	45
68	54
69	57
70	39
71	39
72	42
73	53
74	35
75	40
76	55
77	54
78	28
79	35
80	18
81	25
82	27
83	19
84	7
85	5
86	7
87	7
88	7
89	3
90	5
91	4
92	3
93	1
94	2
end
