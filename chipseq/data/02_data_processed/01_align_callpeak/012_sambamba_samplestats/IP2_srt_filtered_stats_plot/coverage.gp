
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_srt_filtered_stats/coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "IP2_srt_filtered_stats.txt" noenhanced
            set xrange [:6]
            plot '-' with lines notitle
        1	781780210
2	352147449
3	117454119
4	32296366
5	7965744
6	1943212
7	561640
8	229230
9	128112
10	87502
11	65498
12	51215
13	41957
14	34758
15	29071
16	25129
17	21019
18	17926
19	14988
20	12450
21	10585
22	8724
23	7226
24	6025
25	5152
26	4301
27	3632
28	3259
29	2900
30	2710
31	2219
32	2098
33	1866
34	1768
35	1624
36	1541
37	1438
38	1400
39	1312
40	1120
41	998
42	968
43	875
44	869
45	917
46	875
47	818
48	646
49	666
50	601
51	509
52	529
53	454
54	377
55	341
56	352
57	352
58	348
59	371
60	329
61	317
62	321
63	262
64	245
65	239
66	221
67	200
68	198
69	172
70	176
71	182
72	141
73	132
74	145
75	126
76	101
77	70
78	103
79	72
80	54
81	41
82	17
83	15
84	17
85	19
86	26
87	20
88	27
89	13
90	9
92	4
94	1
end
