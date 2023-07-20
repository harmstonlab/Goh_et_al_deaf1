
            set terminal png size 600,400 truecolor
            set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_srt_filtered_stats/coverage.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Number of mapped bases"
            set xlabel "Coverage"
            set log y
            set style fill solid border -1
            set title "IP1_srt_filtered_stats.txt" noenhanced
            set xrange [:5]
            plot '-' with lines notitle
        1	776328081
2	321317720
3	97556370
4	24207694
5	5340990
6	1173091
7	320918
8	131443
9	77941
10	55143
11	42534
12	33947
13	27606
14	22385
15	18273
16	15226
17	12606
18	10655
19	8402
20	6807
21	5449
22	4437
23	3525
24	3166
25	2861
26	2542
27	2357
28	2014
29	1893
30	1719
31	1572
32	1506
33	1451
34	1303
35	1149
36	1080
37	1066
38	962
39	977
40	952
41	819
42	786
43	658
44	603
45	583
46	563
47	402
48	422
49	395
50	409
51	346
52	298
53	299
54	239
55	262
56	232
57	214
58	217
59	183
60	187
61	194
62	172
63	159
64	137
65	97
66	69
67	81
68	58
69	39
70	39
71	26
72	36
73	30
74	32
75	24
76	35
77	36
78	65
79	64
80	56
81	29
82	12
83	14
84	16
85	12
86	5
87	3
88	4
89	7
90	5
91	9
92	1
93	1
end
