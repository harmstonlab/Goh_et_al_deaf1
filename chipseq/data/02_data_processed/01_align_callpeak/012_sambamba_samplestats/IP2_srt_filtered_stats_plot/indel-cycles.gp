
        set terminal png size 600,400 truecolor
        set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_srt_filtered_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP2_srt_filtered_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
1	0
2	318
3	1094
4	2009
5	3204
6	3032
7	2945
8	3069
9	3183
10	3380
11	3435
12	3669
13	3779
14	3964
15	4017
16	3902
17	3953
18	4031
19	4076
20	4178
21	4132
22	4172
23	4178
24	4053
25	4063
26	4070
27	4072
28	4094
29	4214
30	4302
31	4358
32	4415
33	4495
34	4207
35	4363
36	4101
37	4017
38	3814
39	3655
40	3501
41	3577
42	3370
43	3629
44	3939
45	5042
46	3102
47	2504
48	1343
49	447
end
1	1
2	260
3	1548
4	2690
5	2889
6	3021
7	3161
8	3237
9	3462
10	3600
11	3791
12	4103
13	4158
14	4230
15	4274
16	4355
17	4396
18	4516
19	4658
20	4603
21	4487
22	4506
23	4599
24	4553
25	4614
26	4572
27	4564
28	4491
29	4612
30	4606
31	4810
32	4780
33	4769
34	4589
35	4458
36	4306
37	4099
38	3978
39	3847
40	3683
41	3498
42	3344
43	3479
44	3705
45	2473
46	1603
47	466
48	111
49	0
end
