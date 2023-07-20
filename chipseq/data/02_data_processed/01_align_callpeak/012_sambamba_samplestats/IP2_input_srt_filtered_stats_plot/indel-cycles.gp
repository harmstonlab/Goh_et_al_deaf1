
        set terminal png size 600,400 truecolor
        set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP2_input_srt_filtered_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP2_input_srt_filtered_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
1	0
2	305
3	1131
4	2146
5	3434
6	3314
7	3379
8	3301
9	3357
10	3489
11	3719
12	3813
13	4111
14	4202
15	4271
16	4117
17	4167
18	4318
19	4366
20	4242
21	4362
22	4380
23	4326
24	4444
25	4253
26	4426
27	4340
28	4479
29	4336
30	4463
31	4638
32	4576
33	4785
34	4535
35	4471
36	4355
37	4237
38	4100
39	3956
40	3649
41	3530
42	3374
43	3539
44	3615
45	4065
46	2848
47	1967
48	1001
49	320
end
1	1
2	293
3	1668
4	2943
5	3225
6	3257
7	3253
8	3406
9	3570
10	3821
11	4124
12	4375
13	4329
14	4429
15	4442
16	4476
17	4607
18	4592
19	4655
20	4817
21	4836
22	4869
23	4773
24	4806
25	4762
26	4794
27	4810
28	4760
29	4807
30	4753
31	5198
32	5011
33	5100
34	4817
35	4793
36	4540
37	4350
38	4259
39	3989
40	3806
41	3694
42	3501
43	3488
44	3315
45	2320
46	1442
47	421
48	77
49	0
end
