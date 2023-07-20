
        set terminal png size 600,400 truecolor
        set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_input_srt_filtered_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP1_input_srt_filtered_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
2	277
3	1124
4	2106
5	3172
6	3193
7	3175
8	3186
9	3203
10	3356
11	3493
12	3838
13	3960
14	4106
15	4045
16	4066
17	4190
18	4230
19	4169
20	4319
21	4263
22	4094
23	4310
24	4224
25	4145
26	4174
27	4269
28	4321
29	4178
30	4511
31	4581
32	4492
33	4432
34	4465
35	4327
36	4247
37	4103
38	4087
39	3799
40	3665
41	3545
42	3360
43	3611
44	3626
45	4369
46	2925
47	2030
48	980
49	326
end
2	280
3	1516
4	2766
5	2973
6	3196
7	3219
8	3297
9	3517
10	3777
11	3830
12	4102
13	4377
14	4459
15	4525
16	4606
17	4589
18	4533
19	4710
20	4769
21	4723
22	4849
23	4628
24	4793
25	4697
26	4918
27	4537
28	4577
29	4800
30	4671
31	5278
32	5038
33	4916
34	4724
35	4603
36	4556
37	4338
38	3992
39	3977
40	3838
41	3567
42	3471
43	3404
44	3467
45	2438
46	1358
47	435
48	112
49	0
end
