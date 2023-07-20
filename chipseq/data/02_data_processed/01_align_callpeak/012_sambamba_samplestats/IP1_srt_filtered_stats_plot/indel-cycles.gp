
        set terminal png size 600,400 truecolor
        set output "/home/qianhui/deaf1_chipseq/data/02_data_processed/samplestats/IP1_srt_filtered_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP1_srt_filtered_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
2	261
3	1008
4	1865
5	2904
6	2780
7	2778
8	2728
9	2835
10	2996
11	3156
12	3310
13	3582
14	3528
15	3577
16	3660
17	3623
18	3627
19	3649
20	3684
21	3804
22	3748
23	3821
24	3757
25	3787
26	3962
27	3832
28	3816
29	3720
30	3990
31	4022
32	4020
33	4044
34	4028
35	3909
36	3737
37	3654
38	3560
39	3332
40	3185
41	3143
42	2996
43	3168
44	3365
45	4036
46	2534
47	1935
48	997
49	321
end
2	258
3	1360
4	2474
5	2719
6	2821
7	2858
8	2931
9	3122
10	3303
11	3491
12	3729
13	3813
14	3966
15	3981
16	4003
17	4155
18	4091
19	4154
20	4154
21	4185
22	4215
23	4186
24	4204
25	4294
26	4166
27	4126
28	4107
29	4148
30	4165
31	4458
32	4468
33	4388
34	4317
35	4213
36	3967
37	3838
38	3693
39	3551
40	3415
41	3175
42	3101
43	3091
44	3234
45	2034
46	1300
47	413
48	85
49	0
end
