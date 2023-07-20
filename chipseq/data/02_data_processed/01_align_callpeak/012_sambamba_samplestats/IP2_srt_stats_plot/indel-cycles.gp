
        set terminal png size 600,400 truecolor
        set output "IP2_srt_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP2_srt_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
1	0
2	358
3	1256
4	2222
5	3685
6	3543
7	3348
8	3523
9	3614
10	3786
11	3838
12	4106
13	4279
14	4433
15	4546
16	4392
17	4456
18	4434
19	4488
20	4677
21	4618
22	4642
23	4643
24	4485
25	4520
26	4554
27	4549
28	4551
29	4766
30	4927
31	4916
32	4923
33	5054
34	4751
35	4951
36	4685
37	4592
38	4543
39	4265
40	4042
41	4129
42	3951
43	4559
44	4841
45	6246
46	3862
47	3001
48	1619
49	538
end
1	1
2	314
3	1928
4	3748
5	4140
6	4242
7	4422
8	4515
9	5005
10	5165
11	5318
12	5492
13	6041
14	6086
15	5930
16	5629
17	5761
18	5971
19	6046
20	6059
21	5809
22	5958
23	6441
24	6654
25	6169
26	6048
27	5978
28	5954
29	6141
30	6107
31	6490
32	6362
33	6553
34	6458
35	6045
36	6059
37	5844
38	5507
39	5240
40	5333
41	5327
42	4944
43	5066
44	5206
45	3080
46	1953
47	659
48	230
49	0
end
