
        set terminal png size 600,400 truecolor
        set output "IP1_input_srt_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP1_input_srt_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
2	307
3	1263
4	2307
5	3563
6	3628
7	3562
8	3567
9	3545
10	3707
11	3851
12	4226
13	4394
14	4549
15	4461
16	4493
17	4566
18	4573
19	4512
20	4695
21	4683
22	4498
23	4666
24	4599
25	4534
26	4548
27	4660
28	4709
29	4562
30	4993
31	5019
32	4933
33	4884
34	4982
35	4803
36	4667
37	4562
38	4520
39	4268
40	4139
41	4010
42	3871
43	4337
44	4279
45	5295
46	3485
47	2400
48	1194
49	407
end
2	346
3	1789
4	3581
5	3956
6	4139
7	4157
8	4234
9	4668
10	4975
11	5069
12	5179
13	5735
14	5879
15	5900
16	5579
17	5550
18	5649
19	5835
20	5903
21	5710
22	5935
23	6002
24	6315
25	6013
26	6220
27	5655
28	5609
29	5874
30	5807
31	6653
32	6284
33	6186
34	6019
35	5750
36	5951
37	5597
38	5167
39	5085
40	5099
41	4850
42	4664
43	4559
44	4571
45	2958
46	1656
47	583
48	188
49	0
end
