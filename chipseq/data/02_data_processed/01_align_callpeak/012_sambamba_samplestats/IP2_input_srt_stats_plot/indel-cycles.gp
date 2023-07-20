
        set terminal png size 600,400 truecolor
        set output "IP2_input_srt_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP2_input_srt_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
1	0
2	350
3	1267
4	2328
5	3867
6	3783
7	3785
8	3669
9	3740
10	3880
11	4104
12	4175
13	4508
14	4671
15	4668
16	4507
17	4567
18	4699
19	4731
20	4653
21	4780
22	4782
23	4745
24	4869
25	4703
26	4809
27	4804
28	4864
29	4766
30	4895
31	5091
32	5024
33	5236
34	5027
35	4980
36	4834
37	4741
38	4553
39	4452
40	4135
41	3993
42	3852
43	4295
44	4270
45	5024
46	3341
47	2323
48	1212
49	365
end
1	1
2	350
3	1968
4	3832
5	4349
6	4304
7	4337
8	4428
9	4866
10	5109
11	5407
12	5561
13	5830
14	5890
15	6036
16	5513
17	5673
18	5819
19	5905
20	6050
21	5928
22	6121
23	6304
24	6509
25	6151
26	6220
27	6019
28	5965
29	6091
30	6037
31	6647
32	6408
33	6548
34	6270
35	6174
36	6108
37	5872
38	5587
39	5312
40	5214
41	5160
42	4881
43	4863
44	4614
45	2873
46	1752
47	561
48	136
49	0
end
