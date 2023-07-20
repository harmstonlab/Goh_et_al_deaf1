
        set terminal png size 600,400 truecolor
        set output "IP1_srt_stats/indel-cycles.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style line 4 linetype 4  linecolor rgb "blue"
        set style increment user
        set ylabel "Indel count"
        set xlabel "Read Cycle"
        set title "IP1_srt_stats.txt" noenhanced
    plot '-' w l ti 'Insertions', '' w l ti 'Deletions'
2	294
3	1147
4	2085
5	3412
6	3286
7	3202
8	3082
9	3224
10	3387
11	3568
12	3769
13	4040
14	4008
15	4019
16	4085
17	4046
18	3996
19	4035
20	4141
21	4286
22	4198
23	4267
24	4176
25	4240
26	4399
27	4330
28	4287
29	4186
30	4490
31	4549
32	4523
33	4586
34	4580
35	4465
36	4274
37	4198
38	4157
39	3906
40	3698
41	3614
42	3538
43	4102
44	4099
45	5018
46	3128
47	2307
48	1226
49	401
end
2	316
3	1696
4	3567
5	3967
6	4008
7	4077
8	4186
9	4638
10	4818
11	4945
12	5046
13	5559
14	5743
15	5522
16	5268
17	5403
18	5482
19	5525
20	5559
21	5452
22	5676
23	6053
24	6111
25	5911
26	5670
27	5630
28	5509
29	5521
30	5592
31	6234
32	5994
33	6036
34	6013
35	5754
36	5755
37	5545
38	5138
39	4974
40	5049
41	4942
42	4700
43	4614
44	4560
45	2637
46	1601
47	554
48	151
49	0
end
