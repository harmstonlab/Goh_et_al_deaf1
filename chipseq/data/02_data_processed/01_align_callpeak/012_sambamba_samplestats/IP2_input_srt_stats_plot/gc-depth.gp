
            set terminal png size 600,500 truecolor
            set output "IP2_input_srt_stats/gc-depth.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Mapped depth"
            set xlabel "Percentile of mapped sequence ordered by GC content"
            set x2label "GC Content [%]"
            set title "IP2_input_srt_stats.txt" noenhanced
            set x2tics ("30" 0.087,"40" 41.749,"50" 96.034)
            set xtics nomirror
            set xrange [0.1:99.9]

            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#dedede" t '10-90th percentile' , \
                 '-' using 1:2:3 with filledcurve lt 1 lc rgb "#bbdeff" t '25-75th percentile' , \
                 '-' using 1:2 with lines lc rgb "#0084ff" t 'Median'
        0.002	0.000	0.002
0.005	0.002	0.002
0.005	0.002	0.002
0.006	0.002	0.002
0.007	0.002	0.002
0.009	0.002	0.002
0.011	0.002	0.005
0.012	0.002	0.002
0.013	0.005	0.005
0.017	0.002	0.005
0.019	0.005	0.005
0.022	0.002	0.010
0.024	0.005	0.005
0.030	0.002	0.005
0.034	0.005	0.005
0.055	0.002	0.005
0.063	0.005	0.308
0.070	0.002	0.013
0.074	0.005	0.237
0.087	0.002	0.013
0.097	0.007	0.923
0.127	0.002	1.500
0.247	0.463	1.240
0.877	0.535	1.335
2.887	0.570	1.285
6.969	0.600	1.278
13.445	0.615	1.260
21.904	0.608	1.268
31.751	0.603	1.280
41.749	0.598	1.280
51.167	0.600	1.278
59.542	0.620	1.258
67.153	0.640	1.247
73.869	0.647	1.242
79.732	0.650	1.210
84.580	0.645	1.212
88.578	0.630	1.175
91.636	0.608	1.165
94.059	0.608	1.155
96.034	0.575	1.133
97.519	0.620	1.170
98.574	0.620	1.192
99.247	0.632	1.188
99.653	0.592	1.180
99.853	0.637	1.230
99.939	0.572	1.170
99.972	0.565	1.048
99.980	0.538	1.038
99.989	0.665	1.140
99.991	0.002	0.913
99.992	0.068	0.068
99.993	0.002	0.863
99.994	0.973	0.973
99.996	0.038	0.762
99.997	14.170	14.170
99.999	0.002	0.015
100.000	0.007	0.007
end
0.002	0.000	0.002
0.005	0.002	0.002
0.005	0.002	0.002
0.006	0.002	0.002
0.007	0.002	0.002
0.009	0.002	0.002
0.011	0.002	0.005
0.012	0.002	0.002
0.013	0.005	0.005
0.017	0.002	0.002
0.019	0.005	0.005
0.022	0.002	0.002
0.024	0.005	0.005
0.030	0.002	0.002
0.034	0.005	0.005
0.055	0.002	0.002
0.063	0.005	0.007
0.070	0.002	0.005
0.074	0.005	0.013
0.087	0.002	0.005
0.097	0.007	0.543
0.127	0.005	0.855
0.247	0.567	1.090
0.877	0.805	1.130
2.887	0.805	1.092
6.969	0.803	1.085
13.445	0.798	1.080
21.904	0.788	1.082
31.751	0.777	1.100
41.749	0.777	1.112
51.167	0.788	1.128
59.542	0.793	1.122
67.153	0.808	1.120
73.869	0.808	1.115
79.732	0.795	1.092
84.580	0.783	1.077
88.578	0.760	1.053
91.636	0.748	1.025
94.059	0.733	1.005
96.034	0.710	0.978
97.519	0.707	0.988
98.574	0.705	1.005
99.247	0.702	0.978
99.653	0.692	0.993
99.853	0.702	1.015
99.939	0.680	0.998
99.972	0.685	0.980
99.980	0.558	0.993
99.989	0.870	1.107
99.991	0.002	0.913
99.992	0.068	0.068
99.993	0.002	0.863
99.994	0.973	0.973
99.996	0.038	0.762
99.997	14.170	14.170
99.999	0.002	0.015
100.000	0.007	0.007
end
0.002	0.000
0.005	0.002
0.005	0.002
0.006	0.002
0.007	0.002
0.009	0.002
0.011	0.002
0.012	0.002
0.013	0.005
0.017	0.002
0.019	0.005
0.022	0.002
0.024	0.005
0.030	0.002
0.034	0.005
0.055	0.002
0.063	0.005
0.070	0.002
0.074	0.007
0.087	0.002
0.097	0.010
0.127	0.455
0.247	0.923
0.877	0.967
2.887	0.960
6.969	0.950
13.445	0.933
21.904	0.930
31.751	0.933
41.749	0.940
51.167	0.950
59.542	0.950
67.153	0.950
73.869	0.952
79.732	0.935
84.580	0.915
88.578	0.885
91.636	0.868
94.059	0.845
96.034	0.812
97.519	0.810
98.574	0.803
99.247	0.793
99.653	0.783
99.853	0.805
99.939	0.815
99.972	0.808
99.980	0.840
99.989	0.995
99.991	0.002
99.992	0.068
99.993	0.002
99.994	0.973
99.996	0.700
99.997	14.170
99.999	0.002
100.000	0.007
end