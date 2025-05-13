<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[21:02:40] P_table generated
[21:02:40] Range Start: 51 bits
[21:02:40] Range End  : 52 bits
[21:02:40] Block Width: 2^25
[21:02:40] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[21:02:40] Settings written to file
[21:02:40] Creating BloomFile2 with 4 threads
[21:02:40] Creating BloomFile1 with 4 threads
[21:03:41] Writing BloomFile2 to bloom2.bf
[21:03:41] Writing BloomFile1 to bloom1.bf
[21:03:41] Elapsed time: (0)hours (1)minutes (0)seconds


[alexander@alexander-home Point_Search_GMP]$ ./point_search
[21:03:44] S_table generated
[21:03:44] Range Start: 51 bits
[21:03:44] Range End  : 52 bits
[21:03:44] Block Width: 2^25
[21:03:44] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[21:03:44] Loading Bloomfilter bloom1.bf
[21:03:44] Loading Bloomfilter bloom2.bf
[21:03:44] Search in progress...
[21:04:03] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[21:04:03] Privatekey: 000000000000000000000000000000000000000000000000000ad89e2c8e65c3
[21:04:03] Elapsed time: (0)hours (0)minutes (18)seconds

</pre>
