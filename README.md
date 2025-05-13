<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_CPP2]$ ./generate_bloom
[07:48:00] P_table generated
[07:48:00] Range Start: 51 bits
[07:48:00] Range End  : 52 bits
[07:48:00] Block Width: 2^25
[07:48:00] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[07:48:00] Settings written to file
[07:48:01] Creating BloomFile1 with 4 threads
[07:48:01] Creating BloomFile2 with 4 threads
[07:49:04] Writing BloomFile1 to bloom1.bf
[07:49:04] Writing BloomFile2 to bloom2.bf
[07:49:04] Elapsed time: (0)hours (1)minutes (3)seconds

[alexander@alexander-home Point_Search_CPP2]$ ./point_search
[07:59:55] S_table generated
[07:59:55] Range Start: 51 bits
[07:59:55] Range End  : 52 bits
[07:59:55] Block Width: 2^25
[07:59:55] Search Pub : 033195139de0331d7a5cab602c4471f728f2e3fb97ed82f593d49ed30ec3c0ba85
[07:59:55] Loading Bloomfilter bloom1.bf
[07:59:55] Loading Bloomfilter bloom2.bf
[07:59:55] Search in progress...
[08:00:14] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[08:00:14] Privatekey: 000000000000000000000000000000000000000000000000000ad89e2c8e65c3
[08:00:14] Elapsed time: (0)hours (0)minutes (19)seconds

</pre>
