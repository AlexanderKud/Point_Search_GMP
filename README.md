<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[15:56:54] P_table generated
[15:56:54] Range Start: 56 bits
[15:56:54] Range End  : 57 bits
[15:56:54] Block Width: 2^28
[15:56:54] Search Pub : 038188b8729bea9a34ad7b33d44c85f19d93258bab531bf1a7fa9740094a5b1cc0
[15:56:54] Settings written to file
[15:56:54] Creating BloomFile1 with 4 threads
[15:56:54] Creating BloomFile2 with 4 threads
[16:04:56] Writing BloomFile1 to bloom1.bf
[16:04:56] Writing BloomFile2 to bloom2.bf
[16:04:57] Elapsed time: (0)hours (8)minutes (3)seconds


[alexander@alexander-home Point_Search_GMP]$ ./point_search
[16:05:44] S_table generated
[16:05:44] Range Start: 56 bits
[16:05:44] Range End  : 57 bits
[16:05:44] Block Width: 2^28
[16:05:44] Search Pub : 038188b8729bea9a34ad7b33d44c85f19d93258bab531bf1a7fa9740094a5b1cc0
[16:05:44] Loading Bloomfilter bloom1.bf
[16:05:44] Loading Bloomfilter bloom2.bf
[16:05:45] Search in progress...
[16:05:59] BloomFilter Hit bloom1.bf (Even Point) [Higher Range Half]
[16:05:59] Privatekey: 00000000000000000000000000000000000000000000000001e3d553b1d7bc84
[16:05:59] Elapsed time: (0)hours (0)minutes (14)seconds

./generate_bloom uses multiple threads to fill in the bloomfilter binary.
to split the space evenly number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...) divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use.

</pre>
