<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[08:05:27] P_table generated
[08:05:27] Range Start: 54 bits
[08:05:27] Range End  : 55 bits
[08:05:27] Block Width: 2^26
[08:05:27] Search Pub : 03a4e2af6d1c191666070fc618106cee91c44936f0082c832692c2043fd8b67d8c
[08:05:27] Settings written to file
[08:05:27] Creating BloomFile1 with 4 threads
[08:05:27] Creating BloomFile2 with 4 threads
[08:07:27] Writing BloomFile1 to bloom1.bf
[08:07:27] Writing BloomFile2 to bloom2.bf
[08:07:27] Elapsed time: (0)hours (2)minutes (0)seconds


[alexander@alexander-home Point_Search_GMP]$ ./point_search
[08:07:32] S_table generated
[08:07:32] Range Start: 54 bits
[08:07:32] Range End  : 55 bits
[08:07:32] Block Width: 2^26
[08:07:32] Search Pub : 03a4e2af6d1c191666070fc618106cee91c44936f0082c832692c2043fd8b67d8c
[08:07:32] Loading Bloomfilter bloom1.bf
[08:07:32] Loading Bloomfilter bloom2.bf
[08:07:32] Search in progress...
[08:08:22] BloomFilter Hit bloom2.bf (Odd Point) [Lower Range Half]
[08:08:22] Privatekey: 00000000000000000000000000000000000000000000000000546c597195960f
[08:08:22] Elapsed time: (0)hours (0)minutes (49)seconds

./generate_bloom uses multiple threads to fill in the bloomfilter binary.
to split the space evenly number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...) divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use.

</pre>
