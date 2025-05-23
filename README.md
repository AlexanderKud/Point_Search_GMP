<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
  
[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[23:06:16] P_table generated
[23:06:16] Range Start: 54 bits
[23:06:16] Range End  : 55 bits
[23:06:16] Block Width: 2^26
[23:06:16] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[23:06:16] Settings written to file
[23:06:16] Creating BloomFile2 with 4 threads
[23:06:16] Creating BloomFile1 with 4 threads
[23:08:30] Writing BloomFile1 to bloom1.bf
[23:08:31] Writing BloomFile2 to bloom2.bf
[23:08:31] Elapsed time: (0)hours (2)minutes (15)seconds

[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom_batch
[23:10:09] P_table generated
[23:10:09] Range Start: 54 bits
[23:10:09] Range End  : 55 bits
[23:10:09] Block Width: 2^26
[23:10:09] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[23:10:09] Settings written to file
[23:10:09] Creating BloomFile2 with 4 threads
[23:10:09] Creating BloomFile1 with 4 threads
[23:11:24] Writing BloomFile1 to bloom1.bf
[23:11:25] Writing BloomFile2 to bloom2.bf
[23:11:25] Elapsed time: (0)hours (1)minutes (16)seconds

[alexander@alexander-home Point_Search_GMP]$ ./point_search
[23:09:30] S_table generated
[23:09:30] Range Start: 54 bits
[23:09:30] Range End  : 55 bits
[23:09:30] Block Width: 2^26
[23:09:30] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[23:09:30] Loading Bloomfilter bloom1.bf
[23:09:30] Loading Bloomfilter bloom2.bf
[23:09:30] Search in progress...
[23:10:02] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[23:10:02] Privatekey: 0000000000000000000000000000000000000000000000000069fb4a3e8205d5
[23:10:02] Elapsed time: (0)hours (0)minutes (31)seconds


./generate_bloom uses multiple threads and to fill in the bloomfilter binary.
./generate_bloom_batch uses batch inversion.
to split the space evenly, number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use but divided by 2.
because we have two search paths : addition and subtraction.
setting cores beyond hardware concurrency will not yield any additional performance.

</pre>
