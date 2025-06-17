<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>

generate_bloom.cpp
- batch addition
- batch inversion
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)

[alexander@alexander-home Point_Search_GMP]$ ./generate_bloom
[17:13:18] P_table generated
[17:13:18] Range Start: 54 bits
[17:13:18] Range End  : 55 bits
[17:13:18] Block Width: 2^26
[17:13:18] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[17:13:18] Settings written to file
[17:13:18] Creating bloom1 image with 4 threads
[17:13:18] Creating bloom2 image with 4 threads
[17:15:12] Writing bloom1 image to bloom1.bf
[17:15:13] Writing bloom2 image to bloom2.bf
[17:15:13] Elapsed time: (0)hours (1)minutes (44)seconds


[alexander@alexander-home Point_Search_GMP]$ ./point_search
[17:15:28] S_table generated
[17:15:28] Range Start: 54 bits
[17:15:28] Range End  : 55 bits
[17:15:28] Block Width: 2^26
[17:15:28] Search Pub : 03a5f9c69423c70c64fe162af3936014c1346978dccd681fa06a18edaa24e3f7d5
[17:15:28] Loading Bloomfilter bloom1.bf
[17:15:28] Loading Bloomfilter bloom2.bf
[17:15:28] Search in progress...
[17:15:42] BloomFilter Hit bloom2.bf (Odd Point) [Higher Range Half]
[17:15:42] Privatekey: 0000000000000000000000000000000000000000000000000069fb4a3e8205d5
[17:15:42] Elapsed time: (0)hours (0)minutes (12)seconds

./generate_bloom uses multiple threads and batch inversion to fill in the bloomfilter binary.
to split the space evenly, number of cores needs to be some power of two value.
actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...)
divided by 2
actual cores = 8  8 / 2 = 4 cores

./point_search is totally scalable and has no such restriction.
can set any desirable number of cores to use but divided by 2.
because we have two search paths : addition and subtraction.
setting cores beyond hardware concurrency will not yield any additional performance.

</pre>
