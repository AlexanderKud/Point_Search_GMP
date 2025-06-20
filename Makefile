default:
	g++ -O3 -m64 -mssse3 -c -Wno-write-strings secp256k1/secp256k1.cpp -o secp256k1.o
	g++ -O3 -march=native -c util/util.cpp -o util.o
	g++ -O3 -march=native -c -Wno-unused-result generate_bloom.cpp -o generate_bloom.o
	g++ -O3 -march=native -c -Wno-unused-result point_search.cpp -o point_search.o
	g++ -o generate_bloom generate_bloom.o secp256k1.o util.o -lgmpxx -lgmp -fopenmp
	g++ -o point_search point_search.o secp256k1.o util.o -lgmpxx -lgmp
	rm *.o
