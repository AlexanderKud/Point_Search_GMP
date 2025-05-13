#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <vector>
#include <thread>
#include <cmath>
#include <gmpxx.h>
#include <gmp.h>

#include "secp256k1/secp256k1.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
namespace fs = filesystem;

inline uint64_t str_to_uint64(std::string const& value) {
  uint64_t result = 0;
  size_t const length = value.size();
  switch (length) {
    case 20:    result += (value[length - 20] - '0') * 10000000000000000000ULL;
    case 19:    result += (value[length - 19] - '0') * 1000000000000000000ULL;
    case 18:    result += (value[length - 18] - '0') * 100000000000000000ULL;
    case 17:    result += (value[length - 17] - '0') * 10000000000000000ULL;
    case 16:    result += (value[length - 16] - '0') * 1000000000000000ULL;
    case 15:    result += (value[length - 15] - '0') * 100000000000000ULL;
    case 14:    result += (value[length - 14] - '0') * 10000000000000ULL;
    case 13:    result += (value[length - 13] - '0') * 1000000000000ULL;
    case 12:    result += (value[length - 12] - '0') * 100000000000ULL;
    case 11:    result += (value[length - 11] - '0') * 10000000000ULL;
    case 10:    result += (value[length - 10] - '0') * 1000000000ULL;
    case  9:    result += (value[length -  9] - '0') * 100000000ULL;
    case  8:    result += (value[length -  8] - '0') * 10000000ULL;
    case  7:    result += (value[length -  7] - '0') * 1000000ULL;
    case  6:    result += (value[length -  6] - '0') * 100000ULL;
    case  5:    result += (value[length -  5] - '0') * 10000ULL;
    case  4:    result += (value[length -  4] - '0') * 1000ULL;
    case  3:    result += (value[length -  3] - '0') * 100ULL;
    case  2:    result += (value[length -  2] - '0') * 10ULL;
    case  1:    result += (value[length -  1] - '0');
  }
  return result;
}

auto main() -> int {
    
    auto start = std::chrono::high_resolution_clock::now();
    Secp256k1 *secp256k1 = new Secp256k1(); secp256k1->Init();
    
    fs::path current_path = fs::current_path();
    auto file_list = get_files_in_directory(current_path);
    vector<string> targets = {"settings1.txt", "settings2.txt", "bloom1.bf", "bloom2.bf"};
    for (auto i : file_list) {
        for (auto t : targets) {
            if (i == t) { std::remove(t.c_str()); }
        }
    }

    mpz_class pk; pk = 1;
    vector<Point> P_table;
    Point P;
    for (int i = 0; i < 256; i++)
    {
        P = secp256k1->ScalarMultiplication(pk);
        P_table.push_back(P);
        mpz_mul_ui(pk.get_mpz_t(), pk.get_mpz_t(), 2);
    }
    print_time(); cout << "P_table generated" << endl;

    uint64_t range_start, range_end, block_width;
    string temp, search_pub;
    ifstream inFile("settings.txt");
    getline(inFile, temp); range_start = str_to_uint64(temp);
    getline(inFile, temp); range_end = str_to_uint64(temp);
    getline(inFile, temp); block_width = str_to_uint64(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    Point start_point, end_point, puzzle_point, puzzle_point_05, puzzle_point_divide2;
    Point first_point, second_point, P1, P2, Q1, Q2;
    mpz_class stride_sum; stride_sum = 0;
    start_point = P_table[range_start];
    end_point   = P_table[range_end];
    mpz_class div2;
    div2 = 2;
    Point point_05 {mpz_class("00000000000000000000003b78ce563f89a0ed9414f5aa28ad0d96d6795f9c63", 16),
                    mpz_class("c0c686408d517dfd67c2367651380d00d126e4229631fd03f8ff35eef1a61e3c", 16)};
    puzzle_point = secp256k1->ParsePublicKeyHex(search_pub);
    puzzle_point_05 = secp256k1->AddPoints(puzzle_point, point_05);

    puzzle_point_divide2 = secp256k1->PointDivision(puzzle_point, div2);

    first_point  = P_table[range_start - 1];
    second_point = P_table[range_start - 2];

    P1 = secp256k1->SubtractPoints(puzzle_point_divide2, first_point);
    P2 = secp256k1->SubtractPoints(puzzle_point_divide2, second_point);
    Q1 = secp256k1->AddPoints(P1, P2);
    Q2 = secp256k1->AddPoints(puzzle_point_divide2, Q1);

    char strideS[68];
    gmp_snprintf(strideS, 67, "%0.64Zx", stride_sum.get_mpz_t());
    string cpub_w = secp256k1->GetPublicKeyHex(Q2);
    
    ofstream outFile1;
    outFile1.open("settings1.txt", ios::app);
    outFile1 << cpub_w <<'\n';
    outFile1 << strideS << '\n';
    outFile1.close();
    
    ofstream outFile2;
    outFile2.open("settings2.txt", ios::app);
    outFile2 << cpub_w <<'\n';
    outFile2 << strideS << '\n';
    outFile2.close();
    
    print_time(); cout << "Settings written to file" << endl;
    
    using filter = boost::bloom::filter<std::string, 32>;
    uint64_t n_elements = uint64_t(pow(2, block_width) * 1.0);
    double error = 0.0000000001;
    int n_cores = 4;  //actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...) divided by 2
    uint64_t count = uint64_t(pow(2, block_width) / n_cores); // actual cores = 8  8 / 2 = 4 cores for each lambda function
    mpz_class add_key;                                        // should be some power of two to evenly divide the space between threads
    add_key = count;
    Point Add_Point = secp256k1->ScalarMultiplication(add_key);
    
    auto bloom_create1 = [&]() {
        string bloomfile = "bloom1.bf";
        Point P = puzzle_point;
        vector<Point> starting_points;
        for (int i = 0; i < n_cores; i++) {
            starting_points.push_back(P);
            P = secp256k1->AddPoints(P, Add_Point);
        }

        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) {            
            Point current = start_point;
            for (uint64_t i = 0; i < count; i++) {
                bf.insert(secp256k1->GetPublicKeyHex(current));
                current = secp256k1->AddPoints(current, secp256k1->G);
            }
        };
        
        std::thread myThreads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }
    
        print_time(); cout << "Creating BloomFile1 with " << n_cores << " threads" << '\n';

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join();
        }

        print_time(); cout << "Writing BloomFile1 to bloom1.bf" << '\n';
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };

    auto bloom_create2 = [&]() {
        string bloomfile = "bloom2.bf";
        Point P = puzzle_point_05;
        vector<Point> starting_points;
        for (int i = 0; i < n_cores; i++) {
            starting_points.push_back(P);
            P = secp256k1->AddPoints(P, Add_Point);
        }
        
        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) {            
            Point current = start_point;
            for (uint64_t i = 0; i < count; i++) {
                bf.insert(secp256k1->GetPublicKeyHex(current));
                current = secp256k1->AddPoints(current, secp256k1->G);
            }
        };
        
        std::thread myThreads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }
    
        print_time(); cout << "Creating BloomFile2 with " << n_cores << " threads" << '\n';

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join();
        }

        print_time(); cout << "Writing BloomFile2 to bloom2.bf" << '\n'; 
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };

    std::thread thread1(bloom_create1);
    std::thread thread2(bloom_create2);
    
    thread1.join();
    thread2.join();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
}
