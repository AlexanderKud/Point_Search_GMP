#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <vector>
#include <algorithm>
#include <thread>
#include <gmpxx.h>
#include <gmp.h>

#include "secp256k1/secp256k1.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;

auto main() -> int {

    Secp256k1 *secp256k1 = new Secp256k1(); secp256k1->Init();
    int cpuCores = 4; // actual number of processing cores divided by 2
    
    mpz_class pk; pk = 1;
    vector<mpz_class> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        mpz_mul_ui(pk.get_mpz_t(), pk.get_mpz_t(), 2);
    }
    print_time(); cout << "S_table generated" << endl;

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

    mpz_class pre_calc_sum;
    mpz_add(pre_calc_sum.get_mpz_t(), S_table[range_start - 1].get_mpz_t(), S_table[range_start - 2].get_mpz_t());
    
    string bloomfile1 = "bloom1.bf";
    string bloomfile2 = "bloom2.bf";
    using filter = boost::bloom::filter<std::string, 32>;
    
    print_time(); cout << "Loading Bloomfilter bloom1.bf" << endl;
    filter bf1;
    std::ifstream in1(bloomfile1, std::ios::binary);
    std::size_t c1;
    in1.read((char*) &c1, sizeof(c1));
    bf1.reset(c1); // restore capacity
    boost::span<unsigned char> s1 = bf1.array();
    in1.read((char*) s1.data(), s1.size()); // load array
    in1.close();

    
    print_time(); cout << "Loading Bloomfilter bloom2.bf" << endl;
    filter bf2;
    std::ifstream in2(bloomfile2, std::ios::binary);
    std::size_t c2;
    in2.read((char*) &c2, sizeof(c2));
    bf2.reset(c2); // restore capacity
    boost::span<unsigned char> s2 = bf2.array();
    in2.read((char*) s2.data(), s2.size()); // load array
    in2.close();
    
    auto pow10_nums = break_down_to_pow10(uint64_t(pow(2, block_width)));
    vector<Point> pow10_points;
    mpz_class pow_key;
    for (auto n : pow10_nums) {
        pow_key = n;
        pow10_points.push_back(secp256k1->ScalarMultiplication(pow_key));
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto addition_search = [&]() {
        int save_counter = 0;
        string temp;
        Point start_point, stride_point, calc_point;
        mpz_class stride_sum, stride;
        ifstream inFile("settings1.txt");
        getline(inFile, temp);
        start_point = secp256k1->ParsePublicKeyHex(trim(temp));
        getline(inFile, temp);
        stride_sum = mpz_class(trim(temp).data(), 16);
        inFile.close();
        
        stride = pow(2, block_width);
        stride_point = secp256k1->ScalarMultiplication(stride);
        
        int n_cores = cpuCores;
        
        mpz_class offset_Step, vector_Num;
        mpz_fdiv_q_ui(offset_Step.get_mpz_t(), S_table[range_start - 2].get_mpz_t(),  n_cores);
        
        vector_Num = 0;
        vector<mpz_class> offset_Nums;
        for (int i = 0; i < n_cores; i++) {
            offset_Nums.push_back(vector_Num);
            mpz_add(vector_Num.get_mpz_t(), vector_Num.get_mpz_t(), offset_Step.get_mpz_t());
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < n_cores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < n_cores; i++) {
            vector_Point = secp256k1->AddPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }

        auto scalable_addition_search = [&](Point starting_Point, int threadIdx, mpz_class offset, mpz_class stride_Sum) {
            Point starting_point = starting_Point;
            mpz_class stride_sum = stride_Sum;
            string cpub;
            
            while (true) {
                
                cpub = secp256k1->GetPublicKeyHex(starting_point);
                if (bf1.may_contain(cpub)) {
                    print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;
                    Point P = starting_point;
                    vector<uint64_t> privkey_num;
                    int index = 0;
                    string cpub1;
                    for (auto p : pow10_points) {
                        int count = 0;
                        cpub1 = secp256k1->GetPublicKeyHex(P);
                        while (bf1.may_contain(cpub1)) {
                            P = secp256k1->SubtractPoints(P, p);
                            cpub1 = secp256k1->GetPublicKeyHex(P);
                            count += 1;
                        }
                        privkey_num.push_back(pow10_nums[index] * (count - 1));
                        P = secp256k1->AddPoints(P, p);
                        index += 1;
                    }
                    mpz_class Int_steps, Int_temp, privkey;
                    uint64_t steps = 0;
                    for (auto i : privkey_num) { steps += i; }
                    Int_steps = steps;
                    mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), offset.get_mpz_t());
                    mpz_sub(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                    mpz_sub(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                    mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                    calc_point = secp256k1->ScalarMultiplication(privkey);
                    if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                        char privkeyStr[68];
                        gmp_snprintf(privkeyStr, 67, "%0.64Zx", privkey.get_mpz_t());
                        print_time(); cout << "Privatekey: " << privkeyStr << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkeyStr << '\n';
                        outFile.close();
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = end - start;
                        auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                        duration -= hours;
                        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                        duration -= minutes;
                        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                        print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                        exit(0);
                    }
                    print_time(); cout << "False Positive" << endl;
                }
                
                if (bf2.may_contain(cpub)) {
                    print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;
                    Point P = starting_point;
                    vector<uint64_t> privkey_num;
                    int index = 0;
                    string cpub2;
                    for (auto p : pow10_points) {
                        int count = 0;
                        cpub2 = secp256k1->GetPublicKeyHex(P);
                        while (bf2.may_contain(cpub2)) {
                            P = secp256k1->SubtractPoints(P, p);
                            cpub2 = secp256k1->GetPublicKeyHex(P);
                            count += 1;
                        }
                        privkey_num.push_back(pow10_nums[index] * (count - 1));
                        P = secp256k1->AddPoints(P, p);
                        index += 1;
                    }
                    mpz_class Int_steps, Int_temp, privkey;
                    uint64_t steps = 0;
                    for (auto i : privkey_num) { steps += i; }
                    Int_steps = steps;
                    mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), offset.get_mpz_t());
                    mpz_sub(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                    mpz_sub(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                    mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                    mpz_add_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 1);
                    calc_point = secp256k1->ScalarMultiplication(privkey);
                    
                    if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                        char privkeyStr[68];
                        gmp_snprintf(privkeyStr, 67, "%0.64Zx", privkey.get_mpz_t());
                        print_time(); cout << "Privatekey: " << privkeyStr << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkeyStr << '\n';
                        outFile.close();
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = end - start;
                        auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                        duration -= hours;
                        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                        duration -= minutes;
                        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                        print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                        exit(0);
                    }
                    print_time(); cout << "False Positive" << endl;
                }
                
                starting_point = secp256k1->AddPoints(starting_point, stride_point);
                mpz_add(stride_sum.get_mpz_t(), stride_sum.get_mpz_t(), stride.get_mpz_t());
                
                if (threadIdx == 0) {
                    save_counter += 1;
                    if (save_counter % 70000000 == 0) {
                        char strideS[68];
                        gmp_snprintf(strideS, 67, "%0.64Zx", stride_sum.get_mpz_t());
                        cpub = secp256k1->GetPublicKeyHex(starting_point);
                        ofstream outFile;
                        outFile.open("settings1.txt");
                        outFile << cpub <<'\n';
                        outFile << strideS << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings1.txt" << endl;
                    }
                }
            }
        };

        std::thread addition_Threads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            addition_Threads[i] = std::thread(scalable_addition_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < n_cores; i++) {
            addition_Threads[i].join();
        }
    };
    
    auto subtraction_search = [&]() {
        int save_counter = 0;
        string temp;
        Point start_point,stride_point, calc_point;
        mpz_class stride_sum, stride;
        ifstream inFile("settings2.txt");
        getline(inFile, temp);
        start_point = secp256k1->ParsePublicKeyHex(trim(temp));
        getline(inFile, temp);
        stride_sum = mpz_class(trim(temp).data(), 16);
        inFile.close();
        
        stride = pow(2, block_width);
        stride_point = secp256k1->ScalarMultiplication(stride);
        
        int n_cores = cpuCores;
        
        mpz_class offset_Step, vector_Num;
        mpz_fdiv_q_ui(offset_Step.get_mpz_t(), S_table[range_start - 2].get_mpz_t(),  n_cores);
        
        vector_Num = 0;
        vector<mpz_class> offset_Nums;
        for (int i = 0; i < n_cores; i++) {
            offset_Nums.push_back(vector_Num);
            mpz_add(vector_Num.get_mpz_t(), vector_Num.get_mpz_t(), offset_Step.get_mpz_t());
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < n_cores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(offset_Nums[i]));
        }
        
        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);
        
        for (int i = 1; i < n_cores; i++) {
            vector_Point = secp256k1->SubtractPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }
        
        auto scalable_subtraction_search = [&](Point starting_Point, int threadIdx, mpz_class offset, mpz_class stride_Sum) {
            Point starting_point = starting_Point;
            mpz_class stride_sum = stride_Sum;
            string cpub;
            while (true) {
                cpub = secp256k1->GetPublicKeyHex(starting_point);
                if (bf1.may_contain(cpub)) {
                    print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;
                    Point P(starting_point);
                    vector<uint64_t> privkey_num;
                    int index = 0;
                    string cpub1;
                    for (auto p : pow10_points) {
                        int count = 0;
                        cpub1 = secp256k1->GetPublicKeyHex(P);
                        while (bf1.may_contain(cpub1)) {
                            P = secp256k1->SubtractPoints(P, p);
                            cpub1 = secp256k1->GetPublicKeyHex(P);
                            count += 1;
                        }
                        privkey_num.push_back(pow10_nums[index] * (count - 1));
                        P = secp256k1->AddPoints(P, p);
                        index += 1;
                    }
                    mpz_class Int_steps, Int_temp, privkey;
                    uint64_t steps = 0;
                    for (auto i : privkey_num) { steps += i; }
                    Int_steps = steps;
                    mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), Int_steps.get_mpz_t());
                    mpz_add(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                    mpz_add(privkey.get_mpz_t(), privkey.get_mpz_t(), offset.get_mpz_t()); //##############
                    mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                    calc_point = secp256k1->ScalarMultiplication(privkey);
                    if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                        char privkeyStr[68];
                        gmp_snprintf(privkeyStr, 67, "%0.64Zx", privkey.get_mpz_t());
                        print_time(); cout << "Privatekey: " << privkeyStr << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkeyStr << '\n';
                        outFile.close();
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = end - start;
                        auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                        duration -= hours;
                        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                        duration -= minutes;
                        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                        print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                        exit(0);
                    }
                    print_time(); cout << "False Positive" << endl;
                }
                
                if (bf2.may_contain(cpub)) {
                    print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;
                    Point P = starting_point;
                    vector<uint64_t> privkey_num;
                    int index = 0;
                    string cpub2;
                    for (auto p : pow10_points) {
                        int count = 0;
                        cpub2 = secp256k1->GetPublicKeyHex(P);
                        while (bf2.may_contain(cpub2)) {
                            P = secp256k1->SubtractPoints(P, p);
                            cpub2 = secp256k1->GetPublicKeyHex(P);
                            count += 1;
                        }
                        privkey_num.push_back(pow10_nums[index] * (count - 1));
                        P = secp256k1->AddPoints(P, p);
                        index += 1;
                    }
                    mpz_class Int_steps, Int_temp, privkey;
                    uint64_t steps = 0;
                    for (auto i : privkey_num) { steps += i; }
                    Int_steps = steps;
                    mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), Int_steps.get_mpz_t());
                    mpz_add(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                    mpz_add(privkey.get_mpz_t(), privkey.get_mpz_t(), offset.get_mpz_t());
                    mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                    mpz_add_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 1);
                    calc_point = secp256k1->ScalarMultiplication(privkey);
                    if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                        char privkeyStr[68];
                        gmp_snprintf(privkeyStr, 67, "%0.64Zx", privkey.get_mpz_t());
                        print_time(); cout << "Privatekey: " << privkeyStr << endl;
                        ofstream outFile;
                        outFile.open("found.txt", ios::app);
                        outFile << privkeyStr << '\n';
                        outFile.close();
                        auto end = std::chrono::high_resolution_clock::now();
                        auto duration = end - start;
                        auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                        duration -= hours;
                        auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
                        duration -= minutes;
                        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
                        print_time(); cout << "Elapsed time: (" << hours.count() << ")hours (" << minutes.count() << ")minutes (" << seconds.count() << ")seconds\n";
                        exit(0);
                    }
                    print_time(); cout << "False Positive" << endl;
                }
                starting_point = secp256k1->SubtractPoints(starting_point, stride_point);
                mpz_add(stride_sum.get_mpz_t(), stride_sum.get_mpz_t(), stride.get_mpz_t());
                if (threadIdx == 0) {
                    save_counter += 1;
                    if (save_counter % 70000000 == 0) {
                        char strideS[68];
                        gmp_snprintf(strideS, 67, "%0.64Zx", stride_sum.get_mpz_t());
                        cpub = secp256k1->GetPublicKeyHex(starting_point);
                        ofstream outFile;
                        outFile.open("settings2.txt");
                        outFile << cpub <<'\n';
                        outFile << strideS << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings2.txt" << endl;
                    }
                }
            }
        };

        std::thread subtraction_Threads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            subtraction_Threads[i] = std::thread(scalable_subtraction_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < n_cores; i++) {
            subtraction_Threads[i].join();
        }
    };

    print_time(); cout << "Search in progress..." << endl;
    
    std::thread thread1(addition_search);
    std::thread thread2(subtraction_search);
    
    thread1.join();
    thread2.join();
}
