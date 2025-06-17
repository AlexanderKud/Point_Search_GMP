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

static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion(one ModInv for the entire group) using IntGroup class
const mpz_class Fp = mpz_class("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10);

auto main() -> int {

    Secp256k1 *secp256k1 = new Secp256k1(); secp256k1->Init(); // initialize secp256k1 context
    int cpuCores = 4; // actual number of processing cores divided by 2
    
    mpz_class pk; pk = 1; // generating power of two values (2^0..2^256) table
    vector<mpz_class> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        mpz_mul_ui(pk.get_mpz_t(), pk.get_mpz_t(), 2);
    }
    print_time(); cout << "S_table generated" << endl;

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter and a stride size to walk the range
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

    mpz_class pre_calc_sum; // precalculated sum for private key recovering
    mpz_add(pre_calc_sum.get_mpz_t(), S_table[range_start - 1].get_mpz_t(), S_table[range_start - 2].get_mpz_t());
    
    using filter = boost::bloom::filter<std::string, 32>;
    
    string bloomfile1 = "bloom1.bf";
    print_time(); cout << "Loading Bloomfilter bloom1.bf" << endl;
    filter bf1;
    std::ifstream in1(bloomfile1, std::ios::binary);
    std::size_t c1;
    in1.read((char*) &c1, sizeof(c1));
    bf1.reset(c1); // restore capacity
    boost::span<unsigned char> s1 = bf1.array();
    in1.read((char*) s1.data(), s1.size()); // load array
    in1.close();

    string bloomfile2 = "bloom2.bf";
    print_time(); cout << "Loading Bloomfilter bloom2.bf" << endl;
    filter bf2;
    std::ifstream in2(bloomfile2, std::ios::binary);
    std::size_t c2;
    in2.read((char*) &c2, sizeof(c2));
    bf2.reset(c2); // restore capacity
    boost::span<unsigned char> s2 = bf2.array();
    in2.read((char*) s2.data(), s2.size()); // load array
    in2.close();
    
    auto pow10_nums = break_down_to_pow10(uint64_t(pow(2, block_width))); // decomposing the 2^block_width to the power of ten values
    vector<Point> pow10_points;                                           // to get the index of the bloomfilter element fast
    mpz_class pow_key;
    for (auto& n : pow10_nums) { // calculating points corresponding to the decomposition components
        pow_key = n;
        pow10_points.push_back(secp256k1->ScalarMultiplication(pow_key));
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    auto addition_search = [&]() { // addition search for the case when the starting point is behind the target point after calculations
        int save_counter = 0;      // the closer the target point to the center of the range from either side
        string temp;               // the faster collision will happen
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
        
        //start splitting the search according to the chosen number of cpu cores
        mpz_class offset_Step, vector_Num;
        mpz_fdiv_q_ui(offset_Step.get_mpz_t(), S_table[range_start - 2].get_mpz_t(),  cpuCores);
        
        vector_Num = 0;
        vector<mpz_class> offset_Nums;
        for (int i = 0; i < cpuCores; i++) {
            offset_Nums.push_back(vector_Num);
            mpz_add(vector_Num.get_mpz_t(), vector_Num.get_mpz_t(), offset_Step.get_mpz_t());
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < cpuCores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(offset_Nums[i]));
        }

        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);

        for (int i = 1; i < cpuCores; i++) {
            vector_Point = secp256k1->AddPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }
    
        Point addPoints[POINTS_BATCH_SIZE]; // array for batch addition points       
        Point batch_Add = secp256k1->DoublePoint(stride_point);
        addPoints[0] = stride_point;
        addPoints[1] = batch_Add;
        for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in batch addition points array with points
        {
            batch_Add = secp256k1->AddPoints(batch_Add, stride_point);
            addPoints[i] = batch_Add;
        }
        // scalable lambda gets its chunk to search through
        auto scalable_addition_search = [&](Point starting_Point, int threadIdx, mpz_class offset, mpz_class stride_Sum) {

            Point startPoint = starting_Point;
            Point BloomP; // point for insertion of the batch into the bloomfilter
            mpz_class stride_sum = stride_Sum;
            string cpub, xc, xc_sub;
            int index, count;
            vector<uint64_t> privkey_num;
            uint64_t steps;
            mpz_class Int_steps, Int_temp, privkey;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            mpz_class deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            mpz_class pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            mpz_class pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            mpz_class deltaY; // values to store the results of points addition formula
            mpz_class slope[POINTS_BATCH_SIZE];

            mpz_class batch_stride, batch_index;
            mpz_mul_ui(batch_stride.get_mpz_t(), stride.get_mpz_t(), POINTS_BATCH_SIZE);
            
            while (true) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    mpz_sub(deltaX[i].get_mpz_t(), startPoint.x.get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(deltaX[i].get_mpz_t(), deltaX[i].get_mpz_t(), Fp.get_mpz_t()); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv(deltaX); // assign array deltaX to modGroup for batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                    mpz_mul(slope[i].get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                    mpz_mod(slope[i].get_mpz_t(), slope[i].get_mpz_t(), Fp.get_mpz_t());

                    mpz_mul(pointBatchX[i].get_mpz_t(), slope[i].get_mpz_t(), slope[i].get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), startPoint.x.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());

                }
                
                mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                mpz_mul(slope[i].get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                mpz_mod(slope[i].get_mpz_t(), slope[i].get_mpz_t(), Fp.get_mpz_t());

                mpz_mul(pointBatchX[i].get_mpz_t(), slope[i].get_mpz_t(), slope[i].get_mpz_t());
                mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), startPoint.x.get_mpz_t());
                mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());
        
                mpz_sub(pointBatchY[i].get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                mpz_mul(pointBatchY[i].get_mpz_t(), slope[i].get_mpz_t(), pointBatchY[i].get_mpz_t());
                mpz_sub(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), startPoint.y.get_mpz_t());
                mpz_mod(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), Fp.get_mpz_t());
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {

                    xc = secp256k1->GetXHex(pointBatchX[i]);

                    if (bf1.may_contain(xc)) {

                        print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Lower Range Half]" << endl;

                        BloomP.x = pointBatchX[i];
                        mpz_sub(BloomP.y.get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                        mpz_mul(BloomP.y.get_mpz_t(), slope[i].get_mpz_t(), BloomP.y.get_mpz_t());
                        mpz_sub(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), startPoint.y.get_mpz_t());
                        mpz_mod(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), Fp.get_mpz_t());

                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) { // getting the index of the element in the bloomfilter
                            count = 0;
                            xc_sub = secp256k1->GetXHex(BloomP.x);
                            while (bf1.may_contain(xc_sub)) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                xc_sub = secp256k1->GetXHex(BloomP.x);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }

                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; } // we got here the index of the element in the bloomfilter
                        Int_steps = steps; // restoring the private key
                        mpz_mul_ui(batch_index.get_mpz_t(), stride.get_mpz_t(), (i + 1));
                        mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), batch_index.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), offset.get_mpz_t());
                        mpz_sub(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                        mpz_sub(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                        mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2); // we got here the private key
                        calc_point = secp256k1->ScalarMultiplication(privkey);

                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
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
                    
                    if (bf2.may_contain(xc)) {

                        print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Lower Range Half]" << endl;

                        BloomP.x = pointBatchX[i];
                        mpz_sub(BloomP.y.get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                        mpz_mul(BloomP.y.get_mpz_t(), slope[i].get_mpz_t(), BloomP.y.get_mpz_t());
                        mpz_sub(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), startPoint.y.get_mpz_t());
                        mpz_mod(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), Fp.get_mpz_t());

                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            xc_sub = secp256k1->GetXHex(BloomP.x);
                            while (bf2.may_contain(xc_sub)) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                xc_sub = secp256k1->GetXHex(BloomP.x);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }

                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps = steps;
                        mpz_mul_ui(batch_index.get_mpz_t(), stride.get_mpz_t(), (i + 1));
                        mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), batch_index.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), offset.get_mpz_t());
                        mpz_sub(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                        mpz_sub(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                        mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                        mpz_add_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 1);
                        calc_point = secp256k1->ScalarMultiplication(privkey);
                        
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
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
                }
                
                startPoint.x = pointBatchX[POINTS_BATCH_SIZE - 1]; // setting the new startPoint for the next batch iteration
                startPoint.y = pointBatchY[POINTS_BATCH_SIZE - 1];
                
                mpz_add(stride_sum.get_mpz_t(), stride_sum.get_mpz_t(), batch_stride.get_mpz_t());
                    
                if (threadIdx == 0) {
                    save_counter += 1;
                    if (save_counter % 50000 == 0) {
                        char strideS[68];
                        gmp_snprintf(strideS, 67, "%0.64Zx", stride_sum.get_mpz_t());
                        cpub = secp256k1->GetPublicKeyHex(startPoint);
                        ofstream outFile;
                        outFile.open("settings1.txt");
                        outFile << cpub <<'\n';
                        outFile << strideS << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings1.txt" << endl;
                    }
                }
            } // while (true) loop end curly brace
        };

        std::thread addition_Threads[cpuCores];
        for (int i = 0; i < cpuCores; i++) {
            addition_Threads[i] = std::thread(scalable_addition_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < cpuCores; i++) {
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
        
        mpz_class offset_Step, vector_Num;
        mpz_fdiv_q_ui(offset_Step.get_mpz_t(), S_table[range_start - 2].get_mpz_t(),  cpuCores);
        
        vector_Num = 0;
        vector<mpz_class> offset_Nums;
        for (int i = 0; i < cpuCores; i++) {
            offset_Nums.push_back(vector_Num);
            mpz_add(vector_Num.get_mpz_t(), vector_Num.get_mpz_t(), offset_Step.get_mpz_t());
        }
        
        vector<Point> offset_Points;
        offset_Points.push_back(secp256k1->G);
        for (int i = 1; i < cpuCores; i++) {
            offset_Points.push_back(secp256k1->ScalarMultiplication(offset_Nums[i]));
        }
        
        vector<Point> starting_Points;
        Point vector_Point(start_point);
        starting_Points.push_back(vector_Point);
        
        for (int i = 1; i < cpuCores; i++) {
            vector_Point = secp256k1->SubtractPoints(start_point, offset_Points[i]);
            starting_Points.push_back(vector_Point);
        }

        Point addPoints[POINTS_BATCH_SIZE]; // array for batch addition points       
        Point batch_Add = secp256k1->DoublePoint(stride_point);
        addPoints[0] = stride_point;
        addPoints[0].y = secp256k1->ModNeg(addPoints[0].y);
        addPoints[1] = batch_Add;
        addPoints[1].y = secp256k1->ModNeg(addPoints[1].y);
        for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in batch addition points array with points
        {
            batch_Add = secp256k1->AddPoints(batch_Add, stride_point);
            addPoints[i] = batch_Add;
            addPoints[i].y = secp256k1->ModNeg(addPoints[i].y);
        }
        
        auto scalable_subtraction_search = [&](Point starting_Point, int threadIdx, mpz_class offset, mpz_class stride_Sum) {

            Point startPoint = starting_Point;
            Point BloomP; // point for insertion of the batch into the bloomfilter
            mpz_class stride_sum = stride_Sum;
            string cpub, xc, xc_sub;
            int index, count;
            vector<uint64_t> privkey_num;
            uint64_t steps;
            mpz_class Int_steps, Int_temp, privkey;
            
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            mpz_class deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            mpz_class pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            mpz_class pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            mpz_class deltaY; // values to store the results of points addition formula
            mpz_class slope[POINTS_BATCH_SIZE];

            mpz_class batch_stride, batch_index;
            mpz_mul_ui(batch_stride.get_mpz_t(), stride.get_mpz_t(), POINTS_BATCH_SIZE);
            
            while (true) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    mpz_sub(deltaX[i].get_mpz_t(), startPoint.x.get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(deltaX[i].get_mpz_t(), deltaX[i].get_mpz_t(), Fp.get_mpz_t()); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv(deltaX); // assign array deltaX to modGroup for batch inversion
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                    mpz_mul(slope[i].get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                    mpz_mod(slope[i].get_mpz_t(), slope[i].get_mpz_t(), Fp.get_mpz_t());

                    mpz_mul(pointBatchX[i].get_mpz_t(), slope[i].get_mpz_t(), slope[i].get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), startPoint.x.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());

                }
                
                mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                mpz_mul(slope[i].get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                mpz_mod(slope[i].get_mpz_t(), slope[i].get_mpz_t(), Fp.get_mpz_t());

                mpz_mul(pointBatchX[i].get_mpz_t(), slope[i].get_mpz_t(), slope[i].get_mpz_t());
                mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), startPoint.x.get_mpz_t());
                mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());
        
                mpz_sub(pointBatchY[i].get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                mpz_mul(pointBatchY[i].get_mpz_t(), slope[i].get_mpz_t(), pointBatchY[i].get_mpz_t());
                mpz_sub(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), startPoint.y.get_mpz_t());
                mpz_mod(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), Fp.get_mpz_t());
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {

                    xc = secp256k1->GetXHex(pointBatchX[i]);
                     
                    if (bf1.may_contain(xc)) {

                        print_time(); cout << "BloomFilter Hit " << bloomfile1 << " (Even Point) [Higher Range Half]" << endl;

                        BloomP.x = pointBatchX[i];
                        mpz_sub(BloomP.y.get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                        mpz_mul(BloomP.y.get_mpz_t(), slope[i].get_mpz_t(), BloomP.y.get_mpz_t());
                        mpz_sub(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), startPoint.y.get_mpz_t());
                        mpz_mod(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), Fp.get_mpz_t());

                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            xc_sub = secp256k1->GetXHex(BloomP.x);
                            while (bf1.may_contain(xc_sub)) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                xc_sub = secp256k1->GetXHex(BloomP.x);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }

                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps = steps;
                        mpz_mul_ui(batch_index.get_mpz_t(), stride.get_mpz_t(), i + 1);
                        mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), batch_index.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), offset.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                        mpz_add(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                        mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                        calc_point = secp256k1->ScalarMultiplication(privkey);

                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
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
                    
                    if (bf2.may_contain(xc)) {

                        print_time(); cout << "BloomFilter Hit " << bloomfile2 << " (Odd Point) [Higher Range Half]" << endl;

                        BloomP.x = pointBatchX[i];
                        mpz_sub(BloomP.y.get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                        mpz_mul(BloomP.y.get_mpz_t(), slope[i].get_mpz_t(), BloomP.y.get_mpz_t());
                        mpz_sub(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), startPoint.y.get_mpz_t());
                        mpz_mod(BloomP.y.get_mpz_t(), BloomP.y.get_mpz_t(), Fp.get_mpz_t());

                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            xc_sub = secp256k1->GetXHex(BloomP.x);
                            while (bf2.may_contain(xc_sub)) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                xc_sub = secp256k1->GetXHex(BloomP.x);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }

                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps = steps;
                        mpz_mul_ui(batch_index.get_mpz_t(), stride.get_mpz_t(), i + 1);
                        mpz_add(Int_temp.get_mpz_t(), stride_sum.get_mpz_t(), batch_index.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), offset.get_mpz_t());
                        mpz_add(Int_temp.get_mpz_t(), Int_temp.get_mpz_t(), Int_steps.get_mpz_t());
                        mpz_add(privkey.get_mpz_t(), pre_calc_sum.get_mpz_t(), Int_temp.get_mpz_t());
                        mpz_mul_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 2);
                        mpz_add_ui(privkey.get_mpz_t(), privkey.get_mpz_t(), 1);
                        calc_point = secp256k1->ScalarMultiplication(privkey);

                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) { // if cpubs are equal we got it
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
                }
                
                startPoint.x = pointBatchX[POINTS_BATCH_SIZE - 1]; // setting the new startPoint for the next batch iteration
                startPoint.y = pointBatchY[POINTS_BATCH_SIZE - 1];
                
                mpz_add(stride_sum.get_mpz_t(), stride_sum.get_mpz_t(), batch_stride.get_mpz_t());
                    
                if (threadIdx == 0) {
                    save_counter += 1;
                    if (save_counter % 50000 == 0) {
                        char strideS[68];
                        gmp_snprintf(strideS, 67, "%0.64Zx", stride_sum.get_mpz_t());
                        cpub = secp256k1->GetPublicKeyHex(startPoint);
                        ofstream outFile;
                        outFile.open("settings2.txt");
                        outFile << cpub <<'\n';
                        outFile << strideS << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Save Data written to settings2.txt" << endl;
                    }
                }
            } // while (true) loop end curly brace
        };

        std::thread subtraction_Threads[cpuCores];
        for (int i = 0; i < cpuCores; i++) {
            subtraction_Threads[i] = std::thread(scalable_subtraction_search, starting_Points[i], i , offset_Nums[i], stride_sum);
        }

        for (int i = 0; i < cpuCores; i++) {
            subtraction_Threads[i].join();
        }
    };

    print_time(); cout << "Search in progress..." << endl;
    
    std::thread thread1(addition_search);
    std::thread thread2(subtraction_search);
    
    thread1.join();
    thread2.join();
}
