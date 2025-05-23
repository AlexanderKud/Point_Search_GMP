#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>
#include <vector>
#include <thread>
#include <gmpxx.h>
#include <gmp.h>

#include "secp256k1/secp256k1.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
namespace fs = filesystem;

static constexpr int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion(one ModInv for the entire group) using IntGroup class
const mpz_class Fp = mpz_class("115792089237316195423570985008687907853269984665640564039457584007908834671663", 10);

auto main() -> int {
    
    auto start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256k1 *secp256k1 = new Secp256k1(); secp256k1->Init(); // initializing secp256k1 context
    
    fs::path current_path = fs::current_path(); // deleting previous settings and bloom files
    auto file_list = get_files_in_directory(current_path);
    vector<string> targets = {"settings1.txt", "settings2.txt", "bloom1.bf", "bloom2.bf"};
    for (auto i : file_list) {
        for (auto t : targets) {
            if (i == t) { std::remove(t.c_str()); }
        }
    }

    mpz_class pk; pk = 1; // generating power of two points table (2^0..2^256) 
    vector<Point> P_table;
    Point P;
    for (int i = 0; i < 256; i++)
    {
        P = secp256k1->ScalarMultiplication(pk);
        P_table.push_back(P);
        mpz_mul_ui(pk.get_mpz_t(), pk.get_mpz_t(), 2);
    }
    print_time(); cout << "P_table generated" << endl;

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter 
    string temp, search_pub;                      // and a stride size to walk the range
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

    Point puzzle_point, puzzle_point_05, puzzle_point_divide2; // calculation points
    Point first_point, second_point, P1, P2, Q1, Q2;
    mpz_class stride_sum, div2;
    
    stride_sum = 0;
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
    uint64_t n_elements = uint64_t(pow(2, block_width));
    double error = 0.0000000001;
    int n_cores = 4;  //actual number of processing cores but equal to some power of two value(2,4,8,16,32,64,...) divided by 2
    uint64_t count = uint64_t(pow(2, block_width) / n_cores); // actual cores = 8  8 / 2 = 4 cores for each lambda function
    mpz_class add_key;                                        // should be some power of two to evenly divide the space between threads
    add_key = count;                                          // here we have namely 11 threads( 1-main thread 2,3 - lambda functions
    Point Add_Point = secp256k1->ScalarMultiplication(add_key);// and 4 threads inside each lambda for process_chunk
                                                               // execution timings are relative to my pc only yours might be quicker
    Point addPoints[POINTS_BATCH_SIZE]; // array for the batch addition points(1G .. 1024G)
    Point batch_Add = secp256k1->DoublePoint(secp256k1->G); // 2G
    addPoints[0] = secp256k1->G; // 1G
    addPoints[1] = batch_Add;    // 2G
    for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in the batch addition points array with points from(3G .. 1024G)
    {
        batch_Add = secp256k1->AddPoints(batch_Add, secp256k1->G);
        addPoints[i] = batch_Add;
    }
    
    int nbBatch = count / POINTS_BATCH_SIZE; // number of batches for the single thread

    auto bloom_create1 = [&]() {
        string bloomfile = "bloom1.bf";
        Point P = puzzle_point;
        vector<Point> starting_points;
        for (int i = 0; i < n_cores; i++) {
            starting_points.push_back(P);
            P = secp256k1->AddPoints(P, Add_Point);
        }

        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
            mpz_class deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            mpz_class pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            mpz_class pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
                      
            Point startPoint = start_point; // start point
            bf.insert(secp256k1->GetPublicKeyHex(startPoint)); // we insert it instantly into the bloomfilter

            Point BloomP; // point for insertion of the batch into the bloomfilter
            mpz_class deltaY, slope, slopeSquared; // values to store the results of points addition formula
            
            for (int i = 0; i < nbBatch; i++) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    mpz_sub(deltaX[i].get_mpz_t(), startPoint.x.get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(deltaX[i].get_mpz_t(), deltaX[i].get_mpz_t(), Fp.get_mpz_t()); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv(deltaX); // assign array deltaX to modGroup for batch inversion
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // follow points addition formula logic
                    
                    mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                    mpz_mul(slope.get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                    mpz_mod(slope.get_mpz_t(), slope.get_mpz_t(), Fp.get_mpz_t());

                    mpz_mul(slopeSquared.get_mpz_t(), slope.get_mpz_t(), slope.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), slopeSquared.get_mpz_t(), startPoint.x.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());
        
                    mpz_sub(pointBatchY[i].get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                    mpz_mul(pointBatchY[i].get_mpz_t(), slope.get_mpz_t(), pointBatchY[i].get_mpz_t());
                    mpz_sub(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), startPoint.y.get_mpz_t());
                    mpz_mod(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), Fp.get_mpz_t());

                }

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points into the bloomfilter
                    BloomP.x = pointBatchX[i];
                    BloomP.y = pointBatchY[i];
                    bf.insert(secp256k1->GetPublicKeyHex(BloomP));
                }
                
                startPoint.x = pointBatchX[POINTS_BATCH_SIZE - 1]; // setting the new startPoint for the next batch iteration
                startPoint.y = pointBatchY[POINTS_BATCH_SIZE - 1];
            }
        };
        
        std::thread myThreads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }
    
        print_time(); cout << "Creating bloom1 image with " << n_cores << " threads" << '\n';

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join();
        }

        print_time(); cout << "Writing bloom1 image to bloom1.bf" << '\n';
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
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
            mpz_class deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            mpz_class pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
            mpz_class pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
                      
            Point startPoint = start_point; // start point
            bf.insert(secp256k1->GetPublicKeyHex(startPoint)); // we insert it instantly into the bloomfilter

            Point BloomP; // point for insertion of the batch into the bloomfilter
            mpz_class deltaY, slope, slopeSquared; // values to store the results of points addition formula
            
            for (int i = 0; i < nbBatch; i++) {
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    mpz_sub(deltaX[i].get_mpz_t(), startPoint.x.get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(deltaX[i].get_mpz_t(), deltaX[i].get_mpz_t(), Fp.get_mpz_t()); // insert each entry into the deltaX array
                }
    
                modGroup.ModInv(deltaX); // assign array deltaX to modGroup for batch inversion
                
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // follow points addition formula logic
                    
                    mpz_sub(deltaY.get_mpz_t(), startPoint.y.get_mpz_t(), addPoints[i].y.get_mpz_t());
                    mpz_mul(slope.get_mpz_t(), deltaY.get_mpz_t(), deltaX[i].get_mpz_t());
                    mpz_mod(slope.get_mpz_t(), slope.get_mpz_t(), Fp.get_mpz_t());

                    mpz_mul(slopeSquared.get_mpz_t(), slope.get_mpz_t(), slope.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), slopeSquared.get_mpz_t(), startPoint.x.get_mpz_t());
                    mpz_sub(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), addPoints[i].x.get_mpz_t());
                    mpz_mod(pointBatchX[i].get_mpz_t(), pointBatchX[i].get_mpz_t(), Fp.get_mpz_t());
        
                    mpz_sub(pointBatchY[i].get_mpz_t(), startPoint.x.get_mpz_t(), pointBatchX[i].get_mpz_t());
                    mpz_mul(pointBatchY[i].get_mpz_t(), slope.get_mpz_t(), pointBatchY[i].get_mpz_t());
                    mpz_sub(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), startPoint.y.get_mpz_t());
                    mpz_mod(pointBatchY[i].get_mpz_t(), pointBatchY[i].get_mpz_t(), Fp.get_mpz_t());

                }

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points into the bloomfilter
                    BloomP.x = pointBatchX[i];
                    BloomP.y = pointBatchY[i];
                    bf.insert(secp256k1->GetPublicKeyHex(BloomP));
                }
                
                startPoint.x = pointBatchX[POINTS_BATCH_SIZE - 1]; // setting the new startPoint for the next batch iteration
                startPoint.y = pointBatchY[POINTS_BATCH_SIZE - 1];
            }
        };
        
        std::thread myThreads[n_cores];
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }
    
        print_time(); cout << "Creating bloom2 image with " << n_cores << " threads" << '\n';

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join();
        }

        print_time(); cout << "Writing bloom2 image to bloom2.bf" << '\n'; 
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
