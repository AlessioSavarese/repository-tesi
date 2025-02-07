#ifndef TEST_H
#define TEST_H

#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdint>

// Include le classi degli algoritmi
#include "KMinHash.h"
#include "FastSimilaritySketching.h"
#include "LettoreFile.h"  // Per la lettura dei file
#include "JS.h"           // Per la gestione della Jaccard Similarity

#include <chrono>
#include <thread>

class Test {
private:
    // Metodi privati di utilità
    static std::vector<uint64_t> generate_random_set(int n, int m);
    // static void save_result(const std::string& filename, const std::string& algoritmo, int param, int rep, double time, const std::string& strParam);


public:
    // Metodi pubblici
    static void test_time_vs_n(int k_fixed, std::vector<int> n_values, int repetitions, int m);
    static void test_time_vs_k(std::vector<int> k_values, int n_fixed, int repetitions, int m);
    // Metodo statico per il test di qualità
    static void test_quality(int k, int n, int repetitions, int m);
    static void test_separation(size_t t, int m, double gamma, int num_coppie, const std::vector<size_t>& r_values);
};

#endif // TEST_H