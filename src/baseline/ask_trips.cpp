#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "od_baseline.hpp"
#include <chrono>

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    std::string input = "../cds/"+filename+"/baseline.sdsl";
    std::ifstream input_file(input);
    sdsl::od_baseline<uint32_t> k3_od;
    k3_od.load(input_file);    
   
    int matrix_size = k3_od.n_stops();

    std::cout << "Size of the structure: " << sdsl::size_in_bytes(k3_od) << " bytes.\n";
    std::cout << "\nTrip Query: \n";
    int n = 1000000;
    srand(1);
    std::vector<std::pair<int,int>> queries;
    for(int i = 0; i < n; ++i) {
        int o = rand() % matrix_size;
        int d = rand() % matrix_size;
        while(d == o) {
            d = rand() % matrix_size;
        }
        queries.push_back({o,d});

    }
  
    int found = 0;
    std::vector<std::vector<uint16_t>> res; 
    std::chrono::duration<double> elapsed(0.0);

    for(int i = 0; i < n; ++i) {
        int o = queries[i].first;
        int d = queries[i].second;
        for(int k = 0; k < 5; ++k) {
            auto start = std::chrono::steady_clock::now();
            res = k3_od.get_trips(o,d);
            auto end = std::chrono::steady_clock::now();
            elapsed += (end - start)/5;
        }
        found += res.size();
    }
    std::cout << "elapsed: " << elapsed.count() << " seconds found: " << found << '\n';
    std::cout << "per_query: " << (elapsed.count()*1000000.0f)/(n*1.0f) << " microS\n";
    std::cout << "time_per_trip: " << (elapsed.count()*1000000.f)/(found*1.0f) << '\n';
    return 0;
}
