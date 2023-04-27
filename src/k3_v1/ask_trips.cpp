#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "k3_tree_od_v1.hpp"
#include <chrono>

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    std::string input = "../cds/"+filename+"/k3_v1.sdsl";
    std::ifstream input_file(input);
    sdsl::k3_tree_od_v1<2> k3_od;
    k3_od.load(input_file);       
   
    int matrix_size = k3_od.get_matrix_size();

    std::cout << "Size of the structure: " << sdsl::size_in_bytes(k3_od) << " bytes.\n";
    std::cout << "\nTrips Query: \n";
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
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < n; ++i) {
        int o = queries[i].first;
        int d = queries[i].second;

        std::vector<std::vector<int>> res = k3_od.get_trips(o,d);
        found += res.size();

    }
   
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end-start;
    std::cout << "elapsed: " << elapsed.count() << " seconds found: " << found << '\n';
    std::cout << "per_query: " << (elapsed.count()*1000000.0f)/(n*1.0f) << " microS\n";
    std::cout << "time_per_trip: " << (elapsed.count()*1000000.f)/(found*1.0f) << '\n';

    return 0;
}
