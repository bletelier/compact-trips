#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "k3_tree_od_v2.hpp"
#include <chrono>

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    std::string input = "../cds/"+filename+"/k3_v2.sdsl";
    std::ifstream input_file(input);
    sdsl::k3_tree_od_v2<2> k3_od;
    k3_od.load(input_file);    
   
    int matrix_size = k3_od.get_matrix_size();

    std::cout << "Size of the structure: " << sdsl::size_in_bytes(k3_od) << " bytes.\n";

    std::cout << "\nQuantity Query: \n";
    int n = 1000000;
    srand(1);
    std::vector<std::pair<int,int>> queries2;
    for(int i = 0; i < n; ++i) {
        int o = rand() % matrix_size;
        int d = rand() % matrix_size;
        while(d == o) {
            d = rand() % matrix_size;
        }
        queries2.push_back({o,d});
    }
    unsigned long found2 = 0;
    unsigned int res = 0;
    std::chrono::duration<double> elapsed(0.0); 
    for(int i = 0; i < n; ++i) {
        int o = queries2[i].first;
        int d = queries2[i].second;
        for(int t = 0; t < 1; ++t) {
            auto start = std::chrono::steady_clock::now();
            res = k3_od.get_people_quantity(o,d);
            auto end = std::chrono::steady_clock::now();
            elapsed += (end-start)/1;
        }
        found2 += res;
    }

    std::cout << "elapsed: " << elapsed.count() << " seconds, found: " << found2 << '\n';
    std::cout << "per_query: " << (elapsed.count()*1000000.0f)/(n*1.0f) << " microS\n";

    return 0;
}
