#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include "od_rep_v1.hpp"

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    std::string input = "../cds/"+filename+"/nok3_v1.sdsl";
    std::ifstream input_file(input);
    sdsl::od_rep_v1<uint32_t> k3_od;
    k3_od.load(input_file);    
   
    int matrix_size = k3_od.matrix_n();

    std::cout << "Size of the structure: " << sdsl::size_in_bytes(k3_od) << " bytes.\n";
    
    std::cout << "\nQuantity Linestop Query: \n";
    int n = 1000;
    srand(1);
    std::vector<int> queries3;
    for(int i = 0; i < n; ++i) {
        int s = rand() % matrix_size;
        queries3.push_back(s);
    }
    unsigned long found3 = 0;
    unsigned int found_ = 0;
    std::chrono::duration<double> elapsed(0.0);

    for(int i = 0; i < n; ++i) {
        int s = queries3[i];

        for(int k = 0; k < 5; ++k) {
            std::map<std::pair<uint16_t, uint16_t>, uint32_t> res;
            auto start = std::chrono::steady_clock::now();
            found_ = k3_od.get_origin_destinations_linestop(s, res); 
            auto end = std::chrono::steady_clock::now(); 
            elapsed += (end-start)/5;
        }
        found3 += found_;
    }

    std::cout << "elapsed: " << elapsed.count() << " seconds, found: " << found3 << '\n';
    std::cout << "per_query: " << (elapsed.count()*1000000.0f)/(n*1.0f) << " microS\n";
    
    return 0;
}
