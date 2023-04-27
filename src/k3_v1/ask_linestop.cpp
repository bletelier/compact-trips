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

    int n = 1000000;
    srand(1);

     std::cout << "\nQuantity Linestop Query: \n";
    n = 1000;
    srand(1);
    std::vector<int> queries3;
    for(int i = 0; i < n; ++i) {
        int s = rand() % matrix_size;
        queries3.push_back(s);
    }
    unsigned long found3 = 0;
    auto start = std::chrono::steady_clock::now();
    for(int i = 0; i < n; ++i) {
        int s = queries3[i];
        found3 += k3_od.get_people_quantity_on_linestop(s); 
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = end-start;
    std::cout << "elapsed: " << elapsed.count() << " seconds, found: " << found3 << '\n';
    std::cout << "per_query: " << (elapsed.count()*1000000.0f)/(n*1.0f) << " microS\n";
    
    return 0;
}
