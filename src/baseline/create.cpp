#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <stdlib.h>
#include "od_baseline.hpp"

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    int n_trips, n_stops;
    std::cin >> n_trips >> n_stops;
    std::vector<std::vector<uint16_t>> trips;
    for(int i = 0; i < n_trips; ++i) {
        int trip_len; std::cin >> trip_len;
        std::vector<uint16_t> trip(trip_len);
        for(int j = 0; j < trip_len; ++j) {
            int stop; std::cin >> stop;
            if(j >= 1 and j < trip_len - 1)
                trip[j+1] = stop - 1;
            else if (j >= trip_len - 1) trip[1] = stop -1;
            else trip[j] = stop - 1;
        }
        
        trips.push_back(trip);
    }

    sdsl::od_baseline<uint32_t> od_rep (trips, n_stops); 
    std::string out = "../cds/"+filename+"/baseline.sdsl";
    sdsl::store_to_file(od_rep,out);
    std::cout << "Size of the structure: " << sdsl::size_in_bytes(od_rep) << " bytes.\n";
    std::cout << "stored\n";
    return 0;
} 
