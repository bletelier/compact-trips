#include <iostream>
#include <stdio.h>
#include "od_rep_v2.hpp"

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    int n_trips, n_nodes; std::cin >> n_trips >> n_nodes;
     
    std::vector<std::vector<uint32_t>> trips;

    for(int i = 0; i < n_trips; ++i) {
        int trip_len; std::cin >> trip_len;
        std::vector<uint32_t> trip(trip_len, 0);
        for(int j = 0; j < trip_len; ++j) {
            int stop; std::cin >> stop;
            trip[j] = stop - 1;
        }
        trips.push_back(trip);
    }
    sdsl::od_rep_v2<uint32_t> od_rep(trips, n_nodes);

    std::string out= "../cds/"+filename+"/nok3_v2.sdsl";
    sdsl::store_to_file(od_rep, out );
    std::cout << "Size of the structure: " << sdsl::size_in_bytes(od_rep) << " bytes.\n"; 
    std::cout << "stored\n";
   
    return 0;
}
