#include <iostream>
#include <stdio.h>
#include "k3_tree_od_v2.hpp"

int main(int argc, char** argv) {
    if(argc < 2) return 1;
    std::string filename = argv[1];
    int n_trips, n_nodes; std::cin >> n_trips >> n_nodes;
    
    std::vector<std::tuple<long unsigned int, long unsigned int, long unsigned int>> m;

    std::vector<std::vector<int>> trips;
    
    for(int i = 0; i < n_trips; ++i) {
        int trip_len; std::cin >> trip_len;
        std::vector<int> trip(trip_len, 0);
        for(int j = 0; j < trip_len; ++j) {
            int stop; std::cin >> stop;
            trip[j] = stop - 1;
        }
        int o = trip[0];
        int d = trip[trip_len - 1];
        
        if(trip_len == 3) m.push_back({o, d, trip[1]});

        trips.push_back(trip);
    }
    
    sdsl::k3_tree_od_v2<2> k3_od(trips, n_nodes);
    std::string out = "../cds/"+filename+"/k3_v2.sdsl";
    sdsl::store_to_file(k3_od, out);

    std::cout << "Size of the structure: " << sdsl::size_in_bytes(k3_od) << " bytes.\n";
    std::cout << "stored\n";
        
    return 0;
}
