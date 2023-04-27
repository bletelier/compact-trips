#include "utils.hpp"

int main(int argc, char **argv) {
    if(argc != 2) {
        std::cout << "ERROR!, output file name needed\n";
        return 1;
    }
    std::string filename (argv[1]);
    int n_stops, n_lines; std::cin >> n_lines;
    std::vector<std::vector<int>> lines(n_lines);                                                              
    //Read and store the lines input (n_lines lines and each line has size_line stops).
    read_input(n_lines, lines, n_stops);
    //Create the city-map graph using the stops as nodes and the lines as edges. Note that a line has a fully connected subgraph between the stops of that line,
    //since you can do any trip without changing to another line. 
    std::vector<std::vector<std::pair<int,int>>> map(n_stops+1);
    std::vector<std::vector<int>> map_line(n_stops+1, std::vector<int> (n_lines,0));
    std::bitset<BITS> transfer(0x0);
    std::vector<std::bitset<BITS>> bits(n_lines, std::bitset<BITS>(0x0));
    int ns = 0;
    create_map(lines, map, bits, transfer, map_line, ns);

    //show_stops(map, map_line); 
    
    int seed, n_trips; std::cin >> seed >> n_trips;
    simulate_trips(seed, n_trips, map, bits, transfer, map_line, ns, filename);
    return 0;
}
