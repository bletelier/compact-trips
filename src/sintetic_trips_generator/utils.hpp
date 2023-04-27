#include <stdio.h>
#include <vector>
#include <utility>
#include <set>
#include <bitset>
#include <iostream>
#include <stdlib.h>

const int TRESHOLD = 5;
const int BITS = 1200;

void read_input(int n_lines, std::vector<std::vector<int>> &lines, int &n_stops);
void create_map(std::vector<std::vector<int>> &lines, std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, int &ns);
void show_stops(std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::vector<int>> &map_line);
void simulate_trips(int seed, int n_trips, std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, int ns, std::string filename);

bool simulate_t(int start, int len, std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, const char* file);
