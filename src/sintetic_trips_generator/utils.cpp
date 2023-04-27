#include "utils.hpp"
#include <map>

void read_input(int n_lines, std::vector<std::vector<int>> &lines, int &n_stops) {
    std::set<int> used;
    std::map<int,int> idd;
    int id = 1;
    int siz = 0;
    int max_stop = -1;
    for (int i = 0; i < n_lines; ++i) {
        int size_line; std::cin >> size_line;
        siz += size_line;
        for (int j = 0; j < size_line; ++j) {
            int stop; std::cin >> stop;
            if(used.find(stop) == used.end()) idd[stop] = id++; 
            used.insert(stop);
            max_stop = std::max(max_stop, idd[stop]);
            lines[i].push_back(idd[stop]);
        }
    } 
    n_stops = used.size();
}

void create_map(std::vector<std::vector<int>> &lines, std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, int &ns) {
    int c = 0;
    
    for(std::vector<int> line : lines) {
        for(int i = 0; i < line.size(); ++i) {
            bits[c].set(line[i]);
            //std::cout << line.size() << '\n';
            for(int j = 0; j < line.size(); ++j) {
                if(i == j) continue;
                //std::cout << "a\n";
                //std::cout << i << ' ' <<  j << ' '<< map.size() << ' ' << line[i] << ' ' << c << ' ' << line[j] << '\n';
                map[line[i]].push_back({line[j], c});
            }
        }
        //std::cout << "L" << c+1 << ":\t" << bits[c] << '\n';
        c++;
    }
    int n_stops = map.size();
    for(int i = 1; i < map.size(); ++i) {
        
        int choiced_lines = 0;
        std::set<int> used;
        for(std::pair<int,int> line : map[i]) {
            if(used.count(line.second) <= 0)
            {
                used.insert(line.second);
                if(choiced_lines <= 0) map_line[i][line.second] = i;
                else map_line[i][line.second] = n_stops++;
                choiced_lines++;
            }
        }
        if (choiced_lines > 1) transfer.set(i);
        //std::cout <<i << ' ' <<  n_stops << '\n';
    }
    ns = n_stops - 1;
    //std::cout <<"Tr:\t" << transfer << "\n\n";
}

void show_stops(std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::vector<int>> &map_line) {
    for(int i = 1; i < map.size(); ++i) {
        std::cout << i << ": ";
        for(int j = 0; j < map[i].size(); ++j) {
            std::cout << map[i][j].first << "(L" << map[i][j].second+1 << ") ";    
        }
        std::cout << '\n';
    }
    std::cout << "\n\t";
    for(int i = 0, j = 0; i < map_line.size(); ++i) {
        if (i > 0) std::cout << i << '\t';
        for(int ml : map_line[i]) {
            if (i <= 0) std::cout << 'L' << ++j << '\t'; 
            else std::cout << ml << '\t';
        }
        std::cout << '\n';
    }
    //std::cout << '\n';
}

void simulate_trips(int seed, int n_trips, std::vector<std::vector<std::pair<int, int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, int ns, std::string filename) {
    srand(seed);
    std::string sfile = "../trips/"+filename+"_"+std::to_string(n_trips)+"_"+std::to_string(ns)+".txt";
    const char *file = sfile.c_str();
    FILE *fp2 = fopen(file, "w");
    std::cout << sfile << '\n';
    fprintf(fp2, "%d %d\n", n_trips, ns);
    fclose(fp2);
    std::vector<unsigned long> totals (4, 0);
    //std::cout << n_trips << '\n';
    long eqa = 10000;
    for(long i = 0; i < n_trips; ++i) {
        /*int stop = (rand() % (map.size()-1)) + 1;
        while(stop == start) {
            stop = (rand() % (map.size()-1)) + 1;
        }*/
        //std::cout << i+1 << "\nstart: " << start << " stop: " << stop << '\n';
        bool try_again = true;
        //std::cout << i+1 << '\n';
        //else if(rlen <= 100) len = 5;
        int len;
        while(try_again) {
            int rlen = rand() % 100 + 1;
            
            if(rlen <= 50) len = 2;
            else if(rlen <= 85) len = 3;
            else if(rlen <= 98) len = 4;
            else if(rlen <= 100) len = 5;
            int start = (rand() % (map.size()-1)) + 1;
            try_again = simulate_t(start, len, map, bits, transfer, map_line, file);
        }
        totals[len-2]++;
        if ((i+1) == eqa) {std::cout << (eqa/10000)*10 << '\n'; eqa += 10000;}
    }
    unsigned int total = 0;
    for(int i = 0; i < 4; ++i) {
        std::cout << i+2 << ": " << 100.0f*(totals[i]/(n_trips*1.0f)) << "% | " << totals[i] << "\n";
        total+=totals[i]*(i+2);
    }
    std::cout << "avg len: " << total/(n_trips*1.0f) << '\n';
}

bool simulate_t(int start, int len, std::vector<std::vector<std::pair<int,int>>> &map, std::vector<std::bitset<BITS>> &bits, std::bitset<BITS> &transfer, std::vector<std::vector<int>> &map_line, const char* file) {
    int choosed_line;
    std::set<int> used_lines;
    int choosed_stop = start;
    std::set<int> used_stops;
    used_stops.insert(start);
    int last_used = -1;
    std::vector<int> trip;
    std::bitset<BITS> my_transfer = transfer;
    my_transfer[start] = 0;
    trip.push_back(start);
    for(int l = 1; l < len; ++l) {
        
        std::set<int> possible_lines;
        std::vector<std::pair<int,int>> poss = map[choosed_stop];
        //std::cout << "pos: ";
        for(int i = 0; i < poss.size(); ++i) { 
            if(used_lines.find(poss[i].second) == used_lines.end()) {
                possible_lines.insert(poss[i].second);
                //std::cout << poss[i].second << ' ';
            }
        }
        if(possible_lines.size() <= 0) return true;
        //std::cout << '\n';

        int rndl = rand() % possible_lines.size();
        std::set<int>::iterator it = possible_lines.begin();
        for(int i = 0; i < rndl; ++i) {
            it++;
        }
        //std::cout << "the iterator\n";
        //std::cout << choosed_stop <<' ';
        
        //elegir la linea que usare (si no soy transfer solo tengo una opcion)
        //si una linea ya fue usada, no usarla de nuevo.
        if(l == len - 1) {
            //me voy a cualquier lado dentro de la ultima linea que estaba (diferente a inicio)
            
            int rnds = rand() % bits[*it].count() + 1;
            for(int i = 0, j = 0; i < BITS; ++i) {
                if(rnds == j) {
                        if(used_stops.find(i-1) == used_stops.end()) {
                            trip.push_back(*it+1);
                            trip.push_back(i-1);
                            break;
                        }
                        else {
                            rnds = rand() % bits[*it].count() + 1;
                            i = 0;
                            j = 0;
                        }
                }
                if(bits[*it].test(i)) j++;
            }
        }
        else {
            //me voy a cualquier transfer que este en mi linea y me conecte con otra linea que no haya usado
            std::bitset<BITS> and_stop = bits[*it] & my_transfer;
            //std::cout << "aqui\n";
            while(and_stop.count() < 1) {
                possible_lines.erase(it);
                if(possible_lines.size() <= 0) return true;
                rndl = rand() % possible_lines.size();
                it = possible_lines.begin();
                for(int i = 0; i < rndl; ++i) {
                    it++;
                }
                and_stop = bits[*it] & my_transfer;
            }
            //std::cout << "hola" << '\n';
            used_lines.insert(*it);
            //bits[*it][choosed_stop] = 0;   
            /*for(int b = 0; b < 12; ++b) {
                std::cout << and_stop[12-b];
            }
            std::cout << '\n';*/

            int rnds = rand() % and_stop.count() + 1;
            for(int i = 0, j = 0; i < BITS; ++i) {
                if(rnds == j) {
                    trip.push_back(*it+1);
                    trip.push_back(i-1);
                    choosed_stop = i-1;
                    my_transfer[i-1] = 0;
                    used_stops.insert(i-1);
                    //bits[*it][choosed_stop] = 0;   
                    break;
                }
                if(and_stop.test(i)) j++;
            }      
            //std::cout << "choosed: " << rnds << '\n'; 
        }
    }
    
    FILE *fp2 = fopen(file, "a+");
    fprintf(fp2, "%ld\n%d", trip.size() - (trip.size()/2), map_line[trip[0]][trip[1]-1]);
    //fprintf(fp2, "%d %d\n", n_trips, ns);
    for(int i = 2; i < trip.size()-1; i+=2) {
        fprintf(fp2, " %d", map_line[trip[i]][trip[i-1]-1]);
    }

    fprintf(fp2," %d\n",map_line[trip[trip.size()-1]][trip[trip.size()-2]-1]); 
    fclose(fp2);    
    return false;
}
