/*  
 * Created by Benjamin Letelier on 27/04/23.
 *
 * Copyright (C) 2023-current-year, Benjamin Letelier, all rights reserved.
 *
 * 
 * Author's contact: Benjamin Letelier  <bletelier@gmail.com>
 * Facultad de Ciencias de la Ingeniería. Universidad Austral de Chile. Chile.
 * Databases Lab, University of A Coruña. Campus de Elviña s/n. Spain
 *
 * DESCRIPTION
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef LBD_OD_BASELINE_HPP
#define LBD_OD_BASELINE_HPP

#include <sdsl/int_vector.hpp>

//This is an example, create your own namespace.
namespace sdsl {

    template <class t_value>
    class od_baseline {

    public:

        //Definiton of types, useful when you need to know a type inside a class
        typedef uint64_t size_type;
        typedef uint16_t id_type;
        typedef t_value value_type;
        //Typename for obtaining the type inside a class
        typedef int_vector<32> t_iv;
        typedef bit_vector t_bv;
        typedef std::vector<std::vector<uint16_t>> t_vv;

    private:
        t_vv od_trips; //vector<vector<u16>> storing entire trips ordered from o and then d.
        t_iv od_pq;  //people quantity for a trip
        size_type od_n_stops; // == od_pq.size();
        

        void copy(const od_baseline& o)
        {
            od_trips = o.od_trips;
            od_pq = o.od_pq;
            od_n_stops = o.od_n_stops;
        }

        t_value lower_search(int l, int r, id_type val, t_value pos) {
            t_value l_bound = od_trips.size();
            while(l <= r) {
                int m = ((r-l)/2) + l;
                if(od_trips[m][pos] > val) {
                    r = m-1;
                }
                else if(od_trips[m][pos] < val) {
                    l = m+1;
                }
                else {
                    r = m-1;
                    l_bound = m;
                }
            }
            return l_bound;
        }
        
        t_value upper_search(int l, int r, id_type val, t_value pos) {
            t_value r_bound = od_trips.size();
            while(l <= r) {
                int m = ((r-l)/2) + l;
                if(od_trips[m][pos] > val) {
                    r = m-1;
                }
                else if(od_trips[m][pos] < val) {
                    l = m+1;
                }
                else {
                    l = m+1;
                    r_bound = m;
                }
            }
            return r_bound;
        }

    public:

        // If you need access to parts of your class and you don't want to modify their values
        
        const t_iv &pq = od_pq;  //people quantity

		//*******************************************************//
        //******************* CONSTRUCTORS **********************//
        //*******************************************************//

        //Default constructor
        od_baseline() = default;

        //Copy constructor
        od_baseline(const od_baseline& o)
        {
            copy(o);
        }

        //Move constructor
        od_baseline(od_baseline&& o)
        {
            *this = std::move(o);
        }

        /**
         *  Constructor
         *  Container can be any class.
         *  Depending on the implementation it requires some methods.
         *  In this case we need that the container supports iterators, size, and each element a first and second attribute.
         */
        
        od_baseline(std::vector<std::vector<id_type>> &trips, size_type n){
            size_type n_trips = trips.size();
            if(n_trips < 1) {
                throw std::logic_error("There are not enough trips.");
            }
            if(n <= 1) {
                throw std::logic_error("There are not enough linestops.");
            }

            //people quantity that made certain trip.
            std::map<std::vector<id_type>, t_value> trip_ways;
            
            //for a certain (origin, destination) we will store all the trips that exists for that (o,d)
            std::map<std::pair<t_value,t_value>, std::vector<std::pair<id_type, std::vector<id_type>>>> od_ways; 
            std::cout << "storing individual trips\n";
            t_value n_dif_trips = 0;
            for(long i = 0; i < n_trips; ++i) {
                std::vector<id_type> trip = trips[i];
                id_type trip_len = trip.size();
                id_type o = trip[0];
                id_type d = trip[1];
                if(trip_ways.count(trip) <= 0) {
                    trip_ways[trip] = 1;
                    od_ways[{o,d}].push_back({trip_len, trip});
                    n_dif_trips++;
                }
                else trip_ways[trip]++;
            }
            
            std::cout << "done\n";

            std::cout << "creating aggregated and ordered trips and pq vector\n";
            t_vv od_trips_;
            t_iv od_pq_(n_dif_trips, 0);
            t_value cc = 0;
            for(int o = 0; o < n; ++o) {
                for(int d = 0; d < n; ++d) {
                    if(o == d) continue;
                        std::sort(od_ways[{o,d}].begin(), od_ways[{o,d}].end());
                        for(int i = 0; i < od_ways[{o,d}].size(); ++i) {
                            od_trips_.push_back(od_ways[{o,d}][i].second); 
                            od_pq_[cc++] = (trip_ways[od_ways[{o,d}][i].second]);
                    }
                }
            }

            od_trips = t_vv(od_trips_);
            od_pq = t_iv(od_pq_);
            od_n_stops = n;

            
            std::cout << "done\n";
        }


        //*******************************************************//
        //******************** FUNCTIONS ************************//
        //*******************************************************//

        //Number of elements
        inline size_type number_diff_trips()const
        {
            return od_pq.size();
        }

        //Returns if the data structure is empty
        inline bool empty()const
        {
            return od_pq.empty();
        }
        
        //Returns the number of existing linestops
        inline size_type n_stops() const {
            return od_n_stops;
        }


        value_type get_people_quantity(id_type o, id_type d) {
            value_type lid = lower_search(0, od_trips.size() - 1, o, 0);
            if(lid == od_trips.size()) return 0; //There doesn't exists any trip that started with o
            value_type rid = upper_search(lid, od_trips.size() - 1, o, 0);
               
            value_type llid = lower_search(lid, rid, d, 1);
            if(llid == od_trips.size()) return 0; //There doesnt't exists any trip that started with o and ended with d
            value_type rrid = upper_search(llid, rid, d, 1);
           
            value_type total = 0;
        
            for(value_type i = llid; i <= rrid; ++i) {
                total += od_pq[i]; 
            }
            return total;
        }

        t_vv get_trips(id_type o, id_type d) {
            value_type lid = lower_search(0, od_trips.size() - 1, o, 0);
            if(lid == od_trips.size()) return {}; //There doesn't exists any trip that started with o

            value_type rid = upper_search(lid, od_trips.size() - 1, o, 0);

            
            value_type llid = lower_search(lid, rid, d, 1);
            if(llid == od_trips.size()) return {}; //There doesnt't exists any trip that started with o and ended with d
            value_type rrid = upper_search(llid, rid, d, 1);
           
            t_vv res((rrid - llid) + 1);

            for(value_type i = llid, c = 0; i <= rrid; ++i, ++c) {
                res[c].resize(od_trips[i].size());
                for(value_type j = 0; j < od_trips[i].size(); ++j) {
                    if(j == 1) 
                        res[c][od_trips[i].size() - 1] = od_trips[i][j];
                    else if(j > 1)
                        res[c][j-1] = od_trips[i][j];
                    else res[c][j] = od_trips[i][j];
                }
            }
            
            return res;
        }

        value_type get_people_quantity_transfer(id_type s) {
            value_type total = 0;
            for(value_type i = 0; i < od_trips.size(); ++i) {
                if(od_trips[i][0] == s) continue;
                if(od_trips[i][1] == s) continue;
                for(value_type j = 2; j < od_trips[i].size(); ++j) {
                    if(od_trips[i][j] == s) {
                        total += od_pq[i];
                        break;
                    }
                }
            }
            return total;
        }

        bool exist_direct_trip(id_type o, id_type d) {
            value_type lid = lower_search(0, od_trips.size() - 1, o, 0);
            if(lid == od_trips.size()) return false; //There doesn't exists any trip that started with o

            value_type rid = upper_search(lid, od_trips.size() - 1, o, 0);

            
            value_type llid = lower_search(lid, rid, d, 1);
            if(llid == od_trips.size()) return false; //There doesnt't exists any trip that started with o and ended with d

            if(od_trips[llid].size() == 2) return true;
            return false;
        }

        //If len == 0 return all transfers
        //if len >= 1 return just with 'len' transfer
        t_vv get_transfer_trips(id_type o, id_type d, uint16_t len) {
            value_type lid = lower_search(0, od_trips.size() - 1, o, 0);
            if(lid == od_trips.size()) return {}; //There doesn't exists any trip that started with o

            value_type rid = upper_search(lid, od_trips.size() - 1, o, 0);

            
            value_type llid = lower_search(lid, rid, d, 1);
            if(llid == od_trips.size()) return {}; //There doesnt't exists any trip that started with o and ended with d
            value_type rrid = upper_search(llid, rid, d, 1);
           
            t_vv res;
            
            uint16_t len1 = len == 0 ? 1 : len;
            uint16_t len2 = len == 0 ? 253 : len;
            //std::cout << len1 << ' ' << len2 << ' ' << len << '\n';
            for(value_type i = llid, c = 0; i <= rrid; ++i) {
                if((od_trips[i].size() - 2) < len1 or (od_trips[i].size() - 2) > len2) continue;
                //std::cout << "here " << i << ' ' << len1 << ' ' << len2 << ' ' << od_trips[i].size() << '\n';
                res.push_back(std::vector<uint16_t>(od_trips[i].size()));
                //std::cout << res[c].size() << '\n';
                for(value_type j = 0; j < od_trips[i].size(); ++j) {
                    if(j == 1) 
                        res[c][od_trips[i].size() - 1] = od_trips[i][j];
                    else if(j > 1)
                        res[c][j-1] = od_trips[i][j];
                    else res[c][j] = od_trips[i][j];
                }
                c++;
            }
            return res;
        }

        value_type get_origin_destination_linestop(id_type s, std::map<std::pair<id_type, id_type>, value_type> &res) {       
            value_type ress = 0;
            value_type lid = lower_search(0, od_trips.size() - 1, s, 0);
            
            if(lid != od_trips.size()) {
                value_type rid = upper_search(lid, od_trips.size() - 1, s, 0);
                for(value_type i = 0; i < lid; ++i) {
                    if(od_trips[i][1] == s) continue;
                    id_type o = od_trips[i][0] + 1;
                    id_type d = od_trips[i][1] + 1;
                    for(value_type j = 2; j < od_trips[i].size(); ++j) {
                        if(od_trips[i][j] == s) {
                            (res.count({o,d}) == 0) ? res[{o,d}] = 1 : res[{o,d}]++;
                            ress++;
                            break;    
                        }
                    }
                }
                //in order to avoid a if comparison if(i == lid) then i = rid+1, so we can avoid those that start
                //with s and were already calculated
                for(value_type i = rid+1; i < od_trips.size(); ++i) {
                    if(od_trips[i][1] == s) continue;
                    id_type o = od_trips[i][0] + 1;
                    id_type d = od_trips[i][1] + 1;
                    for(value_type j = 2; j < od_trips[i].size(); ++j) {
                        if(od_trips[i][j] == s) { 
                            (res.count({o,d}) == 0) ? res[{o,d}] = 1 : res[{o,d}]++;
                            ress++;
                            break;
                        }
                    }
                } 
            }
            else {
                for(value_type i = 0; i < od_trips.size(); ++i) {
                    if(od_trips[i][1] == s) continue;
                    id_type o = od_trips[i][0] + 1;
                    id_type d = od_trips[i][1] + 1;
                    for(value_type j = 2; j < od_trips[i].size(); ++j) {
                        if(od_trips[i][j] == s) {
                            (res.count({o,d}) == 0) ? res[{o,d}] = 1 : res[{o,d}]++;
                            ress++;
                            break;
                        }
                    }
                }
            }
            return ress;
        }
        
        


        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //Swap method
        void swap(od_baseline& o)
        {
            od_trips.swap(o.od_trips);
            //Swap ranks and selects
            od_pq.swap(o.od_pq);
            
            std::swap(od_n_stops, o.od_n_stops);
        }

        //Assignment Operator
        od_baseline& operator=(const od_baseline& o)
        {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //Assignment Move Operator
        od_baseline& operator=(od_baseline&& o)
        {
            if (this != &o) {
                //Move elements
                od_trips = std::move(o.od_trips);
                od_pq = std::move(o.od_pq);
                
                od_n_stops = std::move(o.od_n_stops);
            }
            return *this;
        }

        //Exist
        inline bool trip_exist(size_type o, size_type d) const{
            return false; //If there is a 0, then it means that at least have 1 trip
        }

        //Access operator
        /*value_type operator[](size_type i) const {
            if(m_bitmap[i]){
                auto n_ones = m_rank_bitmap(i+1);
                auto index = n_ones -1;
                return (value_type) m_data[index];
            }
            return -1;
        }*/
        std::vector<id_type> getV(int i) {
            return od_trips[i];
        }


        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//

        //Serialize to a stream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            //Create a node and store the different elements of our structure in that node
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += od_pq.serialize(out, child, "pq");
            written_bytes += serialize_vector(od_trips, out, child, "trips");
            written_bytes += write_member(od_n_stops, out, child, "n_stops");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //Load from a stream
        void load(std::istream& in)
        {
            od_pq.load(in);
            std::cout << od_pq.size() << '\n';
            od_trips.resize(od_pq.size());
            load_vector(od_trips, in);
            read_member(od_n_stops, in);
        }

    };
}

#endif //LBD_OD_REP_V1_HPP
