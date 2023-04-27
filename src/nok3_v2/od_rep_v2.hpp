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


#ifndef LBD_OD_REP_V2_HPP
#define LBD_OD_REP_V2_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>

//This is an example, create your own namespace.
namespace sdsl {

    template <class t_value>
    class od_rep_v2 {

    public:

        //Definiton of types, useful when you need to know a type inside a class
        typedef uint64_t size_type;
        typedef t_value value_type;
        typedef sdsl::bit_vector t_bv;
        //Typename for obtaining the type inside a class
        typedef typename t_bv::rank_1_type t_rank;
        typedef typename t_bv::rank_0_type t_rank0;
        typedef typename t_bv::select_1_type t_select;
        typedef sdsl::int_vector<> t_iv;
        typedef sdsl::int_vector<32> t_32v;
        typedef sdsl::int_vector<16> t_16v;
        typedef sdsl::int_vector<8> t_8v;
        typedef sdsl::enc_vector<> t_ev;

    private:

        t_bv od_matrix;
        t_rank0 od_matrix_rank0;
        t_select od_matrix_select;
        t_ev od_pq;  //people quantity
        t_32v od_tl;  //trip len
        t_16v od_td;  //transfer description
        
        size_type od_matrix_n;


        void copy(const od_rep_v2& o)
        {
            od_matrix = o.od_matrix;
            od_matrix_rank0 = o.od_matrix_rank0;
            od_matrix_rank0.set_vector(&od_matrix);
            od_matrix_select = o.od_matrix_select;
            od_matrix_select.set_vector(&od_matrix);
            
            od_pq = o.od_pq;  //people quantity
            od_tl = o.od_tl;  //trip len
            od_td = o.od_td;  //transfer description
            
            od_matrix_n = o.od_matrix_n;   
        }

    public:

        // If you need access to parts of your class and you don't want to modify their values
        
        const t_bv &matrix = od_matrix;     
        const t_ev &pq = od_pq;  //people quantity
        const t_ev &tl = od_tl;  //trip len
        const t_16v &td = od_td;  //transfer description



		//*******************************************************//
        //******************* CONSTRUCTORS **********************//
        //*******************************************************//

        //Default constructor
        od_rep_v2() = default;

        //Copy constructor
        od_rep_v2(const od_rep_v2& o)
        {
            copy(o);
        }

        //Move constructor
        od_rep_v2(od_rep_v2&& o)
        {
            *this = std::move(o);
        }

        /**
         *  Constructor
         *  Container can be any class.
         *  Depending on the implementation it requires some methods.
         *  In this case we need that the container supports iterators, size, and each element a first and second attribute.
         */
        
        od_rep_v2(std::vector<std::vector<value_type>> &trips, size_type n){
            size_type n_trips = trips.size();
            if(n_trips < 1) {
                throw std::logic_error("There are not enough trips.");
            }
            if(n <= 1) {
                throw std::logic_error("There are not enough linestops.");
            }

            od_matrix_n = n;

            std::cout << "storing trips\n";
            std::map<std::pair<value_type, value_type>, std::vector<std::pair<value_type, std::vector<value_type>>>> od_ways;
            std::map<std::vector<value_type>, value_type> trip_ways;

            value_type matrix_size = n*n;
            value_type dt_size = 0;

            for(size_type i = 0; i < n_trips; ++i) {
                std::vector<value_type> trip = trips[i];
                value_type trip_len = trip.size();
                value_type o = trip[0];
                value_type d = trip[trip_len - 1];
                if(trip_ways.count(trip) <= 0) {
                    trip_ways[trip] = 1;
                    od_ways[{o,d}].push_back({trip_len, trip});
                    matrix_size++;
                }
                else trip_ways[trip]++;
            }
            std::cout << "done\n";
            
            std::cout << "creating matrix\n";
            t_bv od_matrix_(matrix_size, 0);
            size_type acum = 0;

            for(size_type o = 0; o < n; ++o) {
                for(size_type d = 0; d < n; ++d) {
                    acum += od_ways[{o,d}].size();
                    size_type m_id = ((n * o) + (d+1)) - 1;

                    od_matrix_[m_id + acum] = 1;
                    std::sort(od_ways[{o,d}].begin(), od_ways[{o,d}].end());
                }
            }
            od_matrix = t_bv(od_matrix_);
            std::cout << "done\n";

            std::cout << "creating rank0 and select support\n";
            od_matrix_rank0 = t_rank0(&od_matrix);
            od_matrix_select = t_select(&od_matrix);
            std::cout << "done\n";

            assert(matrix_size - (n*n) == od_matrix_rank0(matrix_size));

            std::cout << "creating pq and tl psums\n";
            t_iv od_pqq_(matrix_size - (n*n), 0);
            t_32v od_tlq_(matrix_size - (n*n), 0);
            acum = 0;
            for(size_type o = 0; o < n; ++o) {
                for(size_type d = 0; d < n; ++d) {
                    acum += od_ways[{o,d}].size();
                    size_type m_id = ((n * o) + (d+1)) - 1;

                    size_type m_id_f = od_matrix_rank0(m_id + acum);

                    for(long t_id = od_ways[{o,d}].size() - 1, c = 1; t_id >= 0; --t_id, ++c) {
                        long pqq_id = m_id_f - c;
                        
                        std::vector<value_type> trip = od_ways[{o,d}][t_id].second;
                        od_pqq_[pqq_id] = trip_ways[trip];
                        od_tlq_[pqq_id] = trip.size() - 2;    
                    }
                }
            }

            t_iv od_pq_(od_pqq_.size() + 1, 0);
            t_32v od_tl_(od_tlq_.size() + 1, 0);
            
            for(size_type i = 0; i < od_pqq_.size(); ++i) {
                od_pq_[i+1] = od_pq_[i] + od_pqq_[i];
                od_tl_[i+1] = od_tl_[i] + od_tlq_[i];
            }

            od_pq = t_ev(od_pq_);
            od_tl = t_32v(od_tl_);
            std::cout << "done\n";
            
            std::cout << "creating td\n";
            t_16v od_td_(od_tl_[od_tlq_.size()], 0);
            for(size_type o = 0; o < n; ++o) {
                for(size_type d = 0; d < n; ++d) {
                    if(o == d) continue;
                    
                    size_type m_id = ((n * o) + (d+1)) - 1;
                    
                    size_type m_id_b = od_matrix_select(m_id) + 1;
                    size_type m_id_e = od_matrix_select(m_id + 1);
                    size_type tl_id_b = od_matrix_rank0(m_id_b);

                    for(size_type t_id = 0; t_id < od_ways[{o,d}].size(); ++t_id) {
                        std::vector<value_type> trip = od_ways[{o,d}][t_id].second;
                        size_type td_id_b = od_tl_[tl_id_b];

                        for(size_type td_id = td_id_b, id = 1; id < trip.size()-1; ++td_id, ++id) {
                            od_td_[td_id] = trip[id]; //TODO: BORRAR + 1
                        }
                        tl_id_b++;   
                    }
                }
            }

            od_td = t_16v(od_td_);
            //util::bit_compress(od_td);
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
        inline size_type matrix_n() const {
            return od_matrix_n;
        }


        value_type get_people_quantity(size_type o, size_type d) {
            if(o == d or o >= od_matrix_n or d >= od_matrix_n) return 0;
            int m_id = ((od_matrix_n * o) + (d+1)) - 1;
            if(m_id <= 0) return 0;
            value_type res = 0;
            int bp_id_s = od_matrix_select(m_id) + 1;
            int bp_id_e = od_matrix_select(m_id+1);
            int ts_id_s = od_matrix_rank0(bp_id_s);
            int ts_id_e = ts_id_s + (bp_id_e - bp_id_s);
            res += (od_pq[ts_id_e] - od_pq[ts_id_s]);
            
            return res;
        }

        std::vector<std::vector<t_value>> get_trips(t_value o, t_value d) {
            if(o == d or o >= od_matrix_n or d >= od_matrix_n) return {};
            int m_id = ((od_matrix_n * o) + (d+1)) - 1;
            if(m_id <= 0) return {};


            int bp_id_s = od_matrix_select(m_id) + 1;
            int bp_id_e = od_matrix_select(m_id+1);
            int ts_id_s = od_matrix_rank0(bp_id_s);
            int ts_id_e = ts_id_s + (bp_id_e - bp_id_s);
            std::vector<std::vector<t_value>> res(ts_id_e - ts_id_s);
            for(size_type id = ts_id_s, c = 0; id < ts_id_e; ++id, ++c) {
                size_type td_id = od_tl[id];
                value_type len = od_tl[id+1] - td_id;
                res[c].resize(len+2);
                res[c][0] = o;
                res[c][len+1] = d;
                for(value_type idd = td_id, cc=1; idd < td_id + len; ++idd,++cc) {
                   res[c][cc] = od_td[idd];
                }
            }
            return res; 
        }

        value_type get_people_quantity_linestop(size_type s) {
            if(s >= od_matrix_n) return 0;
            int m_id_start = ((od_matrix_n * s) + 1) - 1;
            int m_id_end = m_id_start + od_matrix_n;
        
            value_type total = 0;
            //start with s;
           // std::cout << "start: " << s << ' ' << od_matrix.size() - (od_pq.size() - 1) << ' '<<m_id_end <<  "\n";
            int bp_id_s = m_id_start > 0 ? od_matrix_select(m_id_start) + 1:0;
            //std::cout << "ded1\n";
            int bp_id_e = od_matrix_select(m_id_end);
            //std::cout << "ded2\n";
            int ts_id_s = od_matrix_rank0(bp_id_s);
            //std::cout << ts_id_s<<"\n";
            int ts_id_e = od_matrix_rank0(bp_id_e); 
            total += od_pq[ts_id_e] - od_pq[ts_id_s];
            //std::cout << ts_id_e<<"\n";
            //std::cout << od_pq.size()<<"\n";
            int size = od_tl.size() - 1;
            
            for(size_type ts_id = 0; ts_id < ts_id_s; ts_id++) {
                size_type idd = od_tl[ts_id];
                size_type len = od_tl[ts_id + 1] - idd;
                for(size_type td_id = idd; td_id < idd + len; ++td_id) {
                    if(od_td[td_id] == s) {
                        total += od_pq[ts_id+1] - od_pq[ts_id];
                        break;
                    }
                }
            }

            for(size_type ts_id = ts_id_e; ts_id < size; ts_id++) {
                size_type idd = od_tl[ts_id];
                size_type len = od_tl[ts_id + 1] - idd;
                for(size_type td_id = idd; td_id < idd + len; ++td_id) {
                    if(od_td[td_id] == s) {
                        total += od_pq[ts_id+1] - od_pq[ts_id];
                        break;
                    }
                }
            }



            return total;
        
        }
        
        


        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //Swap method
        void swap(od_rep_v2& o)
        {
            od_matrix.swap(o.od_matrix);
            //Swap ranks and selects
            sdsl::util::swap_support(od_matrix_rank0, o.od_matrix_rank0,
                               &od_matrix, &(o.od_matrix));
            sdsl::util::swap_support(od_matrix_select, o.od_matrix_select,
                               &od_matrix, &(o.od_matrix));
            od_pq.swap(o.od_pq);
            od_tl.swap(o.od_tl);
            od_td.swap(o.od_td);
            
            std::swap(od_matrix_n, o.od_matrix_n);
        }

        //Assignment Operator
        od_rep_v2& operator=(const od_rep_v2& o)
        {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //Assignment Move Operator
        od_rep_v2& operator=(od_rep_v2&& o)
        {
            if (this != &o) {
                //Move elements
                od_matrix = std::move(o.od_matrix);
                od_matrix_rank0 = std::move(o.od_matrix_rank0);
                od_matrix_rank0.set_vector(&od_matrix);
                od_matrix_select = std::move(o.od_matrix_select);
                od_matrix_select.set_vector(&od_matrix);
                od_pq = std::move(o.pq);
                od_tl = std::move(o.tl);
                od_td = std::move(o.td);

                od_matrix_n = std::move(o.od_matrix_n);
            }
            return *this;
        }

        //Exist
        inline bool trip_exist(size_type o, size_type d) const{
            if(o == d or o >= od_matrix_n or d >= od_matrix_n) return false;
            size_type m_id = ((od_matrix_n * o) + (d+1)) - 1;
            size_type bv_id_s = od_matrix_select(m_id) + 1;
            return !od_matrix[bv_id_s]; //If there is a 0, then it means that at least have 1 trip
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

        //*******************************************************//
        //********************** FILE ***************************//
        //*******************************************************//

        //Serialize to a stream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            //Create a node and store the different elements of our structure in that node
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += od_matrix.serialize(out, child, "matrix");
            written_bytes += od_matrix_rank0.serialize(out, child, "matrix_rank0");
            written_bytes += od_matrix_select.serialize(out, child, "matrix_select");
            written_bytes += od_pq.serialize(out, child, "pq");
            written_bytes += od_td.serialize(out, child, "td");
            written_bytes += od_tl.serialize(out, child, "tl");
            written_bytes += write_member(od_matrix_n, out, child, "matrix_n");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //Load from a stream
        void load(std::istream& in)
        {
            od_matrix.load(in);
            od_matrix_rank0.load(in, &od_matrix);
            od_matrix_select.load(in, &od_matrix);
            od_pq.load(in);
            od_td.load(in);
            od_tl.load(in);
            read_member(od_matrix_n, in);
        }

    };
}

#endif //LBD_OD_REP_V2_HPP
