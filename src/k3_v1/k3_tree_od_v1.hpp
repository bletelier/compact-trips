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

#ifndef LBD_K3_TREE_OD_V1
#define LBD_K3_TREE_OD_V1

#include "k3_tree_base.hpp"

//! Namespace for the succint data structure library
namespace sdsl {


template<uint8_t t_k1=2, uint8_t t_k2=2,
        typename t_bv=bit_vector,
        typename t_ev=enc_vector<>,
        typename t_iv=int_vector<>,
        typename t_select=typename t_bv::select_1_type,
        typename t_rank=typename t_bv::rank_1_type,
        typename t_rank0=typename t_bv::rank_0_type>
class k3_tree_od_v1: public k3_tree_base<>
{

    static_assert(t_k1>1, "t_k has to be larger than 1.");
    static_assert(t_k2>1, "t_k has to be larger than 1.");

    public:
        typedef k3_tree_base::size_type size_type;
        typedef k3_tree_base::size_type pos_type;
        typedef k3_tree_base::tuple_result tuple_result;
        typedef k3_tree_base::point_type point_type;
        typedef std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> region_result;
    private:
        //! Bit array to store all the bits of the tree, except those in the
        //! last level.
        t_bv        k_t;
        t_rank      k_t_rank;
        //! Bit array to store a bit for each 0 in k_t (1 -> it as at least one point | 0 -> empty node)
        t_bv        k_l;
        t_rank      k_l_rank;

        //k3-tree that has 1 l
        t_bv        k_k3_t;
        t_bv        k_k3_l;
        t_rank      k_k3_t_rank;
        t_rank      k_k3_l_rank;
        t_iv        k_k3_id;

        //not necesary if psum is used
        t_bv        k_nv;


        t_ev        k_nv_psum;
        t_iv        k_id;
        t_iv        k_gr;

        t_bv        k_bp;
        t_rank0     k_bp_rank;
        t_select    k_bp_select;
        t_iv        k_ts;

        t_ev        k_pq;

        uint8_t     k_k1;
        uint8_t     k_k2;
        uint16_t    k_height;
        uint16_t    k_matrix_size;
        // tmp params
        size_type k_k1_3;

    protected:

    void init (size_type size) {
        // Set params
        k_k1 = t_k1;
        k_k1_3 = pow(k_k1, 3);
        k_k2 = t_k2;
        k_height = std::ceil(std::log(size)/std::log(k_k1));
        k_height = k_height > 1 ? k_height : 1; // If size == 0
        k_size = pow(k_k1, k_height);
        k_matrix_size = size;
    }


    //*******************************************************//
    //********************** HELPERS ************************//
    //*******************************************************//

    /*! Get the chunk index ([0, k^3[) of a submatrix point.
     *
     * Gets a point in the global matrix and returns its corresponding chunk
     * in the submatrix specified.
     *
     * \param x x coordinate of the point in the global matrix.
     * \param y y coordinate of the point in the global matrix.
     * \param z z coordinate of the point in the global maxtrix.
     * \param x_0 x offset of the submatix in the global matrix.
     * \param y_0 y offset of the submatix in the global matrix.
     * \param z_0 z offset of the submatix in the global matrix.
     * \param l size of the chunk at the submatrix.
     * \param k the k parameter from the k^3 tree.
     * \returns the index of the chunk containing the point at the submatrix.
     */
    inline uint16_t get_chunk_idx(pos_type x, pos_type y, pos_type z,
                                  pos_type x_0, pos_type y_0, pos_type z_0,
                                  size_type l, uint8_t k)
    {
        return  ((x - x_0) / l) * k * k + ((y - y_0) / l) * k + (z - z_0) / l;
    }

    //*******************************************************//
    //********************* BUILD TREE **********************//
    //*******************************************************//
    //! Build a tree from an point collection
    /*! This method takes a vector of points
     *  and the cube size. And takes linear time over the amount of
     *  points to build the k_3 representation.
     *  \param edges A vector with all the points of the cube, it can
     *               not be empty.
     *  \param size Size of the cube, all the nodes in point must be
     *              within 0 and size ([0, size[).
     */
    void build(std::vector<point_type>& points, const size_type size) {

        // (Initial position, end position, submatrix size, base x coordinate, base y coordinate, base z coordinate)
        typedef std::tuple<pos_type, pos_type, size_type, pos_type, pos_type, pos_type > t_part_tuple;

        // Initialize k3-tree values (ks, height, size, etc)
        init(size);

        // Init bit_vectors
        bit_vector k_t_ = bit_vector(k_k1_3 * k_height * points.size(), 0);
        bit_vector k_l_;

        std::queue<t_part_tuple> q;
        pos_type t = 0, last_level = 0;
        pos_type init_p, end_p, x_0, y_0, z_0, it, t_x, t_y, t_z;
        size_type l = std::pow(k_k1, k_height - 1);
        std::vector<pos_type > pos_by_chunk(k_k1_3 + 1, 0);

        q.push(t_part_tuple(0, points.size(), l, 0, 0, 0));

        while (!q.empty()) {
            std::vector<pos_type> amount_by_chunk(k_k1_3, 0);
            std::tie(init_p, end_p, l, x_0, y_0, z_0) = q.front();

            q.pop();

            // Get size for each chunk
            for (it = init_p; it < end_p; it++) {
                amount_by_chunk[get_chunk_idx(
                        std::get<0>(points[it]), std::get<1>(points[it]), std::get<2>(points[it]),
                x_0, y_0, z_0, l, k_k1)] += 1;
            }
            if (l == 1) {
                if (last_level == 0) {
                    last_level = t;
                    k_l_ = bit_vector(k_t_.size() - last_level, 0);
                    k_t_.resize(last_level);
                    last_level = 1; // if t was 0
                    t = 0; // Restart counter as we're storing at k_l_ now.
                }
                for (it = 0; it < k_k1_3; it++,t++)
                    if (amount_by_chunk[it] != 0)
                        k_l_[t] = 1;
                // At l == 1 we do not put new elements at the queue.
                continue;
            }

            // Set starting position in the vector for each chunk
            pos_by_chunk[0] = init_p;
            for (it = 1; it < k_k1_3; it++)
                pos_by_chunk[it] =
                        pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
            // To handle the last case when it = k_2 - 1
            pos_by_chunk[k_k1_3] = end_p;
            // Push to the queue every non zero elements chunk
            for (it = 0; it < k_k1_3; it++,t++) {
                // If not empty chunk, set bit to 1
                if (amount_by_chunk[it] != 0) {
                    t_x = it / (k_k1 * k_k1);
                    t_y = (it / k_k1) % k_k1;
                    t_z = it % k_k1;
                    k_t_[t] = 1;
                    q.push(t_part_tuple(pos_by_chunk[it],
                                        pos_by_chunk[it + 1],
                                        l / k_k1,
                                        x_0 + t_x * l,
                                        y_0 + t_y * l,
                                        z_0 + t_z * l));
                }
            }
            size_type chunk;

            // Sort edges' vector
            for (unsigned ch = 0; ch < k_k1_3; ch++) {
                size_type be = ch == 0 ? init_p : pos_by_chunk[ch - 1];
                for (it = pos_by_chunk[ch]; it < be + amount_by_chunk[ch];) {
                    chunk = get_chunk_idx(
                            std::get<0>(points[it]), std::get<1>(points[it]), std::get<2>(points[it]),
                            x_0, y_0, z_0, l, k_k1);

                    if (pos_by_chunk[chunk] != it)
                        std::iter_swap(points.begin() + it,
                                       points.begin() + pos_by_chunk[chunk]);
                    else
                        it++;
                    pos_by_chunk[chunk]++;
                }
            }
        }
        k_l_.resize(t);
        //k2_tree_ns::build_template_vector<t_bv>(k_t_, k_l_, k_t, k_l);
        k_t = t_bv(k_t_);
        k_l = t_bv(k_l_);

        k_t_rank = t_rank(&k_t);
        k_l_rank = t_rank(&k_l);
    }
    
    void build2(std::vector<point_type>& points, const size_type size) {

        // (Initial position, end position, submatrix size, base x coordinate, base y coordinate, base z coordinate)
        typedef std::tuple<pos_type, pos_type, size_type, pos_type, pos_type, pos_type > t_part_tuple;

        // Initialize k3-tree values (ks, height, size, etc)
        //init(size);

        // Init bit_vectors
        bit_vector k_t_ = bit_vector(k_k1_3 * k_height * points.size(), 0);
        bit_vector k_l_;

        std::queue<t_part_tuple> q;
        pos_type t = 0, last_level = 0;
        pos_type init_p, end_p, x_0, y_0, z_0, it, t_x, t_y, t_z;
        size_type l = std::pow(k_k1, k_height - 1);
        std::vector<pos_type > pos_by_chunk(k_k1_3 + 1, 0);

        q.push(t_part_tuple(0, points.size(), l, 0, 0, 0));

        while (!q.empty()) {
            std::vector<pos_type> amount_by_chunk(k_k1_3, 0);
            std::tie(init_p, end_p, l, x_0, y_0, z_0) = q.front();

            q.pop();

            // Get size for each chunk
            for (it = init_p; it < end_p; it++) {
                amount_by_chunk[get_chunk_idx(
                        std::get<0>(points[it]), std::get<1>(points[it]), std::get<2>(points[it]),
                x_0, y_0, z_0, l, k_k1)] += 1;
            }
            if (l == 1) {
                if (last_level == 0) {
                    last_level = t;
                    k_l_ = bit_vector(k_t_.size() - last_level, 0);
                    k_t_.resize(last_level);
                    last_level = 1; // if t was 0
                    t = 0; // Restart counter as we're storing at k_l_ now.
                }
                for (it = 0; it < k_k1_3; it++,t++)
                    if (amount_by_chunk[it] != 0)
                        k_l_[t] = 1;
                // At l == 1 we do not put new elements at the queue.
                continue;
            }

            // Set starting position in the vector for each chunk
            pos_by_chunk[0] = init_p;
            for (it = 1; it < k_k1_3; it++)
                pos_by_chunk[it] =
                        pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
            // To handle the last case when it = k_2 - 1
            pos_by_chunk[k_k1_3] = end_p;
            // Push to the queue every non zero elements chunk
            for (it = 0; it < k_k1_3; it++,t++) {
                // If not empty chunk, set bit to 1
                if (amount_by_chunk[it] != 0) {
                    t_x = it / (k_k1 * k_k1);
                    t_y = (it / k_k1) % k_k1;
                    t_z = it % k_k1;
                    k_t_[t] = 1;
                    q.push(t_part_tuple(pos_by_chunk[it],
                                        pos_by_chunk[it + 1],
                                        l / k_k1,
                                        x_0 + t_x * l,
                                        y_0 + t_y * l,
                                        z_0 + t_z * l));
                }
            }
            size_type chunk;

            // Sort edges' vector
            for (unsigned ch = 0; ch < k_k1_3; ch++) {
                size_type be = ch == 0 ? init_p : pos_by_chunk[ch - 1];
                for (it = pos_by_chunk[ch]; it < be + amount_by_chunk[ch];) {
                    chunk = get_chunk_idx(
                            std::get<0>(points[it]), std::get<1>(points[it]), std::get<2>(points[it]),
                            x_0, y_0, z_0, l, k_k1);

                    if (pos_by_chunk[chunk] != it)
                        std::iter_swap(points.begin() + it,
                                       points.begin() + pos_by_chunk[chunk]);
                    else
                        it++;
                    pos_by_chunk[chunk]++;
                }
            }
        }
        k_l_.resize(t);
        //k2_tree_ns::build_template_vector<t_bv>(k_t_, k_l_, k_t, k_l);
        k_k3_t = t_bv(k_t_);
        k_k3_l = t_bv(k_l_);

        k_k3_t_rank = t_rank(&k_k3_t);
        k_k3_l_rank = t_rank(&k_k3_l);
    }



    public:

    //*******************************************************//
    //******************* CONSTRUCTOR ***********************//
    //*******************************************************//
    k3_tree_od_v1() = default;

    k3_tree_od_v1(const k3_tree_od_v1& tr)
    {
        *this = tr;
    }

    k3_tree_od_v1(k3_tree_od_v1&& tr)
    {
        *this = std::move(tr);
    }

    k3_tree_od_v1(std::string filename, size_type size)
    {
        // Initialize parameters
        init(size);

        // Open data file
        std::ifstream data_file(filename);
        assert(data_file.is_open() && data_file.good());

        // Read file and create conceptual tree
        pos_type pos_x, pos_y, pos_z;

        std::vector<point_type> points;
        while (!data_file.eof() && data_file.good()) {

            // Get position (x, y, z)
            read_member(pos_x, data_file);
            read_member(pos_y, data_file);
            read_member(pos_z, data_file);
//            std::cout << "Adding point " << pos_x << ", " << pos_y << ", " << pos_z << ")" << std::endl;

            points.push_back(point_type(pos_x, pos_y, pos_z));
        }

        // Remove duplicates
        std::sort( points.begin(), points.end() );
        points.erase( std::unique( points.begin(), points.end() ), points.end() );

        // Build structures
        build(points, size);
    }

    k3_tree_od_v1(std::vector<std::vector<int>> &trips, size_type size) {
        std::map<std::pair<int,int>, std::vector<std::pair<int, std::vector<int>>>> od_ways;
        std::map<std::vector<int>, int> trip_ways;
        std::map<std::pair<std::pair<int,int>,int>, std::vector<std::pair<int,int>>> transfer_map;
        int bv_size = size*size;
        std::vector<point_type> matrix;
        std::vector<point_type> matrix_otd;
        int size_nv = 0;
        for(int i = 0; i < trips.size(); ++i) {
            std::vector<int> trip = trips[i];
            int trip_len = trip.size();
            int o = trip[0];
            int d = trip[trip_len - 1];
            if(trip_len > 3){
                for(int j = 1; j < trip_len - 1; ++j) {
                    matrix.push_back({o,d,trip[j]});
                }
            }
            else if(trip_len == 3) matrix_otd.push_back({o,d,trip[1]});
            if(trip_ways.count(trip) <= 0) {
                trip_ways[trip] = 1;
                od_ways[{o,d}].push_back({trip_len, trip});
                bv_size++;
            }
            else trip_ways[trip]++;
        }
        
        std::cout << "stored trips\n";

        std::map<std::pair<std::pair<int,int>,int>, int> transfer_1_map;
        
        t_bv k_bp_(bv_size, 0);
        int_vector<> k_ts_(bv_size - (size*size), 0);	
	    t_iv k_pq_(k_ts_.size()+1, 0);
        int acum = 0;
	    for(int x = 0; x < size; ++x) { //TODO: IMPROVE THIS USING EXISTING (X,Y) TRIPS
    	    for(int y = 0; y < size; ++y) {
	        acum += od_ways[{x,y}].size();
                int id = ((size * x) + (y+1)) - 1;
                k_bp_[id+acum] = 1;
                std::sort(od_ways[{x,y}].begin(), od_ways[{x,y}].end());
                int tm_id = 0;
                int tm_1_id = 0;
                int t1_add = 0;
                int t2_add = 0;
                for(int j = 0; j < od_ways[{x,y}].size(); ++j) {
                    int o = od_ways[{x,y}][j].second[0];
                    int d = od_ways[{x,y}][j].second[od_ways[{x,y}][j].second.size() - 1];
                    if(od_ways[{x,y}][j].second.size() <= 2) {
                        t1_add++; 
                        //t2_add++; 
                        continue;
                    }
                    else if(od_ways[{x,y}][j].second.size() == 3) {
                        std::vector<int> tr = od_ways[{x,y}][j].second;
                        transfer_1_map[{{o,d},tr[1]}] = tm_1_id + t1_add;
                        tm_1_id++;
                        //t2_add++;
                        continue;   
                    }
                    for(int k = 1; k < od_ways[{x,y}][j].second.size() - 1; ++k) {
		                if (transfer_map.count({{o,d}, od_ways[{x,y}][j].second[k]}) >= 1) size_nv++;
		                    transfer_map[{{o,d}, od_ways[{x,y}][j].second[k]}].push_back({k-1,tm_id});
		            }
		            tm_id++;	
	            }
	        }
    	}
        std::cout << "created bitmaps\n";
        k_bp = t_bv(k_bp_);
        k_bp_rank = t_rank0(&k_bp);
        acum = 0; 
        t_iv k_pqq_ = t_iv(k_ts_.size(), 0); 
        for(int x = 0; x < size; ++x) {
            for(int y = 0; y < size; ++y) {
                acum += od_ways[{x,y}].size();
                int id = ((size * x) + (y+1)) - 1;
                int idd = k_bp_rank(id+acum);
                for(int ts = od_ways[{x,y}].size() - 1, c = 1; ts >= 0; --ts, ++c) {
                    int t_id = idd - c;
                    std::vector<int> trip = od_ways[{x,y}][ts].second;
                    k_ts_[t_id] = trip.size() - 2;
                    k_pqq_[t_id] = trip_ways[trip]; 
                }

            }
        }
        for(int i = 0; i < k_pqq_.size(); ++i) {
            k_pq_[i+1] = k_pq_[i] + k_pqq_[i];
        }
        k_bp_select = t_select(&k_bp);
        k_ts = int_vector<>(k_ts_);
        k_pq = t_ev(k_pq_);
        std::cout << "created psum for pq\n";
        util::bit_compress(k_ts);
        //util::bit_compress(k_pq);
        
        build(matrix, size);
        std::cout << "builded k3-transfer > 1\n";
    	build2(matrix_otd, size);
        std::cout << "builded k3-transfer = 1\n";
        size_nv += k_l_rank(k_l.size());
        t_bv k_nv_ (size_nv, 0);
        int_vector<> k_qnv(k_l_rank(k_l.size()), 0);
         
        for(int i = 0; i < matrix.size(); ++i) {
            int a = std::get<0>(matrix[i]);
            int b = std::get<1>(matrix[i]);
            int c = std::get<2>(matrix[i]);
            int l_id = get_l_id(a, b, c);
            int l_rank_id = k_l_rank(l_id);
            k_qnv[l_rank_id] = transfer_map[{{a, b}, c}].size();
        }
        std::cout << "created k_qnv\n";
        int_vector<> id(size_nv, 0);
        int_vector<> k3_id(k_k3_l_rank(k_k3_l.size()), 0);
        int_vector<> gr(size_nv, 0);
        int_vector<> qnv_psum(k_qnv.size() + 1, 0);
        for(int i = 0; i < k_qnv.size(); ++i) {
            qnv_psum[i+1] = k_qnv[i] + qnv_psum[i];
            int nv_id = qnv_psum[i+1] - 1;
            k_nv_[nv_id] = 1;
        }
        std::cout << "created k_nv\n";
        for(int i = 0; i < matrix.size(); ++i) {
            int a = std::get<0>(matrix[i]);
            int b = std::get<1>(matrix[i]);
            int c = std::get<2>(matrix[i]);
            int l_id = get_l_id(a, b, c);
            int l_rank_id = k_l_rank(l_id);
            for(int j = qnv_psum[l_rank_id], k = 0; j < qnv_psum[l_rank_id + 1]; ++j, k++) {
                id[j] = transfer_map[{{a,b},c}][k].second;
                gr[j] = transfer_map[{{a,b},c}][k].first;
            }
        }
        std::cout << "created id and gr\n";

        for(int i = 0; i < matrix_otd.size(); ++i) {
            int a = std::get<0>(matrix_otd[i]);
            int b = std::get<1>(matrix_otd[i]);
            int c = std::get<2>(matrix_otd[i]);
            int l_id = get_k3_l_id(a, b, c);
            int l_rank_id = k_k3_l_rank(l_id);
            k3_id[l_rank_id] = transfer_1_map[{{a,b},c}];
        }
        std::cout << "created k3_id\n";

        k_id = int_vector<>(id);
        k_nv_psum = enc_vector<>(qnv_psum);
        k_gr = int_vector<>(gr);
        
        k_k3_id = int_vector<>(k3_id);

        k_nv = t_bv(k_nv_);
   	    sdsl::util::bit_compress(k_id);
   	    sdsl::util::bit_compress(k_gr);
   	    sdsl::util::bit_compress(k_k3_id);
    }
    
    
    k3_tree_od_v1(std::vector<point_type> &points, size_type size)
    {
        build(points, size);
    }

    //*******************************************************//
    //*************** BASIC OPERATIONS **********************//
    //*******************************************************//

    //! Move assignment operator
    k3_tree_od_v1& operator=(k3_tree_od_v1&& tr)
    {
        if (this != &tr) {
            k_t = std::move(tr.k_t);
            k_l = std::move(tr.k_l);
            k_k1 = std::move(tr.k_k1);
            k_k2 = std::move(tr.k_k2);
            k_height = std::move(tr.k_height);
            k_t_rank = std::move(tr.k_t_rank);
            k_t_rank.set_vector(&k_t);
        }
        return *this;
    }

    //! Assignment operator
    k3_tree_od_v1& operator=(k3_tree_od_v1& tr)
    {
        if (this != &tr) {
            k_t = tr.k_t;
            k_l = tr.k_l;
            k_t_rank = tr.k_t_rank;
            k_t_rank.set_vector(&k_t);
            k_k1 = tr.k_k1;
            k_k2 = tr.k_k2;
            k_height = tr.k_height;
        }
        return *this;
    }

    //! Swap operator
    void swap(k3_tree_od_v1& tr)
    {
        if (this != &tr) {
            std::swap(k_t, tr.k_t);
            std::swap(k_l, tr.k_l);
            util::swap_support(k_t_rank, tr.k_t_rank, &k_t, &(tr.k_t));
            std::swap(k_k1, tr.k_k1);
            std::swap(k_k2, tr.k_k2);
            std::swap(k_height, tr.k_height);
        }
    }

    //! Equal operator
    bool operator==(const k3_tree_od_v1& tr) const
    {
        // TODO check the rank support equality?
        if (k_k1 != tr.k_k1 || k_k2 != tr.k_k2 || k_height != tr.k_height)
            return false;
        if (k_t.size() != tr.k_t.size() || k_l.size() != tr.k_l.size())
            return false;
        for (unsigned i = 0; i < k_t.size(); i++)
            if (k_t[i] != tr.k_t[i])
                return false;
        for (unsigned i = 0; i < k_l.size(); i++)
            if (k_l[i] != tr.k_l[i])
                return false;
        return true;
    }

    //*******************************************************//
    //********************** QUERIES ************************//
    //*******************************************************//

    bool get(pos_type pos_x, pos_type pos_y, pos_type pos_z) {

        if (k_t.size() == 0) {
            // Empty matrix
            return false;
        }

        if (pos_x > k_size || pos_y > k_size || pos_z > k_size) {
            return false;
        }

        // Size of the submatrix at the current level
        size_type submatrix_size_l = k_size;
        size_type node_pos;
        size_type children_pos = 0;
        uint8_t k = k_k1;


        for (uint16_t l = 0; l < k_height-1; l++) {
            submatrix_size_l /= k;
            node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
            node_pos += children_pos;

            if (k_t[node_pos] == 0) {
                return false; // Empty submatrix
            } else {
                // Go to next level
                children_pos = k_t_rank (node_pos+1) * k_k1_3;
            }

            // Calculate local position on the current submatrix
            pos_x %= submatrix_size_l;
            pos_y %= submatrix_size_l;
            pos_z %= submatrix_size_l;
        }

        // Last level
        submatrix_size_l /= k;
        node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
        node_pos += children_pos;
        return (k_l[node_pos - k_t.size()] == 1);
    }

    unsigned int get_people_quantity(int o, int d) {
        if(o == d or o >= k_matrix_size or d >= k_matrix_size) return 0;
        int m_id = ((k_matrix_size * o) + (d+1)) - 1;
        if(m_id <= 0) return 0;
        unsigned int res = 0;
        int bp_id_s = k_bp_select(m_id) + 1;
        int bp_id_e = k_bp_select(m_id+1);
        int ts_id_s = k_bp_rank(bp_id_s);
        int ts_id_e = ts_id_s + (bp_id_e - bp_id_s);
        res += (k_pq[ts_id_e] - k_pq[ts_id_s]);
        
        /*for(int i = bp_id_s; i < bp_id_e; ++i) {
            res+=k_pq[ts_id++];
        }
        while(!k_bp[bp_id]) {
            res += k_pq[ts_id];
            bp_id++;
            ts_id++;
        }*/
        return res;
    }
    unsigned int get_origin_destinations_linestop(int s, std::map<std::pair<uint16_t, uint16_t>, uint32_t> &res) {
        unsigned int total = 0;        
        //People who made at least 1 transfer and the transfer is S
        {
            std::vector<tuple_result> points;
            get_region_k3(0,0,s,k_matrix_size-1,k_matrix_size-1,s,points);
            int o, d;
            for(int i = 0; i < points.size(); ++i) {
                o = std::get<0>(points[i]) + 1;
                d = std::get<1>(points[i]) + 1;
                res[{o,d}] = 1;
                total++;
            } 
        }

        //People who made 2 or more transfer and one of those transfer is S
        {
            std::vector<region_result> points;
            get_regionl(0,0,s,k_matrix_size-1,k_matrix_size-1,s,points);
            int o, d;
            for(int i = 0; i < points.size(); ++i) {
                o = std::get<0>(points[i]) + 1;
                d = std::get<1>(points[i]) + 1;
                int l_rank_id = std::get<3>(points[i]);
                int n_t = (k_nv_psum[l_rank_id + 1] - k_nv_psum[l_rank_id]);
                (res.count({o,d}) == 0) ? res[{o,d}] = n_t : res[{o,d}] += n_t;
                total += n_t;
            } 
        }
        return total;
    }

    int get_matrix_size() {
        return k_matrix_size;
    }
    
    std::vector<std::vector<int>> get_trips(int o, int d) {
        if(o == d) return {};
        int m_id = ((k_matrix_size * o) + (d+1)) - 1;
        if(m_id <= 0) return {};
        int bp_id = k_bp_select(m_id) + 1;
        int bp_id_e = k_bp_select(m_id + 1);
        int prev_res = 0;
        int ts_id = k_bp_rank(bp_id);
        int cc = 0;
        if(bp_id_e - bp_id == 0) return {};
        std::vector<std::vector<int>> res(bp_id_e - bp_id);
        bool three_stages_trip = false;
        bool four_plus_stages_trip = false;
        for(int i = bp_id, tsid = ts_id; i < bp_id_e; ++i, ++tsid) {
            int val = k_ts[tsid];
            res[cc].resize(val+2);
            res[cc][0] = o;
            res[cc][val+1] = d;
            cc++;
        }
        if(res[0].size() == 2) prev_res++;
        else if(res[0].size() == 3) three_stages_trip = true;
        if(res.size() > 1 and res[1].size() == 3) three_stages_trip = true;
        if(res[cc-1].size() > 3) four_plus_stages_trip = true;
        //std::cout << '\n';
        if(three_stages_trip) {
            std::vector<tuple_result> points;
            get_region_k3(o,d,0,o,d,k_matrix_size-1,points);
            int t;
            for(int i = 0; i < points.size(); ++i) {
                t = std::get<2>(points[i]);

                res[prev_res++][1] = t;
            }
        }

        if(four_plus_stages_trip) {
            std::vector<region_result> points;
            get_regionl(o,d,0,o,d,k_matrix_size-1,points);
            if(points.size() <= 0) return res;
            
            int t;
            for(int i = 0; i < points.size(); ++i) {
                t = std::get<2>(points[i]);
                int l_rank_id = std::get<3>(points[i]); 
                int nv_id = k_nv_psum[l_rank_id];
                int n_t = (k_nv_psum[l_rank_id+1] - nv_id); //We can try using the k_nv bitvector and 1 select operation, then check further if there is a 1 indicating the end of that (x,y)
                for(int offset = 0; offset < n_t; ++offset) {
                    int gen_id = nv_id + offset;
                    int id = k_id[gen_id];
                    
                    int gr = k_gr[gen_id];
                    
                    
                    res[prev_res + id][gr+1] = t;
                } 
            }
        }
        return res;
        
    }


    int get_k3_l_id(pos_type pos_x, pos_type pos_y, pos_type pos_z) {

        if (k_k3_t.size() == 0) {
            // Empty matrix
            return false;
        }
	
        if (pos_x > k_size || pos_y > k_size || pos_z > k_size) {
            return false;
        }

        // Size of the submatrix at the current level
        size_type submatrix_size_l = k_size;
        size_type node_pos;
        size_type children_pos = 0;
        uint8_t k = k_k1;


        for (uint16_t l = 0; l < k_height-1; l++) {
            submatrix_size_l /= k;
            node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
            node_pos += children_pos;

            if (k_k3_t[node_pos] == 0) {
                return false; // Empty submatrix
            } else {
                // Go to next level
                children_pos = k_k3_t_rank (node_pos+1) * k_k1_3;
            }

            // Calculate local position on the current submatrix
            pos_x %= submatrix_size_l;
            pos_y %= submatrix_size_l;
            pos_z %= submatrix_size_l;
        }


        // Last level
        submatrix_size_l /= k;
        node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
        node_pos += children_pos;
        return (node_pos - k_k3_t.size());
    }

    int get_l_id(pos_type pos_x, pos_type pos_y, pos_type pos_z) {

        if (k_t.size() == 0) {
            // Empty matrix
            return false;
        }
	
        if (pos_x > k_size || pos_y > k_size || pos_z > k_size) {
            return false;
        }

        // Size of the submatrix at the current level
        size_type submatrix_size_l = k_size;
        size_type node_pos;
        size_type children_pos = 0;
        uint8_t k = k_k1;


        for (uint16_t l = 0; l < k_height-1; l++) {
            submatrix_size_l /= k;
            node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
            node_pos += children_pos;

            if (k_t[node_pos] == 0) {
                return false; // Empty submatrix
            } else {
                // Go to next level
                children_pos = k_t_rank (node_pos+1) * k_k1_3;
            }

            // Calculate local position on the current submatrix
            pos_x %= submatrix_size_l;
            pos_y %= submatrix_size_l;
            pos_z %= submatrix_size_l;
        }


        // Last level
        submatrix_size_l /= k;
        node_pos = (pos_x / submatrix_size_l) * k * k + pos_y / submatrix_size_l * k + pos_z / submatrix_size_l;
        node_pos += children_pos;
        return (node_pos - k_t.size());
    }
    size_type get_region(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<tuple_result> &result) {
        if (i_pos_x > k_size || i_pos_y > k_size || i_pos_z > k_size) {
            return 0;
        }

        return get_region_r(i_pos_x, i_pos_y, i_pos_z,
                          std::min(e_pos_x, k_size), std::min(e_pos_y, k_size), std::min(e_pos_z, k_size), result,
                          0, 0, 0, k_size, 0, 0);
    }

    size_type get_region_k3_lr(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<region_result> &result) {
        if (i_pos_x > k_size || i_pos_y > k_size || i_pos_z > k_size) {
            return 0;
        }

        return get_region_k3_rl(i_pos_x, i_pos_y, i_pos_z,
                          std::min(e_pos_x, k_size), std::min(e_pos_y, k_size), std::min(e_pos_z, k_size), result,
                          0, 0, 0, k_size, 0, 0);
    
    
    }
 

    size_type get_region_k3(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<tuple_result> &result) {
        if (i_pos_x > k_size || i_pos_y > k_size || i_pos_z > k_size) {
            return 0;
        }

        return get_region_k3_r(i_pos_x, i_pos_y, i_pos_z,
                          std::min(e_pos_x, k_size), std::min(e_pos_y, k_size), std::min(e_pos_z, k_size), result,
                          0, 0, 0, k_size, 0, 0);
    
    
    }
 
    size_type get_regionl(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<region_result> &result) {
        if (i_pos_x > k_size || i_pos_y > k_size || i_pos_z > k_size) {
            return 0;
        }

        return get_region_rl(i_pos_x, i_pos_y, i_pos_z,
                          std::min(e_pos_x, k_size), std::min(e_pos_y, k_size), std::min(e_pos_z, k_size), result,
                          0, 0, 0, k_size, 0, 0);
    
    
    }
    size_type get_region_2(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<tuple_result> &result) {

        if (i_pos_x > k_size || i_pos_y > k_size || i_pos_z > k_size) {
            return 0;
        }

        return get_region_r(i_pos_x, i_pos_y, i_pos_z,
                          std::min(e_pos_x, k_size), std::min(e_pos_y, k_size), std::min(e_pos_z, k_size), result,
                          0, 0, 0, k_size, 0, 0);
    }


    size_type get_region_l(pos_type i_pos_x, pos_type i_pos_y, pos_type i_pos_z,
                           pos_type e_pos_x, pos_type e_pos_y, pos_type e_pos_z,
                         std::vector<region_result> &result) {

        if (i_pos_x >= k_size || i_pos_y >= k_size || i_pos_z >= k_size) {
            return 0;
        }
	
        typedef std::tuple<pos_type, pos_type, pos_type, pos_type, pos_type, pos_type,
                size_type, pos_type, pos_type, pos_type, size_type, uint16_t> t_part_tuple;
        std::stack<t_part_tuple> q;

        q.push(t_part_tuple(i_pos_x, i_pos_y, i_pos_z, e_pos_x, e_pos_y, e_pos_z,
                            k_size, 0, 0, 0, 0, 0));


        uint8_t k = k_k1;
        pos_type x_i_b, x_e_b, y_i_b, y_e_b, z_i_b, z_e_b;  // First and end child
        pos_type xi, xe, yi, ye, zi, ze;                    // Local positions
        pos_type bx, by, bz;                                // Bases positions
        size_type sub_size;                                 // Size of each submatrix
        size_type pos, pos_children;
        uint16_t l;

        while (!q.empty()) {
            std::tie(xi, yi, zi, xe, ye, ze, sub_size, bx, by, bz, pos_children, l) = q.top();
            q.pop();

            // Calculate submatrix size
            sub_size = sub_size / k;

            // For each child that has cells that overlap with the searched region
            x_i_b = xi/sub_size;
            x_e_b = xe/sub_size;
            y_i_b = yi/sub_size;
            y_e_b = ye/sub_size;
            z_i_b = zi/sub_size;
            z_e_b = ze/sub_size;


            for (size_type x = x_i_b; x <= x_e_b; x++) {
                for (size_type y = y_i_b; y <= y_e_b; y++) {
                    for (size_type z = z_i_b; z <= z_e_b; z++) {
                        pos = pos_children + x * k * k + y * k + z; // Position of the current child

                        if (l < (k_height-1)) {
                            // Internal nodes
                            if (k_t[pos] == 1) {
                                // Continue with the search process
                                q.push(t_part_tuple(x == x_i_b ? xi % sub_size : 0,
                                                    y == y_i_b ? yi % sub_size : 0,
                                                    z == z_i_b ? zi % sub_size : 0,
                                                    x == x_e_b ? xe % sub_size : sub_size - 1,
                                                    y == y_e_b ? ye % sub_size : sub_size - 1,
                                                    z == z_e_b ? ze % sub_size : sub_size - 1,
                                                    sub_size, bx + x * sub_size, by + y * sub_size,  bz + z * sub_size,
                                                    k_t_rank(pos + 1) * k_k1_3, l+1));

                            } // ENF IF t[pos] == 1
                        } else {
                            // Leaves nodes
                            if (k_l[pos - k_t.size()] == 1) {
                                if(bx + x < k_size and by + y < k_size and bz + z < k_size)
                                    result.push_back(region_result(bx + x, by + y, bz + z, k_l_rank(pos - k_t.size())));
                            }
                        } // END IF check level
                    } // END FOR z
                } // END FOR y
            } // END FOR x
        } // END WHILE queue
        return result.size();
    }
    
    size_type get_region_k3_l(pos_type i_pos_x, pos_type i_pos_y, pos_type i_pos_z,
                           pos_type e_pos_x, pos_type e_pos_y, pos_type e_pos_z,
                         std::vector<region_result> &result) {

        if (i_pos_x >= this->k_size || i_pos_y >= this->k_size || i_pos_z >= this->k_size) {
                return 0;
            }

            typedef std::tuple<pos_type, pos_type, pos_type, pos_type, pos_type, pos_type,
                    size_type, pos_type, pos_type, pos_type, size_type, uint16_t> t_part_tuple_r;
            std::stack<t_part_tuple_r> q;

            q.push(t_part_tuple_r(i_pos_x, i_pos_y, i_pos_z, e_pos_x, e_pos_y, e_pos_z,
                                this->k_size, 0, 0, 0, 0, 0));


            uint8_t k = this->k_k1;
            pos_type x_i_b, x_e_b, y_i_b, y_e_b, z_i_b, z_e_b;  // First and end child
            pos_type xi, xe, yi, ye, zi, ze;                    // Local positions
            pos_type bx, by, bz;                                // Bases positions
            size_type sub_size;                                 // Size of each submatrix
            size_type pos, pos_children;
            uint16_t l;

            while (!q.empty()) {
                std::tie(xi, yi, zi, xe, ye, ze, sub_size, bx, by, bz, pos_children, l) = q.top();
                q.pop();

                // Calculate submatrix size
                sub_size = sub_size / k;

                // For each child that has cells that overlap with the searched region
                x_i_b = xi/sub_size;
                x_e_b = xe/sub_size;
                y_i_b = yi/sub_size;
                y_e_b = ye/sub_size;
                z_i_b = zi/sub_size;
                z_e_b = ze/sub_size;


                for (size_type x = x_i_b; x <= x_e_b; x++) {
                    for (size_type y = y_i_b; y <= y_e_b; y++) {
                        for (size_type z = z_i_b; z <= z_e_b; z++) {
                            pos = pos_children + x * k * k + y * k + z; // Position of the current child

                            if (l < (this->k_height-1)) {
                                // Internal nodes
                                if (this->k_k3_t[pos] == 1) {
                                    // Continue with the search process
                                    q.push(t_part_tuple_r(x == x_i_b ? xi % sub_size : 0,
                                                        y == y_i_b ? yi % sub_size : 0,
                                                        z == z_i_b ? zi % sub_size : 0,
                                                        x == x_e_b ? xe % sub_size : sub_size - 1,
                                                        y == y_e_b ? ye % sub_size : sub_size - 1,
                                                        z == z_e_b ? ze % sub_size : sub_size - 1,
                                                        sub_size, bx + x * sub_size, by + y * sub_size,  bz + z * sub_size,
                                                        this->k_k3_t_rank(pos + 1) * this->k_k1_3, l+1));

                                } // ENF IF t[pos] == 1
                            } else {
                                // Leaves nodes
                                if (this->k_k3_l[pos - this->k_k3_t.size()] == 1) {
                                    if(bx + x < k_size and by + y < k_size and bz + z < k_size)
                                        result.push_back({bx + x, by + y, bz + z, k_k3_l_rank(pos - k_k3_t.size())});
                                }
                            } // END IF check level
                        } // END FOR z
                    } // END FOR y
                } // END FOR x
            } // END WHILE queue
            return result.size();
        }


    
   size_type get_region_rl(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<region_result> &result,
                         size_type base_x, size_type base_y, size_type base_z,
                         size_type sub_size, uint16_t level, size_type pos_children) {

       size_type pos;
       uint8_t k = k_k1;
       size_type children_size = sub_size / k;
       size_type b_x, b_y, b_z;
       size_type n_points = 0;


       // For each child that has cells that overlap with the searched region
       for (size_type x = i_pos_x/children_size; x <= (e_pos_x/children_size); x++) {
           b_x = x * children_size;
           for (size_type y = i_pos_y/children_size; y <= (e_pos_y/children_size); y++) {
               b_y = y * children_size;
               for (size_type z = i_pos_z/children_size; z <= (e_pos_z/children_size); z++) {
                   pos = pos_children + x * k * k + y * k + z; // Position of the current child
                   b_z = z * children_size;


                   if (level < (k_height-1)) {
                       // Internal nodes
                       if (k_t[pos] == 1) {
                           // Continue with the search process
                           n_points += get_region_rl(std::max(b_x, i_pos_x) - b_x,
                                      std::max(b_y, i_pos_y) - b_y,
                                      std::max(b_z, i_pos_z) - b_z,
                                      std::min(b_x + children_size-1, e_pos_x) - b_x,
                                      std::min(b_y + children_size-1, e_pos_y) - b_y,
                                      std::min(b_z + children_size-1, e_pos_z) - b_z,
                                      result,
                                      base_x + x * children_size, base_y + y * children_size, base_z + z * children_size,
                                      children_size, level + 1, k_t_rank(pos + 1) * k_k1_3);
                       }
                   } else {
                       // Leaves nodes
                       if (k_l[pos - k_t.size()] == 1) {
                           result.push_back(region_result(base_x + x, base_y + y, base_z + z, k_l_rank(pos - k_t.size())));
                           n_points++;
                       }
                   }
               } // END FOR z
           } // END FOR y
       } // END FOR x
        return n_points;
    };
    size_type get_region_k3_rl(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<region_result> &result,
                         size_type base_x, size_type base_y, size_type base_z,
                         size_type sub_size, uint16_t level, size_type pos_children) {

       size_type pos;
       uint8_t k = k_k1;
       size_type children_size = sub_size / k;
       size_type b_x, b_y, b_z;
       size_type n_points = 0;


       // For each child that has cells that overlap with the searched region
       for (size_type x = i_pos_x/children_size; x <= (e_pos_x/children_size); x++) {
           b_x = x * children_size;
           for (size_type y = i_pos_y/children_size; y <= (e_pos_y/children_size); y++) {
               b_y = y * children_size;
               for (size_type z = i_pos_z/children_size; z <= (e_pos_z/children_size); z++) {
                   pos = pos_children + x * k * k + y * k + z; // Position of the current child
                   b_z = z * children_size;


                   if (level < (k_height-1)) {
                       // Internal nodes
                       if (k_k3_t[pos] == 1) {
                           // Continue with the search process
                           n_points += get_region_k3_rl(std::max(b_x, i_pos_x) - b_x,
                                      std::max(b_y, i_pos_y) - b_y,
                                      std::max(b_z, i_pos_z) - b_z,
                                      std::min(b_x + children_size-1, e_pos_x) - b_x,
                                      std::min(b_y + children_size-1, e_pos_y) - b_y,
                                      std::min(b_z + children_size-1, e_pos_z) - b_z,
                                      result,
                                      base_x + x * children_size, base_y + y * children_size, base_z + z * children_size,
                                      children_size, level + 1, k_k3_t_rank(pos + 1) * k_k1_3);
                       }
                   } else {
                       // Leaves nodes
                       if (k_k3_l[pos - k_k3_t.size()] == 1) {
                           result.push_back(region_result(base_x + x, base_y + y, base_z + z, k_k3_l_rank(pos - k_k3_t.size())));
                           n_points++;
                       }
                   }
               } // END FOR z
           } // END FOR y
       } // END FOR x
        return n_points;
    };


    size_type get_region_k3_r(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<tuple_result> &result,
                         size_type base_x, size_type base_y, size_type base_z,
                         size_type sub_size, uint16_t level, size_type pos_children) {

       size_type pos;
       uint8_t k = k_k1;
       size_type children_size = sub_size / k;
       size_type b_x, b_y, b_z;
       size_type n_points = 0;


       // For each child that has cells that overlap with the searched region
       for (size_type x = i_pos_x/children_size; x <= (e_pos_x/children_size); x++) {
           b_x = x * children_size;
           for (size_type y = i_pos_y/children_size; y <= (e_pos_y/children_size); y++) {
               b_y = y * children_size;
               for (size_type z = i_pos_z/children_size; z <= (e_pos_z/children_size); z++) {
                   pos = pos_children + x * k * k + y * k + z; // Position of the current child
                   b_z = z * children_size;


                   if (level < (k_height-1)) {
                       // Internal nodes
                       if (k_k3_t[pos] == 1) {
                           // Continue with the search process
                           n_points += get_region_k3_r(std::max(b_x, i_pos_x) - b_x,
                                      std::max(b_y, i_pos_y) - b_y,
                                      std::max(b_z, i_pos_z) - b_z,
                                      std::min(b_x + children_size-1, e_pos_x) - b_x,
                                      std::min(b_y + children_size-1, e_pos_y) - b_y,
                                      std::min(b_z + children_size-1, e_pos_z) - b_z,
                                      result,
                                      base_x + x * children_size, base_y + y * children_size, base_z + z * children_size,
                                      children_size, level + 1, k_k3_t_rank(pos + 1) * k_k1_3);
                       }
                   } else {
                       // Leaves nodes
                       if (k_k3_l[pos - k_k3_t.size()] == 1) {
                           result.push_back(tuple_result(base_x + x, base_y + y, base_z + z));
                           n_points++;
                       }
                   }
               } // END FOR z
           } // END FOR y
       } // END FOR x
        return n_points;
    };



    size_type get_region_r(size_type i_pos_x, size_type i_pos_y, size_type i_pos_z,
                         size_type e_pos_x, size_type e_pos_y, size_type e_pos_z,
                         std::vector<tuple_result> &result,
                         size_type base_x, size_type base_y, size_type base_z,
                         size_type sub_size, uint16_t level, size_type pos_children) {

       size_type pos;
       uint8_t k = k_k1;
       size_type children_size = sub_size / k;
       size_type b_x, b_y, b_z;
       size_type n_points = 0;


       // For each child that has cells that overlap with the searched region
       for (size_type x = i_pos_x/children_size; x <= (e_pos_x/children_size); x++) {
           b_x = x * children_size;
           for (size_type y = i_pos_y/children_size; y <= (e_pos_y/children_size); y++) {
               b_y = y * children_size;
               for (size_type z = i_pos_z/children_size; z <= (e_pos_z/children_size); z++) {
                   pos = pos_children + x * k * k + y * k + z; // Position of the current child
                   b_z = z * children_size;


                   if (level < (k_height-1)) {
                       // Internal nodes
                       if (k_t[pos] == 1) {
                           // Continue with the search process
                           n_points += get_region_r(std::max(b_x, i_pos_x) - b_x,
                                      std::max(b_y, i_pos_y) - b_y,
                                      std::max(b_z, i_pos_z) - b_z,
                                      std::min(b_x + children_size-1, e_pos_x) - b_x,
                                      std::min(b_y + children_size-1, e_pos_y) - b_y,
                                      std::min(b_z + children_size-1, e_pos_z) - b_z,
                                      result,
                                      base_x + x * children_size, base_y + y * children_size, base_z + z * children_size,
                                      children_size, level + 1, k_t_rank(pos + 1) * k_k1_3);
                       }
                   } else {
                       // Leaves nodes
                       if (k_l[pos - k_t.size()] == 1) {
                           result.push_back(tuple_result(base_x + x, base_y + y, base_z + z));
                           n_points++;
                       }
                   }
               } // END FOR z
           } // END FOR y
       } // END FOR x
        return n_points;
    };

    //*******************************************************//
    //********************** FILE ***************************//
    //*******************************************************//

    //! Serialize to a stream
    /*! Serialize the k3_tree data structure
     *  \param out Outstream to write the k3_tree.
     *  \param v
     *  \param string_name
     *  \returns The number of written bytes.
     */
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                        std::string name="") const
    {
        structure_tree_node* child = structure_tree::add_child(
                v, name, util::class_name(*this));
        size_type written_bytes = 0;

        written_bytes += write_member(K3_TREE_TYPE_BASIC, out, child, "k3_tree_type");
        written_bytes += k_t.serialize(out, child, "t");
        written_bytes += k_t_rank.serialize(out, child, "t_rank");
        written_bytes += k_l.serialize(out, child, "l");
        written_bytes += k_l_rank.serialize(out, child, "l_rank");
        
        written_bytes += k_k3_t.serialize(out, child, "k3_t");
        written_bytes += k_k3_t_rank.serialize(out, child, "k3_t_rank");
        written_bytes += k_k3_l.serialize(out, child, "k3_l");
        written_bytes += k_k3_l_rank.serialize(out, child, "k3_l_rank");
        //written_bytes += k_k3_id.serialize(out, child, "k3_id");
        

        written_bytes += k_nv_psum.serialize(out, child, "nv_psum");
        written_bytes += k_id.serialize(out, child, "id");
        written_bytes += k_gr.serialize(out, child, "gr");
        
        written_bytes += k_bp.serialize(out, child, "bp");
        written_bytes += k_bp_rank.serialize(out, child, "bp_rank");
        written_bytes += k_bp_select.serialize(out, child, "bp_select");
        written_bytes += k_ts.serialize(out, child, "ts");
        written_bytes += k_pq.serialize(out, child, "pq");
        
        written_bytes += write_member(k_k1, out, child, "k1");
        written_bytes += write_member(k_k2, out, child, "k2");
        written_bytes += write_member(k_height, out, child, "k_height");
        written_bytes += write_member(k_size, out, child, "size");
        written_bytes += write_member(k_matrix_size, out, child, "matrix_size");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

	t_bv get_l() {
		return k_l;
	}
    //! Load from istream
    /*! Serialize the k3_tree from the given istream.
     *  \param istream Stream to load the k3_tree from.
     */
    void load(std::istream& in)
    {
        ushort type;
        read_member(type, in);
        k_t.load(in);
        k_t_rank.load(in);
        k_t_rank.set_vector(&k_t);
        k_l.load(in);
        k_l_rank.load(in);
        k_l_rank.set_vector(&k_l);
        
        k_k3_t.load(in);
        k_k3_t_rank.load(in);
        k_k3_t_rank.set_vector(&k_k3_t);
        k_k3_l.load(in);
        k_k3_l_rank.load(in);
        k_k3_l_rank.set_vector(&k_k3_l);
        //k_k3_id.load(in);
 
        k_nv_psum.load(in);
        k_id.load(in);
        k_gr.load(in);
        
        k_bp.load(in);
        k_bp_rank.load(in);
        k_bp_rank.set_vector(&k_bp);
        k_bp_select.load(in);
        k_bp_select.set_vector(&k_bp);
        k_ts.load(in); 
        k_pq.load(in); 

        read_member(k_k1, in);
        read_member(k_k2, in);
        read_member(k_height, in);
        read_member(k_size, in);
        read_member(k_matrix_size, in);
        k_k1_3 = pow(k_k1, 3);
    }

    void print() {
        std::cout << "k3-tree with size " << k_size << " and height " << k_height << std::endl;
        std::cout << "k1 = " << (uint)k_k1 << " | k2 = " << (uint)k_k2 << std::endl;
        std::cout << "Bitmap t (size " << k_t.size() << ") : ";
//        for (size_type p = 0; p < k_t.size(); p++) {
//            std::cout << k_t[p];
//        }
        std::cout << std::endl;
        std::cout << "Bitmap l (size " << k_l.size() << ") : ";
//        for (size_type p = 0; p < k_l.size(); p++) {
//            std::cout << k_l[p];
//        }
        std::cout << std::endl;
        return;
    }

}; // ENC CLASS k3-tree

} // END NAMESPACE sdsl

#endif // LBD_K3_TREE_OD_V1
