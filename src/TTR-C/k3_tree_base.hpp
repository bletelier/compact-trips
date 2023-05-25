/*  
 * Created by Fernando Silva on 13/06/18.
 *
 * Copyright (C) 2018-current-year, Fernando Silva, all rights reserved.
 *
 * 
 * Author's contact: Fernando Silva  <fernando.silva@udc.es>
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

#ifndef K3_TREE_SDSL_K3_TREE_BASE_H
#define K3_TREE_SDSL_K3_TREE_BASE_H

#include <cstdint>
#include <queue>
#include <tuple>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

//! Namespace for the succint data structure library
namespace sdsl {

    static const ushort K3_TREE_TYPE_BASIC = 0;
    static const ushort K3_TREE_TYPE_POINTS = 1;
    static const ushort K3_TREE_TYPE_LEVEL = 2;
    static const ushort K3_TREE_TYPE_LIDAR_POINTS = 10;


template <typename t_real_pos_type=int_vector<>::size_type>
class k3_tree_base
{

    public:
        typedef int_vector<>::size_type size_type;
        typedef int_vector<>::size_type pos_type;
        typedef std::tuple<size_type, size_type, size_type> tuple_result;
        typedef std::tuple<pos_type, pos_type, pos_type> point_type;


    protected:
    pos_type   k_size;


    public:
    const size_type &ssize = k_size;

    public:
        //*******************************************************//
        //******************* CONSTRUCTOR ***********************//
        //*******************************************************//
        k3_tree_base() {};

        k3_tree_base(const k3_tree_base& tr)
        {
            *this = tr;
        }

        k3_tree_base(k3_tree_base&& tr)
        {
            *this = std::move(tr);
        }

        //*******************************************************//
        //*************** BASIC OPERATIONS **********************//
        //*******************************************************//

        //! Move assignment operator
        k3_tree_base& operator=(k3_tree_base&& tr)
        {
            if (this != &tr) {

            }
            return *this;
        }

        //! Assignment operator
        k3_tree_base& operator=(const k3_tree_base& tr)
        {
            if (this != &tr) {

            }
            return *this;
        }

        //! Swap operator
        void swap(k3_tree_base& tr)
        {
            if (this != &tr) {

            }
        }

        //! Equal operator
        bool operator==(const k3_tree_base& tr) const
        {
            return true;
        }

        //*******************************************************//
        //********************** GETTERS ************************//
        //*******************************************************//

        virtual t_real_pos_type get_min_size_x() {return 0;}
        virtual t_real_pos_type get_min_size_y() {return 0;}
        virtual t_real_pos_type get_min_size_z() {return 0;}
        virtual t_real_pos_type get_max_size_x() {return k_size;}
        virtual t_real_pos_type get_max_size_y() {return k_size;}
        virtual t_real_pos_type get_max_size_z() {return k_size;}

        //*******************************************************//
        //********************** QUERIES ************************//
        //*******************************************************//

        virtual bool get(t_real_pos_type pos_x, t_real_pos_type pos_y, t_real_pos_type pos_z) = 0;
        virtual size_type get_region(t_real_pos_type i_pos_x, t_real_pos_type i_pos_y, t_real_pos_type i_pos_z,
                                     t_real_pos_type e_pos_x, t_real_pos_type e_pos_y, t_real_pos_type e_pos_z,
                         std::vector<tuple_result> &result) = 0;
        virtual size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                                    std::string name="") const = 0;
        virtual void load(std::istream& in) = 0;
        virtual void print() = 0;

        //*******************************************************//
        //*********************** TEST **************************//
        //*******************************************************//
        virtual bool check_values(std::string data_file_name){
            std::ifstream values_file(data_file_name);
            assert(values_file.is_open() && values_file.good());

            bool res = check_values(values_file);
            values_file.close();
            return res;
        }

        virtual bool check_values(std::istream& data_file) {

            // Read file
            size_type pos_x, pos_y, pos_z;
            while (!data_file.eof() && data_file.good()) {

                // Get position (x, y, z)
                read_member(pos_x, data_file);
                read_member(pos_y, data_file);
                read_member(pos_z, data_file);

                if (!get(pos_x, pos_y, pos_z)) {
    #ifndef NDEBUG
                    std::cout << "Failed point (" << pos_x << ", " << pos_y << ", " << pos_z << ") " << std::endl;
    #endif
                    return false;
                }

            }

            return true;
        }

        virtual bool check_values(std::vector<point_type>& points) {

            // Read file
            pos_type pos_x, pos_y, pos_z;
            for (size_type p = 0; p < points.size(); p++) {
                pos_x = std::get<0>(points[p]);
                pos_y = std::get<1>(points[p]);
                pos_z = std::get<2>(points[p]);

                if (!get(pos_x, pos_y, pos_z)) {
    #ifndef NDEBUG
                    std::cout << "Failed point (" << pos_x << ", " << pos_y << ", " << pos_z << ") " << std::endl;
    #endif
                    return false;
                }
            }
            return true;
        }

}; // END CLASS k3_tree_base

} // END NAMESPACE sdsl

#endif //K3_TREE_SDSL_K3_TREE_BASE_H
