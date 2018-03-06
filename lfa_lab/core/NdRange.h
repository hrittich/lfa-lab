/*
  LFA Lab - Library to simplify local Fourier analysis.
  Copyright (C) 2018  Hannah Rittich

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef LFA_NDRANGE_H
#define LFA_NDRANGE_H

#include "Common.h"
#include <iterator>

#include <Eigen/Core>

namespace lfa {

/** Cartesian grid graph.
 *
 * Represents a set of elements of the cartesian product
 * [ 0, 1, ..., n_1 ] x [0, 1, ..., n_2] x ... x [0, ... n_d]
 * */
class NdRange {
    public:
        class iterator;

        NdRange(ArrayFi shape = ArrayFi::Zero(0))
         :  m_shape(shape)
        {}

        ArrayFi shape() const { return m_shape; }
        int size() const { return m_shape.prod(); }

        int dimension() const { return m_shape.rows(); }

        iterator begin();
        iterator end();

        int indexOf(ArrayFi coord) const {
            assert(inRange(coord));

            int num = coord[ dimension() - 1 ];
            for (int d = dimension()-2; d >= 0; d--)
            {
                num = num * m_shape[d] + coord[d];
            }

            return num;
        }

        bool inRange(ArrayFi coord) const {
            assert(coord.rows() == dimension());

            for (int d = 0; d < dimension(); ++d)
            {
                if (coord(d) < 0)
                    return false;

                if (coord(d) >= m_shape(d))
                    return false;
            }

            return true;
        }

        /** The number of elements. */
        size_t size() {
            // compute the number of elements
            if (dimension() < 1)
                return 0;

            size_t s = m_shape[0];
            for (int d = 1; d < dimension(); ++d) {
                s *= m_shape[d];
            }

            return s;
        }

    private:
        ArrayFi m_shape;
};

class NdRange::iterator : public std::iterator<std::input_iterator_tag, ArrayFi>
{
    public:
        iterator()
         :  m_out_of_range(true)
        {}

        iterator(const iterator& rhs)
         :  m_out_of_range(rhs.m_out_of_range),
            m_pos(rhs.m_pos),
            m_grid(rhs.m_grid)
        {
        }

        /** Start at position 0. */
        iterator(const NdRange& grid)
         :  m_out_of_range(false),
            m_pos(ArrayFi::Zero(grid.dimension())),
            m_grid(grid)
        {
        }

        iterator& operator++ ()
        {
            int d = 0;
            m_pos(d)++;
            while (m_pos(d) >= m_grid.m_shape(d) && d+1 < dimension())
            {
                m_pos(d) = 0;
                d+=1;
                m_pos(d)+=1;
            }

            // if m_pos(d) is still >= m_grid->m_shape(d) we are out of range
            if (m_pos(d) >= m_grid.m_shape(d)) {
                m_out_of_range = true;
            }

            return *this;
        }

        bool operator== (const iterator& rhs) const {

            if (m_out_of_range != rhs.m_out_of_range) {
                return false;
            }

            if (m_out_of_range) {
                // both are out_of_range
                return true;
            }

            // non are out of range thus compare position
            return (m_pos == rhs.m_pos).all();
        }

        bool operator!= (const iterator& rhs) const {
            return !(*this == rhs);
        }

        const ArrayFi& operator* () {
            return m_pos;
        }

        const ArrayFi* operator-> () {
            return &m_pos;
        }

        int dimension() { return m_grid.dimension(); }

    private:
        bool m_out_of_range;
        ArrayFi m_pos;
        NdRange m_grid;
};

inline NdRange::iterator NdRange::begin()
{
    return iterator(*this);
}

inline NdRange::iterator NdRange::end()
{
    return iterator();
}

}

#endif
