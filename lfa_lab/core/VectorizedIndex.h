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

#ifndef VECTORIZED_INDEX_H
#define VECTORIZED_INDEX_H

#include "Common.h"

namespace lfa {

/** Given an index on a hypercube returns a corresponding vectorized one. Thus
 * an n-dimensional array can be stored in a 1-dimensional array.
 */
class VectorizedIndex
{
    public:
        class Iterator;

        VectorizedIndex(const VectorizedIndex& rhs)
          : m_cube_length(rhs.m_cube_length),
            m_start(rhs.m_start),
            m_end(rhs.m_end),
            m_offset(rhs.m_offset)
        { }


        VectorizedIndex(
                ArrayFi start = ArrayFi::Zero(1),
                ArrayFi end = ArrayFi::Zero(1) )
        {
            resize(start, end);
        }

        void resize(ArrayFi start, ArrayFi end)
        {
            m_start = start;
            m_end = end;
            m_cube_length = end - start + ArrayFi::Ones(start.size());

            m_offset = - hypercube_index(start);
        }

        /** The dimension of the index set. */
        int dimension() const { return m_cube_length.size(); }

        ArrayFi start() const { return m_start; }
        ArrayFi end() const { return m_end; }

        int operator() (ArrayFi pos) const {
            return index_of(pos);
        }

        /** Compute the index of the element at position pos. */
        int index_of(ArrayFi pos) const {
            assert(inRange(pos));
            return m_offset + hypercube_index(pos);
        }

        /** The number of elements. */
        size_t elements() {
            // compute the number of elements
            size_t s = m_cube_length[0];
            for (int d = 1; d < dimension(); ++d) {
                s *= m_cube_length[d];
            }

            return s;
        }

        bool inRange(ArrayFi pos) const {
            assert(pos.rows() == dimension());

            for (int d = 0; d < dimension(); ++d)
            {
                if (pos(d) < start()(d))
                    return false;

                if (pos(d) > end()(d))
                    return false;
            }

            return true;
        }

        bool operator== (const VectorizedIndex& rhs) const {
            if ( (m_cube_length != rhs.m_cube_length).any() )
                return false;

            if ( (m_start != rhs.m_start).any() )
                return false;

            if ( (m_end != rhs.m_end).any() )
                return false;

            return true;
        }
        bool operator!= (const VectorizedIndex& rhs) const { return !(*this == rhs); }
    private:
        ArrayFi m_cube_length;
        ArrayFi m_start, m_end;
        int m_offset;

        /** Compute the index correspoinding to a hypercube where 0 is the
         * origin. */
        int hypercube_index(ArrayFi pos) const {
            assert(pos.rows() == dimension());
            assert(dimension() > 0);

            int num = pos[ dimension() - 1 ];
            for (int d = dimension()-2; d >= 0; d--)
            {
                num = num * m_cube_length[d] + pos[d];
            }

            return num;
        }
};

class VectorizedIndex::Iterator
{
    public:
        Iterator() {}
        Iterator(const VectorizedIndex& idx)
            : m_idx(&idx)
        {
            m_pos = m_idx->start();
        }

        Iterator(const Iterator& rhs)
         :  m_pos(rhs.m_pos), m_idx(rhs.m_idx)
        {

        }

        Iterator& operator++ ()
        {
            int d = 0;
            m_pos(d)++;
            while (m_pos(d) > m_idx->end()(d) && d+1 < m_idx->dimension())
            {
                m_pos(d) = m_idx->start()(d);
                d+=1;
                m_pos(d)+=1;
            }

            return *this;
        }

        operator bool() const {
            int d = m_idx->dimension()-1;
            return (m_pos(d) <= m_idx->end()(d) );
        }

        int operator* () const {
            return value();
        }

        int value() const {
            return (*m_idx)(m_pos);
        }

        const ArrayFi& pos() { return m_pos; }
    private:
        ArrayFi m_pos;
        const VectorizedIndex* m_idx;
};

}

#endif
