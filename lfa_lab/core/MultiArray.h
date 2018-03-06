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

#ifndef LFA_MULTI_ARRAY_H
#define LFA_MULTI_ARRAY_H

#include "Common.h"
#include "VectorizedIndex.h"

namespace lfa {

  /** Multi-dimensional array. */
  template <typename T>
  class MultiArray
  {
    public:
      template <typename ArrayT, typename ScalarRef> class IteratorBase;
      typedef IteratorBase<MultiArray, T&> Iterator;
      typedef IteratorBase<const MultiArray, T> ConstIterator;
      typedef T ValueType;

      MultiArray(const MultiArray& rhs)
        : index_of(rhs.index_of),
          m_elements(rhs.m_elements)
      { }

      MultiArray(ArrayFi start = ArrayFi::Zero(1),
                 ArrayFi end = ArrayFi::Zero(1) );

      /** Sets the array to the given size and dimension.
       * The contents of the array is undefined. */
      void resize(ArrayFi start, ArrayFi end);

      /** Set all elements of the array to a certain value. */
      void fill(const T& value);

      /** Access the element at position pos. */
      inline T& operator() (const ArrayFi& pos) {
        return m_elements[ index_of(pos) ];
      }

      T operator() (const ArrayFi& pos) const {
        return m_elements[ index_of(pos) ];
      }

      /** The dimension of the array. */
      int dimension() const { return index_of.dimension(); }

      ArrayFi startIndex() const { return index_of.start(); }
      ArrayFi endIndex() const { return index_of.end(); }

      ArrayFi shape() const {
        return endIndex() - startIndex() + ArrayFi::Ones(dimension());
      }

      template <typename Iter>
        void assign(Iter begin, Iter end) {
#ifndef NDEBUG
          vector<double>::iterator it =
#endif
            copy(begin, end, m_elements.begin());

          assert(it == m_elements.end());
        }

      bool operator== (const MultiArray& rhs) const {
        if (index_of != rhs.index_of)
          return false;
        return m_elements == rhs.m_elements;
      }

      bool inIndexRange(ArrayFi idx) { return index_of.inRange(idx); }

      vector<ValueType>& linear_access() { return m_elements; }
    private:
      // computes the index in the array
      VectorizedIndex index_of;
      vector<ValueType> m_elements;
  };

  template <typename T>
  MultiArray<T> :: MultiArray(ArrayFi start, ArrayFi end)
  {
    resize(start, end);
  }

  template <typename T>
  void MultiArray<T> :: resize(ArrayFi start, ArrayFi end)
  {
    index_of.resize(start, end);
    m_elements.resize( index_of.elements() );
  }

  template <typename T>
  void MultiArray<T> :: fill(const T& value)
  {
    std::fill(m_elements.begin(), m_elements.end(), value);
  }

  template <typename T>
  template <typename ArrayT, typename ValueRef>
  class MultiArray<T>::IteratorBase
  {
    public:
      IteratorBase() {}
      IteratorBase(ArrayT& s)
        : m_array(&s), m_it(s.index_of)
      { }


      IteratorBase& operator++ ()
      {
        ++m_it;
        return *this;
      }

      operator bool() {
        return (bool)m_it;
      }

      ValueRef operator* () {
        return value();
      }

      ValueRef value() {
        return array()( m_it.pos() );
      }

      bool beforeCenter() {
        ArrayFi pos = m_it.pos();

        for (int d = array().dimension()-1; d > 0; --d) {
          if (pos(d) > 0) {
            return false;
          }
          if (pos(d) < 0) {
            return true;
          }
        }

        // pos(1) ... pos(d-1) is zero
        return (pos(0) < 0);
      }

      bool atCenter() {
        ArrayFi pos = m_it.pos();
        for (int d = 0; d < array().dimension(); ++d)
        {
          if (pos(d) != 0)
            return false;
        }
        return true;
      }

      bool behindCenter() {
        return !(beforeCenter() || atCenter());
      }

      ArrayT& array() { return *m_array; }

      const ArrayFi& pos() { return m_it.pos(); }
    private:
      ArrayT* m_array;
      VectorizedIndex::Iterator m_it;
  };


  template <typename T>
  std::ostream& operator<< (std::ostream& os, const MultiArray<T>& a)
  {
    for (typename MultiArray<T>::ConstIterator it(a); it; )
    {
      os << it.value() << "  ";
      ++it;
      for (int i = 0; i < it.pos().rows() && it.pos()(i) == a.startIndex()(i); ++i)
        os << std::endl;
    }

    return os;
  }

}

#endif
