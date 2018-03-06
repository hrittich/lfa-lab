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

#ifndef LFA_CART_ITERATOR_H
#define LFA_CART_ITERATOR_H

#include "Common.h"
#include <iterator>

namespace lfa {

  /** Cartesian iterator. Iterates over all combinations of the first and the
   * second iterator.
   */
  template <typename first_iter, typename second_iter>
  class CartIterator
    : public std::iterator<std::forward_iterator_tag,
                           std::pair<first_iter, second_iter> >
  {
    public:

      CartIterator() {
      }

      CartIterator(first_iter first_begin,
                   first_iter first_end,
                   second_iter second_begin,
                   second_iter second_end,
                   bool out_of_range = false)
        : m_first_begin(first_begin),
          m_first_end(first_end),
          m_iter(first_begin, second_begin)
      {
        if (out_of_range) {
          m_iter.second = second_end;
        }
      }

      /** Out of stream iterator. */
      CartIterator(second_iter second_end)
      {

      }

      CartIterator& operator++()
      {
        ++m_iter.first;
        if (m_iter.first == m_first_end)
        {
          m_iter.first = m_first_begin;

          ++m_iter.second;
        }
        return *this;
      }

      CartIterator operator++(int)
      {
        CartIterator old = *this;
        operator++ ();
        return old;
      }

      std::pair<first_iter, second_iter>&
        operator* ()
        {
          return m_iter;
        }

      std::pair<first_iter, second_iter>*
        operator-> ()
        {
          return &m_iter;
        }

      bool operator== (const CartIterator& other) const {
        return m_iter == other.m_iter;
      }

      bool operator!= (const CartIterator& other) const {
        return !(*this == other);
      }


    private:
      first_iter m_first_begin, m_first_end;

      std::pair<first_iter, second_iter> m_iter;

  };

  template <typename first_iter, typename second_iter>
  CartIterator<first_iter, second_iter>
    make_cart_iterator(first_iter first_begin,
                       first_iter first_end,
                       second_iter second_begin,
                       second_iter second_end)
  {
   return CartIterator<first_iter, second_iter> (first_begin, first_end,
       second_begin, second_end);
  }

}

#endif
