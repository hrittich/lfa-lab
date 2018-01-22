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

#ifndef LFA_FREQUENCY_SPACE_H
#define LFA_FREQUENCY_SPACE_H

#include "Common.h"
using namespace Eigen;

namespace lfa {

class FrequencySpace {
    public:
        FrequencySpace(ArrayFi shape, ArrayFd step_size)
         :  m_shape(shape), m_step_size(step_size)
        { }

        ArrayFi shape() const { return m_shape; }
        ArrayFd step_size() const { return m_step_size; }

        int modes() { return m_shape.prod(); }

        bool operator== (const FrequencySpace& rhs) const {
            if ( (shape() != rhs.shape()).any())
                return false;

            return (step_size() == rhs.step_size()).all();
        }

        bool operator!=(const FrequencySpace& rhs) const {
            return !(*this == rhs);
        }

    private:
        ArrayFi m_shape;
        ArrayFd m_step_size;
};

}

#endif
