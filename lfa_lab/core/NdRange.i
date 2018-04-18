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

// vim: set filetype=cpp:

%exception PyNdRangeIterator::__next__ {
    try {
        $action
    }
    catch (stop_iteration) {
        PyErr_SetNone(PyExc_StopIteration);
        SWIG_fail;
    }
}

%{
class stop_iteration : public std::exception { };
%}

%inline{
    class PyNdRangeIterator {
        public:
            PyNdRangeIterator(NdRange grid)
              : m_grid(grid),
                m_pos(grid.begin())
            { }

            ArrayFi __next__() {
                if (m_pos == m_grid.end())
                    throw stop_iteration();

                ArrayFi value = *m_pos;
                ++m_pos;
                return value;
            }

        private:
            NdRange m_grid;
            NdRange::iterator m_pos;
    };
}

class NdRange {
    public:
        NdRange(ArrayFi shape = ArrayFi::Zero(0));

        ArrayFi shape() const;
        int size() const;

        int dimension() const;

        bool inRange(ArrayFi coord);
};

%extend NdRange {
    PyNdRangeIterator __iter__() {
        return PyNdRangeIterator(*$self);
    }
}


