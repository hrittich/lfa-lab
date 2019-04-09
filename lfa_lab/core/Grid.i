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

%feature("autodoc", "Description of an infinite grid.

:param arg1: The number of dimensions or a tuple with spacings per dimension.
:type arg1: int or tuple
:param arg2: The step-sizes per dimension, given as a tuple. (optional)
") Grid;
%feature("autodoc", "The step size.") Grid::step_size;
%feature("autodoc", "The number of dimensions of the grid.") Grid::dimension;
%feature("autodoc", "The coarsening range of `other` with respect to `self`.")
  Grid::factor;
class Grid {
    public:
        Grid(int dimension, ArrayFd step_size = ArrayFd::Zero(0));
        Grid(ArrayFi spacing, ArrayFd step_size = ArrayFd::Zero(0));

        Grid coarse(ArrayFi factor);
        ArrayFi coarsening_factor(const Grid& other);

        ArrayFd step_size();

        int dimension() const;
};

%extend Grid {
    std::string __str__() {
        std::stringstream ss;
        ss << *$self;
        return ss.str();
    }

    bool __eq__(const Grid& other) {
      return (*$self) == other;
    }
}


