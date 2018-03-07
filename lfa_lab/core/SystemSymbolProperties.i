/*
  vim: set filetype=cpp:

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

%import "MatrixContainer.i"

class SystemSymbolProperties {
  public:
    SystemSymbolProperties(int rows,
                           int cols,
                           FoProperties element_properties);

    SystemSymbolProperties operator* (const SystemSymbolProperties& other) const;
    SystemSymbolProperties operator+ (const SystemSymbolProperties& other) const;

    SystemSymbolProperties inverse() const;
    SystemSymbolProperties adjoint() const;

    %extend {
      SystemSymbolProperties __rmul__ (double scalar) {
        return scalar * (*$self);
      }
    }

    int dimension() const;
    ArrayFi adjustResolution(ArrayFi desired_resolution) const;

    Grid outputGrid() const;
    Grid inputGrid() const;

    int rows() const;
    int cols() const;

    FoProperties element_properties();
};

SystemSymbolProperties properties_of_symbol_system(
  MatrixContainer<FoProperties> properties);

%template(FoPropertiesMatrix) MatrixContainer<FoProperties>;

