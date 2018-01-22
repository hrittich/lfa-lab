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

%nodefaultctor FoProperties;
class FoProperties {
    public:
        FoProperties(
                SplitFrequencyDomain output_domain,
                SplitFrequencyDomain input_domain);

        Grid outputGrid();
        Grid inputGrid();

        FoProperties inverse();
        FoProperties adjoint();

        ArrayFi adjustResolution(ArrayFi desired_resolution);

        int dimension() const;

        FoProperties expand(ArrayFi factor) const;
};

%extend FoProperties {
    FoProperties __add__(const FoProperties& other) {
        return *$self + other;
    }

    FoProperties __sub__(const FoProperties& other) {
        return *$self - other;
    }

    FoProperties __mul__(const FoProperties& other) {
        return (*$self) * other;
    }

    FoProperties __rmul__(std::complex<double> scalar) {
        return scalar * (*$self);
    }
}


