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

// There is an adapter class that adds more functionality
// (see lfa/stencil.py).
%rename(_SparseStencil) SparseStencil;
class SparseStencil {
    public:
        SparseStencil();
        SparseStencil(const DenseStencil& other);

        int nonZeros();

        void append(ArrayFi offset, std::complex<double> value);

        %extend {
        void append(ArrayFi offset, double value) {
           $self->append(offset, complex<double> (value));
        }
        }

        int dimension() const;
};

%extend SparseStencil {
    ArrayFi _get_offset_at(int p) { return (*$self)[p].offset; }
    std::complex<double> _get_value_at(int p) { return (*$self)[p].value; }
}

