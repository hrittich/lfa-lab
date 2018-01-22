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

%include "Symbol.i"
%include "NdRange.i"

template<typename T>
class NdArray
{
    public:
        NdArray(ArrayFi shape);

        ArrayFi shape() const;

        bool index_in_range(const ArrayFi& pos) const;

        NdRange indices() const;
};

%extend NdArray {
    T __getitem__(ArrayFi pos) {
        if (!$self->index_in_range(pos)) {
            throw out_of_range("Index out of range.");
        }
        return (*$self)(pos);
    }

    void __setitem__(ArrayFi pos, T& value) {
        if (!$self->index_in_range(pos)) {
            throw out_of_range("Index out of range.");
        }

        (*$self)(pos) = value;
    }
}

%template(SymbolNdArray) NdArray<Symbol>;

