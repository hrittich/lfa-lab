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

class NdRange {
  public:
    NdRange(ArrayFi shape = ArrayFi::Zero(0));

    ArrayFi shape() const;
    int size() const;

    int dimension() const;

    bool inRange(ArrayFi coord);

    NdRangeIterator begin();
    NdRangeIterator end();
};

class NdRangeIterator : public std::iterator<std::input_iterator_tag, ArrayFi>
{
  public:
    bool operator== (const NdRangeIterator& rhs) const;
    bool operator!= (const NdRangeIterator& rhs) const;
};

%extend NdRangeIterator {
  void advance() {
    ++(*$self);
  }
  ArrayFi value() {
    return **$self;
  }
}

%pythoncode %{

def NdRange__iter__(self):
  iter = self.begin()
  while iter != self.end():
    yield iter.value()
    iter.advance()

setattr(NdRange, r'__iter__', NdRange__iter__)

%}


