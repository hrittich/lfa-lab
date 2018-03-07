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

%include "MatrixContainer.i"

class SystemSymbol {
  public:
    explicit SystemSymbol(
        int rows = 0,
        int cols = 0,
        HarmonicClusters output_clusters = HarmonicClusters(),
        HarmonicClusters input_clusters = HarmonicClusters());

    static SystemSymbol Identity(int rows,
        int cols,
        HarmonicClusters output_clusters,
        HarmonicClusters input_clusters);

    void resize(int rows, int cols);

    SystemSymbol operator* (const SystemSymbol& other) const;
    SystemSymbol operator+ (const SystemSymbol& other) const;
    //friend SystemSymbol operator* (double scalar, const SystemSymbol& other);
    SystemSymbol operator- (const SystemSymbol& other) const;

    SystemSymbol inverse() const;

    int rows() const { return m_rows; }
    int cols() const { return m_cols; }

    double spectral_radius() const;
    double spectral_norm() const;

    %extend {
      SystemSymbol __rmul__(double scalar) {
        return scalar * (*$self);
      }

      void __setitem__(ArrayFi index, const Symbol& value) {
        if (index.rows() != 2)
          throw logic_error("Invalid number of indices.");

        (*$self)(index[0], index[1]) = value;
      }

      Symbol __getitem__(ArrayFi index) {
        if (index.rows() != 2)
          throw logic_error("Invalid number of indices.");

        return (*$self)(index[0], index[1]);
      }
    }
};

%template(SymbolMatrix) MatrixContainer<Symbol>;

SystemSymbol combine_symbols_into_system(
  const SystemSymbolProperties& properties,
  const MatrixContainer<Symbol>& symbols); 


