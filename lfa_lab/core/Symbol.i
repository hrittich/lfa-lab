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

%feature("autodoc", "Approximation of a matrix symbols.") Symbol;
%feature("autodoc",
"The norms of the rows of the symbol as an :math:`n`-D array.") Symbol::row_norms;
%feature("autodoc",
"The norms of the columns of the symbol as :math:`n`-D array.") Symbol::col_norms;
%feature("autodoc",
"The norms of the rows of a 1D symbol as a vector.") Symbol::row_norms_1d;
%feature("autodoc",
"The norms of the columns of a 1D symbol as a vector.") Symbol::col_norms_1d;
%feature("autodoc",
"The norms of the rows of a 2D symbol as a matrix.") Symbol::row_norms_2d;
%feature("autodoc",
"The norms of the columns of a 2D symbol as a matrix.") Symbol::col_norms_2d;
%feature("autodoc", "The (spectral) norm of the symbol.") Symbol::spectral_norm;
%feature("autodoc", "The spectral radius of the symbol.") Symbol::spectral_radius;
%feature("autodoc", "The eigenvalues of the symbol as a vector.") Symbol::eigenvalues;
%feature("autodoc", "The dimension of the symbol.") Symbol::dimension;
class Symbol {
    public:
        static Symbol Identity(Grid grid, SamplingProperties conf);
        static Symbol Zero(Grid, SamplingProperties conf);

        NdArray<double> row_norms() const;
        NdArray<double> col_norms() const;

        VectorXcd row_norms_1d() const;
        VectorXcd col_norms_1d() const;
        MatrixXcd row_norms_2d() const;
        MatrixXcd col_norms_2d() const;

        Symbol expand(ArrayFi factor) const;

        Symbol inverse() const;
        Symbol adjoint() const;

        double norm() const;
        double spectral_radius() const;
        double spectral_norm() const;

        VectorXcd eigenvalues() const;

        int dimension() const;

        %pythoncode {
            def symbol(self):
                return self
        }

        const HarmonicClusters& outputClusters() const;
        const HarmonicClusters& inputClusters() const;
};

%extend Symbol {
    Symbol __add__(const Symbol& other) {
        return *$self + other;
    }

    Symbol __sub__(const Symbol& other) {
        return *$self - other;
    }

    Symbol __mul__(const Symbol& other) {
        return (*$self) * other;
    }

    Symbol __rmul__(std::complex<double> scalar) {
        return scalar * (*$self);
    }

    std::string __str__() {
        stringstream s;
        s << *$self;
        return s.str();
    }

    int output_coupling_size() {
        return $self->outputClusters().clusterSize();
    }

    int input_coupling_size() {
        return $self->inputClusters().clusterSize();
    }

    ArrayFi output_coupling_shape() {
        return $self->outputClusters().clusterShape();
    }

    ArrayFi input_coupling_shape() {
        return $self->inputClusters().clusterShape();
    }

    ArrayFi output_shape() {
        return $self->outputClusters().shape();
    }

    ArrayFi input_shape() {
        return $self->inputClusters().shape();
    }
}

