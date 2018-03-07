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

#include <gtest/gtest.h>

#include <lfa_lab/core/lfa.h>

using namespace lfa;

int hash_vector(ArrayFi x) {
    int mult = 1;
    int acc = 0;
    for (int i = 0; i < x.rows(); ++i) {
        acc += x(i) * mult;
        mult *= x.rows();
    }
    return acc;
}

double float_hash(ArrayFi a, ArrayFi b, ArrayFi c, int seed = 0)
{
    vector<int> data;
    for (int i = 0; i < a.rows(); ++i) {
        data.push_back(a(i));
    }
    for (int i = 0; i < b.rows(); ++i) {
        data.push_back(b(i));
    }
    for (int i = 0; i < c.rows(); ++i) {
        data.push_back(c(i));
    }

    return float_hash(data, seed);
}

#if 1
Symbol symbol_hashed_index(int offset,
                           HarmonicClusters out_clusters,
                           HarmonicClusters in_clusters)
{
    Symbol sym(out_clusters, in_clusters);

    for (Symbol::iterator p = sym.begin(); p != sym.end(); ++p)
    {
        double d = ((p.row() == p.col()).all()) ? 2.0 : 0.0;

        *p = float_hash(p.base(), p.row(), p.col(), offset) + d;
    }

    return sym;
}
#else
// WARNING: This generates HIGHLY ill conditioned matrices. DO NOT USE!!!
Symbol symbol_hashed_index(int offset,
                           HarmonicClusters out_clusters,
                           HarmonicClusters in_clusters)
{
    Symbol sym(out_clusters, in_clusters);

    for (Symbol::iterator p = sym.begin(); p != sym.end(); ++p)
    {
        *p = hash_vector(p.base()) * 10000
           + (1+hash_vector(p.row())) * (1+hash_vector(p.col())) * 51
           + offset;
    }

    return sym;
}
#endif

int rotate(int x, int end_val)
{
    return (x + 4) % end_val;
}

ArrayFi rotate(ArrayFi x, ArrayFi end_val)
{
    ArrayFi result(x.rows());
    for (int i = 0; i < result.rows(); ++i)
    {
        result(i) = rotate(x(i), end_val(i));
    }
    return result;
}

double rotated_id(ArrayFi i, ArrayFi j, HarmonicClusters clusters)
{
    ArrayFi pi_i = rotate(i, clusters.clusterShape());

    if ( (pi_i == j).all() ) {
        return 1.0;
    } else {
        return 0.0;
    }
}

TEST(SystemSymbol, norm_identity)
{
    HarmonicClusters clusters(Array2i(2,3), Array2i(3,4));
    SystemSymbol S = SystemSymbol::Identity(2, 2, clusters, clusters);

    double expect = sqrt(clusters.size() * S.rows());

    EXPECT_EQ(expect, S.system_norm());
}


TEST(SystemSymbol, MatMultPermut)
{
    // Permutation matrices are orthogonal, hence we check P * P_T = id

    HarmonicClusters clusters(Array2i(2,3), Array2i(2,5));

    SystemSymbol P(2, 2, clusters, clusters);

    for (int i = 0; i < P.rows(); ++i) {
        for (int j = 0; j < P.cols(); ++j) {
            for (Symbol::iterator p = P(i,j).begin(); p != P(i,j).end(); ++p)
            {
                if (rotate(i, P.rows()) == j) {
                    *p = rotated_id(p.row(), p.col(), clusters);
                } else {
                    *p = 0;
                }
            }
        }
    }

    SystemSymbol P_T(2, 2, clusters, clusters);
    for (int i = 0; i < P_T.rows(); ++i) {
        for (int j = 0; j < P_T.cols(); ++j) {
            for (Symbol::iterator p = P_T(i,j).begin(); p != P_T(i,j).end(); ++p)
            {
                if (rotate(j, P.cols()) == i) {
                    *p = rotated_id(p.col(), p.row(), clusters);
                } else {
                    *p = 0;
                }
            }
        }
    }

    SystemSymbol result = P_T * P;
    SystemSymbol I = SystemSymbol::Identity(P.rows(), P.cols(), clusters, clusters);

    EXPECT_LE( (I - result).system_norm(), 1e-12 );

    EXPECT_LE( (P_T - P.inverse()).system_norm(), 1e-12);
}

TEST(SystemSymbol, ConstDiagBlock)
{
    HarmonicClusters clusters(Array2i(2,3), Array2i(2,1));
    SystemSymbol sym(2, 2, clusters, clusters);

    // | 1 2 |   | -2       1 |
    // | 3 4 | * |  3/2  -1/2 | = id

    Symbol I = Symbol::Identity(clusters, clusters);
    sym(0, 0) = 1 * I;
    sym(0, 1) = 2 * I;
    sym(1, 0) = 3 * I;
    sym(1, 1) = 4 * I;

    // Test inverse
    SystemSymbol expected(2, 2, clusters, clusters);
    expected(0, 0) = -2 * I;
    expected(0, 1) =  1 * I;
    expected(1, 0) =  3.0/2 * I;
    expected(1, 1) = -1.0/2 * I;

    EXPECT_LE( (expected - sym.inverse()).system_norm(), 1e-12 );


    // Test Multiplication
    SystemSymbol sym2(2, 2, clusters, clusters);
    sym2(0, 0) = 5 * I;
    sym2(0, 1) = 6 * I;
    sym2(1, 0) = 7 * I;
    sym2(1, 1) = 8 * I;

    expected(0, 0) = 19 * I;
    expected(0, 1) = 22 * I;
    expected(1, 0) = 43 * I;
    expected(1, 1) = 50 * I;

    EXPECT_LE( (expected - sym * sym2).system_norm(), 1e-12 );
}

TEST(SystemSymbol, inverse)
{
    HarmonicClusters clusters(Array2i(2,3), Array2i(2,1));
    SystemSymbol sym(2, 2, clusters, clusters);

    for (int i = 0; i < sym.rows(); ++i) {
        for (int j = 0; j < sym.cols(); ++j) {
            sym(i, j) = symbol_hashed_index(2*i+j, clusters, clusters);
        }
    }

    SystemSymbol sym_inv = sym.inverse();

    SystemSymbol result = sym_inv * sym;
    SystemSymbol I = SystemSymbol::Identity(sym.rows(), sym.cols(), clusters, clusters);

    EXPECT_LE( (I - result).system_norm(), 1e-12);
}


