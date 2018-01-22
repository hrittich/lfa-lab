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

#include "FrequencyIndices.h"
#include "HarmonicClusters.h"
#include "NdRange.h"
#include "Symbol.h"
#include "Math.h"

using namespace lfa;

TEST(Harmonics, NdRange)
{
    using namespace std;

    NdRange g(Array2i(12,13));

    int p = 0;

    for (NdRange::iterator it = g.begin(); it != g.end(); ++it, ++p) {
        EXPECT_EQ(g.indexOf(*it), p);
        // cout << "No " << p << endl << *it << endl << endl;
    }
}

TEST(Harmonics, DiscreteSet)
{
    FrequencyIndices H(Array2i(6,4));

    vector<ArrayFi> es = H.elements();
    //copy(es.begin(), es.end(), ostream_iterator<ArrayFi>(cout, "\n\n"));

    int x[] = { 0, 2, 4, 0, 2, 4,
                1, 3, 5, 1, 3, 5,
                0, 2, 4, 0, 2, 4,
                1, 3, 5, 1, 3, 5, 99 };
    int y[] = { 0, 0, 0, 2, 2, 2,
                0, 0, 0, 2, 2, 2,
                1, 1, 1, 3, 3, 3,
                1, 1, 1, 3, 3, 3, 99 };

    int cor_inds[] = {
        0, 2, 4, 12, 14, 16,
        1, 3, 5, 13, 15, 17,
        6, 8, 10, 18, 20, 22,
        7, 9, 11, 19, 21, 23, 9999 };

    vector<HarmonicIndices> Hs = H.split(Array2i(2,2));
    int k = 0;
    for (size_t i = 0; i < Hs.size(); ++i) {
        // cout << "====== New Subset ======" << endl;

        vector<ArrayFi> es = Hs[i].elements();
        // copy(es.begin(), es.end(), ostream_iterator<ArrayFi>(cout, "\n\n"));

        vector<int> inds = Hs[i].indices();
        //copy(inds.begin(), inds.end(), ostream_iterator<int>(cout, ", "));
        //cout << endl;

        for (size_t j = 0; j < es.size(); ++j) {
            EXPECT_EQ(es[j](0), x[k]);
            EXPECT_EQ(es[j](1), y[k]);
            EXPECT_EQ(inds[j], cor_inds[k]);

            ++k;
        }
    }

    EXPECT_EQ(k, 24);
}


int harmonic_clusters_test_b1[] =
{ 0, 1,  2, 3,  4, 5,  0, 1,  2, 3,  4, 5 };
int harmonic_clusters_test_c1[] =
{ 0, 0,  0, 0,  0, 0,  1, 1,  1, 1,  1, 1 };

int harmonic_clusters_test_b2[] =
{ 0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3 };
int harmonic_clusters_test_c2[] =
{ 0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2 };

int harmonic_clusters_test_b3[] =
{ 0, 1,  0, 1,  0, 1,  0, 1,  0, 1,  0, 1};
int harmonic_clusters_test_c3[] =
{ 0, 0,  1, 1,  2, 2,  3, 3,  4, 4,  5, 5};


TEST(Harmonics, HarmonicClusters)
{
    ArrayFi H(3), B(3);

    H << 2, 3, 6;
    B << 6, 4, 2;
    HarmonicClusters clusters(B, H);

    NdRange g(H * B);

    // check that global -> local -> global is the identity
    for (NdRange::iterator p = g.begin(); p != g.end(); ++p)
    {
        ArrayFi b = clusters.baseIndex(*p);
        ArrayFi c = clusters.clusterIndex(*p);

        EXPECT_TRUE((clusters.globalIndex(b,c) == *p).all());
    }

    // check with expected data
    for (NdRange::iterator p = g.begin(); p != g.end(); ++p)
    {
        ArrayFi b = clusters.baseIndex(*p);
        ArrayFi c = clusters.clusterIndex(*p);

        EXPECT_EQ(b(0), harmonic_clusters_test_b1[(*p)(0)]);
        EXPECT_EQ(c(0), harmonic_clusters_test_c1[(*p)(0)]);
        EXPECT_EQ(b(1), harmonic_clusters_test_b2[(*p)(1)]);
        EXPECT_EQ(c(1), harmonic_clusters_test_c2[(*p)(1)]);
        EXPECT_EQ(b(2), harmonic_clusters_test_b3[(*p)(2)]);
        EXPECT_EQ(c(2), harmonic_clusters_test_c3[(*p)(2)]);
    }
}



TEST(Harmonics, ClusterSanity)
{
    ArrayFi H1(2), B1(2), H2(2), B2(2);

    H1 << 3, 3;
    B1 << 2, 2;

    H2 << 6, 6;
    B2 << 1, 1;

    HarmonicClusters c1(B1, H1);
    HarmonicClusters c2(B2, H2);

    EXPECT_FALSE(c1.isCompatibleTo(c2));

    EXPECT_EQ(c1.mergeCluster(B1), c2);

    // shape has to be invariant
    EXPECT_TRUE( (c1.mergeCluster(B1).shape() == c1.shape()).all() );
}


TEST(Harmonics, DiagSymbolTest)
{
    ArrayFi H(2), B(2);
    ArrayFi F1(2), F2(2);

    B << 6, 15;
    H << 1, 1;
    F1 << 2, 3;
    F2 << 3, 5;

    HarmonicClusters clusters(B, H);

    Symbol L(clusters, clusters);

    NdRange g(B);

    int n = 0;
    for (NdRange::iterator p = g.begin(); p != g.end(); ++p)
    {
        L.ref(*p, Array2i(0,0), Array2i(0,0)) = ++n;
    }

    // expand in one step
    Symbol E1 = L.expand(B);
    EXPECT_EQ((L.matrix().full() - E1.matrix().full()).norm(), 0);

    // expand in two steps
    Symbol E2 = L.expand(F1).expand(F2);
    EXPECT_EQ((L.matrix().full() - E2.matrix().full()).norm(), 0);
}

TEST(Harmonics, Symbol2Harmonics)
{
    ArrayFi H(2), B(2);
    ArrayFi F(2);

    B << 4, 4;
    H << 2, 2;
    F << 2, 2;

    NdRange base(B);
    NdRange cluster(H);

    HarmonicClusters harmonics(B, H);

    Symbol L(harmonics, harmonics);
    for (NdRange::iterator pb = base.begin();
            pb != base.end(); ++pb)
    {
        for (NdRange::iterator prow = cluster.begin();
                prow != cluster.end(); ++prow)
        {
            for (NdRange::iterator pcol = cluster.begin();
                    pcol != cluster.end(); ++pcol)
            {
                ArrayFi global_row = harmonics.globalIndex(*pb, *prow);
                ArrayFi global_col = harmonics.globalIndex(*pb, *pcol);

                // hash the global coordinates
                L.ref(*pb, *prow, *pcol) =
                    global_row(0)         + global_row(1) * 100 +
                    global_col(0) * 10000 + global_col(1) * 1000000;
            }
        }
    }

    // do the expansion
    L = L.expand(F);

    base = NdRange(B / F);
    cluster = NdRange(H * F);
    HarmonicClusters new_harmonics(B / F, H * F);

    for (NdRange::iterator pb = base.begin();
            pb != base.end(); ++pb)
    {
        for (NdRange::iterator prow = cluster.begin();
                prow != cluster.end(); ++prow)
        {
            for (NdRange::iterator pcol = cluster.begin();
                    pcol != cluster.end(); ++pcol)
            {
                ArrayFi global_row = new_harmonics.globalIndex(*pb, *prow);
                ArrayFi global_col = new_harmonics.globalIndex(*pb, *pcol);

                if ( (harmonics.baseIndex(global_row) !=
                        harmonics.baseIndex(global_col)).any() )
                {
                    // the elements belong to different diagonal blocks in the original
                    // symbol
                    EXPECT_EQ(L.ref(*pb, *prow, *pcol), complex<double>(0));
                } else {
                    // check the global coordinate hash
                    EXPECT_EQ(L.ref(*pb, *prow, *pcol),
                        complex<double> (
                            global_row(0)         + global_row(1) * 100 +
                            global_col(0) * 10000 + global_col(1) * 1000000) );
                }
            }
        }
    }

}

bool similar_arrays(ArrayFd a, ArrayFd b)
{
    ArrayFd ones = ArrayFd::Ones(a.rows());

    return ((a - b).abs() < 1e-12 * ones).all();
}



