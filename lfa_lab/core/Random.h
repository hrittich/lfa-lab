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

#ifndef LFA_MULTIGRID_RANDOM_H
#define LFA_MULTIGRID_RANDOM_H

#include "Common.h"

namespace lfa {

/** Pseudo random number generator which generates integers.
 * Linear congruential generator. */
class LcgPrng
{
    public:
        typedef unsigned int number_t;
        typedef number_t seed_t;

        inline LcgPrng() {
            x = 42;
        }

        inline void seed(seed_t s) {
            x = s & 0x7FFFFFFF;
        }

        inline number_t operator() ()
        {
            // compute (a*x + c) mod 2^31
            x = (a * x + c) & 0x7FFFFFFF;

            // return the first 30 bits
            return x & 0x3FFFFFFF;
        }

        inline static number_t max() {
            return 0x3FFFFFFF;
        }

    private:
        number_t x;

        static const number_t a = 1103515245;
        static const number_t c = 12345;
};

/** Pseudo random number generator which generates integers.
 * Produces "better" numbers than LCG, but is more expensive.  */
class IcgPrng
{
    public:
        typedef unsigned int number_t;
        typedef number_t seed_t;

        inline IcgPrng() {
            m_x = 42;
        }

        inline void seed(seed_t s) {
            m_x = s & 0x7FFFFFFF;
        }

        /* Extended Euclidean algorithm. Solves ax + by = gcd(a, b) and
         * returns gcd(a,  b). */
        static inline number_t extended_gcd(number_t a, number_t b,
                number_t& out_x, number_t& out_y)
        {
            number_t x = 0, y = 1, last_x = 1, last_y = 0;
            number_t quot, rem, aux;

            while (b != 0)
            {
                quot = a / b;
                rem = a % b;
                a = b;
                b = rem;

                aux = last_x - quot*x;
                last_x = x;
                x = aux;

                aux = last_y - quot*y;
                last_y = y;
                y = aux;
            }

            out_x = last_x;
            out_y = last_y;

            return a;
        }

        /** Computes the multiplicative inverse of a mod m. */
        static inline number_t modular_inverse(number_t a, number_t m)
        {
            number_t x, y;
            extended_gcd(a, m, x, y);
            return x;
        }

        /** Compute a pseudorandom number. */
        inline number_t operator() ()
        {
            // note that 0x7FFFFFFF is prime
            number_t x_inv = modular_inverse(m_x, 0x7FFFFFFF);

            m_x = (a * x_inv + c) & 0x7FFFFFFF;

            return m_x & 0x3FFFFFFF;
        }

        /** The largest number that will ever be returned. */
        static inline number_t max() { return 0x3FFFFFFF; }

    private:
        number_t m_x;

        static const number_t a = 1103515245;
        static const number_t c = 12345;
};

/** Generates uniform random numbers in the interval [0,1], given a uniform
 * integer random number generator. */
template <typename T, typename IntegerPrng = IcgPrng>
class UniformPrng
{
    public:
        typedef T number_t;
        typedef typename IntegerPrng::seed_t seed_t;

        inline UniformPrng(IntegerPrng rng = IntegerPrng() )
            : rng(rng)
        { }

        inline void seed(seed_t s) { rng.seed(s); }

        number_t operator() ()
        {
            return static_cast<number_t>(rng()) / IntegerPrng::max();
        }

    private:
        IntegerPrng rng;
};

/** Generates normal distributed random numbers, given a uniform PRNG. */
template <typename T, typename UniPrng = UniformPrng<T> >
class NormalPrng
{
    public:
        typedef T number_t;
        typedef typename UniPrng::seed_t seed_t;

        inline NormalPrng(number_t mean = 0, number_t variance = 1, UniPrng rng = UniPrng() )
            : mean(mean), variance(variance), rng(rng)
        { }

        inline void seed(seed_t s) { rng.seed(s); }

        /** Computes using the Marsaglia polar method. */
        number_t operator() ()
        {
            number_t a1, a2, q, p;

            do {
                a1 = 2.0 * rng() - 1.0;
                a2 = 2.0 * rng() - 1.0;

                q  = (a1 * a1) + (a2 * a2);
            } while ((q == 0.0) || (q >= 1));

            p = sqrt(-2 * log(q) / q);

            number_t z1 = a1 * p; // normal distributed N(0, 1)
            // number_t z2 = a2 * p; // normal distributed N(0, 1) (unneeded)

            return variance * z1 + mean;
        }

    private:
        number_t mean, variance;
        UniPrng rng;
};

/** Compute the mean value of a given sequence. */
template <typename Iterator>
inline typename Iterator::value_type
mean(Iterator begin, Iterator end)
{
    typedef typename Iterator::value_type number_t;
    int n = 0;
    number_t sum = 0;

    while(begin != end)
    {
        sum += *begin;
        ++n;
        ++begin;
    }

    return sum / n;
}

/** compute the variance of a given sequence. */
template <typename Iterator>
inline typename Iterator::value_type
variance(Iterator begin, Iterator end)
{
    typedef typename Iterator::value_type number_t;

    number_t mv = mean(begin, end);
    number_t aux, sum = 0;
    int n;

    while(begin != end)
    {
        aux = *begin - mv;
        sum += aux*aux;
        ++n;
        ++begin;
    }

    return sum / n;
}

}

#endif
