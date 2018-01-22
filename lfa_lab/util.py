# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018  Hannah Rittich
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

from __future__ import division
from future_builtins import zip

import numpy as np
from .core import NdRange

def is_even(x):
    return (x % 2) == 0

def lex_less(a, b):
    """Lexicographical comparison."""

    for p, q in zip(a, b):
        if p != q:
            # the first non-equal entry determines the result
            return p < q;

    # all entries are the same
    return False

def lex_greater(a, b):
    lex_less(b, a)


# Linear interpolation
def lin_int(p0, p1, x):
    x0, y0 = p0
    x1, y1 = p1
    t = (x - x0) / (x1 - x0)
    y = (y1 - y0) * t + y0
    return y



class NdArray:
    """An :math:`n`-D array.

    When constructing an NdArray either `entires` or `shape` must be given.
    When `entries` is given, either `dim` of `bottom_p` must be defined.

    Here is a simple example for the usage of this class::

        x = lfa_lab.util.NdArray([[1,2,3],[4,5,6]], dim=2)
        print(x[1,2]) # prints 6
        x[1,2] = 7
        print(x)      # prints NdArray: [[1, 2, 3], [4, 5, 7]]

    :param entries: The entries of the n-D array as a nested list.
    :type entries: list of (list of ... (list of object) ... )
    :param int dim: The dimension :math:`n` of the array.
    :param bottom_p: A predicate to determine the dimension of the array.
    :type bottom_p: object -> bool
    :param tuple shape: The number of entries per dimension.
    """

    def __init__(self, entries = None, dim = None, bottom_p = None, \
            shape = None):

        assert(entries is None or shape is None)

        self._dim = dim

        if (isinstance(entries, NdArray)):
            # copy the entries out of an NdArray
            entries = entries._entries

        # testing the list type
        if dim is None and bottom_p is None:
            bottom_p = lambda e: type(e) != list

        if bottom_p is not None:
            # find the dimension by probing each level for bottom_p
            self._dim = 0
            aux = entries
            while not bottom_p(aux):
                aux = aux[0]
                self._dim += 1

        if entries is None and shape is not None:
            entries = self._allocate(shape)
            self._dim = len(shape)

        assert(self._dim is not None)

        self._entries = entries

    def _allocate(self, shape):
        if len(shape) == 1:
            return [ None for _ in range(shape[0]) ]
        else:
            return [ self._allocate(shape[1:]) for _ in range(shape[0]) ]

    @property
    def dim(self):
        return self._dim

    def __getitem__(self, pos):
        """Return the element at position `pos`.

        :param tuple pos: The index of the entry.
        """
        def access(entries, pos):
            if len(pos) == 0:
                return entries
            else:
                return access(entries[pos[0]], pos[1:])

        return access(self._entries, pos)

    def __setitem__(self, pos, value):
        assert(len(pos) == self.dim)

        def write(entries, pos):
            if len(pos) == 1:
                entries[pos[0]] = value
            else:
                write(entries[pos[0]], pos[1:])

        write(self._entries, pos)

    def __str__(self):
        return 'NdArray: {}'.format(self._entries)


    def __iter__(self):

        def yield_entries(d, entries):
            if d > 0:
                for e in entries:
                    for ee in yield_entries(d-1, e):
                        yield ee
            else:
                yield entries

        for e in yield_entries(self._dim, self._entries):
            yield e

    @property
    def shape(self):
        s = []
        aux = self._entries
        for i in range(self._dim):
            s.append(len(aux))
            aux = aux[0]
        return s

    def map(self, f):
        # map over a nested list (tree) of depth d
        def map_tree(d, entries):
            if d == 1:
                return map(f, entries)
            else:
                return map(lambda e: map_tree(d-1, e), entries)
        return map_tree(self.dim, self._entries)



def _test():
    a = NdArray([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
    print(a.dim)
    a[1,0,1] = 99
    print(a)

    for i in a:
        print(i)

    b = a.map(lambda x: x * 2)
    print(b)

class ShiftedGrid(object):

    def __init__(self, shift, shape):
        self.shift = np.array(shift)
        self.shape = np.array(shape)

        self._plain = NdRange(self.shape)

    @staticmethod
    def from_corners(start, end):
        start = np.array(start)
        end = np.array(end)
        ones = np.ones_like(start)

        return ShiftedGrid(start, end - start + ones)

    def __iter__(self):
        for p in self._plain:
            yield (self.shift + p)


