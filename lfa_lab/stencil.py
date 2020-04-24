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

from .core.extension import _SparseStencil
from .core import *
from .util import *
from .dag import *
import numpy as np
from six import iteritems

__all__ = [
    'SparseStencil',
    'PeriodicStencil',
    'DenseStencil'
]

class SparseStencil(_SparseStencil):
    r"""Storage for a stencil.

    A (constant) stencil is a map
    :math:`s: \mathbb{Z}^d \to \mathbb{C}`. We assume that only finitely many
    function values are non-zero. Therefore, we can represent the stencil by a
    list containing of pairs consisting of the function argument and the
    corresponding value.

    :param entries:
        The entries of the stencil. This argument is supposed to
        be a list of entries. Each entry is a tuple consisting of the offset
        (a tuple) and the corresponding scalar.
    :type entries: list of (tuple of (tuple, double))
    """

    def __init__(self, entries = None):
        super(SparseStencil, self).__init__()

        if entries is None:
            pass
        elif hasattr(entries, '__iter__'):
            for offset, value in entries:
                self.append(offset, value)
        else:
            raise TypeError('The entries argument is not iterable.')

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        return self.nonZeros()

    def __getitem__(self, i):
        """Return the element at index i."""
        offset = self._get_offset_at(i)
        offset = tuple(map(lambda j: int(j), offset))

        value = self._get_value_at(i)
        if (value.imag == 0):
            value = value.real

        return (offset, value)

    def __str__(self):
        return str(list(self))

    def __repr__(self):
        return '(stencil {})'.format(str(self))

    def _create_empty(self):
        return SparseStencil()

    def map(self, f):
        result = self._create_empty()
        for offset, value in self:
            result.append(*f(offset, value))
        return result

    def filter(self, take_pred):
        result = self._create_empty()
        for offset, value in self:
            if take_pred(offset, value):
                result.append(offset, value)
        return result

    @property
    def dim(self):
        return self.dimension()

    def diag(self):
        import numpy as np

        def is_diag(offset, value):
            return (np.array(offset, np.int) == 0).all()

        return self.filter(is_diag)

    def lower(self):
        import numpy as np
        zero = np.zeros(self.dim)

        return self.filter(lambda o, v: lex_less(o, zero))

    def upper(self):
        import numpy as np
        zero = np.zeros(self.dim)

        return self.filter(lambda o, v: lex_less(zero, o))

    def transpose(self):
        return self.map(lambda o, v: (-np.array(o), v))

    def conjugate(self):
        return self.map(lambda o, v: (o, v.conjugate()))

    def adjoint(self):
        return self.transpose().conjugate()

    def __rmul__(self, other):
        """Scalar multiplication"""
        return self.map(lambda o, v: (o, other*v))

    def __mul__(self, other):
        """Stencil composition."""

        entries = {}

        for (o1, v1) in self:
            for (o2, v2) in other:
                o1 = np.array(o1)
                o2 = np.array(o2)
                d = tuple(o1 + o2)

                if d in entries:
                    entries[d] += v1 * v2
                else:
                    entries[d] = v1 * v2

        return SparseStencil(iteritems(entries))



class DenseStencil(object):

    def __init__(self, origin, entries):
        """
            entries = nested list of entries
        """

        self._origin = np.array(origin)
        self._entries = NdArray(entries)

    def __iter__(self):
        indices = NdRange(self._entries.shape)
        for i in indices:
            yield (np.array(i) + self._origin, \
                   self._entries[i])

    def __getitem__(self, i):
        return self._entries[np.array(i) - self._origin]

    def __setitem__(self, i, v):
        self._entries[np.array(i) - self._origin] = v

    def to_sparse(self):
        return SparseStencil(entries=self)


class PeriodicStencil(object):
    """Storage for a periodic stencil.

    A periodic stencil is a family of stencils
    :math:`\{ s_\mathbf{x} \}_\mathbf{x \in \mathbb{Z}^d}` such that
    :math:`s_\mathbf{x} = s_{\mathbf{x}'}` for
    :math:`(\mathbf{x} - \mathbf{x}') \in \mathbf{p} \mathbb{Z}^d`.
    We call :math:`\mathbf{p}` the period of the stencil.
    We store a periodic stencil by storing the stencils for
    :math:`0 \le \mathbf{x} < \mathbf{p}` in an instance of
    :py:class:`lfa_tool.util.NdArray`.

    :param entries:
        The entries as an :py:class:`lfa_tool.util.NdArray` or as a proper
        argument to the constructor of :py:class:`lfa_tool.util.NdArray`.
    """

    def __init__(self, entries):
        self.entries = NdArray(entries)

        # ToDo: Check that they are all defined on the same grid

    @property
    def dim(self):
        return self.entries.dim

    def map(self, f):
        return PeriodicStencil(self.entries.map(f))

    def diag(self):
        return self.map(lambda s: s.diag())

    def lower(self):
        return self.map(lambda s: s.lower())

    def upper(self):
        return self.map(lambda s: s.upper())

