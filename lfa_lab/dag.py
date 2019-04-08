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

from .core import *
from .util import NdArray
import numpy as np
from abc import ABCMeta, abstractmethod
from six import with_metaclass

__all__ = [
    'StencilNode',
    'IdentityNode',
    'BlockNode',
    'PeriodicStencilNode',
    'FlatInterpolationNode',
    'FlatRestrictionNode',
    'ZeroNode',
    'HpFilterNode',
    'SystemNode'
]

default_resolution = 32

class Node(object):
    """This node represents general operators whose symbols can be computed.

    The addition, multiplication and power operators exist for this class.

    This class should not be instanciated directly.
    """

    def __init__(self):
        self._marked = False
        self._ref_count = 0

    def __add__(self, other):
        """The sum of `self` and `other`."""
        return NodeAdd(self, other)

    def __sub__(self, other):
        """The difference of `self` and `other`."""
        return self + (-1) * other

    def __mul__(self, other):
        """The composition of `self` and `other`."""
        return NodeMul(self, other)

    def __rmul__(self, other):
        """Scalar multiplication of `self` and `other`."""
        return NodeScalarMul(other, self)

    def __pow__(self, p):
        """Compute the power of an operator."""
        assert(p >= 1)
        aux = self
        for i in range(1, p):
            aux = aux * self
        return aux

    def inverse(self):
        """The inverse of the operator.

        :rtype: Node
        """
        return NodeInverse(self)

    def adjoint(self):
        """The adjoint of the operator.

        :rtype: Node
        """
        return NodeAdjoint(self)

    def __getitem__(self, index):
        """Access the element of a system.

        :rtype: Node
        """
        return NodeSubscript(self, index)

    @property
    def output_grid(self):
        return self.properties.outputGrid()

    @property
    def input_grid(self):
        return self.properties.inputGrid()

    @property
    def dim(self):
        return self.properties.dimension()

    def inc_ref(self):
        self._ref_count += 1

    def dec_ref(self):
        assert(self._ref_count > 0)
        self._ref_count -= 1

        if self._ref_count == 0:
            del self._symbol

    # Graph algorithms
    def _unmark_all(self):
        if self._marked:
            self._marked = False
            for d in self.dependencies:
                d._unmark_all()

    def _walk_dependencies_first(self, f):
        """Walks the DAG and calls f on every node. The function f is called
        on the dependencies first."""

        def visit(node):
            if not node._marked:
                node._marked = True
                for d in node.dependencies:
                    visit(d)
                f(node)

        visit(self)
        self._unmark_all()

    def symbol(self, desired_resolution = None):
        """The symbol of the operator.

        :rtype: Symbol
        """
        global default_resolution

        # ensure that we are not deleted
        self.inc_ref()

        d = self.properties.dimension()
        if desired_resolution is None:
            desired_resolution = np.ones(d) * default_resolution

        # compute the resolution
        resolution = self.properties.adjustResolution(desired_resolution)
        conf = SamplingProperties(resolution, self.properties.inputGrid())

        def set_configuration(n):
            n.configuration = conf

        # increase the refcount for all dependencies
        def ref_dependencies(node):
            for dep in node.dependencies:
                dep.inc_ref()

        def compute_and_deref(node):
            node.compute_symbol()

            for dep in node.dependencies:
                dep.dec_ref()

        self._walk_dependencies_first(ref_dependencies)
        self._walk_dependencies_first(set_configuration)
        self._walk_dependencies_first(compute_and_deref)

        symbol = self._symbol
        self.dec_ref() # free memory

        return symbol

class Splitable(with_metaclass(ABCMeta)):

    @abstractmethod
    def matching_identity(self):
        pass

    @abstractmethod
    def matching_zero(self):
        pass

    @abstractmethod
    def diag(self):
        pass

    @abstractmethod
    def upper(self):
        pass

    @abstractmethod
    def lower(self):
        pass

class IdentityNode(Node, Splitable):
    """The identity on a certain grid."""

    def __init__(self, grid):
        super(IdentityNode, self).__init__()
        self.grid = grid
        self.dependencies = []

        domain = SplitFrequencyDomain(grid)
        self.properties = FoProperties(domain, domain)

    def compute_symbol(self):
        self._symbol = Symbol.Identity(self.grid, self.configuration)

    def matching_identity(self):
        return self

    def matching_zero(self):
        return ZeroNode(self.grid)

    def diag(self):
        return self

    def upper(self):
        return ZeroNode(self.grid)

    def lower(self):
        return ZeroNode(self.grid)

class ZeroNode(Node):

    def __init__(self, grid):
        super(ZeroNode, self).__init__()
        self.grid = grid
        self.dependencies = []

        domain = SplitFrequencyDomain(grid)
        self.properties = FoProperties(domain, domain)

    def compute_symbol(self):
        self._symbol = Symbol.Zero(self.grid, self.configuration)

    def matching_identity(self):
        return IdentityNode(self.grid)

    def matching_zero(self):
        self

    def diag(self):
        return self

    def upper(self):
        return self

    def lower(self):
        return self

class GeneratorNode(Node):
    """
        A node whose symbol is created by a generator.
    """
    def __init__(self, generator, dependencies = []):
        super(GeneratorNode, self).__init__()

        self._generator = generator
        self.properties = self._generator.properties()
        self.dependencies = dependencies

    def compute_symbol(self):
        self._symbol = self._generator.generate(self.configuration)

class StencilNode(GeneratorNode):
    """An Operator given by a stencil.

    To construct a stencil node, use the
    :py:func:`lfa_lab.operator.from_stencil` method.

    :param SparseStencil stencil:
      The stencil that should be turned into an operator.
    :param Grid grid: The corresponding grid.

    :ivar stencil: The stencil of the operator.
    :ivar grid: The grid of the operator.
    """
    def __init__(self, stencil, grid):
        gen = FoStencil(stencil, grid)
        self.stencil = stencil
        self.grid = grid

        super(StencilNode, self).__init__(gen)

    def matching_identity(self):
        """An identity operator that has the same output and input grid as the
        current one."""
        return IdentityNode(self.grid)

    def matching_zero(self):
        """A zero operator that has the same output and input grid as the
        current one."""
        return ZeroNode(self.grid)

    def diag(self):
        """A stencil operator that was constructed using the diagonal entries
        of the stencil of this operator."""
        return StencilNode(self.stencil.diag(), self.grid)

    def upper(self):
        """A stencil operator that was constructed using the strictly upper
        triangular part of the stencil of this operator."""
        return StencilNode(self.stencil.upper(), self.grid)

    def lower(self):
        """A stencil operator that was constructed using the strictly lower
        triangular part of the stencil of this operator."""
        return StencilNode(self.stencil.lower(), self.grid)


class FlatRestrictionNode(GeneratorNode):
    __slots__ = 'output_grid'

    def __init__(self, output_grid, input_grid):
        gen = flat_restriction_sb(output_grid, input_grid)
        self.output_grid = output_grid

        super(FlatRestrictionNode, self).__init__(gen)

    def matching_identity(self):
        raise Exception('No matching identity available.')

    def matching_zero(self):
        # ToDo: Create a proper ZeroNode between two grids
        return ZeroNode(self.output_grid) * self

class FlatInterpolationNode(GeneratorNode):
    __slots__ = 'output_grid'

    def __init__(self, output_grid, input_grid):
        gen = flat_interpolation_sb(output_grid, input_grid)
        self.output_grid = output_grid

        super(FlatInterpolationNode, self).__init__(gen)

    def matching_identity(self):
        raise Exception('No matching identity available.')

    def matching_zero(self):
        # ToDo: Create a proper ZeroNode between two grids 
        return ZeroNode(self.output_grid) * self

class NodeAdd(Node):

    def __init__(self, a, b):
        super(NodeAdd, self).__init__()

        self._a = a
        self._b = b
        self.dependencies = [a, b]
        self.properties = a.properties + b.properties

    def compute_symbol(self):
        self._symbol = self._a._symbol + self._b._symbol

class NodeMul(Node):

    def __init__(self, a, b):
        super(NodeMul, self).__init__()

        self._a = a
        self._b = b
        self.dependencies = [a, b]
        self.properties = a.properties * b.properties

    def compute_symbol(self):
        self._symbol = self._a._symbol * self._b._symbol

    def matching_identity(self):
        # ToDo: matching_identitys can be constructed from the properties.
        # Hence a single implementation is sufficient. The same is true for
        # matching_zero.
        if isinstance(self.properties, SystemSymbolProperties):
            return SystemNode.make_identity(self.properties.outputGrid(),
                                            self.properties.rows(),
                                            self.properties.cols())
        else:
            return IdentityNode(self.properties.outputGrid())

    def matching_zero(self):
        return self._a.matching_zero() * self._b.matching_zero()

class NodeScalarMul(Node):

    def __init__(self, a, b):
        super(NodeScalarMul, self).__init__()

        self._a = a
        self._b = b
        self.dependencies = [b]
        self.properties = a * b.properties

    def compute_symbol(self):
        self._symbol = self._a * self._b._symbol


class NodeInverse(Node):

    def __init__(self, other):
        super(NodeInverse, self).__init__()

        self._other = other
        self.dependencies = [other]
        self.properties = other.properties.inverse()

    def compute_symbol(self):
        self._symbol = self._other._symbol.inverse()

class NodeAdjoint(Node):
    def __init__(self, other):
        super(NodeAdjoint, self).__init__()

        self._other = other
        self.dependencies = [other]
        self.properties = other.properties.adjoint()

    def compute_symbol(self):
        self._symbol = self._other._symbol.adjoint()

class NodeSubscript(Node):
    def __init__(self, container_node, index):
        super(NodeSubscript, self).__init__()

        self._container_node = container_node
        self._index = index
        self.dependencies = [container_node]
        self.properties = container_node.properties.element_properties()

    def compute_symbol(self):
        self._symbol = self._container_node._symbol[self._index]

class BlockNode(Node):

    def __init__(self, scalars):
        """
        scalars is a nested list containing the scalar nodes
        """
        super(BlockNode, self).__init__()

        self._scalars = NdArray(scalars, bottom_p = lambda e: hasattr(e, 'properties'))

        first_scalar = self._scalars[np.zeros(self._scalars.dim, np.int)]

        grid = first_scalar.properties.outputGrid()
        period = self._scalars.shape

        self._generator = BlockSb(grid, period)

        self.dependencies = list(self._scalars)
        self.properties = self._generator.properties()

    def compute_symbol(self):

        # Gather the symbols
        symbols = SymbolNdArray(self._scalars.shape)
        for i in symbols.indices():
            symbols[i] = self._scalars[i]._symbol

        # Generate the combined symbol
        self._generator.scalarSymbols(symbols)
        self._symbol = self._generator.generate(self.configuration)


# ToDo: Remove this class.
class PeriodicStencilNode(BlockNode, Splitable):

    def __init__(self, stencils, grid):
        # convert all stencils to symbols
        self.grid = grid
        self.stencils = stencils
        ops = stencils.entries.map(lambda s: StencilNode(s, grid))

        super(PeriodicStencilNode, self).__init__(ops)

    def matching_identity(self):
        return IdentityNode(self.grid)

    def matching_zero(self):
        return ZeroNode(self.grid)

    def diag(self):
        """Diagonal part of the stencil."""
        return PeriodicStencilNode(self.stencils.diag(), self.grid)

    def upper(self):
        """Upper part of the stencil."""
        return PeriodicStencilNode(self.stencils.upper(), self.grid)

    def lower(self):
        """Lower part of the stencil."""
        return PeriodicStencilNode(self.stencils.lower(), self.grid)

class HpFilterNode(GeneratorNode):
    """High pass filter symbol."""

    def __init__(self, fine_grid, coarse_grid):

        super(HpFilterNode, self).__init__(
                HpFilterSb(fine_grid, coarse_grid))

class SystemNode(Node,Splitable):
    """System of symbols.

    :param entries: The entries as a list of lists.
    """

    def __init__(self, entries):
        super(SystemNode, self).__init__()

        ps = FoPropertiesMatrix(len(entries), len(entries[0]))
        self.dependencies = []
        for i in range(ps.rows()):
            assert(len(entries[i]) == ps.cols())
            for j in range(ps.cols()):
                ps[i,j] = entries[i][j].properties
                self.dependencies.append(entries[i][j])

        self._entries = entries
        self.properties = properties_of_symbol_system(ps)

    def compute_symbol(self):
        s = SymbolMatrix(self.properties.rows(), self.properties.cols())
        for i in range(self.properties.rows()):
            for j in range(self.properties.cols()):
                s[i,j] = self._entries[i][j]._symbol

        self._symbol = combine_symbols_into_system(self.properties, s)

    @property
    def grid(self):
        return self._entries[0][0].grid

    @staticmethod
    def make_identity(grid, rows, cols):
        I = IdentityNode(grid)
        Z = ZeroNode(grid)

        def entry(i,j):
            if i == j:
                return I
            else:
                return Z

        def row(i):
            return [ entry(i,j) for j in range(cols) ]

        new_entries = [ row(i) for i in range(rows) ]
        return SystemNode(new_entries)

    def matching_identity(self):
        # ToDo: assert equal input and output
        return SystemNode.make_identity(self.properties.inputGrid(),
                                        self.properties.rows(),
                                        self.properties.cols())

    def matching_zero(self):
        Z = self._entries[0][0].matching_zero()

        new_entries = [ [ Z for j in range(self.properties.rows()) ]
                        for j in range(self.properties.cols()) ]
        return SystemNode(new_entries)

    def diag(self):
        new_entries = []
        for i, row in enumerate(self._entries):
            new_row = []
            for j, entry in enumerate(row):
                if i == j:
                    new_row.append(entry.diag())
                else:
                    new_row.append(ZeroNode(entry.grid))
            new_entries.append(new_row)
        return SystemNode(new_entries)

    def elementwise_diag(self):
        new_entries = \
            list(map(lambda l: list(map(lambda e: e.diag(), l)), self._entries))

        return SystemNode(new_entries)

    def upper(self):
        pass

    def lower(self):
        pass


