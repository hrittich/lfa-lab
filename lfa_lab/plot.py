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

"""Plotting commands for Fourier symbols.
"""

import numpy as np

__all__ = [
    'plot_2d',
    'plot_1d',
]

try:
    import matplotlib.pyplot as plt
    # from mpl_toolkits.mplot3d import Axes3D
    import mpl_toolkits.mplot3d

    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False

default_options = {
    'norm_type': 'columns',
    'style_2d': 'image',
    'style_3d': 'isosurface',
    'vmin': None,
    'vmax': None,
    'interpolation': 'bilinear',
    'colormap': 'RdBu_r',
    'new_figure': True,
}

def mesh(x, y, z):
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d

    fig = plt.gcf()
    fig.clear()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, z)
    return fig

def contour(x, y, z):
    import matplotlib.pyplot as plt

    fig = plt.gca()
    fig.clear()
    cs = plt.contour(x, y, z)
    plt.clabel(cs)
    return fig

def image(x, y, z, options):
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    import matplotlib.ticker as ticker

    fig = plt.gca()
    fig.clear()
    #fig.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    fig.xaxis.set_major_formatter(ticker.FormatStrFormatter(
        r'%.1f $\frac{2\pi}{h_1}$'))
    #fig.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    fig.yaxis.set_major_formatter(ticker.FormatStrFormatter(
        r'%.1f $\frac{2\pi}{h_2}$'))
    cs = plt.imshow(z,
                    interpolation=options['interpolation'],
                    extent=(0,1,0,1),
                    vmin=options['vmin'], vmax=options['vmax'],
                    cmap=get_cmap(options['colormap']))
    plt.colorbar()
    return fig

def plot_2d(sym, **kwargs):
    """Plot the sampling of a symbol.

    :param Symbol sym: The symbol that should be plotted.
    :param str norm_type: The type of the norm. Possible values: 'rows',
        'output', 'columns', and 'input'.
    :param str style_2d: The plotting style for 2D symbols. It can be set to
        'mesh', 'image', or 'contour'.
    """
    import matplotlib.pyplot as plt

    smpl = sym.symbol()

    options = default_options.copy()
    options.update(kwargs)

    if options['new_figure']: plt.figure()

    norm_type = options['norm_type']
    if norm_type == 'rows' or norm_type == 'output':
        M = np.abs(smpl.row_norms_2d())
    elif norm_type == 'columns' or norm_type == 'input':
        M = np.abs(smpl.col_norms_2d())
    else:
        raise Exception('Unknown norm type "{}".'.format(norm_type))

    (x,y) = np.meshgrid(np.linspace(0, 1, M.shape[0]),
                        np.linspace(0, 1, M.shape[1]))

    style_2d = options['style_2d']
    if style_2d == 'mesh':
        return mesh(x, y, M)
    elif style_2d == 'image':
        return image(x, y, M, options)
    elif style_2d == 'contour':
        return contour(x, y, M)
    else:
        raise Exception('Unknown style "{}".'.format(style_2d))


def plot_1d(op, **kwargs):
    """Plot a 1d symbol.

    :param str norm_type: The type of the norm. Possible values: 'rows',
        'output', 'columns', and 'input'.
    """
    import matplotlib.pyplot as plt

    options = default_options.copy()
    options.update(kwargs)

    if options['new_figure']:
        fig = plt.figure()
    else:
        fig = plt.gcf()

    smpl = op.symbol()
    if options['norm_type'] == 'rows' or options['norm_type'] == 'output':
        v = np.abs(smpl.row_norms_1d())
    elif options['norm_type'] == 'columns' or options['norm_type'] == 'input':
        v = np.abs(smpl.col_norms_1d())
    else:
        raise Exception('Unknown norm type "{}".'.format(options['norm_type']))

    t = np.linspace(0, 1, len(v))

    plt.plot(t, v)
    return fig

