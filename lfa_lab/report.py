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

""" Produces a report for an operator """

from . import plot
import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    plt = None

def print_report(E, title=''):
    """Print a report about the operator E on the screen.

    This opens a figure. You may have to call
    matplotlib.pyplot.show() to show it.
    """

    smpl = E.symbol()

    if len(title) > 0:
        print('')
        print(title)
        print(('=' * len(title)))
        print('')
    print(('r(E)       = {0}'.format(smpl.spectral_radius())))
    print(('|| E ||    = {0}'.format(smpl.spectral_norm())))
    print(('Coupling   = {} x {}'.format(
        np.array(smpl.output_coupling_shape()),
        np.array(smpl.input_coupling_shape()))))
    print(('Resolution = {} x {}'.format(
        np.array(smpl.output_shape()),
        np.array(smpl.input_shape()))))


    if smpl.dimension() == 1:
        plot.plot_1d(smpl)
    elif smpl.dimension() == 2:
        plt.figure()
        plt.subplot(1,2, 1)

        plot.plot_2d(smpl, norm_type = 'rows', new_figure=False)
        plt.title('{} (output)'.format(title))

        plt.subplot(1,2, 2)
        plot.plot_2d(smpl, norm_type = 'columns', new_figure=False)
        plt.title('{} (input)'.format(title))

    else:
        print('I am sorry. I can not plot {}-dimensional data.'
                .format(smpl.dimension()))

def save_report(E, file_name_prefix, title='', standalone=False):
    r"""Store a report about an operator as a LaTeX file.

    :param str file_name_prefix: The file name without file extension.
    :param str title: A short description of the operator.
    :param bool standalone:
        Whether the TeX file contains the preamble or is suitable for
        embedding into another document.
    """

    template = r'''
\begin{{minipage}}{{0.38\textwidth}}
\begin{{tabular}}{{@{{}}lp{{11em}}@{{}}}}
  \toprule
  \multicolumn{{2}}{{c}}{{Properties}} \\
  \midrule
  Info      & {title} \\
  $r(E)$    & ${spectral_radius}$ \\
  $\| E \|$ & ${spectral_norm}$ \\
  Coupling  & \verb\{output_coupling_shape}\ $\times$
              \verb\{input_coupling_shape}\ \\
            & ${output_coupling} \times {input_coupling}$ \\
  Resol.    & \verb\{output_res}\ $\times$
              \verb\{input_res}\ \\
  \bottomrule
\end{{tabular}}
\end{{minipage}}
\hfill
\begin{{minipage}}{{0.60\textwidth}}
  \centering
  \includegraphics[width=\textwidth]{{{file_name}_in_spec}}
  Input Spectrum
  \includegraphics[width=\textwidth]{{{file_name}_out_spec}}
  Output Spectrum
\end{{minipage}}
'''

    template_doc = r'''
\documentclass[DIV=15]{{scrartcl}}
\usepackage{{amsmath}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\begin{{document}}
{contents}
\end{{document}}'''

    smpl = E.symbol()

    contents = template.format(
        file_name = file_name_prefix,
        spectral_radius = smpl.spectral_radius(),
        spectral_norm = smpl.spectral_norm(),
        title=title,
        output_coupling=smpl.output_coupling_size(),
        input_coupling=smpl.input_coupling_size(),
        output_coupling_shape=np.array(smpl.output_coupling_shape()),
        input_coupling_shape=np.array(smpl.input_coupling_shape()),
        output_res=np.array(smpl.output_shape()),
        input_res=np.array(smpl.input_shape()))

    if standalone:
        contents = template_doc.format(contents=contents)

    with open(file_name_prefix + '.tex', 'w') as tex_fp:
        tex_fp.write(contents)

    #plt.figure()

    plot.plot_2d(smpl, norm_type = 'columns')
    plt.axis('scaled')
    plt.savefig(file_name_prefix + '_in_spec.eps', format='eps')
    plt.close()

    plot.plot_2d(smpl, norm_type = 'rows')
    plt.axis('scaled')
    plt.savefig(file_name_prefix + '_out_spec.eps', format='eps')
    plt.close()



