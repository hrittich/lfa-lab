###################
Homepage of LFA Lab
###################

Introduction
============

Local Fourier analysis uses the discrete Time Fourier transform. The transform
represents grid functions using the integral

.. math::

   f(x) =
   \int_{\Theta_\mathbf{h}} \hat{f}(\theta) e^{\mathrm{i} \langle \mathbf{x},
    \theta \rangle} \,\mathrm{d}\theta

where :math:`\Theta_\mathbf{h} := [ 0, \tfrac{2\pi}{h_1} ) \times \cdots \times [ 0, \tfrac{2\pi}{h_d} )`.

Sponsors
========

The development of LFA Lab was partly supported by the following institutions.

.. image:: logo-dfg.png
   :target: http://www.dfg.de/
   :alt: DFG, Deutsche Forschungsgemeinschaft
   :height: 2.5em

.. image:: logo-sppexa.png
   :target: http://www.sppexa.de/
   :alt: SPPEXA, Software for Exascale Computing
   :height: 2.5em

.. image:: logo-exastencils.png
   :target: http://www.exastencils.org/
   :alt: ExaStencils, Advances Stencil-Code Engineering
   :height: 2.5em


Contents
========

.. toctree::
  :maxdepth: 2

  install
  tutorial
  api
  citation
  download
  license
  glossary

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

