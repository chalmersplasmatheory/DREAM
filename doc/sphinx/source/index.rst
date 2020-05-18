.. DREAM documentation master file, created by
   sphinx-quickstart on Mon Mar  9 17:21:34 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DREAM's documentation!
=================================

DREAM (for *Disruption and Runaway Electron Avoidance Model*) is a scientific
code which solves a system of equations describing the evolution of a tokamak
plasma during a disruption. The tool is primarily designed for studying the
evolution of runaway electrons in the disruption.

The code is separated into two parts: a C++ library, which constitutes the
computational kernel of the code, and a Python interface for easy access and
configuration of simulations. This documentation describes both parts of the
code, with the documentation for the former being mostly intended for code
developers, while the latter is intended for code users.

For more details on the physical models implemented in DREAM, please consult
the physics manuals available in the DREAM GitHub repository.

The DREAM code is developed and maintained by the
`Plasma Theory group <https://ft.nephy.chalmers.se/>`_ at Chalmers University of
Technology, Gothenburg, Sweden.

.. toctree::
   :maxdepth: 4

   intro
   compiling
   kernel/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
