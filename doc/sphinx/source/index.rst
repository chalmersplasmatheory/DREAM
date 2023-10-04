.. DREAM documentation master file, created by
   sphinx-quickstart on Mon Mar  9 17:21:34 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: _static/figs/logo1.svg
   :width: 100%
   :alt: DREAM

Welcome to the documentation for DREAM!
=======================================

DREAM (for *Disruption and Runaway Electron Analysis Model*) is a scientific
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

.. note::

    The full code is freely available on
    `GitHub <https://github.com/chalmersplasmatheory/DREAM>`_.


Table of contents
*****************
.. toctree::
   :maxdepth: 2

   compiling
   kernel/index
   frontend/index
   atomicdata
   start/index
   imas/index
   todo

.. raw:: html

   <div class="admonition note" style="width:49%; float:left">
     <a href="compiling.html"><p class="admonition-title">Compiling</p></a>
     <p>Step-by-step guide for setting up DREAM</p>
   </div>
   <div class="admonition note" style="width:49%; float:right">
     <a href="start/index.html"><p class="admonition-title">Getting started</p></a>
     <p>Examples of how to work with DREAM</p>
   </div>
   <div style="clear:both"></div>
   <div class="admonition note" style="width:49%; float:left">
     <a href="frontend/api/Settings/DREAMSettings.html">
       <p class="admonition-title">Settings documentation</p>
     </a>
     <p>Browse all settings available in DREAM</p>
   </div>
   <div class="admonition note" style="width:49%; float:right">
     <a href="frontend/api/DREAMOutput.html">
       <p class="admonition-title">Output documentation</p>
     </a>
     <p>Browse the Python output interface documentation</p>
   </div>
   <div style="clear:both"></div>
   <div class="admonition note" style="width:49%; float:left">
     <a href="atomicdata.html">
       <p class="admonition-title">Atomic data</p>
     </a>
     <p>Learn how to add support for more elements</p>
   </div>
   <div class="admonition note" style="width:49%; float:right">
     <a href="start/examples/index.html">
       <p class="admonition-title">Examples</p>
     </a>
     <p>View examples of DREAM simulations</p>
   </div>
   <div class="admonition note" style="width:49%; float:left">
     <a href="imas/index.html">
       <p class="admonition-title">IMAS interface</p>
     </a>
     <p>See the details of the IMAS interfacing</p>
   </div>
   <div style="clear:both"></div>

Introduction
************
The **Disruption Runaway Electron Analysis Model** DREAM is designed to
simulate runaway electron generation during a tokamak disruption.

For users
---------
If you would like to run simulations with DREAM, you should have a look at the
documentation of the :ref:`python-frontend`. In particular, options for setting 
up the equation can be found on the :ref:`ds-eqsys` settings page.

For developers
--------------
If you are trying to understand how DREAM works on a deeper level, you are most
likely interested in the documentation for the :ref:`cpp-kernel`.

For physicists
--------------
If you are curious about the physical model implemented in DREAM, you should
have a look at the physics notes located under
`doc/notes/ <https://github.com/chalmersplasmatheory/DREAM/tree/master/doc/notes>`_
in the DREAM git repository.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
