.. _dream-eqget:

Equilibrium tools
=================
DREAM provides some Python classes and graphical tools for working with magnetic
equilibria in various formats. As described on the page about the
:ref:`_radialgrid-numerical`, DREAM uses a LUKE format for equilibrium data.
We however also provide a number of classes for converting data from other
formats to the LUKE format. Currently, data can be exported from ASDEX Upgrade
shotfiles and EQDSK files.

.. contents:: Page overview
   :local:
   :depth: 3

eqget GUI
---------
The ``eqget`` GUI can be launched by running the Python script
``tools/eqget/eqget.py``. It allows you to load and visualize magnetic
equilibrium data in a number of supported formats (listed above).

.. figure:: _static/figs/eqget-GUI.png
   :width: 60%
   :align: center
   :alt: Screenshot of the DREAM equilibrium visualization tool

Python classes
--------------
A number of scripts are provided for working with magnetic equilibria for
DREAM. At the time of writing, the following scripts are available:

- ``AUG.py`` -- loading data from ASDEX Upgrade shotfiles
- ``EQDSK.py`` -- class for loading data from EQDSK files
- ``EqFile.py`` -- loading data from LUKE files (standard format used in DREAM)

Examples
*******
**ASDEX Upgrade shotfiles**

.. code-block:: python

   import AUG

   ...
   eq = AUG.get_LUKE(shot=shot, time=time, npsi=80, ntheta=80, filename=None)

   # eq.keys() = ['id','Rp','psi_apRp','theta','ptx','pty','ptBx','ptBy','ptBPHI']

**EQDSK files**

.. code-block:: python

   from EQDSK import EQDSK

   ...
   eq = EQDSK(filename, cocos=1, process=True, override_psilim=False)
   luke_eq = eq.get_LUKE(npsi=80, ntheta=80)
   # or
   eq.save_LUKE(filename, npsi=80, ntheta=80)

.. note::

   The parameter ``override_psilim`` can be used to work around poorly resolved
   EQDSK equilibria, for which the poloidal flux of the last closed flux surface
   (LCFS) is inaccurately specified in the file. If ``override_psilim=True``,
   the flux surface with ``psin=1-2e-4`` (``psin`` being the normalized poloidal
   flux, with ``psin=1`` on the LCFS) is used as the LCFS. If a numerical value
   is specified for ``override_psilim``, the flux surface with
   ``psin=1-override_psilim`` is used as the LCFS.


