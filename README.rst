This repository contains the accompanying python code for the paper

 "The Order of the Monster Finite Simple Group"

by Gerald Hoehn and Martin Seysen.


Requirements
------------


The python code can be executed on a 64-bit Linux, Windows, or macOS
system with python 3.8 or higher.

The following external software packages are required.

The python packages numpy, sympy, and mmgroup are required.
The standard way to install these packages is:

.. code-block::

    python3 -m pip install --upgrade numpy sympy mmgroup

This command also updates outdated versions of these packages,
which is highly recommended.
In a Windows system the string "python3" should be replaced by "python".

For details of the installation of the mmgroup package we refer to:

https://mmgroup.readthedocs.io/en/latest/api.html#installation-and-test


The GAP computer algebra system is also required. 
For installing GAP we refer to:

https://www.gap-system.org/install/

https://www.math.rwth-aachen.de/~Frank.Luebeck/GAPrsync/index.html


Quick installation and test
---------------------------

After installing the required external packages you may simply
clone this repo to you local computer and run:

.. code-block::

   cd monster_order
   cd axis_orbits
   python3 axis.py --all 
   python3 baby_axis.py --all 

In a Windows system the string "python3" should be replaced by "python".

Then you may follow the detailed instructions in file

monster_order/axis_orbits/readme.txt


