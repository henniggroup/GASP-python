GASP is a genetic algorithm for structure and phase prediction written in Python and interfaced to GULP_, LAMMPS_ and VASP_. It can search for the structures of clusters, 2D materials, wires, and bulk materials and do both fixed-composition and phase diagram searches.

.. _VASP: http://www.vasp.at/
.. _LAMMPS: http://lammps.sandia.gov/
.. _GULP: https://gulp.curtin.edu.au/gulp/ 


Getting GASP
============
It is easiest to install GASP and all its dependencies into a conda_ environment. GASP makes extensive use of pymatgen_, an open source Python library for materials analysis, and these instructions have been adapted from pymatgen. More details on installing pymatgen can be found at http://pymatgen.org/installation.html.

If pymatgen is already installed, steps 1-3 may be skipped.

.. _conda: http://conda.pydata.org/docs/index.html 
.. _pymatgen: http://pymatgen.org/

1. Install conda 
----------------

Download and install the version of conda for your operating system from http://conda.pydata.org/miniconda.html. Although GASP is compatible with both Python 2.7 and 3.5, pymatgen recommends using Python 3.5. 

For Windows, make sure you have the Miniconda3 installer, and simply double-click the .exe file. 

For Mac or Linux, run the bash script::  

    # if Mac
    bash Miniconda3-latest-MacOSX-x86_64.sh

    # if Linux
    bash Miniconda3-latest-Linux-x86_64.sh

The installer will ask you to approve the terms of the license and then tell you where the installation will be located. The default location is fine, so press Enter. Finally, the installer will ask if you would like it to add the install location to your PATH by prepending a line to your .bash_profile (for Mac) or .bashrc (for Linux) file. You do, so type 'yes' and press Enter.

After completing the installation, create a new terminal in order for the environment variables added by conda to be loaded.


2. Create a conda environment
-----------------------------

To create a new conda environment named 'my_gasp', type::

    conda create --name my_gasp python

When conda asks you::

    proceed ([y]/n)?

Type 'y' and press Enter.

Now activate the environment so that packages can be installed into it::

    # if Mac or Linux
    source activate my_gasp

    # if Windows
    activate my_gasp


3. Install pymatgen and its dependencies 
----------------------------------------

pymatgen recommends using the gcc compiler. To do so, type::

    export CC=gcc 

Install numpy, scipy, matplotlib and pymatgen with pip::

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install pymatgen

For Mac, an additional step is needed in order to use the scripts included with GASP for making plots. These scripts depend on the matplotlib_ library, which requires a framework build of Python to run properly on Mac OS X. However, a regular Python build comes with conda by default. To install a framework build in your conda environment, type::

    conda install python.app  

See the 'Visualizing output' section of the the `usage file`_ for more information on making plots.

.. _matplotlib: http://matplotlib.org/index.html 


4. Install GASP-python
----------------------

Clone the repository from github::

    git clone https://github.com/henniggroup/GASP-python.git

Move into the 'GASP-python' directory and run the setup.py script::

    cd GASP-python
    python setup.py develop


Using GASP
==========

See the `usage file`_.

.. _usage file: docs/usage.md


License
=======

GASP-python is released under the MIT License::

    Copyright (c) 2016-2017 Henniggroup Cornell/University of Florida

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Contributing
============

We try to follow the PEP8 coding style used by pymatgen: http://pymatgen.org/contributing.html#coding-guidelines

Authors
=======

Benjamin Revard

Richard G. Hennig    
