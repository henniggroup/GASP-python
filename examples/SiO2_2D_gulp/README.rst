Files to run a genetic algorithm structure search for 2D SiO2 with a maximum layer thickness of 5 Angstroms, using the potential of Vashishta *et al.* that ships with GULP.  


Description of files
====================

ga_input.yaml 

	Main input file for GASP. See the `usage file`_ for descriptions of each keyword and option. 

.. _usage file: ../../docs/usage.md


header_file 

	Header for GULP input files created by GASP.


potential_file 

	GULP file containing the empirical potential.


callgulp 

	Script called by GASP to run a GULP calculation.


Running the search
==================

Instructions for running this search on your system are given below. 

Note that we assume GASP-python and GULP are already installed on your system. If this is not the case, see the `main README file`_ and the `GULP documentation`_ for instructions on how to get GASP-python and GULP, respectively. 

.. _main README file: ../../README.rst
.. _GULP documentation: https://nanochemistry.curtin.edu.au/gulp/request.cfm?rel=download

1. Copy all files to the location on your computer where you would like to run the search.

2. Modify the location of the GULP header file given in ga_input.yaml (on line 6). Replace it with the location of header_file on your computer.  

3. Modify the location of the GULP potential file given in ga_input.yaml (on line 7). Replace it with the location of potential_file on your computer.  

4. Modify the location of the GULP binary given in callgulp (on line 15). Replace it with the location of the GULP binary on your computer. 

5. Move callgulp to a location in your system's PATH and make it executable.  

6. To start the search and save the output to a file called ga_output, move into the folder containing ga_input.yaml and type::

	run.py ga_input.yaml 2>&1 | tee ga_output
