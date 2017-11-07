Files to run a genetic algorithm search for the Al-ZrAl-CuAl partial phase diagram, using the EAM potentials of Cheng *et al.* and LAMMPS for the energy calculations.  


Description of files
====================

ga_input.yaml 

	Main input file for GASP. See the `usage file`_ for descriptions of each keyword and option. 

.. _usage file: ../../docs/usage.md


in.min 

	LAMMPS input script.


ZrCuAl.eam.alloy 

	LAMMPS file containing the EAM potentials.


ref_states/

	Directory containing the files for the structures at the endpoints of the composition space (Al, ZrAl, CuAl) in POSCAR format.


calllammps 

	Script called by GASP to run a LAMMPS calculation.


Running the search
==================

Instructions for running this search on your system are given below. 

Note that we assume GASP-python and LAMMPS are already installed on your system. If this is not the case, see the `main README file`_ and the `LAMMPS documentation`_ for instructions on how to get GASP-python and LAMMPS, respectively. 

.. _main README file: ../../README.rst
.. _LAMMPS documentation: http://lammps.sandia.gov/download.html 

1. Copy all files to the location on your computer where you would like to run the search.

2. Modify the location of the LAMMPS input script given in ga_input.yaml (on line 8). Replace it with the location of the in.min file on your computer.  

3. Modify the location of the ref_states directory given in ga_input.yaml (on line 12). Replace it with the location of the ref_states directory on your computer.

4. Modify the location of the LAMMPS potential file given in in.min (on line 8). Replace it with the location of the ZrCuAl.eam.alloy file on your computer. Note that 'Zr Al Cu' must appear at the end of the line after the location. 

5. Modify the location of the LAMMPS binary given in calllammps (on line 15). Replace it with the location of the LAMMPS binary on your computer. 

6. Move calllammps to a location in your system's PATH and make it executable.  

7. To start the search and save the output to a file called ga_output, move into the folder containing ga_input.yaml and type::

	run.py ga_input.yaml 2>&1 | tee ga_output

Shown below is the partial phase diagram GASP found when we ran this search.

.. image:: partial_phase_diagram.png
	:height: 283 px
	:width: 387 px
	:scale: 100 %
	:align: center
