# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import yaml


"""
Energy Calculators module:

This module contains the classes used to compute the energies of structures
with external energy codes. All energy calculator classes must implement a
do_energy_calculation() method.

1. VaspEnergyCalculator: for using VASP to compute energies

2. LammpsEnergyCalculator: for using LAMMSP to compute energies

3. GulpEnergyCalculator: for using GULP to compute energies

"""

from gasp.general import Cell

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.io.lammps.data import LammpsData, LammpsBox, ForceField, Topology
import pymatgen.command_line.gulp_caller as gulp_caller

import shutil
import subprocess
import os
import collections

from pymatgen.io.ase import AseAtomsAdaptor
from ase.calculators.calculator import kptdensity2monkhorstpack as kdens2mp
from ase.calculators.espresso import Espresso
from ase.io import read
import threading

import time

class QEEnergyCalculator(threading.Thread):
    """
    Calculates the energy of an organism using VASP.
    """

    def __init__(self, qe_calc_setting_yaml, slurm_script_template, ncore, geometry):
        '''
        Makes a QuantumEspressoEnergyCalculator.

        Args:
            pwi_file: the path to the pwi file

            geometry: the Geometry of the search
        '''
        super(QEEnergyCalculator, self).__init__()
        self.name = 'quantum_espresso'
        self.slurm_script_template = slurm_script_template
        self.ncore_base = ncore

        # paths to the INCAR, KPOINTS and POTCARs files
        self.AAA = AseAtomsAdaptor()
        self.qe_dft_calc = qe_calc_setting_yaml
        with open(self.qe_dft_calc,'r') as f:
            self.calc_config=yaml.safe_load(f)
        self.job_id = None

    def do_energy_calculation(self, organism, dictionary, key,
                              composition_space):
        """
        Calculates the energy of an organism using Quantum Espresso, and stores the relaxed
        organism in the provided dictionary at the provided key. If the
        calculation fails, stores None in the dictionary instead.

        Args:
            organism: the Organism whose energy we want to calculate

            dictionary: a dictionary in which to store the relaxed Organism

            key: the key specifying where to store the relaxed Organism in the
                dictionary

            composition_space: the CompositionSpace of the search

        Precondition: the garun directory and temp subdirectory exist, and we
            are currently located inside the garun directory

        TODO: maybe use the custodian package for error handling
        """
        self._stop_event = threading.Event()
        # make the job directory
        job_dir_path = str(os.getcwd()) + '/temp/' + str(organism.id)
        os.mkdir(job_dir_path)
        
        pwi_path = job_dir_path + '/' + str(organism.id)
        self.write_input_file(organism,pwi_path)

        shutil.copy(self.slurm_script_template, job_dir_path)
        slurm_script_name = os.path.basename(self.slurm_script_template)

        self.create_slurm_script(slurm_script_name,pwi_path)
        os.makedirs(job_dir_path+'/'+'outputs_slurm')

        # write out the unrelaxed structure to a poscar file
        organism.cell.to(fmt='poscar', filename=job_dir_path + '/POSCAR.' +
                         str(organism.id) + '_unrelaxed')

        # run 'callqe' script as a subprocess to run Quantum Espresso
        print('Starting quantum espresso calculation on organism {} '.format(organism.id))
        # devnull = open(os.devnull, 'w')
        # try:
        #     self.submit_slurm_job()
        #     subprocess.run(['callqe', pwi_path+'.pwi'],
        #                     )
        # except:
        #     print('Error submitting Quantum Espresso on organism {} '.format(organism.id))
        #     dictionary[key] = None
        #     return
        # parse the relaxed structure from the CONTCAR file
        #submit slurm job
        # current_dir = os.getcwd()
        # os.chdir(job_dir_path)
        result=subprocess.run(['sbatch', 'run.sh'],stdout=subprocess.PIPE,cwd=job_dir_path)
        self.job_id = result.stdout.decode().split()[-1]
        # os.chdir(current_dir)
        # running = True 
        # while running:
        #     err_file_path = job_dir_path+'/'+'outputs_slurm'+'/'+'err'
        #     if os.path.exists(err_file_path):
        #         with open(err_file_path,'r') as f:
        #             for line in f:
        #                 if 'Job' in line and 'Ended' in line:
        #                     running = False
        # Check the status regularly
        while not self._stop_event.is_set():
            status = self.check_job_status()
            if status != "RUNNING":
                break
            time.sleep(10)
            #need to manual cancel job from slurm if terminate run.py
        try:
            print(job_dir_path + '/' + str(organism.id) + '.pwo')
            ase_relaxed_cell = read(job_dir_path + '/' + str(organism.id) + '.pwo')
            relaxed_cell = self.AAA.get_structure(ase_relaxed_cell)
            #relaxed_cell = Cell.from_file(job_dir_path + '/CONTCAR')
        except:
            print('Error reading structure of organism {} from pwo file '.format(organism.id))
            dictionary[key] = None
            return

        # check if the qe calculation converged
        converged = False
        with open(job_dir_path + '/' + str(organism.id) + '.pwo') as f:
            for line in f:
                if 'JOB' in line and 'DONE.' in line:
                    converged = True
        if not converged:
            print('QE relaxation of organism {} did not converge '.format(
                organism.id))
            dictionary[key] = None
            return
        enthalpy = ase_relaxed_cell.get_potential_energy()
        organism.cell = relaxed_cell
        organism.total_energy = enthalpy
        organism.epa = enthalpy/organism.cell.num_sites
        print('Setting energy of organism {} to {} '
              'eV/atom '.format(organism.id, organism.epa))
        dictionary[key] = organism

    def check_job_status(self):
        # Check the status of the job using squeue
        result = subprocess.run(['squeue', '-j', self.job_id], stdout=subprocess.PIPE)
        lines = result.stdout.decode().split('\n')
        if len(lines) > 2: # Job is still in the queue
            return "RUNNING"
        else:
            return "FINISHED"

    def stop(self):
        self._stop_event.set()

    def create_slurm_script(self,slurm_script_name,pwi_path):
        dir_path = os.path.dirname(pwi_path)
        base_name = os.path.basename(pwi_path)
        
        if self.natoms > 0 and self.natoms <= 16:
            self.ncore = self.ncore_base * 1
        elif self.natoms > 16 and self.natoms < 48:
            self.ncore = self.ncore_base * 4
        else:
            self.ncore = self.ncore_base * 8
        # elif self.natoms > 64:
        #     self.ncore = self.ncore_base * 8

        # Read the content of the file
        with open(dir_path+'/'+slurm_script_name, 'r') as f:
            content = f.read()
        content = content.replace('100',str(self.ncore))

        # Write the modified content back to the file
        with open(dir_path+'/'+slurm_script_name, 'w') as f:
            f.write(content)
            f.write(f"cd {dir_path}\n\n")
            f.write(f"mpirun -np {self.ncore} $qe7_pw < {base_name}.pwi > {base_name}.pwo\n\n")
            f.write("echo ' '\n")
            f.write("echo 'Job Ended'")


    def write_input_file(self,organism,pwi_path):
        ase_atoms = self.AAA.get_atoms(organism.cell)
        self.natoms=len(ase_atoms)
        # calc data
        input_data = self.calc_config["input_data"]
        pseudopotentials = self.calc_config["pseudopotentials"]
        if self.calc_config["kpts"]:
            kpts = tuple(self.calc_config["kpts"])
        else:
            kpts = kdens2mp(ase_atoms, kptdensity=3.5, even=True)
        calc_obj = Espresso(
            input_data=input_data,
            pseudopotentials=pseudopotentials,
            kpts=kpts,
            label=str(pwi_path),
        )
        calc_obj.write_input(ase_atoms)

class VaspEnergyCalculator(object):
    """
    Calculates the energy of an organism using VASP.
    """

    def __init__(self, incar_file, kpoints_file, potcar_files, geometry):
        '''
        Makes a VaspEnergyCalculator.

        Args:
            incar_file: the path to the INCAR file

            kpoints_file: the path to the KPOINTS file

            potcar_files: a dictionary containing the paths to the POTCAR
                files, with the element symbols as keys

            geometry: the Geometry of the search
        '''

        self.name = 'vasp'

        # paths to the INCAR, KPOINTS and POTCARs files
        self.incar_file = incar_file
        self.kpoints_file = kpoints_file
        self.potcar_files = potcar_files

    def do_energy_calculation(self, organism, dictionary, key,
                              composition_space):
        """
        Calculates the energy of an organism using VASP, and stores the relaxed
        organism in the provided dictionary at the provided key. If the
        calculation fails, stores None in the dictionary instead.

        Args:
            organism: the Organism whose energy we want to calculate

            dictionary: a dictionary in which to store the relaxed Organism

            key: the key specifying where to store the relaxed Organism in the
                dictionary

            composition_space: the CompositionSpace of the search

        Precondition: the garun directory and temp subdirectory exist, and we
            are currently located inside the garun directory

        TODO: maybe use the custodian package for error handling
        """

        # make the job directory
        job_dir_path = str(os.getcwd()) + '/temp/' + str(organism.id)
        os.mkdir(job_dir_path)

        # copy the INCAR and KPOINTS files to the job directory
        shutil.copy(self.incar_file, job_dir_path)
        shutil.copy(self.kpoints_file, job_dir_path)

        # sort the organism's cell and write to POSCAR file
        organism.cell.sort()
        organism.cell.to(fmt='poscar', filename=job_dir_path + '/POSCAR')

        # get a list of the element symbols in the sorted order
        symbols = []
        for site in organism.cell.sites:
            if site.specie.symbol not in symbols:
                symbols.append(site.specie.symbol)

        # write the POTCAR file by concatenating the appropriate elemental
        # POTCAR files
        total_potcar_path = job_dir_path + '/POTCAR'
        with open(total_potcar_path, 'w') as total_potcar_file:
            for symbol in symbols:
                with open(self.potcar_files[symbol], 'r') as potcar_file:
                    for line in potcar_file:
                        total_potcar_file.write(line)

        # run 'callvasp' script as a subprocess to run VASP
        print('Starting VASP calculation on organism {} '.format(organism.id))
        devnull = open(os.devnull, 'w')
        try:
            subprocess.call(['callvasp', job_dir_path], stdout=devnull,
                            stderr=devnull)
        except:
            print('Error running VASP on organism {} '.format(organism.id))
            dictionary[key] = None
            return

        # parse the relaxed structure from the CONTCAR file
        try:
            relaxed_cell = Cell.from_file(job_dir_path + '/CONTCAR')
        except:
            print('Error reading structure of organism {} from CONTCAR '
                  'file '.format(organism.id))
            dictionary[key] = None
            return

        # check if the VASP calculation converged
        converged = False
        with open(job_dir_path + '/OUTCAR') as f:
            for line in f:
                if 'reached' in line and 'required' in line and \
                        'accuracy' in line:
                    converged = True
        if not converged:
            print('VASP relaxation of organism {} did not converge '.format(
                organism.id))
            dictionary[key] = None
            return

        # parse the internal energy and pV (if needed) and compute the enthalpy
        pv = 0
        with open(job_dir_path + '/OUTCAR') as f:
            for line in f:
                if 'energy(sigma->0)' in line:
                    u = float(line.split()[-1])
                elif 'enthalpy' in line:
                    pv = float(line.split()[-1])
        enthalpy = u + pv

        organism.cell = relaxed_cell
        organism.total_energy = enthalpy
        organism.epa = enthalpy/organism.cell.num_sites
        print('Setting energy of organism {} to {} '
              'eV/atom '.format(organism.id, organism.epa))
        dictionary[key] = organism


class LammpsEnergyCalculator(object):
    """
    Calculates the energy of an organism using LAMMPS.
    """

    def __init__(self, input_script, geometry):
        """
        Makes a LammpsEnergyCalculator.

        Args:
            input_script: the path to the lammps input script

            geometry: the Geometry of the search

        Precondition: the input script exists and is valid
        """

        self.name = 'lammps'

        # the path to the lammps input script
        self.input_script = input_script

    def do_energy_calculation(self, organism, dictionary, key,
                              composition_space):
        """
        Calculates the energy of an organism using LAMMPS, and stores the
        relaxed organism in the provided dictionary at the provided key. If the
        calculation fails, stores None in the dictionary instead.

        Args:
            organism: the Organism whose energy we want to calculate

            dictionary: a dictionary in which to store the relaxed Organism

            key: the key specifying where to store the relaxed Organism in the
                dictionary

            composition_space: the CompositionSpace of the search

        Precondition: the garun directory and temp subdirectory exist, and we
            are currently located inside the garun directory
        """

        # make the job directory
        job_dir_path = str(os.getcwd()) + '/temp/' + str(organism.id)
        os.mkdir(job_dir_path)

        # copy the lammps input script to the job directory and get its path
        shutil.copy(self.input_script, job_dir_path)
        script_name = os.path.basename(self.input_script)
        input_script_path = job_dir_path + '/' + str(script_name)

        # write the in.data file
        self.conform_to_lammps(organism.cell)
        self.write_data_file(organism, job_dir_path, composition_space)

        # write out the unrelaxed structure to a poscar file
        organism.cell.to(fmt='poscar', filename=job_dir_path + '/POSCAR.' +
                         str(organism.id) + '_unrelaxed')

        # run 'calllammps' script as a subprocess to run LAMMPS
        print('Starting LAMMPS calculation on organism {} '.format(
            organism.id))
        try:
            lammps_output = subprocess.check_output(
                ['calllammps', input_script_path], stderr=subprocess.STDOUT)
            # convert from bytes to string (for Python 3)
            lammps_output = lammps_output.decode('utf-8')
        except subprocess.CalledProcessError as e:
            # write the output of a bad LAMMPS call to for the user's reference
            with open(job_dir_path + '/log.lammps', 'w') as log_file:
                log_file.write(e.output.decode('utf-8'))
            print('Error running LAMMPS on organism {} '.format(organism.id))
            dictionary[key] = None
            return

        # write the LAMMPS output
        with open(job_dir_path + '/log.lammps', 'w') as log_file:
            log_file.write(lammps_output)

        # parse the relaxed structure from the atom.dump file
        symbols = []
        all_elements = composition_space.get_all_elements()
        for element in all_elements:
            symbols.append(element.symbol)
        try:
            relaxed_cell = self.get_relaxed_cell(
                job_dir_path + '/dump.atom', job_dir_path + '/in.data',
                symbols)
        except:
            print('Error reading structure of organism {} from LAMMPS '
                  'output '.format(organism.id))
            dictionary[key] = None
            return

        # parse the total energy from the log.lammps file
        try:
            total_energy = self.get_energy(job_dir_path + '/log.lammps')
        except:
            print('Error reading energy of organism {} from LAMMPS '
                  'output '.format(organism.id))
            dictionary[key] = None
            return

        # check that the total energy isn't unphysically large
        # (can be a problem for empirical potentials)
        epa = total_energy/organism.cell.num_sites
        if epa < -50:
            print('Discarding organism {} due to unphysically large energy: '
                  '{} eV/atom.'.format(organism.id, str(epa)))
            dictionary[key] = None
            return

        organism.cell = relaxed_cell
        organism.total_energy = total_energy
        organism.epa = epa
        print('Setting energy of organism {} to {} eV/atom '.format(
            organism.id, organism.epa))
        dictionary[key] = organism

    def conform_to_lammps(self, cell):
        """
        Modifies a cell to satisfy the requirements of lammps, which are:

            1. the lattice vectors lie in the principal directions

            2. the maximum extent in the Cartesian x-direction is achieved by
                lattice vector a

            3. the maximum extent in the Cartesian y-direction is achieved by
                lattice vector b

            4. the maximum extent in the Cartesian z-direction is achieved by
                lattice vector c

        by taking supercells along lattice vectors when needed.

        Args:
            cell: the Cell to modify
        """

        cell.rotate_to_principal_directions()
        lattice_coords = cell.lattice.matrix
        ax = lattice_coords[0][0]
        bx = lattice_coords[1][0]
        cx = lattice_coords[2][0]
        by = lattice_coords[1][1]
        cy = lattice_coords[2][1]
        if ax < bx or ax < cx:
            cell.make_supercell([2, 1, 1])
            self.conform_to_lammps(cell)
        elif by < cy:
            cell.make_supercell([1, 2, 1])
            self.conform_to_lammps(cell)

    def write_data_file(self, organism, job_dir_path, composition_space):
        """
        Writes the file (called in.data) containing the structure that LAMMPS
        reads.

        Args:
            organism: the Organism whose structure to write

            job_dir_path: path to the job directory (as a string) where the
                file will be written

            composition_space: the CompositionSpace of the search
        """

        # get xhi, yhi and zhi from the lattice vectors
        lattice_coords = organism.cell.lattice.matrix
        xhi = lattice_coords[0][0]
        yhi = lattice_coords[1][1]
        zhi = lattice_coords[2][2]
        box_bounds = [[0.0, xhi], [0.0, yhi], [0.0, zhi]]

        # get xy, xz and yz from the lattice vectors
        xy = lattice_coords[1][0]
        xz = lattice_coords[2][0]
        yz = lattice_coords[2][1]
        box_tilts = [xy, xz, yz]

        # make a LammpsBox object from the bounds and tilts
        lammps_box = LammpsBox(box_bounds, tilt=box_tilts)

        # parse the element symbols and atom_style from the lammps input
        # script, preserving the order in which the element symbols appear
        # TODO: not all formats give the element symbols at the end of the line
        #       containing the pair_coeff keyword. Find a better way.
        elements_dict = collections.OrderedDict()
        num_elements = len(composition_space.get_all_elements())

        is_single_element = (num_elements == 1)
        if is_single_element:
            single_element = composition_space.get_all_elements()
            elements_dict[single_element[0].symbol] = single_element[0]

        with open(self.input_script, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'atom_style' in line:
                    atom_style_in_script = line.split()[1]
                elif not is_single_element and 'pair_coeff' in line:
                    element_symbols = line.split()[-1*num_elements:]

        if not is_single_element:
            for symbol in element_symbols:
                elements_dict[symbol] = Element(symbol)

        # make a LammpsData object and use it write the in.data file
        force_field = ForceField(elements_dict.items())
        topology = Topology(organism.cell.sites)
        lammps_data = LammpsData.from_ff_and_topologies(
            lammps_box, force_field, [topology],
            atom_style=atom_style_in_script)
        lammps_data.write_file(job_dir_path + '/in.data')

    def get_relaxed_cell(self, atom_dump_path, data_in_path, element_symbols):
        """
        Parses the relaxed cell from the dump.atom file.

        Returns the relaxed cell as a Cell object.

        Args:
            atom_dump_path: the path (as a string) to the dump.atom file

            in_data_path: the path (as a string) to the in.data file

            element_symbols: a tuple containing the set of chemical symbols of
                all the elements in the compositions space
        """

        # read the dump.atom file as a list of strings
        with open(atom_dump_path, 'r') as atom_dump:
            lines = atom_dump.readlines()

        # get the lattice vectors
        a_data = lines[5].split()
        b_data = lines[6].split()
        c_data = lines[7].split()

        # parse the tilt factors
        xy = float(a_data[2])
        xz = float(b_data[2])
        yz = float(c_data[2])

        # parse the bounds
        xlo_bound = float(a_data[0])
        xhi_bound = float(a_data[1])
        ylo_bound = float(b_data[0])
        yhi_bound = float(b_data[1])
        zlo_bound = float(c_data[0])
        zhi_bound = float(c_data[1])

        # compute xlo, xhi, ylo, yhi, zlo and zhi according to the conversion
        # given by LAMMPS
        # http://lammps.sandia.gov/doc/Section_howto.html#howto-12
        xlo = xlo_bound - min([0.0, xy, xz, xy + xz])
        xhi = xhi_bound - max([0.0, xy, xz, xy + xz])
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max([0.0, yz])
        zlo = zlo_bound
        zhi = zhi_bound

        # construct a Lattice object from the lo's and hi's and tilts
        a = [xhi - xlo, 0.0, 0.0]
        b = [xy, yhi - ylo, 0.0]
        c = [xz, yz, zhi - zlo]
        relaxed_lattice = Lattice([a, b, c])

        # get the number of atoms
        num_atoms = int(lines[3])

        # get the atom types and their Cartesian coordinates
        types = []
        relaxed_cart_coords = []
        for i in range(num_atoms):
            atom_info = lines[9 + i].split()
            types.append(int(atom_info[1]))
            relaxed_cart_coords.append([float(atom_info[2]) - xlo,
                                        float(atom_info[3]) - ylo,
                                        float(atom_info[4]) - zlo])

        # read the atom types and corresponding atomic masses from in.data
        with open(data_in_path, 'r') as data_in:
            lines = data_in.readlines()
        types_masses = {}
        for i in range(len(lines)):
            if 'Masses' in lines[i]:
                for j in range(len(element_symbols)):
                    types_masses[int(lines[i + j + 2].split()[0])] = float(
                        lines[i + j + 2].split()[1])

        # map the atom types to chemical symbols
        types_symbols = {}
        for symbol in element_symbols:
            for atom_type in types_masses:
                # round the atomic masses to one decimal point for comparison
                if format(float(Element(symbol).atomic_mass), '.1f') == format(
                        types_masses[atom_type], '.1f'):
                    types_symbols[atom_type] = symbol

        # make a list of chemical symbols (one for each site)
        relaxed_symbols = []
        for atom_type in types:
            relaxed_symbols.append(types_symbols[atom_type])

        return Cell(relaxed_lattice, relaxed_symbols, relaxed_cart_coords,
                    coords_are_cartesian=True)

    def get_energy(self, lammps_log_path):
        """
        Parses the final energy from the log.lammps file written by LAMMPS.

        Returns the total energy as a float.

        Args:
            lammps_log_path: the path (as a string) to the log.lammps file
        """

        # read the log.lammps file as a list of strings
        with open(lammps_log_path, 'r') as f:
            lines = f.readlines()

        # get the last line with the keywords (where the final energy is)
        match_strings = ['Step', 'Temp', 'E_pair', 'E_mol', 'TotEng']
        for i in range(len(lines)):
            if all(match in lines[i] for match in match_strings):
                energy = float(lines[i + 2].split()[4])
        return energy


class GulpEnergyCalculator(object):
    """
    Calculates the energy of an organism using GULP.
    """

    def __init__(self, header_file, potential_file, geometry):
        """
        Makes a GulpEnergyCalculator.

        Args:
            header_file: the path to the gulp header file

            potential_file: the path to the gulp potential file

            geometry: the Geometry of the search

        Precondition: the header and potential files exist and are valid
        """

        self.name = 'gulp'

        # the paths to the header and potential files
        self.header_path = header_file
        self.potential_path = potential_file

        # read the gulp header and potential files
        with open(header_file, 'r') as gulp_header_file:
            self.header = gulp_header_file.readlines()
        with open(potential_file, 'r') as gulp_potential_file:
            self.potential = gulp_potential_file.readlines()

        # for processing gulp input and output
        self.gulp_io = gulp_caller.GulpIO()

        # whether the anions and cations are polarizable in the gulp potential
        self.anions_shell, self.cations_shell = self.get_shells()

        # determine which lattice parameters should be relaxed
        # and make the corresponding flags for the input file
        #
        # relax a, b, c, alpha, beta, gamma
        if geometry.shape == 'bulk':
            self.lattice_flags = None
        # relax a, b and gamma but not c, alpha and beta
        elif geometry.shape == 'sheet':
            self.lattice_flags = '1 1 0 0 0 1'
        # relax c, but not a, b, alpha, beta and gamma
        elif geometry.shape == 'wire':
            self.lattice_flags = '0 0 1 0 0 0'
        # don't relax any of the lattice parameters
        elif geometry.shape == 'cluster':
            self.lattice_flags = '0 0 0 0 0 0'

    def get_shells(self):
        """
        Determines whether the anions and cations have shells by looking at the
        potential file.

        Returns two booleans indicating whether the anions and cations have
        shells, respectively.
        """

        # get the symbols of the elements with shells
        shells = []
        for line in self.potential:
            if 'shel' in line:
                line_parts = line.split()
                shells.append(str(line_parts[line_parts.index('shel') - 1]))
        shells = list(set(shells))

        # determine whether the elements with shells are anions and/or cations
        anions_shell = False
        cations_shell = False
        for symbol in shells:
            element = Element(symbol)
            if element in gulp_caller._anions:
                anions_shell = True
            elif element in gulp_caller._cations:
                cations_shell = True
        return anions_shell, cations_shell

    def do_energy_calculation(self, organism, dictionary, key,
                              composition_space):
        """
        Calculates the energy of an organism using GULP, and stores the relaxed
        organism in the provided dictionary at the provided key. If the
        calculation fails, stores None in the dictionary instead.

        Args:
            organism: the Organism whose energy we want to calculate

            dictionary: a dictionary in which to store the relaxed Organism

            key: the key specifying where to store the relaxed Organism in the
                dictionary

            composition_space: the CompositionSpace of the search

        Precondition: the garun directory and temp subdirectory exist, and we
            are currently located inside the garun directory

        TODO: maybe use the custodian package for error handling
        """

        # make the job directory
        job_dir_path = str(os.getcwd()) + '/temp/' + str(organism.id)
        os.mkdir(job_dir_path)

        # just for testing, write out the unrelaxed structure to a poscar file
        # organism.cell.to(fmt='poscar', filename= job_dir_path +
        #    '/POSCAR.' + str(organism.id) + '_unrelaxed')

        # write the GULP input file
        gin_path = job_dir_path + '/' + str(organism.id) + '.gin'
        self.write_input_file(organism, gin_path)

        # run 'calllgulp' script as a subprocess to run GULP
        print('Starting GULP calculation on organism {} '.format(organism.id))
        try:
            gulp_output = subprocess.check_output(['callgulp', gin_path],
                                                  stderr=subprocess.STDOUT)
            # convert from bytes to string (for Python 3)
            gulp_output = gulp_output.decode('utf-8')
        except subprocess.CalledProcessError as e:
            # write the output of a bad GULP call to for the user's reference
            with open(job_dir_path + '/' + str(organism.id) + '.gout',
                      'w') as gout_file:
                gout_file.write(e.output.decode('utf-8'))
            print('Error running GULP on organism {} '.format(organism.id))
            dictionary[key] = None
            return

        # write the GULP output for the user's reference
        with open(job_dir_path + '/' + str(organism.id) + '.gout',
                  'w') as gout_file:
            gout_file.write(gulp_output)

        # check if not converged (part of this is copied from pymatgen)
        conv_err_string = 'Conditions for a minimum have not been satisfied'
        gradient_norm = self.get_grad_norm(gulp_output)
        if conv_err_string in gulp_output and gradient_norm > 0.1:
            print('The GULP calculation on organism {} did not '
                  'converge '.format(organism.id))
            dictionary[key] = None
            return

        # parse the relaxed structure from the gulp output
        try:
            # TODO: change this line if pymatgen fixes the gulp parser
            relaxed_cell = self.get_relaxed_cell(gulp_output)
        except:
            print('Error reading structure of organism {} from GULP '
                  'output '.format(organism.id))
            dictionary[key] = None
            return

        # parse the total energy from the gulp output
        try:
            total_energy = self.get_energy(gulp_output)
        except:
            print('Error reading energy of organism {} from GULP '
                  'output '.format(organism.id))
            dictionary[key] = None
            return

        # sometimes gulp takes a supercell
        num_atoms = self.get_num_atoms(gulp_output)

        organism.cell = relaxed_cell
        organism.epa = total_energy/num_atoms
        organism.total_energy = organism.epa*organism.cell.num_sites
        print('Setting energy of organism {} to {} eV/atom '.format(
            organism.id, organism.epa))
        dictionary[key] = organism

    def write_input_file(self, organism, gin_path):
        """
        Writes the gulp input file.

        Args:
            organism: the Organism whose energy we want to calculate

            gin_path: the path to the GULP input file
        """

        # get the structure lines
        structure_lines = self.gulp_io.structure_lines(
            organism.cell, anion_shell_flg=self.anions_shell,
            cation_shell_flg=self.cations_shell, symm_flg=False)
        structure_lines = structure_lines.split('\n')
        del structure_lines[-1]  # remove empty line that gets added

        # GULP errors out if too many decimal places in lattice parameters
        lattice_parameters = structure_lines[1].split()
        for i in range(len(lattice_parameters)):
            lattice_parameters[i] = round(float(lattice_parameters[i]), 10)

        rounded_lattice_parameters = ""
        for lattice_parameter in lattice_parameters:
            rounded_lattice_parameters += str(lattice_parameter)
            rounded_lattice_parameters += " "
        structure_lines[1] = rounded_lattice_parameters

        # add flags for relaxing lattice parameters and ion positions
        if self.lattice_flags is not None:
            structure_lines[1] = structure_lines[1] + self.lattice_flags
            for i in range(3, len(structure_lines)):
                structure_lines[i] = structure_lines[i] + ' 1 1 1'

        # add newline characters to the end of each of the structure lines
        for i in range(len(structure_lines)):
            structure_lines[i] = structure_lines[i] + '\n'

        # construct complete input
        gulp_input = self.header + structure_lines + self.potential

        # print gulp input to a file
        with open(gin_path, 'w') as gin_file:
            for line in gulp_input:
                gin_file.write(line)

    def get_grad_norm(self, gout):
        """
        Parses the final gradient norm from the GULP output.

        Args:
            gout: the GULP output, as a string
        """

        output_lines = gout.split('\n')
        for line in output_lines:
            if 'Final Gnorm' in line:
                line_parts = line.split()
                return float(line_parts[3])

    def get_energy(self, gout):
        """
        Parses the final energy from the GULP output.

        Args:
            gout: the GULP output, as a string
        """

        output_lines = gout.split('\n')
        for line in output_lines:
            if 'Final energy' in line:
                return float(line.split()[3])

    def get_num_atoms(self, gout):
        """
        Parses the number of atoms used by GULP in the calculation.

        Args:
            gout: the GULP output, as a string
        """

        output_lines = gout.split('\n')
        for line in output_lines:
            if 'Total number atoms' in line:
                line_parts = line.split()
                return int(line_parts[-1])

    # This method is copied from GulpIO.get_relaxed_structure, and I modified
    # it slightly to get it to work.
    # TODO: if pymatgen fixes this method, then I can delete this.
    # Alternatively, could submit a pull request with my fix
    def get_relaxed_cell(self, gout):
        # Find the structure lines
        structure_lines = []
        cell_param_lines = []
        output_lines = gout.split("\n")
        no_lines = len(output_lines)
        i = 0
        # Compute the input lattice parameters
        while i < no_lines:
            line = output_lines[i]
            if "Full cell parameters" in line:
                i += 2
                line = output_lines[i]
                a = float(line.split()[8])
                alpha = float(line.split()[11])
                line = output_lines[i + 1]
                b = float(line.split()[8])
                beta = float(line.split()[11])
                line = output_lines[i + 2]
                c = float(line.split()[8])
                gamma = float(line.split()[11])
                i += 3
                break
            elif "Cell parameters" in line:
                i += 2
                line = output_lines[i]
                a = float(line.split()[2])
                alpha = float(line.split()[5])
                line = output_lines[i + 1]
                b = float(line.split()[2])
                beta = float(line.split()[5])
                line = output_lines[i + 2]
                c = float(line.split()[2])
                gamma = float(line.split()[5])
                i += 3
                break
            else:
                i += 1

        while i < no_lines:
            line = output_lines[i]
            if "Final fractional coordinates of atoms" in line or \
                    "Final asymmetric unit coordinates" in line:  # Ben's add
                # read the site coordinates in the following lines
                i += 6
                line = output_lines[i]
                while line[0:2] != '--':
                    structure_lines.append(line)
                    i += 1
                    line = output_lines[i]
                    # read the cell parameters
                i += 9
                line = output_lines[i]
                if "Final cell parameters" in line:
                    i += 3
                    for del_i in range(6):
                        line = output_lines[i + del_i]
                        cell_param_lines.append(line)
                break
            else:
                i += 1

        # Process the structure lines
        if structure_lines:
            sp = []
            coords = []
            for line in structure_lines:
                fields = line.split()
                if fields[2] == 'c':
                    sp.append(fields[1])
                    coords.append(list(float(x) for x in fields[3:6]))
        else:
            raise IOError("No structure found")

        if cell_param_lines:
            a = float(cell_param_lines[0].split()[1])
            b = float(cell_param_lines[1].split()[1])
            c = float(cell_param_lines[2].split()[1])
            alpha = float(cell_param_lines[3].split()[1])
            beta = float(cell_param_lines[4].split()[1])
            gamma = float(cell_param_lines[5].split()[1])
        latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        return Cell(latt, sp, coords)
