# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Organism Creators module:

This module contains the classes used to create organisms for the initial
population.

1. RandomOrganismCreator: creates random organisms

2. FileOrganismCreator: creates organisms by reading their structures from
        files

"""
from gasp.general import Organism, Cell

from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition

from fractions import Fraction
import warnings
import os
import math


class RandomOrganismCreator(object):
    """
    Creates random organisms for the initial population.
    """

    def __init__(self, random_org_parameters, composition_space, constraints):
        """
        Makes a RandomOrganismCreator, and sets default parameter values if
        necessary.

        Args:
            random_org_parameters: the parameters for generating random
                organisms

            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search
        """

        self.name = 'random organism creator'

        # defaults
        #
        # number of random organisms to make (only used for epa searches)
        self.default_number = 28
        # max number of atoms
        self.default_max_num_atoms = min(constraints.min_num_atoms + 6,
                                         constraints.max_num_atoms)
        # volume scaling behavior
        # default volumes per atom of elemental ground state structures
        # computed from structures on materials project (materialsproject.org)
        self.all_default_vpas = {'H': 13.89, 'He': 15.79, 'Li': 20.12,
                                 'Be': 7.94, 'C': 10.58, 'N': 42.73,
                                 'O': 13.46, 'F': 16.00, 'Ne': 19.93,
                                 'Na': 37.12, 'Mg': 23.04, 'Al': 16.47,
                                 'Si': 20.44, 'P': 23.93, 'S': 36.03,
                                 'Cl': 34.90, 'Ar': 44.87, 'K': 73.51,
                                 'Ca': 42.42, 'Sc': 24.64, 'Ti': 17.11,
                                 'V': 13.41, 'Cr': 11.57, 'Mn': 11.04,
                                 'Fe': 11.55, 'Co': 10.92, 'Ni': 10.79,
                                 'Cu': 11.82, 'Zn': 15.56, 'Ga': 20.34,
                                 'Ge': 23.92, 'As': 22.45, 'Se': 38.13,
                                 'Br': 37.53, 'Kr': 65.09, 'Rb': 90.44,
                                 'Sr': 54.88, 'Y': 32.85, 'Zr': 23.50,
                                 'Nb': 18.31, 'Mo': 15.89, 'Tc': 14.59,
                                 'Ru': 13.94, 'Rh': 14.25, 'Pd': 15.45,
                                 'Ag': 18.00, 'Cd': 23.28, 'In': 27.56,
                                 'Sn': 36.70, 'Sb': 31.78, 'Te': 35.03,
                                 'I': 50.34, 'Xe': 83.51, 'Cs': 116.17,
                                 'Ba': 63.64, 'Hf': 22.50, 'Ta': 18.25,
                                 'W': 16.19, 'Re': 15.06, 'Os': 14.36,
                                 'Ir': 14.55, 'Pt': 15.72, 'Au': 18.14,
                                 'Hg': 31.45, 'Tl': 31.13, 'Pb': 32.30,
                                 'Bi': 36.60, 'La': 37.15, 'Ce': 26.30,
                                 'Pr': 36.47, 'Nd': 35.44, 'Pm': 34.58,
                                 'Sm': 33.88, 'Eu': 46.28, 'Gd': 33.33,
                                 'Tb': 32.09, 'Dy': 31.57, 'Ho': 31.45,
                                 'Er': 30.90, 'Tm': 30.30, 'Yb': 40.45,
                                 'Lu': 29.43, 'Ac': 45.52, 'Th': 32.03,
                                 'Pa': 25.21, 'U': 19.98, 'Np': 18.43,
                                 'Pu': 18.34}

        self.default_vpas = self.get_default_vpas(composition_space)

        # set to defaults
        if random_org_parameters in (None, 'default'):
            self.number = self.default_number
            self.max_num_atoms = self.default_max_num_atoms
            self.vpas = self.default_vpas
        # parse the parameters and set to defaults if necessary
        else:
            # the number to make
            if 'number' not in random_org_parameters:
                self.number = self.default_number
            elif random_org_parameters['number'] in (None, 'default'):
                self.number = self.default_number
            else:
                self.number = random_org_parameters['number']

            # the max number of atoms
            if 'max_num_atoms' not in random_org_parameters:
                self.max_num_atoms = self.default_max_num_atoms
            elif random_org_parameters['max_num_atoms'] in (None, 'default'):
                self.max_num_atoms = self.default_max_num_atoms
            elif random_org_parameters['max_num_atoms'] > \
                    constraints.max_num_atoms:
                print('The value passed to the "max_num_atoms" keyword in the '
                      'InitialPopulation block may not exceed the value passed'
                      ' to the "max_num_atoms" keyword in the Constraints '
                      'block.')
                print('Quitting...')
                quit()
            elif random_org_parameters['max_num_atoms'] < \
                    constraints.min_num_atoms:
                print('The value passed to the "max_num_atoms" keyword in the '
                      'InitialPopulation block may not be smaller than the '
                      'value passed to the "min_num_atoms" keyword in the '
                      'Constraints block.')
                print('Quitting...')
                quit()
            else:
                self.max_num_atoms = random_org_parameters['max_num_atoms']

            # volume scaling
            self.vpas = self.default_vpas
            if 'volumes_per_atom' not in random_org_parameters:
                pass
            elif random_org_parameters['volumes_per_atom'] in (None,
                                                               'default'):
                pass
            else:
                # replace the specified volumes per atom with the given values
                for symbol in random_org_parameters['volumes_per_atom']:
                    self.vpas[symbol] = random_org_parameters[
                        'volumes_per_atom'][symbol]

        self.num_made = 0  # number added to initial population
        self.is_successes_based = True  # it's based on number added
        self.is_finished = False

    def get_default_vpas(self, composition_space):
        """
        Returns a dictionary containing the default volumes per atom for all
        the elements in the composition space.

        Args:
            composition_space: the CompositionSpace of the search
        """

        default_vpas = {}
        for element in composition_space.get_all_elements():
            default_vpas[element.symbol] = self.all_default_vpas[
                element.symbol]
        return default_vpas

    def create_organism(self, id_generator, composition_space, constraints,
                        random):
        """
        Creates a random organism for the initial population.

        Returns a random organism, or None if an error was encountered during
        volume scaling.

        Note: for phase diagram searches, this is will not create structures
            with compositions equivalent to the endpoints of the composition
            space. Reference structures at those compositions should be
            provided with the FileOrganismCreator.

        Args:
            id_generator: the IDGenerator used to assign id numbers to all
                organisms

            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG
        """

        # make a random lattice
        random_lattice = self.make_random_lattice(constraints, random)

        # get a list of species for the random organism
        species = self.get_species_list(composition_space, constraints, random)
        if species is None:  # could happen for pd searches...
            return None

        # for each specie, generate a set of random fractional coordinates
        random_coordinates = []
        for _ in range(len(species)):
            random_coordinates.append([random.random(), random.random(),
                                       random.random()])

        # make a random cell
        random_cell = Cell(random_lattice, species, random_coordinates)

        # optionally scale the volume of the random structure
        if not self.scale_volume(random_cell):
            return None  # sometimes pymatgen's scaling algorithm crashes

        # make the random organism
        random_org = Organism(random_cell, id_generator, self.name,
                              composition_space)
        print('Random organism creator making organism {} '.format(
            random_org.id))
        return random_org

    def make_random_lattice(self, constraints, random):
        """
        Returns a random lattice that satisfies the constraints on maximum and
        minimum lengths and angles.

        Args:
            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG
        """

        # make three random lattice vectors that satisfy the length constraints
        a = constraints.min_lattice_length + random.random()*(
            constraints.max_lattice_length - constraints.min_lattice_length)
        b = constraints.min_lattice_length + random.random()*(
            constraints.max_lattice_length - constraints.min_lattice_length)
        c = constraints.min_lattice_length + random.random()*(
            constraints.max_lattice_length - constraints.min_lattice_length)

        # make three random lattice angles that satisfy the angle constraints
        alpha = constraints.min_lattice_angle + random.random()*(
            constraints.max_lattice_angle - constraints.min_lattice_angle)
        beta = constraints.min_lattice_angle + random.random()*(
            constraints.max_lattice_angle - constraints.min_lattice_angle)
        gamma = constraints.min_lattice_angle + random.random()*(
            constraints.max_lattice_angle - constraints.min_lattice_angle)

        # build the random lattice
        return Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    def get_species_list(self, composition_space, constraints, random):
        """
        Returns a list containing the species in the random organism.

        Args:
            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG
        """
        if composition_space.objective_function == 'epa':
            return self.get_epa_species_list(composition_space, constraints,
                                             random)
        elif composition_space.objective_function == 'pd':
            return self.get_pd_species_list(composition_space, constraints,
                                            random)

    def get_epa_species_list(self, composition_space, constraints, random):
        """
        Returns a list containing the species in the random organism.

        Precondition: the composition space contains only one endpoint
            (it's a fixed-composition search)

        Args:
            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG

        Description:

            1. Computes the minimum and maximum number of formula units from
                the minimum (constraints.min_num_atoms) and maximum
                (self.max_num_atoms) number of atoms and the number of atoms
                per formula unit.

            2. Gets a random number of formula units within the range allowed
                by the minimum and maximum number of formula units.

            3. Computes the number of atoms of each species from the random
                number of formula units.
        """
        # get random number of formula units and resulting number of atoms
        reduced_formula = composition_space.endpoints[0].reduced_composition
        num_atoms_in_formula = reduced_formula.num_atoms
        max_num_formulas = int(math.floor(
            self.max_num_atoms/num_atoms_in_formula))
        min_num_formulas = int(math.ceil(
            constraints.min_num_atoms/num_atoms_in_formula))
        # round up the next formula unit if necessary
        if max_num_formulas < min_num_formulas:
            max_num_formulas += 1
        random_num_formulas = random.randint(min_num_formulas,
                                             max_num_formulas)

        # add the right number of each specie
        species = []
        for specie in reduced_formula:
            for _ in range(random_num_formulas*int(reduced_formula[specie])):
                species.append(specie)
        return species

    def get_pd_species_list(self, composition_space, constraints, random):
        """
        Returns a list containing the species in the random organism.

        Precondition: the composition space contains multiple endpoints
            (it's a fixed-composition search)

        Args:
            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG

        Description:

            1. Gets a random fraction of each composition space endpoint such
                that the fractions sum to 1.

            2. Computes the fraction of each specie from the fraction of each
                endpoint and the amount of each specie within each endpoint.

            3. Approximates the fraction of each specie as a rational number
                with a maximum possible denominator of self.max_num_atoms.

            4. Takes the product of the denominators of all the species'
                rational fractions, and then multiplies each specie's rational
                fraction by this product to obtain the number of atoms of that
                species.

            5. Reduces the composition obtained in 4. to the smallest possible
                number of atoms.

            6. Checks that the resulting number of atoms satisfies the maximum
                (self.max_num_atoms) and minimum (constraints.min_num_atoms)
                number of atoms constraints, and that the resulting composition
                is not equivalent to one of the endpoint compositions.
        """

        # get random fractions for each endpoint that sum to 1 (i.e., a
        # random location in the composition space)
        random.shuffle(composition_space.endpoints)  # to remove bias
        frac_sum = 0
        endpoint_fracs = {}
        for i in range(len(composition_space.endpoints) - 1):
            next_frac = random.uniform(0, 1.0 - frac_sum)
            endpoint_fracs[composition_space.endpoints[i]] = next_frac
            frac_sum += next_frac
        endpoint_fracs[composition_space.endpoints[-1]] = 1.0 - frac_sum

        # compute amount of each element from amount of each endpoint
        all_elements = composition_space.get_all_elements()
        element_amounts = {}
        for element in all_elements:
            element_amounts[element] = 0
        for formula in endpoint_fracs:
            for element in formula:
                element_amounts[element] += endpoint_fracs[
                    formula]*formula[element]

        # normalize the amounts of the elements
        amounts_sum = 0
        for element in element_amounts:
            amounts_sum += element_amounts[element]
        for element in element_amounts:
            element_amounts[element] = element_amounts[element]/amounts_sum

        # approximate the decimal amount of each element as a fraction
        # (rational number)
        rational_amounts = {}
        for element in element_amounts:
            rational_amounts[element] = Fraction(
                element_amounts[element]).limit_denominator(
                    self.max_num_atoms)

        # multiply the denominators together, then multiply each fraction
        # by this result to get the number of atoms of each element
        denom_product = 1.0
        for element in rational_amounts:
            denom_product *= rational_amounts[element].denominator
        for element in rational_amounts:
            element_amounts[element] = round(float(
                denom_product)*rational_amounts[element])

        # make a Composition object from the amounts of each element and reduce
        # it if possible
        random_composition = Composition(element_amounts)
        reduced_composition = random_composition.reduced_composition
        num_atoms = int(reduced_composition.num_atoms)

        # check the min and max number of atoms constraints
        if num_atoms > self.max_num_atoms or num_atoms < \
                constraints.min_num_atoms:
            return None

        # check the composition - don't want endpoints
        for endpoint in composition_space.endpoints:
            if endpoint.almost_equals(reduced_composition):
                return None

        # save the element objects
        species = []
        for specie in reduced_composition:
            for _ in range(int(reduced_composition[specie])):
                species.append(specie)
        return species

    def scale_volume(self, random_cell):
        """
        Scales the volume of the random cell according the values in
        self.vpas.

        Returns a boolean indicating whether volume scaling was completed
        without errors.

        Args:
            random_cell: the random Cell whose volume to possibly scale
        """

        # compute the volume to scale to
        composition = random_cell.composition
        total_volume = 0
        for specie in composition:
            total_volume += composition[specie]*self.vpas[specie.symbol]

        # scale the volume
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            random_cell.scale_lattice(total_volume)
            if str(random_cell.lattice.a) == 'nan' or \
                    random_cell.lattice.a > 100:
                return False
            else:
                return True

    def update_status(self):
        '''
        Increments num_made, and if necessary, updates is_finished.
        '''
        self.num_made = self.num_made + 1
        print('Organisms left for {}: {} '.format(
            self.name, self.number - self.num_made))
        if self.num_made == self.number:
            self.is_finished = True


class FileOrganismCreator(object):
    """
    Creates organisms from files (poscar or cif) for the initial population.
    """

    def __init__(self, path_to_folder):
        """
        Makes a FileOrganismCreator.

        Args:
            path_to_folder: the path to the folder containing the files from
                which to make organisms

        Precondition: the folder exists and contains files
        """

        self.name = 'file organism creator'
        self.path_to_folder = path_to_folder
        self.files = [f for f in os.listdir(self.path_to_folder) if
                      os.path.isfile(os.path.join(self.path_to_folder, f))]
        self.number = len(self.files)
        self.num_made = 0  # number of attempts (usually number of files given)
        self.is_successes_based = False  # it's based on number attempted
        self.is_finished = False

    def create_organism(self, id_generator, composition_space, constraints,
                        random):
        """
        Creates an organism for the initial population from a poscar or cif
        file.

        Returns an organism, or None if one could not be created.

        Args:
            id_generator: the IDGenerator used to assign id numbers to all
                organisms

            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            random: a copy of Python's built in PRNG

        TODO: the last three arguments are never actually used in this method,
            but I included them so the method has the same arguments as
            RandomOrganismCreator.create_organism() to allow the
            create_organism method to be called on both RandomOrganismCreator
            and FileOrganismCreator without having to know in advance which one
            it is. Maybe there's a better way to deal with this...
        """

        if self.files[self.num_made - 1].endswith('.cif') or self.files[
                self.num_made - 1].startswith('POSCAR'):
            try:
                new_cell = Cell.from_file(
                    str(self.path_to_folder) + "/" + str(
                        self.files[self.num_made - 1]))
                new_org = Organism(new_cell, id_generator, self.name,
                                   composition_space)
                print('Making organism {} from file: {} '.format(
                    new_org.id, self.files[self.num_made - 1]))
                self.update_status()
                return new_org
            except:
                print('Error reading structure from file: {} '.format(
                    self.files[self.num_made - 1]))
                self.update_status()
                return None
        else:
            print('File {} has invalid extension - file must end with .cif or '
                  'begin with POSCAR '.format(self.files[self.num_made - 1]))
            self.update_status()
            return None

    def get_cells(self):
        """
        Creates cells from the files and puts them in a list.

        Returns the list of Cell objects.

        Used for checking if all the composition space endpoint are included
        for phase diagram searches.
        """

        file_cells = []
        for cell_file in self.files:
            if cell_file.endswith('.cif') or cell_file.startswith(
                    'POSCAR'):
                try:
                    new_cell = Cell.from_file(
                        str(self.path_to_folder) + "/" + str(cell_file))
                    file_cells.append(new_cell)
                except:
                    pass
        return file_cells

    def update_status(self):
        """
        Increments num_made, and if necessary, updates is_finished.
        """

        self.num_made = self.num_made + 1
        print('Organisms left for {}: {} '.format(
            self.name, self.number - self.num_made))
        if self.num_made == len(self.files):
            self.is_finished = True
