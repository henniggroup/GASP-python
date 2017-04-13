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
import numpy as np


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
        # allow structure with compositions at the endpoints (for pd searches)
        self.default_allow_endpoints = True
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
            self.allow_endpoints = self.default_allow_endpoints
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

            # allowing composition space endpoints (only used for pd searches)
            if 'allow_endpoints' not in random_org_parameters:
                self.allow_endpoints = self.default_allow_endpoints
            elif random_org_parameters['allow_endpoints'] in (None, 'default'):
                self.allow_endpoints = self.default_allow_endpoints
            else:
                self.allow_endpoints = random_org_parameters['allow_endpoints']

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

            5. Checks if the total number of atoms exceeds self.max_num_atoms.
                If so, reduce the amount of each atom with a multiplicative
                factor.

            6. Reduces the resulting composition (i.e., find the smallest
                number of atoms needed to describe the composition).

            7. Optionally increases the number of atoms (w/o changing the
                composition) such that the min num atoms constraint is
                satisfied if possible.

            8. Checks that the resulting number of atoms satisfies the maximum
                (self.max_num_atoms) number of atoms constraint, and optionally
                checks that the resulting composition is not equivalent to one
                of the endpoint compositions.
        """

        # get random fractions for each endpoint that sum to one 1 (i.e., a
        # random location in the composition space
        fracs = self.get_random_endpoint_fractions(composition_space, random)
        composition_space.endpoints.sort()
        endpoint_fracs = {}
        for i in range(len(fracs)):
            endpoint_fracs[composition_space.endpoints[i]] = fracs[i]

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

        # see how many total atoms we have
        num_atoms = 0
        for element in element_amounts:
            num_atoms += element_amounts[element]

        # reduce the number of atoms of each element if needed
        if num_atoms > self.max_num_atoms:
            numerator = random.randint(
                int(round(0.5*(constraints.min_num_atoms +
                               self.max_num_atoms))), self.max_num_atoms)
            factor = numerator/num_atoms
            for element in element_amounts:
                element_amounts[element] = round(
                    factor*element_amounts[element])

        # make a Composition object from the amounts of each element
        random_composition = Composition(element_amounts)
        random_composition = random_composition.reduced_composition

        # possibly increase the number of atoms by a random (allowed) amount
        min_multiple = int(
            math.ceil(constraints.min_num_atoms/random_composition.num_atoms))
        max_multiple = int(
            math.floor(self.max_num_atoms/random_composition.num_atoms))
        if max_multiple > min_multiple:
            random_multiple = random.randint(min_multiple, max_multiple)
            bigger_composition = {}
            for element in random_composition:
                bigger_composition[element] = \
                    random_multiple*random_composition[element]
            random_composition = Composition(bigger_composition)

        # check the max number of atoms constraints (should be ok)
        if int(random_composition.num_atoms) > self.max_num_atoms:
            return None

        # check the composition - only allow endpoints if specified
        if not self.allow_endpoints:
            for endpoint in composition_space.endpoints:
                if endpoint.almost_equals(
                        random_composition.reduced_composition):
                    return None

        # save the element objects
        species = []
        for specie in random_composition:
            for _ in range(int(random_composition[specie])):
                species.append(specie)
        return species

    def get_random_endpoint_fractions(self, composition_space, random):
        """
        Uniformly samples the composition space. Returns a list containing the
        fractions of each endpoint composition (that sum to 1).

        Args:
            composition_space: the CompositionSpace of the search

            random: a copy of Python's built-in PRNG

        Description:

            1. Computes vectors that span the normalized composition space
                (e.g., the triangular facet for a ternary system) by
                subtracting the first composition fraction unit vector from the
                others.

            2. Takes a random linear combination of these vectors by
                multiplying each one by a uniform random number and then taking
                their sum.

            3. Adds the first composition unit vector to the result from step 2
                to obtain a vector with random fractions of each endpoint
                composition.

            4. Checks that the vector from step 3 lies in the portion of the
                plane that corresponds to normalized amounts. This is done be
                checking that amount of the first endpoint composition is
                non-negative. If it's negative, calls itself recursively until
                a valid solution is found.
        """

        # compute the vectors corresponding to the needed binary edges of the
        # phase diagram (w.r.t. to the first endpoint of the composition space)
        num_endpoints = len(composition_space.endpoints)
        bindary_edges = []
        for i in range(1, num_endpoints):
            edge = [-1]
            for j in range(1, num_endpoints):
                if j == i:
                    edge.append(1)
                else:
                    edge.append(0)
            bindary_edges.append(np.array(edge))

        # take a linear combination of the edge vectors, where the weight of
        # each vector is drawn from a uniform distribution
        weighted_average = random.random()*bindary_edges[0]
        for i in range(1, len(bindary_edges)):
            weighted_average = np.add(weighted_average,
                                      random.random()*bindary_edges[i])

        # add the first unit vector to the weighted average of the edge
        # vectors to obtain the fractions of each endpoint
        endpoint_fracs = weighted_average.tolist()
        endpoint_fracs[0] = endpoint_fracs[0] + 1

        # check that the computed fraction of the first endpoint is not less
        # than zero. If it is, try again.
        if endpoint_fracs[0] < 0:
            return self.get_random_endpoint_fractions(composition_space,
                                                      random)
        else:
            return endpoint_fracs

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
