# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Variations module:

This module contains the classes used to create organisms for the initial
population.

1. RandomOrganismCreator: creates random organisms

2. FileOrganismCreator: creates organisms by reading their structures from
        files

"""
from general import Organism

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition

from fractions import Fraction
import warnings
import numpy as np
import os


class RandomOrganismCreator(object):
    """
    Creates random organisms for the initial population
    """

    def __init__(self, random_org_parameters, composition_space):
        """
        Creates a RandomOrganismCreator.

        Args:
            random_org_parameters: the parameters for generating random
                organisms

            composition_space: a CompositionSpace object

            random: Python's PRNG
        """

        # the name of this creator
        self.name = 'random organism creator'
        # the default number of random organisms to make (only used for epa
        # searches)
        self.default_number = 28
        # the default volume scaling behavior
        self.default_volume = 'from_atomic_radii'

        # if entire random_org_parameters is None or 'default', then set to
        # defaults
        if random_org_parameters in (None, 'default'):
            self.number = self.default_number
            self.volume = self.default_volume

        # otherwise, parse the parameters and set to defaults if necessary
        else:
            # the number to make
            if 'number' not in random_org_parameters:
                # use the default value if the flag hasn't been used
                self.number = self.default_number
            elif random_org_parameters['number'] in (None, 'default'):
                # use the default value if the flag was left blank or set to
                # 'default'
                self.number = self.default_number
            else:
                # otherwise, parse the value from the parameters
                self.number = random_org_parameters['number']

            # the volume to scale them to
            if 'volume' not in random_org_parameters:
                # use the default value if the flag hasn't been used
                self.volume = self.default_volume
            elif random_org_parameters['volume'] in (None, 'default'):
                # use the default value if the flag was left blank or set to
                # 'default'
                self.volume = self.default_volume
            else:
                # otherwise, parse the value from the parameters
                self.volume = random_org_parameters['volume']

        # to keep track of how many have been made
        # for random organism creator, it's the number of organisms that have
        # been added to the initial population
        self.num_made = 0
        # whether it's based on number created or number added
        self.is_successes_based = True
        # whether it's finished
        self.is_finished = False

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
            id_generator: an IDGenerator object

            composition_space: a CompositionSpace object

            constraints: a Constraints object

            random: copy of Python's built in PRNG
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
        random_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # get a list of elements for the random organism
        if composition_space.objective_function == 'epa':
            reduced_formula = composition_space.endpoints[
                0].reduced_composition
            num_atoms_in_formula = reduced_formula.num_atoms
            max_num_formulas = int(
                constraints.max_num_atoms/num_atoms_in_formula)
            min_num_formulas = int(
                constraints.min_num_atoms/num_atoms_in_formula)
            if min_num_formulas == 0:
                min_num_formulas = 1
            # get a random number of formula units, and the resulting random
            # number of atoms
            random_num_formulas = random.randint(min_num_formulas,
                                                 max_num_formulas)
            num_atoms = int(random_num_formulas*num_atoms_in_formula)
            # add the right number of each element
            elements = []
            for element in reduced_formula:
                # for some reason, reduced_formula[element] is a float for
                # elemental structures, so have to cast it to an int below
                for _ in range(random_num_formulas*int(
                        reduced_formula[element])):
                    elements.append(element)

        elif composition_space.objective_function == 'pd':
            # randomize the order of the endpoints in the list so we don't bias
            # toward the first one
            random.shuffle(composition_space.endpoints)
            # get random fractions for each endpoint (i.e., a random location
            # in the composition space) that sum to 1
            frac_sum = 0
            endpoint_fracs = {}
            for i in range(len(composition_space.endpoints) - 1):
                next_frac = random.uniform(0, 1.0 - frac_sum)
                endpoint_fracs[composition_space.endpoints[i]] = next_frac
                frac_sum += next_frac
            # assign the last one whatever is left
            endpoint_fracs[composition_space.endpoints[-1]] = 1.0 - frac_sum
            # make a dictionary with all the elements as keys and initialize
            # the values to zero
            all_elements = composition_space.get_all_elements()
            element_amounts = {}
            for element in all_elements:
                element_amounts[element] = 0
            # compute the amounts of each element from the amounts of each
            # endpoint
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
            # set the maximum denominator to the max number of atoms  divided
            # by the number of endpoints
            max_denom = int(constraints.max_num_atoms/len(
                composition_space.endpoints))
            # approximate the decimal amount of each element as a fraction
            # (rational number)
            rational_amounts = {}
            for element in element_amounts:
                rational_amounts[element] = Fraction(
                    element_amounts[element]).limit_denominator(
                        random.randint(int(0.9*max_denom), max_denom))
            # multiply all the denominators together
            denom_product = 1.0
            for element in rational_amounts:
                denom_product *= rational_amounts[element].denominator
            # multiply each fraction by this big denominator to get the amounts
            for element in rational_amounts:
                element_amounts[element] = float(
                    denom_product)*rational_amounts[element]
            # now round all the element amounts to the nearest integer
            # (although they should already be integers...)
            for element in element_amounts:
                element_amounts[element] = round(element_amounts[element])
            # make a Composition object from the amounts of each element, and
            # get the reduced composition
            random_composition = Composition(element_amounts)
            reduced_composition = random_composition.reduced_composition
            # set the number of atoms
            num_atoms = int(reduced_composition.num_atoms)
            # check the min and max number of atoms constraints
            if num_atoms > constraints.max_num_atoms or num_atoms < \
                    constraints.min_num_atoms:
                return None
            # check that the composition isn't at one of the composition space
            # endpoints
            for endpoint in composition_space.endpoints:
                if endpoint.almost_equals(reduced_composition):
                    return None
            # put the element objects in a list, with one entry for each atom
            elements = []
            for element in reduced_composition:
                for _ in range(int(reduced_composition[element])):
                    elements.append(element)

        # for each element, generate a set of random fractional coordinates
        # TODO: this doesn't ensure the structure will satisfy the per-species
        #     mids, and in fact most won't. It's ok because they'll just fail
        #     development, but there might be a better way...
        # TODO: also, this doesn't ensure the structure will satisfy the max
        #     size constraint in Geometry. There could be a way to do this by
        #     applying more stringent constraints on how the random lattice
        #     vector lengths and angles are chosen when we have a non-bulk
        #     geometry...
        random_coordinates = []
        for _ in range(num_atoms):
            random_coordinates.append([random.random(), random.random(),
                                       random.random()])

        # make a random structure from the random lattice, random species, and
        # random coordinates
        random_structure = Structure(random_lattice, elements,
                                     random_coordinates)

        # optionally scale the volume of the random structure
        # The actual volume per atom of most elements is about two times the
        # volume per atom computed from the atomic radii. We use that
        # relationship to compute to what volume the the cell should be scaled.
        # Since the per-species mids are just fractions of the atomic radii,
        # this helps remove bias introduced by the per-species mids.
        if self.volume == 'from_atomic_radii':
            # compute volumes per atom (in Angstrom^3) of each element in the
            # random organism

            # the sum of the volumes of all the atoms in the structure, where
            # the volume of each atom is computed from its atomic radius
            total_atomic_volume = 0
            for element in random_structure.composition:
                # compute the volume of each type of atom from it's atomic
                # radius
                volume_per_atom = (4.0/3.0)*np.pi*np.power(
                    element.atomic_radius, 3)
                # get how many atoms of this type (element) there are in the
                # structure
                num_atoms = random_structure.composition[element]
                # increment the total atomic volume of the structure
                total_atomic_volume += volume_per_atom*num_atoms

            # the empirical ratio of the actual volume per atom and that
            # obtained by summing atomic volumes (computed from atomic radii)
            scale_factor = 2

            # Scale the volume of the random organism to satisfy the computed
            # mean volume per atom. The structure.scale_lattice method
            # sometimes throws divide-by-zero runtime warnings. To keep these
            # from getting printed to the output, we're temporarily suppressing
            # warnings here. When the divide-by-zero error happens, the lattice
            # either gets scaled to a huge volume, or else scaling fails and
            # the lattice vectors are set to NaN. After doing the scaling, the
            # first lattice vector is checked to see if either of these error
            # occured.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                # check if the volume scaling worked
                random_structure.scale_lattice(
                    scale_factor*total_atomic_volume)
                if str(random_structure.lattice.a) == 'nan' or \
                        random_structure.lattice.a > 100:
                    return None

        elif self.volume == 'random':
            # no volume scaling
            pass

        else:
            # scale to the given volume per atom
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                random_structure.scale_lattice(self.volume*len(
                    random_structure.sites))
                # check if the volume scaling worked
                if str(random_structure.lattice.a) == 'nan' or \
                        random_structure.lattice.a > 100:
                    return None

        # return a random organism with the scaled random structure
        random_org = Organism(random_structure, id_generator, self.name)
        print('Random organism creator making organism {} '.format(
            random_org.id))

        # for testing
        # print('Random organism composition: {}'.format(
        #    random_org.composition.formula.replace(' ', '')))

        return random_org

    def update_status(self):
        '''
        Increments num_made, and if necessary, updates is_finished
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
        Creates a FileOrganismCreator.

        Args:
            path_to_folder: the path to the folder containing the files from
                which to make organisms

        Precondition: the folder exists and contains files
        """

        # the name of this creator
        self.name = 'file organism creator'
        # all the files in the given folder
        self.path_to_folder = path_to_folder
        self.files = [f for f in os.listdir(self.path_to_folder) if
                      os.path.isfile(os.path.join(self.path_to_folder, f))]
        self.number = len(self.files)

        # to keep track of how many have been made
        # for file organism creator, it's the number of attempts to make
        # organisms from files (usually the number of files provided)
        self.num_made = 0
        # whether it's based on number created or number added
        self.is_successes_based = False
        # whether it's finished
        self.is_finished = False

    def create_organism(self, id_generator, composition_space, constraints,
                        random):
        """
        Creates an organism for the initial population from a poscar or cif
        file.

        Returns an organism, or None if one could not be created.

        Args:
            id_generator: an IDGenerator object

            composition_space: a CompositionSpace object

            constraints: a Constraints object

            random: Python's built in PRNG

        TODO: the last three arguments are never actually used in this method,
            but I included them so the method has the same arguments as
            RandomOrganismCreator.creatorOrganism() to allow the
            create_organism method to be called on both RandomOrganismCreator
            and FileOrganismCreator without having to know in advance which one
            it is. Maybe there's a better way to deal with this...
        """

        if self.files[self.num_made - 1].endswith('.cif') or self.files[
                self.num_made - 1].startswith('POSCAR'):
            try:
                new_struct = Structure.from_file(
                    str(self.path_to_folder) + "/" + str(
                        self.files[self.num_made - 1]))
                new_org = Organism(new_struct, id_generator, self.name)
                print('Making organism {} from file: {} '.format(
                    new_org.id, self.files[self.num_made - 1]))
                # update status each time the method is called, since this is
                # an attempts-based creator
                self.update_status()
                return new_org
            # return None if a structure couldn't be read from a file
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

    def get_structures(self):
        """
        Creates a structures from the files and puts them in a list.

        Returns the list of Structure objects.

        Used for checking if all the composition space endpoint are included
        for phase diagram searches.
        """

        file_structures = []
        for structure_file in self.files:
            if structure_file.endswith('.cif') or structure_file.startswith(
                    'POSCAR'):
                try:
                    new_struct = Structure.from_file(
                        str(self.path_to_folder) + "/" + str(structure_file))
                    file_structures.append(new_struct)
                except:
                    pass
        return file_structures

    def update_status(self):
        """
        Increments num_made, and if necessary, updates is_finished
        """

        self.num_made = self.num_made + 1
        print('Organisms left for {}: {} '.format(
            self.name, self.number - self.num_made))
        if self.num_made == len(self.files):
            self.is_finished = True
