# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Variations module:

This module contains the classes used to create offspring organisms from parent
organisms.

1. Mating: creates an offspring organism by combining two parent organisms

2. StructureMut: creates an offspring organism by mutating the structure of a
        parent organism

3. NumStoichsMut: creates an offspring organism by adding or removing atoms
        from a parent organism

4. Permutation: creates an offspring organism by swapping atoms in a parent
        organism

"""

from general import Organism

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, Specie

import copy
import numpy as np


class Mating(object):
    """
    An operator that creates an offspring organism by combining the structures
    of two parent organisms.
    """

    def __init__(self, mating_params):
        """
        Creates a Mating operator

        Args:
            mating_params: The parameters for doing the mating operation, as a
                dictionary.

        Precondition: the 'fraction' parameter in mating_params is not
            optional, and it is assumed that mating_params contains this
            parameter.
        """

        # the name of this variation
        self.name = 'mating'

        # parse the fraction from the parameters. This argument is not
        # optional, so it doesn't have a default value here
        self.fraction = mating_params['fraction']

        # default values
        # the average (fractional) location along the randomly chosen lattice
        # vector to make the cut
        self.default_mu_cut_loc = 0.5
        # the standard deviation of the (fractional) location along the
        # randomly chosen lattice vector to make the cut
        self.default_sigma_cut_loc = 0.5
        # the probability of randomly shifting the atoms along the lattice
        # vector of the cut
        self.default_shift_prob = 1.0
        # the probability that one of the parents will be doubled before doing
        # the variation
        self.default_doubling_prob = 0.1
        # whether or not to grow the smaller parent (by taking a supercell) to
        # the approximate size the larger parent before doing the variation
        self.default_grow_parents = True
        # the cutoff distance (as fraction of atomic radius) below which to
        # merge sites with the same element in the offspring structure
        self.default_merge_cutoff = 1.0

        # the mean of the cut location
        if 'mu_cut_loc' not in mating_params:
            # use the default value if the flag hasn't been used
            self.mu_cut_loc = self.default_mu_cut_loc
        elif mating_params['mu_cut_loc'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.mu_cut_loc = self.default_mu_cut_loc
        else:
            # otherwise, parse the value from the parameters
            self.mu_cut_loc = mating_params['mu_cut_loc']

        # the standard deviation of the cut location
        if 'sigma_cut_loc' not in mating_params:
            # use the default value if the flag hasn't been used
            self.sigma_cut_loc = self.default_sigma_cut_loc
        elif mating_params['sigma_cut_loc'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.sigma_cut_loc = self.default_sigma_cut_loc
        else:
            # otherwise, parse the value from the parameters
            self.sigma_cut_loc = mating_params['sigma_cut_loc']

        # the probability of shifting the atoms
        if 'shift_prob' not in mating_params:
            # use the default value if the flag hasn't been used
            self.shift_prob = self.default_shift_prob
        elif mating_params['shift_prob'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.shift_prob = self.default_shift_prob
        else:
            # otherwise, parse the value from the parameters
            self.shift_prob = mating_params['shift_prob']

        # the probability of doubling one of the parents before doing the
        # variation
        if 'doubling_prob' not in mating_params:
            # use the default value if the flag hasn't been used
            self.doubling_prob = self.default_doubling_prob
        elif mating_params['doubling_prob'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.doubling_prob = self.default_doubling_prob
        else:
            # otherwise, parse the value from the parameters
            self.doubling_prob = mating_params['doubling_prob']

        # whether to grow the smaller of the parents before doing the variation
        if 'grow_parents' not in mating_params:
            # use the default value if the flag hasn't been used
            self.grow_parents = self.default_grow_parents
        elif mating_params['grow_parents'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.grow_parents = self.default_grow_parents
        else:
            # otherwise, parse the value from the parameters
            self.grow_parents = mating_params['grow_parents']

        # the cutoff distance (as fraction of atomic radius) below which to
        # merge sites in the offspring structure
        if 'merge_cutoff' not in mating_params:
            # use the default value if the flag hasn't been used
            self.merge_cutoff = self.default_merge_cutoff
        elif mating_params['merge_cutoff'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.merge_cutoff = self.default_merge_cutoff
        else:
            # otherwise, parse the value from the parameters
            self.merge_cutoff = mating_params['merge_cutoff']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the mating operation.

        Returns the resulting offspring as an Organism.

        Args:
            pool: the Pool of Organisms

            random: Python's built in PRNG

            geometry: the Geometry object

            id_generator: the IDGenerator

        Description:

            Creates an offspring organism by combining chunks cut from two
            parent organisms.

                1. Selects two organisms from the pool to act as parents, and
                    makes copies of them.

                2. Optionally doubles one of the parents. This occurs with
                    probability self.doubling_prob, and if it happens, the
                    parent with the smallest cell volume is doubled.

                3. Optionally grows one the parents to the approximate size of
                    the other parent, if self.grow_parents is True. The number
                    of times to double the smaller parent is determined by the
                    ratio of the cell volumes of the parents.

                4. Randomly selects one of the three lattice vectors to slice.

                5. Determines the cut location (in fractional coordinates) by
                    drawing from a Gaussian with mean self.mu_cut_loc and
                    standard deviation self.sigma_cut_loc.

                6. In each parent, optionally shift the atoms (in fractional
                    space, with probability self.shift_prob) by an amount drawn
                    from a uniform distribution along the direction of the
                    lattice vector to cut. For non-bulk geometries, shift only
                    occurs if the lattice vector to cut is along a periodic
                    direction.

                7. Copy the sites from the first parent organism with
                    fractional coordinate less than the randomly chosen cut
                    location along the randomly chosen lattice vector to the
                    offspring organism, and do the same for the second parent
                    organism, except copy the sites with fractional coordinate
                    greater than the cut location.

                8. Merge sites in offspring structure that have the same
                    element and are closer than self.merge_sites times the
                    atomic radius of the element.
        """

        # select two parent organisms from the pool
        parent_orgs = pool.select_organisms(2, random)

        # make deep copies of the parent organisms
        parent_1 = copy.deepcopy(parent_orgs[0])
        parent_2 = copy.deepcopy(parent_orgs[1])

        # optionally double one of the parents
        if random.random() < self.doubling_prob:
            # pick the smallest parent (based on cell volume)
            vol_1 = parent_1.structure.lattice.volume
            vol_2 = parent_2.structure.lattice.volume
            if vol_1 < vol_2:
                self.double_parent(parent_1, geometry, random)
            else:
                self.double_parent(parent_2, geometry, random)

        # grow the smaller parent if specified
        if self.grow_parents:
            # pick the smallest parent (based on cell volume)
            vol_1 = parent_1.structure.lattice.volume
            vol_2 = parent_2.structure.lattice.volume
            if vol_1 < vol_2:
                volume_ratio = vol_2/vol_1
                parent_to_grow = parent_1
            else:
                volume_ratio = vol_1/vol_2
                parent_to_grow = parent_2
            # compute how many times to double the smaller parent based on
            # ratio of cell volumes
            num_doubles = self.get_num_doubles(volume_ratio)
            # double the smaller parent the computed number of times
            for _ in range(num_doubles):
                self.double_parent(parent_to_grow, geometry, random)

        # lists to hold the species of the sites contributed from each parent
        species_from_parent_1 = []
        species_from_parent_2 = []

        # loop needed here because sometimes no atoms get contributed from one
        # of the parents, so have to try again
        while len(species_from_parent_1) == 0 or len(
                species_from_parent_2) == 0:
            # lists to hold the species and fractional coordinates of the sites
            # contributed from each parent.
            species_from_parent_1 = []
            frac_coords_from_parent_1 = []
            species_from_parent_2 = []
            frac_coords_from_parent_2 = []

            # randomly select the lattice vector to cut
            cut_vector_index = random.randint(0, 2)

            # draw the random cut location from a Gaussian, and make sure it's
            # between 0 and 1
            cut_location = random.gauss(self.mu_cut_loc, self.sigma_cut_loc)
            while cut_location > 1 or cut_location < 0:
                cut_location = random.gauss(self.mu_cut_loc,
                                            self.sigma_cut_loc)

            # possibly shift the atoms in the first parent along the cut vector
            if random.random() < self.shift_prob:
                shifted_parent_1 = self.do_random_shift(
                    parent_1, cut_vector_index, geometry, random)

            # possibly shift the atoms in the second parent along the cut
            # vector
            if random.random() < self.shift_prob:
                shifted_parent_2 = self.do_random_shift(
                    parent_2, cut_vector_index, geometry, random)

            # get the species and fractional coordinates of each site in parent
            # 1 with fractional coordinate along the cut vector less than the
            # cut location
            for site in shifted_parent_1.structure.sites:
                if site.frac_coords[cut_vector_index] < cut_location:
                    species_from_parent_1.append(site.species_and_occu)
                    frac_coords_from_parent_1.append(site.frac_coords)

            # get the species and fractional coordinates of each site in parent
            # 2 with fractional coordinate along the cut vector greater than
            # the cut location
            for site in shifted_parent_2.structure.sites:
                if site.frac_coords[cut_vector_index] > cut_location:
                    species_from_parent_2.append(site.species_and_occu)
                    frac_coords_from_parent_2.append(site.frac_coords)

        # combine the information for the sites contributed by each parent
        offspring_species = species_from_parent_1 + species_from_parent_2
        offspring_frac_coords = frac_coords_from_parent_1 + \
            frac_coords_from_parent_2

        # compute the lattice vectors of the offspring by taking the average
        # of the lattice vectors of the parents
        offspring_lengths = 0.5*(np.array(
            shifted_parent_1.structure.lattice.abc) + np.array(
                shifted_parent_2.structure.lattice.abc))
        offspring_angles = 0.5*(np.array(
            shifted_parent_1.structure.lattice.angles) + np.array(
                shifted_parent_2.structure.lattice.angles))
        offspring_lattice = Lattice.from_lengths_and_angles(offspring_lengths,
                                                            offspring_angles)

        # make the offspring structure from the offspring lattice, species and
        # fractional coordinates
        offspring_structure = Structure(offspring_lattice, offspring_species,
                                        offspring_frac_coords)

        # merge sites in the offspring structure
        offspring_structure = self.merge_sites(offspring_structure)

        # make the offspring organism from the offspring structure
        offspring = Organism(offspring_structure, id_generator, self.name)

        # print out a message
        print('Creating offspring organism {} from parent organisms {} and {} '
              'with the mating variation '.format(offspring.id,
                                                  parent_orgs[0].id,
                                                  parent_orgs[1].id))
        return offspring

    def get_num_doubles(self, volume_ratio):
        """
        Returns the number of times to double a cell based on the given volume
        ratio. Essentially maps the volume ratio to a step function that
        approximates the base 2 logarithm.

        Args:
            volume_ratio: the ratio of the volumes of the two parent organisms.
        """

        if volume_ratio < 1.5:
            return 0
        elif volume_ratio >= 1.5 and volume_ratio < 3:
            return 1
        elif volume_ratio >= 3 and volume_ratio < 6:
            return 2
        elif volume_ratio >= 6 and volume_ratio < 12:
            return 3
        elif volume_ratio >= 12 and volume_ratio < 24:
            return 4
        elif volume_ratio >= 24 and volume_ratio < 48:
            return 5
        elif volume_ratio >= 48 and volume_ratio < 96:
            return 6
        else:
            return 7

    def double_parent(self, organism, geometry, random):
        """
        Takes a supercell of the organism. For bulk geometries, the supercell
        is taken in the direction of the organism's shortest lattice vector.
        For non-bulk geometries, the supercell is taken in the direction of a
        randomly chosen lattice vector.

        Modifies the structure of the organism.

        Args:
            organism: the Organism to take the supercell of

            geometry: a Geometry object
        """

        if geometry.shape == 'bulk':
            # get the index of the smallest lattice vector of the organism
            lattice_lengths = organism.structure.lattice.abc
            smallest_vector = min(lattice_lengths)
            doubling_index = lattice_lengths.index(smallest_vector)
        else:
            # if geometry is not bulk, then  pick a random lattice vector to
            # double
            doubling_index = random.choice([0, 1, 2])
        # take a supercell of the organism along the smallest lattice vector
        scaling_factors = [1, 1, 1]
        scaling_factors[doubling_index] = 2
        organism.structure.make_supercell(scaling_factors)

    def do_random_shift(self, organism, lattice_vector_index, geometry,
                        random):
        """
        Makes a copy of the organism, and shifts all the atoms in the copy by a
        random amount (drawn from uniform distribution) along the lattice
        vector specified by the given index. After shifting the atoms by the
        random amount along the specified lattice vector, checks if any atoms
        lie outside the cell. If so, they are replaced with their periodic
        images inside the cell.

        Note: this method checks the geometry, and will not do the shift if the
        specified lattice vector is not in a periodic direction because that
        could destroy some of the local structure of a non-bulk structure.

        Args:
            organism: the Organism whose Structure to change by shifting the
                atoms

            lattice_vector_index: the index (0, 1 or 2) of the lattice vector
                along which to shift the atoms

            geometry: the Geometry object

            random: Python's built in PRNG
        """

        # if shape is cluster, then no shift
        if geometry.shape == 'cluster':
            return organism
        # if shape is wire, then don't shift if not c lattice vector
        elif geometry.shape == 'wire' and lattice_vector_index != 2:
            return organism
        # if shape is sheet, then don't shift if c lattice vector
        elif geometry.shape == 'sheet' and lattice_vector_index == 2:
            return organism
        # otherwise, shift all the atoms along the specified lattice vector by
        # a random amount
        else:
            # make a copy of the organism
            shifted_org = copy.deepcopy(organism)

            # randomly select the fractional amount to shift the atoms along
            # the lattice vector
            shift_vector = random.random(
                )*shifted_org.structure.lattice.matrix[lattice_vector_index]
            site_indices = [i for i in range(len(shifted_org.structure.sites))]
            shifted_org.structure.translate_sites(site_indices, shift_vector,
                                                  frac_coords=False,
                                                  to_unit_cell=False)

            # translate the atoms back into the cell in case the shifts moved
            # some of them outside of it
            # dictionary to map site indices to translation vectors
            translation_vectors = {}

            # compute the translation vector for each site (combination of the
            # lattice vectors)
            for site in shifted_org.structure.sites:
                translation_vector = []
                for coord in site.frac_coords:
                    if coord < 0.0:
                        translation_vector.append(1.0)
                    elif coord > 1.0:
                        translation_vector.append(-1.0)
                    else:
                        translation_vector.append(0.0)
                # save the computed translation vector in the dictionary
                translation_vectors[
                    shifted_org.structure.sites.index(
                        site)] = translation_vector

            # move the sites back into the cell
            for key in translation_vectors:
                shifted_org.structure.translate_sites(key,
                                                      translation_vectors[key],
                                                      frac_coords=True,
                                                      to_unit_cell=False)
            return shifted_org

    def merge_sites(self, structure):
        """
        Merges sites in the structure that have the same element and are closer
        that self.merge_sites times the atomic radius. Merging means replacing
        two sites in the structure with one site at the mean of the two sites
        positions.

        Returns a new structure with the sites merged

        Args:
            structure: the Structure whose sites to merge
        """

        # list to hold the data for the new sites (produced from merging)
        species = []
        frac_coords = []
        # list to hold the indices of the sites that have been merged (so they
        # won't be in the new structure)
        merged_indices = []

        # go through the structure site by site and merge pairs of sites if
        # needed
        for site in structure.sites:
            # check that the site hasn't already been merged
            if structure.sites.index(site) not in merged_indices:
                symbol = site.specie.symbol
                element = Element(site.specie.symbol)
                a_radius = element.atomic_radius
                for other_site in structure.sites:
                    # check that the other site hasn't already been merged
                    if structure.sites.index(other_site) not in merged_indices:
                        # check that the other site is not the site, and that
                        # it has the same symbol as the site
                        if other_site != site and other_site.specie.symbol == \
                                    symbol:
                            # check the distance between the sites
                            if site.distance(
                                    other_site) < a_radius*self.merge_cutoff:
                                # record that both of these sites have been
                                # merged
                                merged_indices.append(
                                    structure.sites.index(other_site))
                                merged_indices.append(
                                    structure.sites.index(site))
                                # compute the fractional coordinates of the new
                                # site
                                new_frac_coords = np.add(
                                    site.frac_coords, other_site.frac_coords)/2
                                species.append(site.specie)
                                frac_coords.append(new_frac_coords)

        # get the data for the sites that were NOT merged
        for site in structure.sites:
            if structure.sites.index(site) not in merged_indices:
                species.append(site.specie)
                frac_coords.append(site.frac_coords)

        # make a new structure
        new_structure = Structure(structure.lattice, species, frac_coords)

        # if any merges were done, call merge_sites recursively on the new
        # structure
        if len(new_structure.sites) < len(structure.sites):
            return self.merge_sites(new_structure)
        else:
            return new_structure


class StructureMut(object):
    """
    An operator that creates an offspring organism by mutating the structure
    of a parent organism.
    """

    def __init__(self, structure_mut_params):
        """
        Creates a Mutation operator.

        Args:
            structure_mut_params: The parameters for doing the mutation
                operation, as a dictionary

        Precondition: the 'fraction' parameter in structure_mut_params is not
            optional, and it is assumed that structure_mut_params contains this
            parameter
        """

        # the name of this variation
        self.name = 'structure mutation'

        # parse the fraction from the parameters. This argument is not
        # optional, so it doesn't have a default value here
        self.fraction = structure_mut_params['fraction']

        # the default values

        # the fraction of atoms to perturb (on average)
        self.default_frac_atoms_perturbed = 1.0
        # the standard deviation of the perturbation of an atomic coordinate
        # (in Angstroms)
        self.default_sigma_atomic_coord_perturbation = 0.5
        # the maximum allowed perturbation of an atomic coordinate
        # (in Angstroms)
        self.default_max_atomic_coord_perturbation = 5.0
        # the standard deviation of the non-identity components of the elements
        # of the strain matrix
        # note: the non-identity components of the strain matrix elements are
        # constrained to be between -1 and 1
        self.default_sigma_strain_matrix_element = 0.2

        # the fraction of atoms to perturb, on average
        if 'frac_atoms_perturbed' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        elif structure_mut_params['frac_atoms_perturbed'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        else:
            # otherwise, parse the value from the parameters
            self.frac_atoms_perturbed = structure_mut_params[
                'frac_atoms_perturbed']

        # the standard deviation of the perturbation of each atomic coordinate
        if 'sigma_atomic_coord_perturbation' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_atomic_coord_perturbation = \
                self.default_sigma_atomic_coord_perturbation
        elif structure_mut_params['sigma_atomic_coord_perturbation'] in (
                None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.sigma_atomic_coord_perturbation = \
                self.default_sigma_atomic_coord_perturbation
        else:
            # otherwise, parse the value from the parameters
            self.sigma_atomic_coord_perturbation = structure_mut_params[
                'sigma_atomic_coord_perturbation']

        # the maximum allowed magnitude of perturbation of each atomic coord
        if 'max_atomic_coord_perturbation' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.max_atomic_coord_perturbation = \
                self.default_max_atomic_coord_perturbation
        elif structure_mut_params['max_atomic_coord_perturbation'] in (
                None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.max_atomic_coord_perturbation = \
                self.default_max_atomic_coord_perturbation
        else:
            # otherwise, parse the value from the parameters
            self.max_atomic_coord_perturbation = structure_mut_params[
                'max_atomic_coord_perturbation']

        # the standard deviation of the magnitude of the non-identity
        # components of the elements in the strain matrix
        if 'sigma_strain_matrix_element' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_strain_matrix_element = \
                self.default_sigma_strain_matrix_element
        elif structure_mut_params['sigma_strain_matrix_element'] in (
                None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.sigma_strain_matrix_element = \
                self.default_sigma_strain_matrix_element
        else:
            # otherwise, parse the value from the parameters
            self.sigma_strain_matrix_element = structure_mut_params[
                'sigma_strain_matrix_element']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the structural mutation operation.

        Returns the resulting offspring as an Organism.

         Args:
            pool: the Pool of Organisms

            random: Python's built in PRNG

            geometry: the Geometry object

            id_generator: the IDGenerator

        Description:

            Creates an offspring organism by perturbing the atomic positions
            and lattice vectors of the parent structure.

                1. Selects a parent organism from the pool and makes a copy of
                    it

                2. Perturbs the atomic coordinates of each site with
                    probability self.frac_atoms_perturbed. The perturbation of
                    each atomic coordinate is drawn from a Gaussian with mean
                    zero and standard deviation
                    self.sigma_atomic_coord_perturbation. The magnitude of each
                    atomic coordinate perturbation is constrained to not exceed
                    self.max_atomic_coord_perturbation.

                3. The lattice vectors are perturbed by taking the product of
                    each lattice vector with a strain matrix. The strain matrix
                    is defined as

                        I + E

                    where I is the 3x3 identity matrix and E is a 3x3
                    perturbation matrix whose elements are distinct and drawn
                    from a Gaussian with mean zero and standard deviation
                    self.sigma_strain_matrix_element and are constrained to lie
                    between -1 and 1
        """

        # select a parent organism from the pool
        parent_org = pool.select_organisms(1, random)

        # make a deep copy of the structure of the parent organism, so that the
        # subsequent mutation doesn't affect the structure of the parent
        structure = copy.deepcopy(parent_org[0].structure)

        # for each site in the structure, determine whether to perturb it
        for site in structure.sites:
            if random.random() < self.frac_atoms_perturbed:
                # generate three random perturbation magnitudes (in Cartesian
                # coordinates), one for each atomic coordinate, and make sure
                # they aren't larger than the max allowed perturbation

                # the first one
                nudge_x = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_x) > self.max_atomic_coord_perturbation:
                    nudge_x = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # the second one
                nudge_y = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_y) > self.max_atomic_coord_perturbation:
                    nudge_y = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # the third one
                nudge_z = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_z) > self.max_atomic_coord_perturbation:
                    nudge_z = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # make the perturbation vector
                perturbation_vector = [nudge_x, nudge_y, nudge_z]
                # translate the site be the computed perturbation vector
                structure.translate_sites(
                    structure.sites.index(site), perturbation_vector,
                    frac_coords=False, to_unit_cell=False)

        # compute the random non-identity components of the nine elements of
        # the strain matrix
        epsilons = []
        for _ in range(9):
            epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            # make sure it's in [-1, 1]
            while np.abs(epsilon) > 1:
                epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            epsilons.append(epsilon)

        # construct the strain matrix (I + epsilon_ij), and randomly assign
        # positive or negative directions to the non-identity components
        row_1 = [1 + epsilons[0], epsilons[1], epsilons[2]]
        row_2 = [epsilons[3], 1 + epsilons[4], epsilons[5]]
        row_3 = [epsilons[6], epsilons[7], 1 + epsilons[8]]
        strain_matrix = np.array([row_1, row_2, row_3])

        # compute new lattice vectors by applying the strain matrix to the
        # lattice vectors
        new_a = strain_matrix.dot(structure.lattice.matrix[0])
        new_b = strain_matrix.dot(structure.lattice.matrix[1])
        new_c = strain_matrix.dot(structure.lattice.matrix[2])
        new_lattice = Lattice([new_a, new_b, new_c])

        # assign the new lattice to the structure (this doesn't change the
        # sites' fractional coordinates, but it does change their Cartesian
        # coordinates)
        structure.modify_lattice(new_lattice)

        # make a new offspring organism
        offspring = Organism(structure, id_generator, self.name)

        # make sure all the site are within the cell (some of them could have
        # been pushed outside by the atomic coordinate perturbations)
        offspring.translate_atoms_into_cell()

        # print out a message
        print('Creating offspring organism {} from parent organism {} with '
              'the structure mutation variation '.format(offspring.id,
                                                         parent_org[0].id))
        return offspring


class NumStoichsMut(object):
    """
    An operator that creates an offspring organism by mutating the number of
    stoichiometries' worth of atoms in the parent organism.
    """

    def __init__(self, num_stoichs_mut_params):
        """
        Creates a NumStoichsMut operator

        Args:
            num_stoichs_mut_params: The parameters for doing the NumStoichsMut
                operation, as a dictionary

        Precondition: the 'fraction' parameter in num_stoichs_mut_params is not
            optional, and it is assumed that num_stoichs_mut_params contains
            this parameter
        """

        # the name of this variation
        self.name = 'number of stoichiometries mutation'

        # parse the fraction from the parameters. This argument is not
        # optional, so it doesn't have a default value here
        self.fraction = num_stoichs_mut_params['fraction']

        # the default values

        # the average number of stoichimetries to add
        self.default_mu_num_adds = 0
        # the standard deviation of the number of stoichiometries to add
        self.default_sigma_num_adds = 1
        # whether to scale the volume of the offspring to equal that of the
        # parent
        self.default_scale_volume = True

        # the average number of stoichiometries to add
        if 'mu_num_adds' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.mu_num_adds = self.default_mu_num_adds
        elif num_stoichs_mut_params['mu_num_adds'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.mu_num_adds = self.default_mu_num_adds
        else:
            # otherwise, parse the value from the parameters
            self.mu_num_adds = num_stoichs_mut_params['mu_num_adds']

        # the standard deviation of the number of stoichiometries
        if 'sigma_num_adds' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_num_adds = self.default_sigma_num_adds
        elif num_stoichs_mut_params['sigma_num_adds'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.sigma_num_adds = self.default_sigma_num_adds
        else:
            # otherwise, parse the value from the parameters
            self.sigma_num_adds = num_stoichs_mut_params['sigma_num_adds']

        # whether to scale the volume of the offspring
        if 'scale_volume' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.scale_volume = self.default_scale_volume
        elif num_stoichs_mut_params['scale_volume'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.scale_volume = self.default_scale_volume
        else:
            # otherwise, parse the value from the parameters
            self.scale_volume = num_stoichs_mut_params['scale_volume']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the number of stoichiometries mutation operation

        Returns the resulting offspring as an Organism.

        Args:
            pool: the Pool of Organisms

            random: Python's built in PRNG

            geometry: the Geometry object

            id_generator: the IDGenerator

        Description:

            Creates an offspring organism by adding or removing a random number
            of stoichiometries' worth of atoms to or from the parent structure.

                1. Selects a parent organism from the pool and makes a copy of
                    it

                2. Computes the number of stoichiometries to add or remove by
                    drawing from a Gaussian with mean self.mu_num_adds and
                    standard deviation self.sigma_num_adds and rounding the
                    result to the nearest integer

                3. Computes the number of atoms of each type to add or remove,
                    and does the additions or removals

                4. If self.scale_volume is True, scales the new structure to
                    have the same volume per atom as the parent
        """

        # select a parent organism from the pool
        parent_org = pool.select_organisms(1, random)

        # make a deep copy of the structure of the parent organism
        structure = copy.deepcopy(parent_org[0].structure)
        # get the total number of atoms in the parent
        parent_num_atoms = len(structure.sites)
        # get the reduced composition of the parent
        reduced_composition = structure.composition.reduced_composition
        # get the volume per atom of the parent
        vol_per_atom = structure.lattice.volume/len(structure.sites)

        # loop needed here in case no atoms got added or removed, in which case
        # need to try again. This happens if the randomly chosen number of
        # atoms to remove exceeds the number of atoms in the cell
        while len(structure.sites) == parent_num_atoms:

            # compute the number of stoichiometries to add or remove, and make
            # sure it's not zero
            num_add = int(round(random.gauss(self.mu_num_adds,
                                             self.sigma_num_adds)))
            while num_add == 0:
                num_add = int(round(random.gauss(self.mu_num_adds,
                                                 self.sigma_num_adds)))

            # compute the number of each type of atom to add (or remove), and
            # also the total number of atoms to add or remove
            amounts_to_add = {}
            total_add = 0
            for key in reduced_composition:
                amounts_to_add[key] = int(num_add*reduced_composition[key])
                total_add = total_add + int(num_add*reduced_composition[key])

            # if num_add is positive, put the new atoms in the cell at random
            # locations
            if num_add > 0:
                for key in amounts_to_add:
                    for _ in range(amounts_to_add[key]):
                        frac_coords = [random.random(), random.random(),
                                       random.random()]
                        structure.append(Specie(key, 0), frac_coords)
                # to remove the oxidation state of 0 we had to specify above
                # in the structure.append() method
                structure.remove_oxidation_states()
                # to make sure all the atoms of the same element are listed
                # consecutively, so they get printed consecutively in the
                # poscar file
                structure.sort()

            # if num_adds is negative and the structure contains enough atoms,
            # randomly remove the right number of atoms of each element
            elif num_add < 0 and -1*total_add < len(structure.sites):
                site_indices_to_remove = []
                for key in amounts_to_add:
                    for _ in range(0, -1*amounts_to_add[key]):
                        # pick a random site in the structure that has the
                        # right element and that hasn't already been picked
                        random_site = random.choice(structure.sites)
                        while str(random_site.specie.symbol) != str(
                                key) or structure.sites.index(
                                    random_site) in site_indices_to_remove:
                            random_site = random.choice(structure.sites)
                        # record the index to of the site that has been
                        # designated for removal
                        site_indices_to_remove.append(
                            structure.sites.index(random_site))
                # remove the chosen sites
                structure.remove_sites(site_indices_to_remove)

        # optionally scale the volume after atoms have been added or removed
        if self.scale_volume:
            structure.scale_lattice(vol_per_atom*len(structure.sites))

        # create a new organism from the structure and return it
        offspring = Organism(structure, id_generator, self.name)

        print('Creating offspring organism {} from parent organism {} with '
              'the number of stoichiometries mutation variation '.format(
                  offspring.id, parent_org[0].id))

        return offspring


class Permutation(object):
    """
    An operator that creates an offspring organism by swapping atomic species
    in a parent organism.
    """

    def __init__(self, permutation_params, composition_space):
        """
        Creates a Permutation operator.

        Args:
            permutation_params: The parameters for doing the permutation
                operation, as a dictionary

            composition_space: the CompositionSpace object

        Precondition: the 'fraction' parameter in permutation_params is not
            optional, and it is assumed that permutation_params contains this
            parameter
        """

        # the name of this variation
        self.name = 'permutation'

        # parse the fraction from the parameters. This argument is not
        # optional, so it doesn't have a default value here
        self.fraction = permutation_params['fraction']

        # the max number of times to try getting a parent organism from the
        # pool that can undergo at least one swap. This could be a problem if,
        # e.g., the composition space is a single pure element but the
        # Permutation variation has non-zero fraction
        self.max_num_selects = 1000

        # the default values

        # the average number of pairs to swap
        self.default_mu_num_swaps = 2
        # the standard deviation of pairs to swap
        self.default_sigma_num_swaps = 1
        # which atomic pairs to swap
        self.default_pairs_to_swap = composition_space.get_all_swappable_pairs(
            )

        # the average number of swaps
        if 'mu_num_swaps' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.mu_num_swaps = self.default_mu_num_swaps
        elif permutation_params['mu_num_swaps'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.mu_num_swaps = self.default_mu_num_swaps
        else:
            # otherwise, parse the value from the parameters
            self.mu_num_swaps = permutation_params['mu_num_swaps']

        # the standard deviation of the number of swaps
        if 'sigma_num_swaps' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.sigma_num_swaps = self.default_sigma_num_swaps
        elif permutation_params['sigma_num_swaps'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.sigma_num_swaps = self.default_sigma_num_swaps
        else:
            # otherwise, parse the value from the parameters
            self.sigma_num_swaps = permutation_params['sigma_num_swaps']

        # which pairs of atoms to swap
        if 'pairs_to_swap' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.pairs_to_swap = self.default_pairs_to_swap
        elif permutation_params['pairs_to_swap'] in (None, 'default'):
            # use the default value if the flag was left blank or set to
            # 'default'
            self.pairs_to_swap = self.default_pairs_to_swap
        else:
            # otherwise, parse the value from the parameters
            self.pairs_to_swap = permutation_params['pairs_to_swap']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the permutation operation.

        Returns the resulting offspring as an Organism, or None if no offspring
        could be created.

        Args:
            pool: the Pool of Organisms

            random: Python's built in PRNG

            geometry: the Geometry object

            id_generator: the IDGenerator

        Description:

            Creates and offspring organism by swapping the elements of some of
            the sites in the parent structure.

                1. Selects a parent organism from the pool that is able to have
                    at least one of the allowed swaps done on it

                2. Computes the number of swaps to try to do by drawing from a
                    Gaussian with mean self.mu_num_swaps and standard
                    deviation self.sigma_num_swaps and rounding to the nearest
                    integer

                3. Tries to do the computed number of allowed swaps by randomly
                    electing an allowed pair to swap and then randomly
                    selecting sites in the structure with elements of the
                    allowed pair. This is repeated until either the computed
                    number of swaps have been done or no more swaps are
                    possible with the parent structure
        """

        # select a parent organism from the pool
        parent_org = pool.select_organisms(1, random)
        # make a deep copy of the structure of the parent organism
        structure = copy.deepcopy(parent_org[0].structure)
        # keep trying until we get a parent that has at least one possible swap
        possible_swaps = self.get_possible_swaps(structure)
        num_selects = 0
        while len(possible_swaps) == 0 and num_selects < self.max_num_selects:
            # select a parent organism from the pool
            parent_org = pool.select_organisms(1, random)
            # make a deep copy of the structure of the parent organism
            structure = copy.deepcopy(parent_org[0].structure)
            possible_swaps = self.get_possible_swaps(structure)
            num_selects = num_selects + 1

        # if the maximum number of selections have been made, then this isn't
        # working and it's time to stop
        if num_selects >= self.max_num_selects:
            return None

        # compute a positive random number of swaps to do
        num_swaps = int(round(random.gauss(self.mu_num_swaps,
                                           self.sigma_num_swaps)))
        while num_swaps <= 0:
            num_swaps = int(round(random.gauss(self.mu_num_swaps,
                                               self.sigma_num_swaps)))

        # try to select the computed number of swaps
        num_swaps_selected = 0
        pair_indices = []
        structure_to_check = copy.deepcopy(structure)
        # keep getting more swaps until either we've got enough or no more
        # swaps are possible
        while num_swaps_selected < num_swaps and len(possible_swaps) > 0:
            # pick a random pair to swap that we know is possible
            swap = random.choice(possible_swaps)
            symbols = swap.split()
            # keep trying until we find a site with the first element in the
            # pair
            site_1 = random.choice(structure_to_check.sites)
            while str(site_1.specie.symbol) != symbols[0]:
                site_1 = random.choice(structure_to_check.sites)
            # keep trying until we find a site with the second element in the
            # pair
            site_2 = random.choice(structure_to_check.sites)
            while str(site_2.specie.symbol) != symbols[1]:
                site_2 = random.choice(structure_to_check.sites)
            # record the indices (w.r.t. to the unchanged structure) for this
            # pair of sites
            pair_index = [structure.sites.index(site_1),
                          structure.sites.index(site_2)]
            pair_indices.append(pair_index)
            num_swaps_selected = num_swaps_selected + 1
            # remove these two sites from the structure to check
            structure_to_check.remove_sites(
                [structure_to_check.sites.index(site_1),
                 structure_to_check.sites.index(site_2)])
            # update the possible swaps
            possible_swaps = self.get_possible_swaps(structure_to_check)

        # do the swaps with the selected pairs
        for pair_index in pair_indices:
            species_1 = structure.sites[pair_index[0]].specie
            species_2 = structure.sites[pair_index[1]].specie
            structure.replace(pair_index[0], species_2)
            structure.replace(pair_index[1], species_1)

        # make a new organism from the structure and return it
        offspring = Organism(structure, id_generator, self.name)

        # print out a message
        print('Creating offspring organism {} from parent organism {} with '
              'the permutation variation '.format(offspring.id,
                                                  parent_org[0].id))
        return offspring

    def get_possible_swaps(self, structure):
        """
        Returns a list of swaps that are possible to do, based on what atoms
        are in the cell and which pairs are in self.pairs_to_swap. The returned
        list is a sublist of self.pairs_to_swap. Does not change the structure.

        Args:
            structure: the Structure object to check
        """

        possible_pairs = []
        for pair in self.pairs_to_swap:
            symbols = pair.split()
            # check if the structure contains the elements in the chosen pair
            has_element_1 = False
            has_element_2 = False
            for site in structure.sites:
                if str(site.specie.symbol) == symbols[0]:
                    has_element_1 = True
                if str(site.specie.symbol) == symbols[1]:
                    has_element_2 = True
            # append the pair if both elements were found in the structure
            if has_element_1 and has_element_2:
                possible_pairs.append(pair)
        return possible_pairs
