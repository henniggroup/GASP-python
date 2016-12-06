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
import warnings


class Mating(object):
    """
    An operator that creates an offspring organism by combining the structures
    of two parent organisms.
    """

    def __init__(self, mating_params):
        """
        Makes a Mating operator, and sets default parameter values if
        necessary.

        Args:
            mating_params: The parameters for doing the mating operation, as a
                dictionary.

        Precondition: the 'fraction' parameter in mating_params is not
            optional, and it is assumed that mating_params contains this
            parameter.
        """

        self.name = 'mating'
        self.fraction = mating_params['fraction']  # not optional

        # default values
        #
        # the average (fractional) location along the randomly chosen lattice
        # vector to make the cut
        self.default_mu_cut_loc = 0.5
        # the standard deviation of the (fractional) location along the
        # randomly chosen lattice vector to make the cut
        self.default_sigma_cut_loc = 0.5
        # the probability of randomly shifting the atoms along the lattice
        # vector of the cut before making the cut
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
            self.mu_cut_loc = self.default_mu_cut_loc
        elif mating_params['mu_cut_loc'] in (None, 'default'):
            self.mu_cut_loc = self.default_mu_cut_loc
        else:
            self.mu_cut_loc = mating_params['mu_cut_loc']

        # the standard deviation of the cut location
        if 'sigma_cut_loc' not in mating_params:
            self.sigma_cut_loc = self.default_sigma_cut_loc
        elif mating_params['sigma_cut_loc'] in (None, 'default'):
            self.sigma_cut_loc = self.default_sigma_cut_loc
        else:
            self.sigma_cut_loc = mating_params['sigma_cut_loc']

        # the probability of shifting the atoms in the parents
        if 'shift_prob' not in mating_params:
            self.shift_prob = self.default_shift_prob
        elif mating_params['shift_prob'] in (None, 'default'):
            self.shift_prob = self.default_shift_prob
        else:
            self.shift_prob = mating_params['shift_prob']

        # the probability of doubling one of the parents
        if 'doubling_prob' not in mating_params:
            self.doubling_prob = self.default_doubling_prob
        elif mating_params['doubling_prob'] in (None, 'default'):
            self.doubling_prob = self.default_doubling_prob
        else:
            self.doubling_prob = mating_params['doubling_prob']

        # whether to grow the smaller of the parents before doing the variation
        if 'grow_parents' not in mating_params:
            self.grow_parents = self.default_grow_parents
        elif mating_params['grow_parents'] in (None, 'default'):
            self.grow_parents = self.default_grow_parents
        else:
            self.grow_parents = mating_params['grow_parents']

        # the cutoff distance (as fraction of atomic radius) below which to
        # merge sites in the offspring structure
        if 'merge_cutoff' not in mating_params:
            self.merge_cutoff = self.default_merge_cutoff
        elif mating_params['merge_cutoff'] in (None, 'default'):
            self.merge_cutoff = self.default_merge_cutoff
        else:
            self.merge_cutoff = mating_params['merge_cutoff']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the mating operation.

        Returns the resulting offspring Organism.

        Args:
            pool: the Pool of Organisms

            random: copy of Python's built in PRNG

            geometry: the Geometry of the search

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

        Description:

            Creates an offspring organism by combining chunks cut from two
            parent organisms.

                1. Selects two organisms from the pool to act as parents, and
                    makes copies of their structures.

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
                    lattice vector to cut. For non-bulk geometries, the shift
                    only occurs if the lattice vector to cut is along a
                    periodic direction.

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

        # select two parent organisms from the pool and get their structures
        parent_orgs = pool.select_organisms(2, random)
        structure_1 = copy.deepcopy(parent_orgs[0].structure)
        structure_2 = copy.deepcopy(parent_orgs[1].structure)

        # optionally double one of the parents
        if random.random() < self.doubling_prob:
            vol_1 = structure_1.lattice.volume
            vol_2 = structure_2.lattice.volume
            if vol_1 < vol_2:
                self.double_parent(structure_1, geometry, random)
            else:
                self.double_parent(structure_2, geometry, random)

        # grow the smaller parent if specified
        if self.grow_parents:
            vol_1 = structure_1.lattice.volume
            vol_2 = structure_2.lattice.volume
            if vol_1 < vol_2:
                volume_ratio = vol_2/vol_1
                parent_to_grow = structure_1
            else:
                volume_ratio = vol_1/vol_2
                parent_to_grow = structure_2
            num_doubles = self.get_num_doubles(volume_ratio)
            for _ in range(num_doubles):
                self.double_parent(parent_to_grow, geometry, random)

        species_from_parent_1 = []
        species_from_parent_2 = []

        # make the cut
        # loop needed here because sometimes no atoms get contributed from one
        # of the parents, so have to try again
        while len(species_from_parent_1) == 0 or len(
                species_from_parent_2) == 0:
            species_from_parent_1 = []
            frac_coords_from_parent_1 = []
            species_from_parent_2 = []
            frac_coords_from_parent_2 = []

            # get the lattice vector to cut and the cut location
            cut_vector_index = random.randint(0, 2)
            cut_location = random.gauss(self.mu_cut_loc, self.sigma_cut_loc)
            while cut_location > 1 or cut_location < 0:
                cut_location = random.gauss(self.mu_cut_loc,
                                            self.sigma_cut_loc)

            # possibly shift the atoms in each parent along the cut vector
            if random.random() < self.shift_prob:
                self.do_random_shift(structure_1, cut_vector_index, geometry,
                                     random)
            if random.random() < self.shift_prob:
                self.do_random_shift(structure_2, cut_vector_index, geometry,
                                     random)

            # get the site contributions of each parent
            for site in structure_1.sites:
                if site.frac_coords[cut_vector_index] < cut_location:
                    species_from_parent_1.append(site.species_and_occu)
                    frac_coords_from_parent_1.append(site.frac_coords)
            for site in structure_2.sites:
                if site.frac_coords[cut_vector_index] > cut_location:
                    species_from_parent_2.append(site.species_and_occu)
                    frac_coords_from_parent_2.append(site.frac_coords)

        # combine the information for the sites contributed by each parent
        offspring_species = species_from_parent_1 + species_from_parent_2
        offspring_frac_coords = frac_coords_from_parent_1 + \
            frac_coords_from_parent_2

        # compute the lattice vectors of the offspring
        offspring_lengths = 0.5*(np.array(structure_1.lattice.abc) +
                                 np.array(structure_2.lattice.abc))
        offspring_angles = 0.5*(np.array(structure_1.lattice.angles) +
                                np.array(structure_2.lattice.angles))
        offspring_lattice = Lattice.from_lengths_and_angles(offspring_lengths,
                                                            offspring_angles)

        # make the offspring structure and merge close sites
        offspring_structure = Structure(offspring_lattice, offspring_species,
                                        offspring_frac_coords)
        offspring_structure = self.merge_sites(offspring_structure)

        # make the offspring organism from the offspring structure
        offspring = Organism(offspring_structure, id_generator, self.name)
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

    def double_parent(self, structure, geometry, random):
        """
        Modifies a structure by taking a supercell. For bulk geometries, the
        supercell is taken in the direction of the shortest lattice vector. For
        non-bulk geometries, the supercell is taken in the direction of a
        randomly chosen lattice vector.

        Args:
            structure: the Structure to take the supercell of

            geometry: the Geometry of the search

            random: a copy of Python's PRNG
        """

        if geometry.shape == 'bulk':
            lattice_lengths = structure.lattice.abc
            smallest_vector = min(lattice_lengths)
            doubling_index = lattice_lengths.index(smallest_vector)
        else:
            doubling_index = random.choice([0, 1, 2])
        scaling_factors = [1, 1, 1]
        scaling_factors[doubling_index] = 2
        structure.make_supercell(scaling_factors)

    def do_random_shift(self, structure, lattice_vector_index, geometry,
                        random):
        """
        Modifies a structure by shifting the atoms along one of the lattice
        vectors. Makes sure all the atoms lie inside the cell after the shift
        by replacing them with their periodic images if necessary.

        Note: this method checks the geometry, and will not do the shift if the
        specified lattice vector is not in a periodic direction because that
        could destroy some of the local structure.

        Args:
            structure: the Structure to shift

            lattice_vector_index: the index (0, 1 or 2) of the lattice vector
                along which to shift the atoms

            geometry: the Geometry of the search

            random: copy of Python's built in PRNG
        """

        if geometry.shape == 'cluster':
            pass
        elif geometry.shape == 'wire' and lattice_vector_index != 2:
            pass
        elif geometry.shape == 'sheet' and lattice_vector_index == 2:
            pass
        else:
            # do the shift
            shift_vector = random.random()*structure.lattice.matrix[
                lattice_vector_index]
            site_indices = [i for i in range(len(structure.sites))]
            structure.translate_sites(site_indices, shift_vector,
                                      frac_coords=False, to_unit_cell=False)

            # translate the sites back into the cell if needed
            translation_vectors = {}
            for site in structure.sites:
                translation_vector = []
                for coord in site.frac_coords:
                    if coord < 0.0:
                        translation_vector.append(1.0)
                    elif coord > 1.0:
                        translation_vector.append(-1.0)
                    else:
                        translation_vector.append(0.0)
                translation_vectors[
                    structure.sites.index(site)] = translation_vector
            for key in translation_vectors:
                structure.translate_sites(key, translation_vectors[key],
                                          frac_coords=True, to_unit_cell=False)

    def merge_sites(self, structure):
        """
        Merges sites in the structure that have the same element and are closer
        than self.merge_sites times the atomic radius to each other. Merging
        means replacing two sites in the structure with one site at the mean of
        the two sites (Cartesian) positions.

        Returns a new structure with the sites merged.

        Args:
            structure: the Structure whose sites to merge
        """

        species = []
        frac_coords = []
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
                            # make a new site if the two are close enough
                            if site.distance(
                                    other_site) < a_radius*self.merge_cutoff:
                                merged_indices.append(
                                    structure.sites.index(other_site))
                                merged_indices.append(
                                    structure.sites.index(site))
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
        Makes a Mutation operator, and sets default parameter values if
        necessary.

        Args:
            structure_mut_params: The parameters for doing the mutation
                operation, as a dictionary

        Precondition: the 'fraction' parameter in structure_mut_params is not
            optional, and it is assumed that structure_mut_params contains this
            parameter
        """

        self.name = 'structure mutation'
        self.fraction = structure_mut_params['fraction']  # not optional

        # the default values
        #
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
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        elif structure_mut_params['frac_atoms_perturbed'] in (None, 'default'):
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        else:
            self.frac_atoms_perturbed = structure_mut_params[
                'frac_atoms_perturbed']

        # the standard deviation of the perturbation of each atomic coordinate
        if 'sigma_atomic_coord_perturbation' not in structure_mut_params:
            self.sigma_atomic_coord_perturbation = \
                self.default_sigma_atomic_coord_perturbation
        elif structure_mut_params['sigma_atomic_coord_perturbation'] in (
                None, 'default'):
            self.sigma_atomic_coord_perturbation = \
                self.default_sigma_atomic_coord_perturbation
        else:
            self.sigma_atomic_coord_perturbation = structure_mut_params[
                'sigma_atomic_coord_perturbation']

        # the maximum allowed magnitude of perturbation of each atomic coord
        if 'max_atomic_coord_perturbation' not in structure_mut_params:
            self.max_atomic_coord_perturbation = \
                self.default_max_atomic_coord_perturbation
        elif structure_mut_params['max_atomic_coord_perturbation'] in (
                None, 'default'):
            self.max_atomic_coord_perturbation = \
                self.default_max_atomic_coord_perturbation
        else:
            self.max_atomic_coord_perturbation = structure_mut_params[
                'max_atomic_coord_perturbation']

        # the standard deviation of the magnitude of the non-identity
        # components of the elements in the strain matrix
        if 'sigma_strain_matrix_element' not in structure_mut_params:
            self.sigma_strain_matrix_element = \
                self.default_sigma_strain_matrix_element
        elif structure_mut_params['sigma_strain_matrix_element'] in (
                None, 'default'):
            self.sigma_strain_matrix_element = \
                self.default_sigma_strain_matrix_element
        else:
            self.sigma_strain_matrix_element = structure_mut_params[
                'sigma_strain_matrix_element']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the structure mutation operation.

        Returns the resulting offspring as an Organism.

         Args:
            pool: the Pool of Organisms

            random: a copy of Python's built in PRNG

            geometry: the Geometry of the search

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

        Description:

            Creates an offspring organism by perturbing the atomic positions
            and lattice vectors of the parent structure.

                1. Selects a parent organism from the pool and makes a copy of
                    its structure.

                2. Perturbs the atomic coordinates of each site with
                    probability self.frac_atoms_perturbed. The perturbation of
                    each atomic coordinate is drawn from a Gaussian with mean
                    zero and standard deviation
                    self.sigma_atomic_coord_perturbation. The magnitude of each
                    atomic coordinate perturbation is constrained to not exceed
                    self.max_atomic_coord_perturbation.

                3. Perturbs the lattice vectors by taking the product of each
                    lattice vector with a strain matrix. The strain matrix is

                        I + E

                    where I is the 3x3 identity matrix and E is a 3x3
                    perturbation matrix whose elements are distinct and drawn
                    from a Gaussian with mean zero and standard deviation
                    self.sigma_strain_matrix_element and constrained to lie
                    between -1 and 1.
        """

        # select a parent organism from the pool and get its structure
        parent_org = pool.select_organisms(1, random)
        structure = copy.deepcopy(parent_org[0].structure)

        # perturb the site coordinates
        self.perturb_atomic_coords(structure, random)

        # perturb the lattice vectors
        self.perturb_lattice_vectors(structure, random)

        # create a new organism from the perturbed structure
        offspring = Organism(structure, id_generator, self.name)
        offspring.translate_atoms_into_cell()
        print('Creating offspring organism {} from parent organism {} with '
              'the structure mutation variation '.format(offspring.id,
                                                         parent_org[0].id))
        return offspring

    def perturb_atomic_coords(self, structure, random):
        """
        Modifies a structure by perturbing the coordinates of its sites. The
        probability that each site is perturbed is self.frac_atoms_perturbed,
        and the perturbations along each Cartesian coordinate are drawn from a
        Gaussian with mean zero and standard deviation
        self.sigma_atomic_coord_perturbation. The magnitude of each atomic
        coordinate perturbation is constrained to not exceed
        self.max_atomic_coord_perturbation.

        Args:
            structure: the Structure whose site coordinates are perturbed

            random: a copy of Python's built in PRNG
        """

        # for each site in the structure, possibly randomly perturb it
        for site in structure.sites:
            if random.random() < self.frac_atoms_perturbed:
                # perturbation along x-coordinate
                nudge_x = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_x) > self.max_atomic_coord_perturbation:
                    nudge_x = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # perturbation along y-coordinate
                nudge_y = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_y) > self.max_atomic_coord_perturbation:
                    nudge_y = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # perturbation along z-coordinate
                nudge_z = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_z) > self.max_atomic_coord_perturbation:
                    nudge_z = random.gauss(
                        0, self.sigma_atomic_coord_perturbation)
                # translate the site by the random coordinate perturbations
                perturbation_vector = [nudge_x, nudge_y, nudge_z]
                structure.translate_sites(
                    structure.sites.index(site), perturbation_vector,
                    frac_coords=False, to_unit_cell=False)

    def perturb_lattice_vectors(self, structure, random):
        """
        Modifies a structure by perturbing its lattice vectors. Each lattice
        vector is  multiplied by a strain matrix:

                        I + E

        where I is the 3x3 identity matrix and E is a 3x3 perturbation matrix
        whose elements are distinct and drawn from a Gaussian with mean zero
        and standard deviation self.sigma_strain_matrix_element and constrained
        to lie between -1 and 1.

        Args:
            structure: the Structure whose site coordinates are perturbed

            random: a copy of Python's built in PRNG
        """

        # compute the random non-identity components of the nine elements of
        # the strain matrix, and make sure they're in [-1, 1]
        epsilons = []
        for _ in range(9):
            epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            while np.abs(epsilon) > 1:
                epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            epsilons.append(epsilon)

        # construct the strain matrix I + epsilon_ij
        row_1 = [1 + epsilons[0], epsilons[1], epsilons[2]]
        row_2 = [epsilons[3], 1 + epsilons[4], epsilons[5]]
        row_3 = [epsilons[6], epsilons[7], 1 + epsilons[8]]
        strain_matrix = np.array([row_1, row_2, row_3])

        # apply the strain matrix to the lattice vectors
        new_a = strain_matrix.dot(structure.lattice.matrix[0])
        new_b = strain_matrix.dot(structure.lattice.matrix[1])
        new_c = strain_matrix.dot(structure.lattice.matrix[2])
        new_lattice = Lattice([new_a, new_b, new_c])
        structure.modify_lattice(new_lattice)


class NumStoichsMut(object):
    """
    An operator that creates an offspring organism by mutating the number of
    stoichiometries' worth of atoms in the parent organism.
    """

    def __init__(self, num_stoichs_mut_params):
        """
        Makes a NumStoichsMut operator, and sets default parameter values if
        necessary.

        Args:
            num_stoichs_mut_params: The parameters for doing the NumStoichsMut
                operation, as a dictionary

        Precondition: the 'fraction' parameter in num_stoichs_mut_params is not
            optional, and it is assumed that num_stoichs_mut_params contains
            this parameter
        """

        self.name = 'number of stoichiometries mutation'
        self.fraction = num_stoichs_mut_params['fraction']  # not optional

        # the default values
        #
        # the average number of stoichimetries to add
        self.default_mu_num_adds = 0
        # the standard deviation of the number of stoichiometries to add
        self.default_sigma_num_adds = 1
        # whether to scale the volume of the offspring to equal that of the
        # parent
        self.default_scale_volume = True

        # the average number of stoichiometries to add
        if 'mu_num_adds' not in num_stoichs_mut_params:
            self.mu_num_adds = self.default_mu_num_adds
        elif num_stoichs_mut_params['mu_num_adds'] in (None, 'default'):
            self.mu_num_adds = self.default_mu_num_adds
        else:
            self.mu_num_adds = num_stoichs_mut_params['mu_num_adds']

        # the standard deviation of the number of stoichiometries
        if 'sigma_num_adds' not in num_stoichs_mut_params:
            self.sigma_num_adds = self.default_sigma_num_adds
        elif num_stoichs_mut_params['sigma_num_adds'] in (None, 'default'):
            self.sigma_num_adds = self.default_sigma_num_adds
        else:
            self.sigma_num_adds = num_stoichs_mut_params['sigma_num_adds']

        # whether to scale the volume of the offspring
        if 'scale_volume' not in num_stoichs_mut_params:
            self.scale_volume = self.default_scale_volume
        elif num_stoichs_mut_params['scale_volume'] in (None, 'default'):
            self.scale_volume = self.default_scale_volume
        else:
            self.scale_volume = num_stoichs_mut_params['scale_volume']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the number of stoichiometries mutation operation.

        Returns the resulting offspring as an Organism.

        Args:
            pool: the Pool of Organisms

            random: a copy of Python's built in PRNG

            geometry: the Geometry of the search

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

        Description:

            Creates an offspring organism by adding or removing a random number
            of stoichiometries' worth of atoms to or from the parent structure.

                1. Selects a parent organism from the pool and makes a copy of
                    its structure.

                2. Computes the number of stoichiometries to add or remove by
                    drawing from a Gaussian with mean self.mu_num_adds and
                    standard deviation self.sigma_num_adds and rounding the
                    result to the nearest integer.

                3. Computes the number of atoms of each type to add or remove,
                    and does the additions or removals.

                4. If self.scale_volume is True, scales the new structure to
                    have the same volume per atom as the parent.
        """

        # select a parent organism from the pool and get its structure
        parent_org = pool.select_organisms(1, random)
        structure = copy.deepcopy(parent_org[0].structure)
        parent_num_atoms = len(structure.sites)
        reduced_composition = structure.composition.reduced_composition
        vol_per_atom = structure.lattice.volume/len(structure.sites)

        # loop needed here in case no atoms got added or removed, in which case
        # need to try again. This happens if the randomly chosen number of
        # atoms to remove exceeds the number of atoms in the cell
        while len(structure.sites) == parent_num_atoms:
            # compute a non-zero number of stoichiometries to add or remove
            num_add = int(round(random.gauss(self.mu_num_adds,
                                             self.sigma_num_adds)))
            while num_add == 0:
                num_add = int(round(random.gauss(self.mu_num_adds,
                                                 self.sigma_num_adds)))

            # compute the number of each type of atom to add (or remove)
            amounts_to_add = {}
            total_add = 0
            for key in reduced_composition:
                amounts_to_add[key] = int(num_add*reduced_composition[key])
                total_add = total_add + int(num_add*reduced_composition[key])

            # if adding, put the new atoms in the cell at random locations
            if num_add > 0:
                for key in amounts_to_add:
                    for _ in range(amounts_to_add[key]):
                        frac_coords = [random.random(), random.random(),
                                       random.random()]
                        structure.append(Specie(key, 0), frac_coords)
                structure.remove_oxidation_states()
                structure.sort()

            # if removing, take out random atoms (but right number of each)
            elif num_add < 0 and -1*total_add < len(structure.sites):
                site_indices_to_remove = []
                for key in amounts_to_add:
                    for _ in range(0, -1*amounts_to_add[key]):
                        random_site = random.choice(structure.sites)
                        while str(random_site.specie.symbol) != str(
                                key) or structure.sites.index(
                                    random_site) in site_indices_to_remove:
                            random_site = random.choice(structure.sites)
                        site_indices_to_remove.append(
                            structure.sites.index(random_site))
                structure.remove_sites(site_indices_to_remove)

        # optionally scale the volume after atoms have been added or removed
        if self.scale_volume:
            # this is to suppress the warnings produced if the scale_lattice
            # method fails
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                structure.scale_lattice(vol_per_atom*len(structure.sites))
                if str(structure.lattice.a) == 'nan' or structure.lattice.a > \
                        100:
                    return self.do_variation(pool, random, geometry,
                                             id_generator)

        # create a new organism from the structure
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
        Makes a Permutation operator, and sets default parameter values if
        necessary.

        Args:
            permutation_params: The parameters for doing the permutation
                operation, as a dictionary

            composition_space: the CompositionSpace of the search

        Precondition: the 'fraction' parameter in permutation_params is not
            optional, and it is assumed that permutation_params contains this
            parameter
        """

        self.name = 'permutation'
        self.fraction = permutation_params['fraction']  # not optional

        # The max number of times to try getting a parent organism from the
        # pool that can undergo at least one swap. This is needed in case,
        # e.g., the composition space is a single pure element but the
        # Permutation variation has non-zero fraction.
        self.max_num_selects = 1000

        # the default values
        #
        # the average number of pairs to swap
        self.default_mu_num_swaps = 2
        # the standard deviation of pairs to swap
        self.default_sigma_num_swaps = 1
        # which atomic pairs to swap
        self.default_pairs_to_swap = composition_space.get_all_swappable_pairs(
            )

        # the average number of swaps
        if 'mu_num_swaps' not in permutation_params:
            self.mu_num_swaps = self.default_mu_num_swaps
        elif permutation_params['mu_num_swaps'] in (None, 'default'):
            self.mu_num_swaps = self.default_mu_num_swaps
        else:
            self.mu_num_swaps = permutation_params['mu_num_swaps']

        # the standard deviation of the number of swaps
        if 'sigma_num_swaps' not in permutation_params:
            self.sigma_num_swaps = self.default_sigma_num_swaps
        elif permutation_params['sigma_num_swaps'] in (None, 'default'):
            self.sigma_num_swaps = self.default_sigma_num_swaps
        else:
            self.sigma_num_swaps = permutation_params['sigma_num_swaps']

        # which pairs of atoms to swap
        if 'pairs_to_swap' not in permutation_params:
            self.pairs_to_swap = self.default_pairs_to_swap
        elif permutation_params['pairs_to_swap'] in (None, 'default'):
            self.pairs_to_swap = self.default_pairs_to_swap
        else:
            self.pairs_to_swap = permutation_params['pairs_to_swap']

    def do_variation(self, pool, random, geometry, id_generator):
        """
        Performs the permutation operation.

        Returns the resulting offspring as an Organism, or None if no offspring
        could be created.

        Args:
            pool: the Pool of Organisms

            random: a copy of Python's built in PRNG

            geometry: the Geometry of the search

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

        Description:

            Creates and offspring organism by swapping the elements of some of
            the sites in the parent structure.

                1. Selects a parent organism from the pool that is able to have
                    at least one of the allowed swaps done on it.

                2. Computes the number of swaps to try to do by drawing from a
                    Gaussian with mean self.mu_num_swaps and standard
                    deviation self.sigma_num_swaps and rounding to the nearest
                    integer.

                3. Tries to do the computed number of allowed swaps by randomly
                    electing an allowed pair to swap and then randomly
                    selecting sites in the structure with elements of the
                    allowed pair. This is repeated until either the computed
                    number of swaps have been done or no more swaps are
                    possible with the parent structure.
        """

        # select a parent organism from the pool and get its structure
        parent_org = pool.select_organisms(1, random)
        structure = copy.deepcopy(parent_org[0].structure)

        # keep trying until we get a parent that has at least one possible swap
        possible_swaps = self.get_possible_swaps(structure)
        num_selects = 0
        while len(possible_swaps) == 0 and num_selects < self.max_num_selects:
            parent_org = pool.select_organisms(1, random)
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

        # try to select the computed number of swaps - keep getting more until
        # either we've got enough or no more swaps are possible
        num_swaps_selected = 0
        pair_indices = []
        structure_to_check = copy.deepcopy(structure)
        while num_swaps_selected < num_swaps and len(possible_swaps) > 0:
            # pick a random pair to swap that we know is possible
            swap = random.choice(possible_swaps)
            symbols = swap.split()
            # find sites with the elements to swap
            site_1 = random.choice(structure_to_check.sites)
            while str(site_1.specie.symbol) != symbols[0]:
                site_1 = random.choice(structure_to_check.sites)
            site_2 = random.choice(structure_to_check.sites)
            while str(site_2.specie.symbol) != symbols[1]:
                site_2 = random.choice(structure_to_check.sites)
            # record the indices (w.r.t. to the unchanged structure) for this
            # pair of sites
            pair_index = [structure.sites.index(site_1),
                          structure.sites.index(site_2)]
            pair_indices.append(pair_index)
            num_swaps_selected += 1
            # remove these two sites from the structure to check
            structure_to_check.remove_sites(
                [structure_to_check.sites.index(site_1),
                 structure_to_check.sites.index(site_2)])
            possible_swaps = self.get_possible_swaps(structure_to_check)

        # do the swaps with the selected pairs
        for pair_index in pair_indices:
            species_1 = structure.sites[pair_index[0]].specie
            species_2 = structure.sites[pair_index[1]].specie
            structure.replace(pair_index[0], species_2)
            structure.replace(pair_index[1], species_1)

        # make a new organism from the structure
        offspring = Organism(structure, id_generator, self.name)
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
            has_element_1 = False
            has_element_2 = False
            for site in structure.sites:
                if str(site.specie.symbol) == symbols[0]:
                    has_element_1 = True
                if str(site.specie.symbol) == symbols[1]:
                    has_element_2 = True
            if has_element_1 and has_element_2:
                possible_pairs.append(pair)
        return possible_pairs
