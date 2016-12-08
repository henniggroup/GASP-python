# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Development module:

This module contains the classes used when developing an organism
before and after it is submitted for an energy calculation.

1. Constraints: contains the constraints placed on structures

2. Developer: develops organisms

3. RedundancyGuard: checks if an organism is redundant

4. Geometry: contains geometry-specific data for non-bulk searches, including
    any additional constraints

"""

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element, DummySpecie
from pymatgen.core.sites import Site
from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_matcher import ElementComparator

import warnings
import numpy as np


class Constraints(object):
    '''
    Represents the general constraints imposed on structures considered by the
    algorithm.
    '''

    def __init__(self, constraints_parameters, composition_space):
        '''
        Makes a Constraints, and sets default parameter values if necessary.

        Args:
            constraints_parameters: a dictionary of parameters

            composition_space: the CompositionSpace of the search
        '''

        # default values
        self.default_min_num_atoms = 2
        self.default_max_num_atoms = 50
        self.default_min_lattice_length = 0.5
        self.default_max_lattice_length = 20
        self.default_min_lattice_angle = 40
        self.default_max_lattice_angle = 140
        self.default_mid_factor = 0.7

        # set to defaults
        if constraints_parameters in (None, 'default'):
            self.set_all_to_defaults(composition_space)
        # parse the parameters and set to defaults if necessary
        else:
            # min number of atoms
            if 'min_num_atoms' not in constraints_parameters:
                self.min_num_atoms = self.default_min_num_atoms
            elif constraints_parameters['min_num_atoms'] in (None, 'default'):
                self.min_num_atoms = self.default_min_num_atoms
            else:
                self.min_num_atoms = constraints_parameters['min_num_atoms']

            # max number of atoms
            if 'max_num_atoms' not in constraints_parameters:
                self.max_num_atoms = self.default_max_num_atoms
            elif constraints_parameters['max_num_atoms'] in (None, 'default'):
                self.max_num_atoms = self.default_max_num_atoms
            else:
                self.max_num_atoms = constraints_parameters['max_num_atoms']

            # min lattice length
            if 'min_lattice_length' not in constraints_parameters:
                self.min_lattice_length = self.default_min_lattice_length
            elif constraints_parameters['min_lattice_length'] in (None,
                                                                  'default'):
                self.min_lattice_length = self.default_min_lattice_length
            else:
                self.min_lattice_length = constraints_parameters[
                    'min_lattice_length']

            # max lattice length
            if 'max_lattice_length' not in constraints_parameters:
                self.max_lattice_length = self.default_max_lattice_length
            elif constraints_parameters['max_lattice_length'] in (None,
                                                                  'default'):
                self.max_lattice_length = self.default_max_lattice_length
            else:
                self.max_lattice_length = \
                    constraints_parameters['max_lattice_length']

            # min lattice angle
            if 'min_lattice_angle' not in constraints_parameters:
                self.min_lattice_angle = self.default_min_lattice_angle
            elif constraints_parameters['min_lattice_angle'] in (None,
                                                                 'default'):
                self.min_lattice_angle = self.default_min_lattice_angle
            else:
                self.min_lattice_angle = \
                    constraints_parameters['min_lattice_angle']

            # max lattice angle
            if 'max_lattice_angle' not in constraints_parameters:
                self.max_lattice_angle = self.default_max_lattice_angle
            elif constraints_parameters['max_lattice_angle'] in (None,
                                                                 'default'):
                self.max_lattice_angle = self.default_max_lattice_angle
            else:
                self.max_lattice_angle = constraints_parameters[
                    'max_lattice_angle']

            # the per-species minimum interatomic distances
            if 'per_species_mids' not in constraints_parameters:
                self.set_all_mids_to_defaults(composition_space)
            elif constraints_parameters['per_species_mids'] in (None,
                                                                'default'):
                self.set_all_mids_to_defaults(composition_space)
            else:
                self.per_species_mids = constraints_parameters[
                    'per_species_mids']
                # check if any of the specified pairs needs a default mid
                for key in self.per_species_mids:
                    if self.per_species_mids[key] in (None, 'default'):
                        elements = key.split()
                        radius1 = Element(elements[0]).atomic_radius
                        radius2 = Element(elements[1]).atomic_radius
                        self.per_species_mids[key] = self.default_mid_factor*(
                            radius1 + radius2)
                # check for missing pairs, and set default mids for them
                self.set_some_mids_to_defaults(composition_space)

    def set_all_to_defaults(self, composition_space):
        '''
        Sets all general constraints (those in Constraints block of input file)
        to default values.

        Args:
            composition_space: the CompositionSpace of the search
        '''

        self.min_num_atoms = self.default_min_num_atoms
        self.max_num_atoms = self.default_max_num_atoms
        self.min_lattice_length = self.default_min_lattice_length
        self.max_lattice_length = self.default_max_lattice_length
        self.min_lattice_angle = self.default_min_lattice_angle
        self.max_lattice_angle = self.default_max_lattice_angle
        self.set_all_mids_to_defaults(composition_space)

    def set_all_mids_to_defaults(self, composition_space):
        '''
        Sets all the per-species mids to default values based on atomic radii.

        Args:
            composition_space: the CompositionSpace of the search
        '''

        elements = composition_space.get_all_elements()
        self.per_species_mids = {}
        for i in range(0, len(elements)):
            for j in range(i, len(elements)):
                self.per_species_mids[
                    str(elements[i].symbol + " " + elements[j].symbol)
                    ] = self.default_mid_factor*(elements[i].atomic_radius +
                                                 elements[j].atomic_radius)

    def set_some_mids_to_defaults(self, composition_space):
        '''
        Compares all the possible pairs of elements to what is contained in
        self.per_species_mids. If any pairs are missing, adds them with default
        values.

        Args:
            composition_space: the CompositionSpace of the search
        '''

        # find all the missing pairs
        elements = composition_space.get_all_elements()
        missing_pairs = []
        for i in range(0, len(elements)):
            for j in range(i, len(elements)):
                # check both possible orders
                test_key1 = elements[i].symbol + " " + elements[j].symbol
                test_key2 = elements[j].symbol + " " + elements[i].symbol
                if test_key1 not in self.per_species_mids and test_key2 \
                        not in self.per_species_mids:
                    missing_pairs.append(test_key1)

        # calculate the per species mids for all the missing pairs
        for pair in missing_pairs:
            p = pair.split()
            self.per_species_mids[str(pair)] = self.default_mid_factor*(
                Element(p[0]).atomic_radius + Element(p[1]).atomic_radius)

    def get_max_mid(self):
        '''
        Returns largest per-species minimum interatomic distance constraint.
        '''

        max_mid = 0
        for key in self.per_species_mids:
            if self.per_species_mids[key] > max_mid:
                max_mid = self.per_species_mids[key]
        return max_mid


class Developer(object):
    '''
    A Developer object is used to develop an organism before evaluating its
    energy or adding it to the pool or initial population. Doesn't do
    redundancy checking.
    '''

    def __init__(self, developer_parameters, geometry):
        '''
        Makes a Developer, and sets default parameter values if necessary.

        Args:
            niggli: a boolean indicating whether or not to do Niggli cell
                reduction

            scale_density: a boolean indicating whether or not to scale the
                density

            geometry: the Geometry of the search
        '''

        # defaults
        self.default_niggli = True
        if geometry.shape == 'bulk':
            self.default_scale_density = True
        else:
            self.default_scale_density = False

        # set to defaults
        if developer_parameters in (None, 'default'):
            self.niggli = self.default_niggli
            self.scale_density = self.default_scale_density
        # parse the parameters and set to defaults if necessary
        else:
            # niggli
            if 'niggli' not in developer_parameters:
                self.niggli = self.default_niggli
            elif developer_parameters['niggli'] in (None, 'default'):
                self.niggli = self.default_niggli
            else:
                self.niggli = developer_parameters['niggli']

            # scale density
            if 'scale_density' not in developer_parameters:
                self.scale_density = self.default_scale_density
            elif developer_parameters['scale_density'] in (None, 'default'):
                self.scale_density = self.default_scale_density
            else:
                self.scale_density = developer_parameters['scale_density']

    def develop(self, organism, composition_space, constraints, geometry,
                pool):
        '''
        Develops an organism. Can modify the structure of an organism through
        Niggli cell reduction and volume scaling.

        Returns a boolean indicating whether the organism survived development.

        Args:
            organism: the Organism to develop

            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

            geometry: the Geometry of the search

            pool: the Pool
        '''

        # check the constraints on the number of atoms
        if not self.satisfies_num_atoms_constraints(organism, constraints):
            return False

        # check if the organism is is the composition space
        if not self.is_in_composition_space(organism, composition_space, pool):
            return False

        # optionally do Niggli cell reduction
        if self.niggli:
            if not self.niggli_reduction(organism, geometry):
                return False

        # optionally scale the volume per atom if the organism is unrelaxed
        if self.scale_density and len(
                pool.promotion_set) > 0 and organism.epa is None:
            if not self.scale_volume(organism, composition_space, pool):
                return False

        # check the lattice length and angle constraints
        if not self.satisfies_lattice_constraints(organism, constraints):
            return False

        # check the per-species minimum interatomic distance constraints
        if not self.satisfies_mids_constraints(organism, constraints):
            return False

        # check any geometry-specific constraints
        if not self.satisfies_geometry_constraints(organism, geometry):
            return False

        return True

    def satisfies_num_atoms_constraints(self, organism, constraints):
        """
        Returns a boolean indicating whether the organism satisfies the
        constraints on the number of atoms.

        Args:
            organism: the Organism to check

            constraints: the Constraints of the search
        """

        # check max num atoms constraint
        if len(organism.structure.sites) > constraints.max_num_atoms:
            print("Organism {} failed max number of atoms constraint ".format(
                organism.id))
            return False

        # check min num atoms constraint
        if len(organism.structure.sites) < constraints.min_num_atoms:
            print("Organism {} failed min number of atoms constraint ".format(
                organism.id))
            return False
        return True

    def is_in_composition_space(self, organism, composition_space, pool):
        """
        Returns a boolean indicating whether the organism is in the composition
        space.

        Args:
            organism: the Organism to check

            composition_space: the CompositionSpace of the search

            pool: the Pool
        """

        # for epa searches
        if composition_space.objective_function == 'epa':
            return self.is_in_composition_space_epa(organism,
                                                    composition_space)
        # for pd searches
        elif composition_space.objective_function == "pd":
            return self.is_in_composition_space_pd(organism, composition_space,
                                                   pool)

    def is_in_composition_space_epa(self, organism, composition_space):
        """
        Returns a boolean indicating whether the organism has the required
        composition (for fixed-compsition searches).

        Args:
            organism: the Organism to check

            composition_space: the CompositionSpace of the search
        """

        reduced_composition = composition_space.endpoints[
            0].reduced_composition
        org_reduced_composition = organism.composition.reduced_composition
        if not reduced_composition.almost_equals(org_reduced_composition):
            print("Organism {} has incorrect composition ".format(organism.id))
            return False
        return True

    def is_in_composition_space_pd(self, organism, composition_space, pool):
        """
        Returns a boolean indicating whether the organism is in the composition
        space. If the initial population is finished, then returns False for
        organisms with endpoint compositions.

        Args:
            organism: the Organism to check

            composition_space: the CompositionSpace of the search

            pool: the Pool
        """

        # cast the endpoints to PDEntries (just make up some energies)
        pdentries = []
        for endpoint in composition_space.endpoints:
            pdentries.append(PDEntry(endpoint, -10))
        pdentries.append(PDEntry(organism.composition, -10))

        # make a CompoundPhaseDiagram and use it to check if the organism
        # is in the composition space from how many entries it returns
        composition_checker = CompoundPhaseDiagram(
            pdentries, composition_space.endpoints)
        if len(composition_checker.transform_entries(
                    pdentries, composition_space.endpoints)[0]) == len(
                        composition_space.endpoints):
            print("Organism {} lies outside the composition space ".format(
                organism.id))
            return False

        # check the endpoints if we're not making the initial population
        elif len(pool.to_list()) > 0:
            for endpoint in composition_space.endpoints:
                if endpoint.almost_equals(
                        organism.composition.reduced_composition):
                    print('Organism {} is at a composition endpoint '.format(
                        organism.id))
                    return False
        return True

    def niggli_reduction(self, organism, geometry):
        """
        Returns a boolean indicating whether Niggli cell reduction did not
        fail.

        Args:
            organism: the Organism whose structure to Niggli reduce

            geometry: the Geometry of the search
        """

        if geometry.shape == 'bulk':
            try:  # sometimes pymatgen's reduction routine fails
                organism.structure = \
                    organism.structure.get_reduced_structure()
                organism.rotate_to_principal_directions()
            except:
                print('Niggli cell reduction failed on organism {} during '
                      'development '.format(organism.id))
                return False
        elif geometry.shape == 'sheet':
            try:
                organism.reduce_sheet_cell()
                organism.rotate_to_principal_directions()
            except:
                print('2D Niggli cell reduction failed on organism {} '
                      'during development '.format(organism.id))
                return False
            # TODO: call special cell reduction for other geometries here if
            # needed (doesn't makes sense for wires or clusters)
        return True

    def scale_volume(self, organism, composition_space, pool):
        """
        Returns a boolean indicating whether volume scaling did not fail.

        Args:
            organism: the Organism whose volume to scale

            composition_space: the CompositionSpace of the search

            pool: the Pool
        """

        if composition_space.objective_function == 'epa':
            return self.scale_volume_epa(organism, pool)
        elif composition_space.objective_function == 'pd':
            return self.scale_volume_pd(organism, composition_space, pool)

    def scale_volume_epa(self, organism, pool):
        """
        Returns a boolean indicating whether volume scaling did not fail.

        Args:
            organism: the Organism whose volume to scale

            pool: the Pool

        Description:
            Scales the volume per atom of the organism to the average volume
            per atom of the organisms in the promotion set, plus 10%.
        """

        # compute the voluem to scale to
        vpa_sum = 0
        for org in pool.promotion_set:
            vpa_sum += org.structure.volume/len(org.structure.sites)
        vpa_mean = 1.1*(vpa_sum/len(pool.promotion_set))
        num_atoms = len(organism.structure.sites)
        new_vol = vpa_mean*num_atoms
        # this is to suppress the warnings produced if the
        # scale_lattice method fails
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            organism.structure.scale_lattice(new_vol)
            if str(organism.structure.lattice.a) == 'nan' or \
                    organism.structure.lattice.a > 100:
                print('Volume scaling failed on organism {} during '
                      'development '.format(organism.id))
                return False
        return True

    def scale_volume_pd(self, organism, composition_space, pool):
        """
        Returns a boolean indicating whether volume scaling did not fail.

        Args:
            organism: the Organism whose volume to scale

            composition_space: the CompositionSpace of the search

            pool: the Pool

        Description:

            1. Computes the decomposition of the organism - that is, which
                structures on the convex hull (and their relative amounts) that
                the organism would decompose to, based on its composition.

            2. Computes the weighted average volume per atom of the structures
                in the composition.

            3. Scales the volume per atom of the organism to the computed
                value, plus 10%.
        """

        # make CompoundPhaseDiagram and PDAnalyzer objects
        pdentries = []
        for org in pool.promotion_set:
            pdentries.append(PDEntry(org.composition, org.total_energy))
        compound_pd = CompoundPhaseDiagram(pdentries,
                                           composition_space.endpoints)
        pdanalyzer = PDAnalyzer(compound_pd)

        # transform the organism's composition
        transformed_entry = compound_pd.transform_entries(
            [PDEntry(organism.composition, 10)], composition_space.endpoints)

        # get the transformed species and amounts
        transformed_list = str(transformed_entry[0][0]).split()
        del transformed_list[0]
        popped = ''
        while popped != 'with':
            popped = transformed_list.pop()

        # separate the dummy species symbols from the amounts
        symbols = []
        amounts = []
        for entry in transformed_list:
            split_entry = entry.split('0+')
            symbols.append(split_entry[0])
            amounts.append(float(split_entry[1]))

        # make a dictionary mapping dummy species to amounts
        dummy_species_amounts = {}
        for i in range(len(symbols)):
            dummy_species_amounts[DummySpecie(symbol=symbols[i])] = amounts[i]

        # make Composition object with dummy species, get decomposition
        dummy_comp = Composition(dummy_species_amounts)
        decomp = pdanalyzer.get_decomposition(dummy_comp)

        # get original compositions and amounts from the decomposition
        fractions = []
        comps = []
        for item in decomp:
            fractions.append(decomp[item])
            first_split = str(item).split(',')
            second_split = first_split[0].split()
            while second_split[0] != 'composition':
                del second_split[0]
            del second_split[0]
            # build the composition string
            comp_string = ''
            for symbol in second_split:
                comp_string += str(symbol)
            comps.append(Composition(comp_string))

        # get weighted average volume per atom of the organisms in the
        # decomposition
        vpa_mean = 0
        for i in range(len(comps)):
            for org in pool.promotion_set:
                if (comps[i].reduced_composition).almost_equals(
                        org.composition.reduced_composition):
                    vpa_mean += (org.structure.volume/len(
                        org.structure.sites))*fractions[i]

        # compute the new volume and scale to it
        vpa_mean = 1.1*vpa_mean
        num_atoms = len(organism.structure.sites)
        new_vol = vpa_mean*num_atoms
        # this is to suppress the warnings produced if the
        # scale_lattice method fails
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            organism.structure.scale_lattice(new_vol)
            if str(organism.structure.lattice.a) == 'nan' or \
                    organism.structure.lattice.a > 100:
                print('Volume scaling failed on organism {} during '
                      'development '.format(organism.id))
                return False
        return True

    def satisfies_lattice_constraints(self, organism, constraints):
        """
        Returns a boolean indicating whether the organism satisfies the
        constraints on the lattice vector lengths and angles.

        Args:
            organism: the Organism to check

            constraints: the Constraints of the search
        """

        # check the max and min lattice length constraints
        lengths = organism.structure.lattice.abc
        for length in lengths:
            if length > constraints.max_lattice_length:
                print('Organism {} failed max lattice length '
                      'constraint '.format(organism.id))
                return False
            elif length < constraints.min_lattice_length:
                print('Organism {} failed min lattice length '
                      'constraint '.format(organism.id))
                return False

        # check the max and min lattice angle constraints
        angles = organism.structure.lattice.angles
        for angle in angles:
            if angle > constraints.max_lattice_angle:
                print('Organism {} failed max lattice angle '
                      'constraint '.format(organism.id))
                return False
            elif angle < constraints.min_lattice_angle:
                print('Organism {} failed min lattice angle '
                      'constraint '.format(organism.id))
                return False
        return True

    def satisfies_mids_constraints(self, organism, constraints):
        """
        Returns a boolean indicating whether the organism satisfies the
        per-species minimum interatomic distance constraints.

        Args:
            organism: the Organism to check

            constraints: the Constraints of the search
        """

        # check the per-species minimum interatomic distance constraints
        species_symbols = organism.structure.symbol_set
        for site in organism.structure.sites:
            for species_symbol in species_symbols:
                # We don't know the ordering in per_species_mids, so try both
                test_key1 = species_symbol + " " + site.specie.symbol
                test_key2 = site.specie.symbol + " " + species_symbol
                if test_key1 in constraints.per_species_mids:
                    mid = constraints.per_species_mids[test_key1]
                elif test_key2 in constraints.per_species_mids:
                    mid = constraints.per_species_mids[test_key2]
                # get all the sites within a sphere of radius mid centered on
                # the current site
                neighbors = organism.structure.get_neighbors(site, mid)
                # check each neighbor in the sphere to see if it has the
                # forbidden type
                for neighbor in neighbors:
                    if neighbor[0].specie.symbol == species_symbol:
                        print('Organism {} failed per-species minimum '
                              'interatomic distance constraint '.format(
                                  organism.id))
                        return False
        return True

    def satisfies_geometry_constraints(self, organism, geometry):
        """
        Returns a boolean indicating whether the organism satisfies the
        constraints associated with the geometry (max and min size, etc.).

        Args:
            organism: the Organism to check

            geometry: the Geometry of the search
        """

        # check the max size constraint (can only fail for non-bulk geometries)
        if geometry.get_size(organism) > geometry.max_size:
            print("Organism {} failed max size constraint ".format(
                organism.id))
            return False

        # check the min size constraint (can only fail for non-bulk geometries)
        if geometry.get_size(organism) < geometry.min_size:
            print("Organism {} failed min size constraint ".format(
                organism.id))
            return False
        return True


class RedundancyGuard(object):
    '''
    A RedundancyGuard object is used to check if an Organism is redundant with
    other organisms already seen by the alrogithm.
    '''

    def __init__(self, redundancy_parameters):
        '''
        Makes a RedundancyGuard, and sets default parameter values if
        necessary.

        TODO: pymatgen's structure matcher assumes periodic boundary
            conditions, but this doesn't make sense for all geometries...

        Args:
            redundancy parameters: a dictionary of parameters
        '''

        # defaults
        #
        # lattice length tolerance, in fractional coordinates
        self.default_lattice_length_tol = 0.05
        # lattice angle tolerance, in degrees
        self.default_lattice_angle_tol = 2
        # site tolerance, in fraction of average free length per atom
        self.default_site_tol = 0.1
        # whether to transform to primitive cells before comparing
        self.default_use_primitive_cell = True
        # whether to check if structures are equal to supercells of each other
        self.default_attempt_supercell = True
        # the epa difference interval
        self.default_epa_diff = 0.0

        # set to defaults
        if redundancy_parameters in (None, 'default'):
            self.set_all_to_defaults()
        # parse the parameters, and set to defaults if necessary
        else:
            # lattice length tolerance
            if 'lattice_length_tol' not in redundancy_parameters:
                self.lattice_length_tol = self.default_lattice_length_tol
            elif redundancy_parameters['lattice_length_tol'] in (None,
                                                                 'default'):
                self.lattice_length_tol = self.default_lattice_length_tol
            else:
                self.lattice_length_tol = redundancy_parameters[
                    'lattice_length_tol']

            # lattice angle tolerance
            if 'lattice_angle_tol' not in redundancy_parameters:
                self.lattice_angle_tol = self.default_lattice_angle_tol
            elif redundancy_parameters['lattice_angle_tol'] in (None,
                                                                'default'):
                self.lattice_angle_tol = self.default_lattice_angle_tol
            else:
                self.lattice_angle_tol = redundancy_parameters[
                    'lattice_angle_tol']

            # site tolerance
            if 'site_tol' not in redundancy_parameters:
                self.site_tol = self.default_site_tol
            elif redundancy_parameters['site_tol'] in (None, 'default'):
                self.site_tol = self.default_site_tol
            else:
                self.site_tol = redundancy_parameters['site_tol']

            # whether to use primitive cells
            if 'use_primitive_cell' not in redundancy_parameters:
                self.use_primitive_cell = self.default_use_primitive_cell
            elif redundancy_parameters['use_primitive_cell'] in (None,
                                                                 'default'):
                self.use_primitive_cell = self.default_use_primitive_cell
            else:
                self.use_primitive_cell = redundancy_parameters[
                    'use_primitive_cell']

            # whether to try matching supercells
            if 'attempt_supercell' not in redundancy_parameters:
                self.attempt_supercell = self.default_attempt_supercell
            elif redundancy_parameters['attempt_supercell'] in (None,
                                                                'default'):
                self.attempt_supercell = self.default_attempt_supercell
            else:
                self.attempt_supercell = redundancy_parameters[
                    'attempt_supercell']

            # epa difference
            if 'epa_diff' not in redundancy_parameters:
                self.epa_diff = self.default_epa_diff
            elif redundancy_parameters['epa_diff'] in (None, 'default'):
                self.epa_diff = self.default_epa_diff
            else:
                self.epa_diff = redundancy_parameters['epa_diff']

        # make the StructureMatcher object
        #
        # the first False is to prevent the matcher from scaling the volumes,
        # and the second False is to prevent subset matching
        self.structure_matcher = StructureMatcher(
            self.lattice_length_tol, self.site_tol, self.lattice_angle_tol,
            self.use_primitive_cell, False, self.attempt_supercell, False,
            ElementComparator())

    def set_all_to_defaults(self):
        '''
        Sets all the redundancy parameters to default values.
        '''

        self.lattice_length_tol = self.default_lattice_length_tol
        self.lattice_angle_tol = self.default_lattice_angle_tol
        self.site_tol = self.default_site_tol
        self.use_primitive_cell = self.default_use_primitive_cell
        self.attempt_supercell = self.default_attempt_supercell
        self.epa_diff = self.default_epa_diff

    def check_redundancy(self, new_organism, orgs_list):
        '''
        Checks for redundancy, both structural and if specified, epa (d-value).

        Returns the organism with which new_organism is redundant, or None if
        no redundancy.

        Args:
            new_organism: the Organism to check for redundancy

            orgs_list: the list containing all Organisms to check against
        '''

        # if new_organism isn't relaxed, then just check structures
        if new_organism.epa is None:
            for organism in orgs_list:
                if new_organism.id != organism.id:  # just in case
                    # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure,
                                                  organism.structure):
                        print('Organism {} failed structural redundancy - '
                              'looks like organism {} '.format(new_organism.id,
                                                               organism.id))
                        return organism

        # if new_organism is relaxed, only check against relaxed organisms
        else:
            for organism in orgs_list:
                if new_organism.id != organism.id and organism.epa is not None:
                    # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure,
                                                  organism.structure):
                        print('Organism {} failed structural redundancy - '
                              'looks like organism {} '.format(new_organism.id,
                                                               organism.id))
                        return organism
                    # if specified, check how close their epa's are
                    if self.epa_diff > 0:
                        if abs(new_organism.epa -
                               organism.epa) < self.epa_diff:
                            print('Organism {} failed energy per atom '
                                  'redundancy - looks like organism '
                                  '{} '.format(new_organism.id, organism.id))
                            return organism
        return None


class Geometry(object):
    '''
    Represents the geometry data, including any geometry-specific constraints
    (max_size, etc.).
    '''

    def __init__(self, geometry_parameters):
        '''
        Makes a Geometry, and sets default parameter values if necessary.

        Args:
            geometry_parameters: a dictionary of parameters
        '''

        # default values
        self.default_shape = 'bulk'
        self.default_max_size = np.inf
        self.default_min_size = -np.inf
        self.default_padding = 10  # this is only used for non-bulk shapes

        # set to defaults
        if geometry_parameters in (None, 'default'):
            self.shape = self.default_shape
            self.max_size = self.default_max_size
            self.min_size = self.default_min_size
            self.padding = None
        # parse the parameters, and set to defaults if necessary
        else:
            if 'shape' not in geometry_parameters:
                self.shape = self.default_shape
                self.max_size = self.default_max_size
                self.min_size = self.default_min_size
                self.padding = None
            elif geometry_parameters['shape'] in (None, 'default'):
                self.shape = self.default_shape
                self.max_size = self.default_max_size
                self.min_size = self.default_min_size
                self.padding = None
            else:
                self.shape = geometry_parameters['shape']

                # max size
                if 'max_size' not in geometry_parameters:
                    self.max_size = self.default_max_size
                elif geometry_parameters['max_size'] in (None, 'default'):
                    self.max_size = self.default_max_size
                else:
                    self.max_size = geometry_parameters['max_size']

                # min size
                if 'min_size' not in geometry_parameters:
                    self.min_size = self.default_min_size
                elif geometry_parameters['min_size'] in (None, 'default'):
                    self.min_size = self.default_min_size
                else:
                    self.min_size = geometry_parameters['min_size']

                # padding
                if 'padding' not in geometry_parameters:
                    self.padding = self.default_padding
                elif geometry_parameters['padding'] in (None, 'default'):
                    self.padding = self.default_padding
                else:
                    self.padding = geometry_parameters['padding']

    def pad(self, organism):
        '''
        Modifies the structure of an organism to make it conform to the
        required shape. For sheet, wire and cluster geometries, this means
        adding vacuum padding to the cell. For bulk, the structure is
        unchanged. Used to pad a structure prior to an energy calculation.

        Args:
            organism: the Organism to pad
        '''

        # call appropriate padding algorithm
        if self.shape == 'sheet':
            self.pad_sheet(organism)
        elif self.shape == 'wire':
            self.pad_wire(organism)
        elif self.shape == 'cluster':
            self.pad_cluster(organism)

    def pad_sheet(self, organism):
        '''
        Modifies the structure of an organism by adding vertical vacuum padding
        and making the c-lattice vector normal to the plane of the sheet. The
        atoms are shifted to the center of the padded sheet.

        Args:
            organism: the Organism to pad
        '''

        # make the padded structure
        organism.rotate_to_principal_directions()
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        cart_bounds = organism.get_bounding_box(cart_coords=True)
        minz = cart_bounds[2][0]
        maxz = cart_bounds[2][1]
        layer_thickness = maxz - minz
        ax = organism.structure.lattice.matrix[0][0]
        bx = organism.structure.lattice.matrix[1][0]
        by = organism.structure.lattice.matrix[1][1]
        padded_lattice = Lattice([[ax, 0.0, 0.0], [bx, by, 0.0],
                                  [0.0, 0.0, layer_thickness + self.padding]])
        padded_structure = Structure(padded_lattice, species, cartesian_coords,
                                     coords_are_cartesian=True)
        organism.structure = padded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the vertical center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] -
                                        frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def pad_wire(self, organism):
        '''
        Modifies the structure of an organism by making the c lattice vector
        parallel to z-axis, and adds vacuum padding around the structure in the
        x and y directions by replacing a and b lattice vectors with padded
        vectors along the x and y axes, respectively. The atoms are shifted to
        the center of the padded cell.

        Args:
            organism: the Organism to pad
        '''

        # make the padded structure
        organism.rotate_c_parallel_to_z()
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        cart_bounds = organism.get_bounding_box(cart_coords=True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        cz = organism.structure.lattice.matrix[2][2]
        padded_lattice = Lattice([[x_extent + self.padding, 0, 0],
                                  [0, y_extent + self.padding, 0], [0, 0, cz]])
        padded_structure = Structure(padded_lattice, species, cartesian_coords,
                                     coords_are_cartesian=True)
        organism.structure = padded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the horizontal center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] -
                                        frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] -
                                        frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def pad_cluster(self, organism):
        '''
        Modifies the structure of an organism by replacing the three lattice
        vectors with ones along the three Cartesian directions and adding
        vacuum padding to each one.

        Args:
            organism: the Organism to pad
        '''

        # make the padded structure
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        cart_bounds = organism.get_bounding_box(cart_coords=True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        z_min = cart_bounds[2][0]
        z_max = cart_bounds[2][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min
        padded_lattice = Lattice([[x_extent + self.padding, 0, 0],
                                  [0, y_extent + self.padding, 0],
                                  [0, 0, z_extent + self.padding]])
        padded_structure = Structure(padded_lattice, species, cartesian_coords,
                                     coords_are_cartesian=True)
        organism.structure = padded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] -
                                        frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] -
                                        frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] -
                                        frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def unpad(self, organism, constraints):
        '''
        Modifies the structure of an organism by removing vacuum padding, which
        returns an organism's structure to a form used by the variations. Does
        nothing if shape is bulk.

        Args:
            organism: the Organism to unpad

            constraints: the Constraints of the search
        '''

        # call the appropriate unpadding algorithm
        if self.shape == 'sheet':
            self.unpad_sheet(organism, constraints)
        elif self.shape == 'wire':
            self.unpad_wire(organism, constraints)
        elif self.shape == 'cluster':
            self.unpad_cluster(organism, constraints)

    def unpad_sheet(self, organism, constraints):
        '''
        Modifies the structure of an organism by removing vertical vacuum
        padding, leaving only enough to satisfy the per-species MID
        constraints, and makes the c-lattice vector normal to the plane of the
        sheet (if it isn't already).

        Args:
            organism: the Organism to unpad

            constraints: the Constraints of the search
        '''

        # make the unpadded structure
        organism.rotate_to_principal_directions()
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        layer_thickness = self.get_layer_thickness(organism)
        max_mid = constraints.get_max_mid() + 0.01  # just to be safe...
        ax = organism.structure.lattice.matrix[0][0]
        bx = organism.structure.lattice.matrix[1][0]
        by = organism.structure.lattice.matrix[1][1]
        unpadded_lattice = Lattice([[ax, 0.0, 0.0], [bx, by, 0.0],
                                    [0.0, 0.0, layer_thickness + max_mid]])
        unpadded_structure = Structure(unpadded_lattice, species,
                                       cartesian_coords,
                                       coords_are_cartesian=True)
        organism.structure = unpadded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the vertical center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] -
                                        frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def unpad_wire(self, organism, constraints):
        '''
        Modifies the structure of an organism by removing horizontal vacuum
        padding around a wire, leaving only enough to satisfy the per-species
        MID constraints, and makes the three lattice vectors lie along the
        three Cartesian directions.


        Args:
            organism: the Organism to unpad

            constraints: the Constraints of the search
        '''

        # make the unpadded structure
        organism.rotate_c_parallel_to_z()
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        cart_bounds = organism.get_bounding_box(cart_coords=True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        cz = organism.structure.lattice.matrix[2][2]
        max_mid = constraints.get_max_mid() + 0.01  # just to be safe...
        unpadded_lattice = Lattice([[x_extent + max_mid, 0.0, 0.0],
                                    [0, y_extent + max_mid, 0.0],
                                    [0.0, 0.0, cz]])
        unpadded_structure = Structure(unpadded_lattice, species,
                                       cartesian_coords,
                                       coords_are_cartesian=True)
        organism.structure = unpadded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the horizontal center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] -
                                        frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] -
                                        frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def unpad_cluster(self, organism, constraints):
        '''
        Modifies the structure of an organism by removing vacuum padding in
        every direction, leaving only enough to satisfy the per-species MID
        constraints, and makes the three lattice vectors lie along the three
        Cartesian directions.

        Args:
            organism: the Organism to unpad

            constraints: the Constraints of the search
        '''

        # make the unpadded structure
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        cart_bounds = organism.get_bounding_box(cart_coords=True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        z_min = cart_bounds[2][0]
        z_max = cart_bounds[2][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min
        max_mid = constraints.get_max_mid() + 0.01  # just to be safe...
        unpadded_lattice = Lattice([[x_extent + max_mid, 0.0, 0.0],
                                    [0, y_extent + max_mid, 0.0],
                                    [0.0, 0.0, z_extent + max_mid]])
        unpadded_structure = Structure(unpadded_lattice, species,
                                       cartesian_coords,
                                       coords_are_cartesian=True)
        organism.structure = unpadded_structure

        # translate the atoms back into the cell if needed, and shift them to
        # the center
        organism.translate_atoms_into_cell()
        frac_bounds = organism.get_bounding_box(cart_coords=False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] -
                                        frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] -
                                        frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] -
                                        frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector,
                                           frac_coords=True,
                                           to_unit_cell=False)

    def get_size(self, organism):
        '''
        Returns the size of an organism with a non-bulk shape. Returns 0 if the
        shape is bulk.

        Args:
            organism: the Organism whose size to get
        '''

        # call the appropriate method
        if self.shape == 'bulk':
            return 0
        elif self.shape == 'sheet':
            return self.get_layer_thickness(organism)
        elif self.shape == 'wire':
            return self.get_wire_diameter(organism)
        elif self.shape == 'cluster':
            return self.get_cluster_diameter(organism)

    def get_layer_thickness(self, organism):
        '''
        Returns the layer thickness of a sheet structure, which is the maximum
        vertical distance between atoms in the cell.

        Precondition: the organism has already been put into sheet format (c
            lattice vector parallel to the z-axis and a and b lattice vectors
            in the x-y plane)
        '''

        cart_bounds = organism.get_bounding_box(cart_coords=True)
        layer_thickness = cart_bounds[2][1] - cart_bounds[2][0]
        return layer_thickness

    def get_wire_diameter(self, organism):
        '''
        Returns the diameter of a wire structure, defined as the maximum
        distance between atoms projected to the x-y plane.

        Precondition: the organism has already been put into wire format (c
            lattice vector is parallel to z-axis and a and b lattice vectors in
            the x-y plane), and all sites are located inside the cell (i.e.,
            have fractional coordinates between 0 and 1).
        '''

        max_distance = 0
        for site_i in organism.structure.sites:
            # make Site versions of each PeriodicSite so that the computed
            # distance won't include periodic images
            non_periodic_site_i = Site(site_i.species_and_occu,
                                       [site_i.coords[0], site_i.coords[1],
                                        0.0])
            for site_j in organism.structure.sites:
                non_periodic_site_j = Site(site_j.species_and_occu,
                                           [site_j.coords[0], site_j.coords[1],
                                            0.0])
                distance = non_periodic_site_i.distance(non_periodic_site_j)
                if distance > max_distance:
                    max_distance = distance
        return max_distance

    def get_cluster_diameter(self, organism):
        '''
        Returns the diameter of a cluster structure, defined as the maximum
        distance between atoms in the cell.

        Precondition: all sites are located inside the cell (i.e., have
            fractional coordinates between 0 and 1)
        '''

        max_distance = 0
        for site_i in organism.structure.sites:
            # make Site versions of each PeriodicSite so that the computed
            # distance won't include periodic images
            non_periodic_site_i = Site(site_i.species_and_occu, site_i.coords)
            for site_j in organism.structure.sites:
                non_periodic_site_j = Site(site_j.species_and_occu,
                                           site_j.coords)
                distance = non_periodic_site_i.distance(non_periodic_site_j)
                if distance > max_distance:
                    max_distance = distance
        return max_distance
