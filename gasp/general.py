# coding: utf-8
# Copyright(c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
General module:

This module contains several classes central to the algorithm.

1. IDGenerator: generates consecutive integer ID numbers for organisms

2. Organism: represents an organism, which consists primarily of a
         structure and an energy

3. Cell: represents a structure, and extends pymatgen.core.structure.Structure

4. InitialPopulation: represents the initial population of organisms,
         produced from one or more of the organism creators

5. Pool: represents the population of organisms, after the initial population

6. OffspringGenerator: high-level class for making offspring organisms

7. SelectionProbDist: specifies the distribution from which selection
        probabilities are drawn

8. CompositionSpace: specifies the composition or composition range

9. StoppingCriteria: specifies when a search should stop

"""

from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.transformations.standard_transformations import \
    RotationTransformation

import os
import copy
import math
import numpy as np
from collections import deque
from scipy.spatial.qhull import ConvexHull


class IDGenerator(object):
    """
    Generates successive integer ID numbers, starting from 1.
    """

    def __init__(self):
        """
        Makes an IDGenerator.
        """

        self.id = 0

    def make_id(self):
        """
        Returns the next id number.
        """

        self.id += 1
        return self.id


class Organism(object):
    """
    An organism, consisting primarily of a cell and an energy, as well as
    several derived quantities.
    """

    def __init__(self, cell, id_generator, maker):
        """
        Makes an Organism.

        Args:
            cell: the cell of this organism, as a Cell object

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

            maker: the name of algorithm that made the organism, as a string.
                Either a creator or a variation
        """

        self.cell = cell
        self.composition = self.cell.composition
        self.total_energy = None
        # objective function value
        self.value = None
        # energy per atom
        self.epa = None
        self.fitness = None
        # selection probability of this organism
        self.selection_prob = None
        # starting point of interval on number line of size selection_prob
        self.selection_loc = None
        # whether the organism is part of initial population or pool
        self.is_active = False
        # unique id number
        self._id = id_generator.make_id()
        # the name of the algorithm that created this organism
        self.made_by = maker

    # This keeps the id (sort of) immutable by causing an exception to be
    # raised if the id is attempted to be set with org.id = some_id.
    @property
    def id(self):
        return self._id


class Cell(Structure):
    """
    Represents a cell. Provides additional functionality to the
    pymatgen.core.structure.Structure class.
    """

    def rotate_to_principal_directions(self):
        """
        Rotates the cell into the principal directions. That is, lattice vector
        a is parallel to the Cartesian x-axis, lattice vector b lies in the
        Cartesian x-y plane and the z-component of lattice vector c is
        positive.

        Note: this method doesn't change the fractional coordinates of the
        sites. However, the Cartesian coordinates may be changed.
        """

        # rotate about the z-axis to align a vertically with the x-axis
        rotation = RotationTransformation(
            [0, 0, 1], 180 - (180/np.pi)*np.arctan2(
                    self.lattice.matrix[0][1],
                    self.lattice.matrix[0][0]))
        new_structure = rotation.apply_transformation(self)
        self.modify_lattice(new_structure.lattice)

        # rotate about the y-axis to make a parallel to the x-axis
        rotation = RotationTransformation(
                [0, 1, 0], (180/np.pi)*np.arctan2(
                    self.lattice.matrix[0][2],
                    self.lattice.matrix[0][0]))
        new_structure = rotation.apply_transformation(self)
        self.modify_lattice(new_structure.lattice)

        # rotate about the x-axis to make b lie in the x-y plane
        rotation = RotationTransformation(
                [1, 0, 0], 180 - (180/np.pi)*np.arctan2(
                    self.lattice.matrix[1][2],
                    self.lattice.matrix[1][1]))
        new_structure = rotation.apply_transformation(self)
        self.modify_lattice(new_structure.lattice)

        # make sure they are all pointing in positive directions
        if self.lattice.matrix[0][0] < 0:
            # rotate about y-axis to make a positive
            rotation = RotationTransformation([0, 1, 0], 180)
            new_structure = rotation.apply_transformation(self)
            self.modify_lattice(new_structure.lattice)
        if self.lattice.matrix[1][1] < 0:
            # rotate about x-axis to make b positive
            rotation = RotationTransformation([1, 0, 0], 180)
            new_structure = rotation.apply_transformation(self)
            self.modify_lattice(new_structure.lattice)
        if self.lattice.matrix[2][2] < 0:
            # mirror c across the x-y plane to make it positive
            # a and b
            a = self.lattice.matrix[0]
            b = self.lattice.matrix[1]
            # the components of c
            cx = self.lattice.matrix[2][0]
            cy = self.lattice.matrix[2][1]
            cz = -1*self.lattice.matrix[2][2]
            self.modify_lattice(Lattice([a, b, [cx, cy, cz]]))

    def rotate_c_parallel_to_z(self):
        """
        Rotates the cell such that lattice vector c is parallel to the
        Cartesian z-axis.

        Note: this method doesn't change the fractional coordinates of the
        sites. However, the Cartesian coordinates may be changed.
        """

        # rotate about the z-axis until c lies in the x-z plane
        rotation = RotationTransformation(
                [0, 0, 1], 180 - (180/np.pi)*np.arctan2(
                    self.lattice.matrix[2][1],
                    self.lattice.matrix[2][0]))
        new_structure = rotation.apply_transformation(self)
        self.modify_lattice(new_structure.lattice)

        # rotate about the y-axis to make c parallel to the z-axis
        rotation = RotationTransformation(
            [0, 1, 0], 180 - (180/np.pi)*np.arctan2(
                self.lattice.matrix[2][0],
                self.lattice.matrix[2][2]))
        new_structure = rotation.apply_transformation(self)
        self.modify_lattice(new_structure.lattice)

        # make sure c is pointing along the positive z-axis
        if self.lattice.matrix[2][2] < 0:
            # rotate 180 degrees about the x-axis
            rotation = RotationTransformation([1, 0, 0], 180)
            new_structure = rotation.apply_transformation(self)
            self.modify_lattice(new_structure.lattice)

    def translate_atoms_into_cell(self):
        """
        Translates all the atoms into the cell, so that their fractional
        coordinates are between 0 and 1.

        Precondition: all the atoms can fit inside the structure's lattice,
            and there is at least 0.001 of (fractional) extra space along
            each lattice vector.
        """

        # get the bounding box of the atoms, in fractional coordinates
        bounding_box = self.get_bounding_box(cart_coords=False)
        translation_vector = []

        # determine the needed shift along each lattice vector
        for i in range(3):
            if bounding_box[i][0] < 0.0:
                translation_vector.append(-1*(bounding_box[i][0]) + 0.001)
            elif bounding_box[i][1] >= 1.0:
                translation_vector.append(-1*(bounding_box[i][1] + 0.001) +
                                          1.0)
            else:
                translation_vector.append(0.0)

        # do the shift
        site_indices = [i for i in range(len(self.sites))]
        self.translate_sites(site_indices, translation_vector,
                             frac_coords=True, to_unit_cell=False)

    def reduce_cell(self):
        """
        Applies Niggli cell reduction to cell.

        Returns a boolean indicated whether cell reduction was successful.

        Note: pymatgen's Niggli cell reduction algorithm sometimes moves the
            atoms' relative positions a little (I've seen up to 0.5 A...)
        """
        # get the reduced structure
        try:  # sometimes pymatgen's reduction routine fails
            reduced_structure = self.get_reduced_structure()
        except:
            return False

        # modify the cell to correspond to the reduced structure
        rcartesian_coords = reduced_structure.cart_coords
        rspecies = reduced_structure.species
        self.modify_lattice(reduced_structure.lattice)
        site_indices = []
        for i in range(len(self.sites)):
            site_indices.append(i)
        self.remove_sites(site_indices)
        for i in range(len(rcartesian_coords)):
            self.append(rspecies[i], rcartesian_coords[i],
                        coords_are_cartesian=True)

        # rotate the cell into the principal directions
        self.rotate_to_principal_directions()
        return True

    def reduce_sheet_cell(self, geometry, constraints):
        """
        Applies Niggli cell reduction to a sheet cell. The idea is to make
        lattice vector c vertical and add lots of vertical vacuum so that the
        standard reduction algorithm only changes the a and b lattice vectors.

        Returns a boolean indicating whether cell reduction was successful.

        Args:
            geometry: the Geometry of the search

            constraints: the Constraints of the search
        """
        # pad the cell with 100 Angstroms of vertical vacuum
        geometry.pad(self, padding=100)

        # reduce the padded cell
        successful_reduction = self.reduce_cell()

        # unpad the reduced cell
        geometry.unpad(self, constraints)

        return successful_reduction

    def get_bounding_box(self, cart_coords=True):
        """
        Returns the smallest and largest coordinates in each dimension of all
        the sites in a cell, as a list of three lists:

            [[min, max], [min, max], [min, max]]

        where the first inner list contains data for the Cartesian x-coordinate
        or lattice vector a, the second inner list contains data for the
        Cartesian y-coordinate or lattice vector b, and the third inner list
        contains data for the Cartesian z-coordinate or lattice vector c.

        Args:
            cart_coords: whether to give the result in Cartesian or fractional
                coordinates
        """

        if cart_coords:
            coords = self.cart_coords
        else:
            coords = self.frac_coords

        # find the largest and smallest coordinates in each dimension
        minx = np.inf
        maxx = -np.inf
        miny = np.inf
        maxy = -np.inf
        minz = np.inf
        maxz = -np.inf
        for coord in coords:
            if coord[0] < minx:
                minx = coord[0]
            if coord[0] > maxx:
                maxx = coord[0]
            if coord[1] < miny:
                miny = coord[1]
            if coord[1] > maxy:
                maxy = coord[1]
            if coord[2] < minz:
                minz = coord[2]
            if coord[2] > maxz:
                maxz = coord[2]
        return [[minx, maxx], [miny, maxy], [minz, maxz]]


class InitialPopulation():
    """
    The initial population of organisms.
    """

    def __init__(self, run_dir_name):
        """
        Makes an InitialPopulation.

        Args:
            run_dir_name: the name (not path) of the garun directory where the
                search will be done
        """

        self.run_dir_name = run_dir_name
        self.initial_population = []

    def add_organism(self, organism_to_add, composition_space):
        """
        Adds a relaxed organism to the initial population.

        Args:
            organism_to_add: the Organism to add to the initial population
        """

        organism_to_add.cell.sort()
        organism_to_add.cell.to('poscar', os.getcwd() + '/POSCAR.' +
                                str(organism_to_add.id))
        print('Adding organism {} to the initial population'.format(
            organism_to_add.id))
        self.initial_population.append(organism_to_add)
        organism_to_add.is_active = True

    def replace_organism(self, old_org, new_org, composition_space):
        """
        Replaces an organism in the initial population with a new organism.

        Precondition: old_org is a current member of the initial population

        Args:
            old_org: the Organism in the initial population to replace

            new_org: the new Organism to replace the old one
        """

        new_org.cell.sort()
        new_org.cell.to('poscar', os.getcwd() + '/POSCAR.' +
                        str(new_org.id))
        print('Replacing organism {} with organism {} in the initial '
              'population'.format(old_org.id, new_org.id))

        # remove the redundant organism and add the new one
        for org in self.initial_population:
            if org.id == old_org.id:
                self.initial_population.remove(org)
                org.is_active = False
        self.initial_population.append(new_org)
        new_org.is_active = True

    def get_progress(self, composition_space):
        """
        Returns either the best energy per atom (for fixed-composition
        search) or the area/volume of the convex hull (for phase diagram
        search).

        Args:
            composition_space: the CompositionSpace of the search
        """

        if composition_space.objective_function == 'epa':
            return self.get_best_epa()
        elif composition_space.objective_function == 'pd':
            return self.get_convex_hull_area(composition_space)

    def get_best_epa(self):
        """
        Returns the epa of the best organism in the initial population.
        """

        sorted_list = copy.deepcopy(self.initial_population)
        sorted_list.sort(key=lambda x: x.epa, reverse=False)
        return sorted_list[0].epa

    def get_convex_hull_area(self, composition_space):
        """
        Returns the area/volume of the current convex hull.

        Args:
            composition_space: the CompositionSpace of the search
        """

        # check if the initial population contains organisms at all the
        # endpoints of the composition space
        if self.has_endpoints(composition_space) and \
                self.has_non_endpoint(composition_space):

            # compute and print the area or volume of the convex hull
            pdentries = []
            for organism in self.initial_population:
                pdentries.append(PDEntry(organism.composition,
                                         organism.total_energy))
            compound_pd = CompoundPhaseDiagram(pdentries,
                                               composition_space.endpoints)

            # get the data for the convex hull
            qhull_data = compound_pd.qhull_data
            # for some reason, the last point is positive, so remove it
            hull_data = np.delete(qhull_data, -1, 0)

            # make a ConvexHull object from the hull data
            # Sometime this fails, saying that only two points were given to
            # construct the convex hull, even though the if statement above
            # checks that there are enough points.
            try:
                convex_hull = ConvexHull(hull_data)
            except:
                return None
            if len(composition_space.endpoints) == 2:
                return convex_hull.area
            else:
                return convex_hull.volume

    def has_endpoints(self, composition_space):
        """
        Checks if there are organisms in the initial population at all of
        endpoints of the composition space.

        Returns a boolean indicating whether or not there are organisms at all
        the endpoints.

        Args:
            composition_space: the CompositionSpace of the search
        """

        for endpoint in composition_space.endpoints:
            has_endpoint = False
            for organism in self.initial_population:
                if endpoint.almost_equals(
                        organism.composition.reduced_composition):
                    has_endpoint = True
            if not has_endpoint:
                return False
        return True

    def has_non_endpoint(self, composition_space):
        """
        Checks that the initial population contains at least one organism not
        at one of the endpoint compositions.

        Returns a boolean indicating whether or not there is an organism not
        at one of the endpoints.

        Args:
            composition_space: the CompositionSpace of the search
        """

        for organism in self.initial_population:
            not_endpoint = True
            for endpoint in composition_space.endpoints:
                if endpoint.almost_equals(
                        organism.composition.reduced_composition):
                    not_endpoint = False
            if not_endpoint:
                return True
        return False


class Pool(object):
    """
    Pool to hold all the organisms that are currently candidates for selection
    to become parents.

    Composed of two parts: a promotion set containing the best few organisms,
    and a queue containing the rest of the organisms in the pool.
    """
    def __init__(self, pool_params, composition_space, run_dir_name):
        """
        Makes a Pool, and sets default parameter values if necessary.

        Args:
            pool_params: the parameters of the pool, as a dictionary

            composition_space: the CompositionSpace of the search

            run_dir_name: the name (not path) of the garun directory
        """

        # the default values
        if len(composition_space.endpoints) == 1:
            # then we're doing a fixed composition search
            self.default_size = 20
        elif len(composition_space.endpoints) == 2:
            # then we're doing a binary phase diagram search search
            self.default_size = 25
        elif len(composition_space.endpoints) == 3:
            # then we're doing a ternary phase diagram search
            self.default_size = 75
        elif len(composition_space.endpoints) == 4:
            # then we're doing a quaternary phase diagram search
            self.default_size = 150

        # this only gets used for fixed composition searches
        self.default_num_promoted = 3

        # if entire Pool block was set to default or left blank
        if pool_params in (None, 'default'):
            self.size = self.default_size
            self.num_promoted = self.default_num_promoted
        else:
            # check if size parameter was used
            if 'size' not in pool_params:
                self.size = self.default_size
            elif pool_params['size'] in (None, 'default'):
                self.size = self.default_size
            else:
                self.size = pool_params['size']

            # check if num_promoted parameter was used
            if 'num_promoted' not in pool_params:
                self.num_promoted = self.default_num_promoted
            elif pool_params['num_promoted'] in (None, 'default'):
                self.num_promoted = self.default_num_promoted
            else:
                self.num_promoted = pool_params['num_promoted']

        # the best few organisms in the pool
        self.promotion_set = []
        # the rest of the organisms in the pool
        self.queue = deque()
        # the parameters for the selection distribution
        self.selection = []
        # the number of organisms added to the pool (excluding the initial
        # population)
        self.num_adds = 0
        # the name (not  path) of the garun directory
        self.run_dir_name = run_dir_name

    def add_initial_population(self, initial_population, composition_space):
        """
        Adds the organisms of the initial population to the pool.

        Args:
            initial_population: the InitialPopulation containing the relaxed
                organisms of the initial population

            composition_space: the CompositionSpace of the search

        Precondition: the pool is empty (no organisms have been previously
            added to it)
        """

        print('Populating the pool with the initial population...')
        organisms_list = initial_population.initial_population

        # populate the promotion set and queue
        if composition_space.objective_function == 'epa':
            for organism in organisms_list:
                organism.value = organism.epa
            organisms_list.sort(key=lambda x: x.value, reverse=False)

            # assign them to either the promotion set or the queue
            for i in range(len(organisms_list)):
                if i < self.num_promoted:
                    self.promotion_set.append(organisms_list[i])
                else:
                    self.queue.appendleft(organisms_list[i])

        elif composition_space.objective_function == 'pd':
            try:
                self.compute_pd_values(organisms_list, composition_space)
            except:
                print('Error: could not construct the phase diagram because '
                      'there is not a structure at one or more of the '
                      'composition space endpoints.')
                print('One of the provided reference structures probably '
                      'failed development, or else there was an error when '
                      'calculating its energy. See previous output.')
                print('Quitting...')
                quit()

            # assign the ones on the convex hull to the promotion set and the
            # rest to the queue
            for organism in organisms_list:
                if organism.value < 0.000000001:
                    self.promotion_set.append(organism)
                else:
                    self.queue.appendleft(organism)

        # compute the fitnesses and selection probabilities of the organisms in
        # the initial population
        self.compute_fitnesses()
        self.compute_selection_probs()

        organisms_list.sort(key=lambda x: x.fitness, reverse=True)
        print('Summary of the initial population: ')
        for organism in organisms_list:
            print('Organism {} has value {} and fitness {} and selection '
                  'probability {} '.format(organism.id, organism.value,
                                           organism.fitness,
                                           organism.selection_prob))

    def add_organism(self, organism_to_add, composition_space):
        """
        Adds a new organism to the pool.

        For fixed-composition searches, if the new organism is better than one
        of the organisms currently in the promotion set, then it is added to
        the promotion set, and the worst organism in the promotion set is moved
        to the back of the queue. Otherwise, the new organism is just appended
        to the back of the queue.

        For phase diagram searches, if the new organism is on the convex hull,
        then it is added to the promotion set and any organisms in the
        promotion set that are no longer on the convex hull are moved the
        queue. Otherwise, the new organism is added to the back of the queue.

        Args:
            organism: the Organism to add to the pool

            composition_space: the CompositionSpace of the search
        """

        print('Adding organism {} to the pool '.format(organism_to_add.id))

        self.num_adds = self.num_adds + 1
        organism_to_add.cell.sort()
        organism_to_add.cell.to('poscar', os.getcwd() + '/POSCAR.' +
                                str(organism_to_add.id))
        organism_to_add.is_active = True

        if composition_space.objective_function == 'epa':
            organism_to_add.value = organism_to_add.epa
            worst_organism = self.get_worst_in_promotion_set()
            if organism_to_add.value < worst_organism.value:
                self.promotion_set.append(organism_to_add)
                self.promotion_set.remove(worst_organism)
                self.queue.appendleft(worst_organism)
            else:
                self.queue.appendleft(organism_to_add)

        elif composition_space.objective_function == 'pd':
            organisms_list = self.to_list()
            organisms_list.append(organism_to_add)
            self.compute_pd_values(organisms_list, composition_space)
            self.check_promotion_set_pd()
            if organism_to_add.value == 0.0:
                self.promotion_set.append(organism_to_add)
            else:
                self.queue.appendleft(organism_to_add)

    def get_worst_in_promotion_set(self):
        """
        Returns the organism in the promotion set with the worst (largest)
        objective function value.
        """

        worst_organism = self.promotion_set[0]
        for organism in self.promotion_set:
            if organism.value > worst_organism.value:
                worst_organism = organism
        return worst_organism

    def check_promotion_set_pd(self):
        """
        For phase diagram searches, checks whether promotion set and queue
        memberships are correct (all organisms with value 0 in promotion set,
        all organisms with value > 0 in queue), and moves organisms from the
        promotion set to the queue (and vice versa) if needed.

        Precondition: a phase diagram search is being done, and all organisms
            in the pool have up-to-date values (set by calling
            compute_pd_values)
        """

        # get organisms in the promotion set that are no longer on the convex
        # hull
        organisms_to_demote = []
        for organism in self.promotion_set:
            if organism.value > 0.000000001:
                organisms_to_demote.append(organism)

        # move all organism that are no longer on the convex hull from the
        # promotion set to the queue
        for organism_to_demote in organisms_to_demote:
            self.promotion_set.remove(organism_to_demote)
            self.queue.appendleft(organism_to_demote)

        # get organisms in the queue that are now on the convex hull
        organisms_to_promote = []
        for organism in self.queue:
            if organism.value < 0.000000001:
                organisms_to_promote.append(organism)

        # move all the organisms on the convex hull to the promotion set
        for organism_to_promote in organisms_to_promote:
            self.queue.remove(organism_to_promote)
            self.promotion_set.append(organism_to_promote)

    def replace_organism(self, old_org, new_org, composition_space):
        """
        Replaces an organism in the pool with a new organism. The new organism
        has the same location in the pool as the old one.

        Precondition: old_org is a member of the current pool.

        Args:
            old_org: the Organism in the pool to replace

            new_org: the new Organism to replace the old one
        """

        print('Replacing organism {} with organism {} in the pool '.format(
            old_org.id, new_org.id))

        new_org.cell.sort()
        new_org.cell.to('poscar', os.getcwd() + '/POSCAR.' + str(new_org.id))

        # set new objective function value
        if composition_space.objective_function == 'epa':
            new_org.value = new_org.epa
        elif composition_space.objective_function == 'pd':
            organisms_list = self.to_list()
            organisms_list.append(new_org)
            self.compute_pd_values(organisms_list, composition_space)

        # add new organism and remove old one
        if old_org in self.promotion_set:
            self.promotion_set.remove(old_org)
            self.promotion_set.append(new_org)
        elif old_org in self.queue:
            queue_list = list(self.queue)
            queue_list.insert(queue_list.index(old_org), new_org)
            queue_list.remove(old_org)
            self.queue = deque(queue_list)
        old_org.is_active = False
        new_org.is_active = True

        # for pd searches, check promotion set and queue memberships
        if composition_space.objective_function == 'pd':
            self.check_promotion_set_pd()

    def compute_pd_values(self, organisms_list, composition_space):
        """
        Constructs a convex hull from the provided organisms and sets the
        organisms' values to their distances from the convex hull.

        Returns the CompoundPhaseDiagram object computed from the organisms
        in organisms_list.

        Args:
            organisms_list: a list of Organisms whose values we need to compute

            composition_space: the CompositionSpace of the search
        """

        # create a PDEntry object for each organism in the list of organisms
        pdentries = {}
        for organism in organisms_list:
            pdentries[organism.id] = PDEntry(organism.composition,
                                             organism.total_energy)

        # put the pdentries in a list
        pdentries_list = []
        for organism_id in pdentries:
            pdentries_list.append(pdentries[organism_id])

        # create a compound phase diagram object from the list of pdentries
        compound_pd = CompoundPhaseDiagram(pdentries_list,
                                           composition_space.endpoints)

        # create a pd analyzer object from the compound phase diagram
        pd_analyzer = PDAnalyzer(compound_pd)

        # transform the pdentries and put them in a dictionary, with the
        # organism id's as the keys
        transformed_pdentries = {}
        for org_id in pdentries:
            transformed_pdentries[org_id] = compound_pd.transform_entries(
                    [pdentries[org_id]], composition_space.endpoints)[0][0]

        # put the values in a dictionary, with the organism id's as the keys
        values = {}
        for org_id in pdentries:
            values[org_id] = pd_analyzer.get_e_above_hull(
                transformed_pdentries[org_id])

        # assign values to the organisms
        for organism in organisms_list:
            organism.value = values[organism.id]

        return compound_pd

    def compute_fitnesses(self):
        """
        Calculates and assigns fitnesses to all organisms in the pool. If
        selection.num_parents is less than pool.size, then the fitnesses are
        computed relative to the selection.num_parents + 1 best organisms in
        the pool.

        Precondition: the organisms in the pool all have up-to-date values
        """

        # get the best organisms that can be parents
        if self.selection.num_parents < self.size:
            best_organisms = self.get_n_best_organisms(
                self.selection.num_parents + 1)
        else:
            best_organisms = self.get_n_best_organisms(self.size)

        # get the best and worst values in this subset of organisms
        best_value = best_organisms[0].value
        worst_value = best_organisms[-1].value

        # compute the fitnesses of the organisms in this subset
        for organism in best_organisms:
            organism.fitness = (organism.value - worst_value)/(best_value -
                                                               worst_value)

        # assign fitnesses of zero to organisms not in this subset
        for organism in self.to_list():
            if organism not in best_organisms:
                organism.fitness = 0.0

    def get_n_best_organisms(self, n):
        """
        Returns a list containing the n best organisms in the pool, sorted in
        order of increasing value.

        Precondition: all the organisms in the pool have up-to-date values

        Args:
            n: the number of best organisms to get
        """

        pool_list = self.to_list()
        pool_list.sort(key=lambda x: x.value, reverse=False)
        return pool_list[:n]

    def compute_selection_probs(self):
        """
        Calculates and assigns selection probabilities to all the organisms in
        the pool.

        Precondition: the organisms in the pool all have up-to-date fitnesses
        """

        # get a list of the organisms in the pool, sorted by increasing fitness
        organisms = self.to_list()
        organisms.sort(key=lambda x: x.fitness, reverse=False)

        # get the denominator of the expression for selection probability
        denom = 0
        for organism in organisms:
            denom = denom + math.pow(organism.fitness, self.selection.power)

        # compute the selection probability and interval location for each
        # organism in best_organisms. The interval location is defined here as
        # the starting point of the interval
        selection_loc = 0
        for organism in organisms:
            organism.selection_prob = math.pow(organism.fitness,
                                               self.selection.power)/denom
            organism.selection_loc = selection_loc
            selection_loc = selection_loc + organism.selection_prob

    def select_organisms(self, n, random):
        """
        Randomly selects n distinct organisms from the pool based on their
        selection probabilities.

        Returns a list containing n organisms.

        Args:
            n: how many organisms to select from the pool

            random: a copy of Python's built in PRNG

        Precondition: all the organisms in the pool have been assigned
            selection probabilities.
        """

        selected_orgs = []
        pool_list = self.to_list()

        # keep going until we've selected enough
        while True:
            rand = random.random()
            for organism in pool_list:
                if rand >= organism.selection_loc and rand < (
                        organism.selection_loc + organism.selection_prob):
                    # check that we haven't already selected this one
                    if organism not in selected_orgs:
                        selected_orgs.append(organism)
                        # check if we've selected enough
                        if len(selected_orgs) == n:
                            return selected_orgs

    def print_summary(self, composition_space):
        """
        Prints out a summary of the organisms in the pool.

        Args:
            composition_space: the CompositionSpace of the search
        """

        pool_list = self.to_list()
        pool_list.sort(key=lambda x: x.value, reverse=False)
        print('Summary of the pool:')
        for organism in pool_list:
            print('Organism {} has value {} and fitness {} and selection '
                  'probability {} '.format(organism.id, organism.value,
                                           organism.fitness,
                                           organism.selection_prob))

    def get_progress(self, composition_space):
        """
        Returns either the best energy per atom (for fixed-composition
        search) or the area/volume of the convex hull (for phase diagram
        search).

        Args:
            composition_space: the CompositionSpace of the search
        """

        if composition_space.objective_function == 'epa':
            return self.get_best_epa()
        elif composition_space.objective_function == 'pd':
            return self.get_convex_hull_area(composition_space)

    def get_best_epa(self):
        """
        Prints out the id and value of the best organism in the pool.
        """

        pool_list = self.to_list()
        pool_list.sort(key=lambda x: x.fitness, reverse=True)
        return pool_list[0].epa

    def get_convex_hull_area(self, composition_space):
        """
        Prints out the area or volume of the convex hull defined by the
        organisms in the promotion set.
        """

        # make a phase diagram of just the organisms in the promotion set
        # (on the lower convex hull)
        pdentries = []
        for organism in self.promotion_set:
            pdentries.append(PDEntry(organism.composition,
                                     organism.total_energy))
        compound_pd = CompoundPhaseDiagram(pdentries,
                                           composition_space.endpoints)

        # get the data for the convex hull
        qhull_data = compound_pd.qhull_data
        # for some reason, the last point is positive, so remove it
        hull_data = np.delete(qhull_data, -1, 0)

        # make a ConvexHull object from the hull data
        try:
            convex_hull = ConvexHull(hull_data)
        except:
            return None
        if len(composition_space.endpoints) == 2:
            return convex_hull.area
        else:
            return convex_hull.volume

    def to_list(self):
        """
        Returns a list containing all the organisms in the pool.
        """

        return self.promotion_set + list(self.queue)


class OffspringGenerator(object):
    """
    This class handles generating offspring organisms from the pool and the
    variations. Basically a place for the make_offspring_organism method to
    live.
    """

    def make_offspring_organism(self, random, pool, variations, geometry,
                                id_generator, whole_pop, developer,
                                redundancy_guard, composition_space,
                                constraints):
        """
        Returns a developed, non-redundant, unrelaxed offspring organism
        generated with one of the variations.

        Args:
            random: a copy of Python's PRNG

            pool: the Pool

            variations: list of all the Variations used in this search

            geometry: the Geometry of the search

            id_generator: the IDGenerator used to assign id numbers to all
                organisms

            whole_pop: list containing copies of the organisms to check for
                redundancy

            developer: the Developer of the search

            redundancy_guard: the RedundancyGuard of the search

            composition_space: the CompositionSpace of the search

            constraints: the Constraints of the search

        Description:

            1. Randomly selects one of the Variations (based on the fraction
                values of the Variations) with which to make an offspring
                Organism.

            2. Tries to generate a developed, non-redundant, unrelaxed
                offspring Organism with the selected Variation.

            3. If no success after 1000 attempts, randomly selects a different
                Variation and tries to make an offspring with that one.

            4. If no success after trying all the Variations, cycles through
                them again.
        """

        # holds the variations that have been tried
        tried_variations = []
        # the maximum number of times to try making a valid offspring with a
        # given variation before giving up and trying a different variation
        max_num_tries = 1000

        while True:
            variation = self.select_variation(random, tried_variations,
                                              variations)
            num_tries = 0
            while num_tries < max_num_tries:
                offspring = variation.do_variation(pool, random, geometry,
                                                   constraints, id_generator)
                if developer.develop(
                        offspring, composition_space, constraints, geometry,
                        pool) and (redundancy_guard.check_redundancy(
                            offspring, whole_pop, geometry) is None):
                    return offspring
                else:
                    num_tries = num_tries + 1
            tried_variations.append(variation)

            # if we've tried all the variations without success, then cycle
            # through them again
            if len(tried_variations) == len(variations):
                tried_variations = []

    def select_variation(self, random, tried_variations, variations):
        """
        Selects a variation that hasn't been tried yet based on their selection
        probabilities.

        Args:
            random: a copy of Python's PRNG

            tried_variations: list of Variations that have already been
                unsuccessfully tried

            variations: list of all the Variations used in this search
        """

        while True:
            rand = random.random()
            lower_bound = 0
            for i in range(len(variations)):
                if rand > lower_bound and rand <= (
                        lower_bound + variations[i].fraction) and \
                        variations[i] not in tried_variations:
                    return variations[i]
                else:
                    lower_bound = lower_bound + variations[i].fraction


class SelectionProbDist(object):
    """
    Defines the probability distribution over the fitnesses of the organisms
    in the pool that determines their selection probabilities.

    Comprised of two numbers: the number of organisms to select from, and the
    selection power.
    """

    def __init__(self, selection_params, pool_size):
        """
        Makes a SelectionProbDist, and sets default parameter values if
        necessary.

        Args:
            selection_params: the parameters defining the distribution, as a
                dictionary

            pool_size: the size of the Pool (how many organisms it contains)
        """

        # default values
        self.default_num_parents = pool_size
        self.default_power = 1

        # if entire Selection block was left blank or set to default
        if selection_params in (None, 'default'):
            self.num_parents = self.default_num_parents
            self.power = self.default_power
        else:
            # check the num_parents parameter
            if 'num_parents' not in selection_params:
                self.num_parents = self.default_num_parents
            elif selection_params['num_parents'] in (None, 'default'):
                self.num_parents = self.default_num_parents
            # check that the given number isn't larger than the pool size
            elif selection_params['num_parents'] > pool_size:
                self.num_parents = self.default_num_parents
            else:
                self.num_parents = selection_params['num_parents']

            # check the selection_power parameter
            if 'power' not in selection_params:
                self.power = self.default_power
            elif selection_params['power'] in (None, 'default'):
                self.power = self.default_power
            else:
                self.power = selection_params['power']


class CompositionSpace(object):
    """
    Represents the composition space to be searched by the algorithm.
    """

    def __init__(self, endpoints):
        """
        Makes a CompositionSpace, which is list of
        pymatgen.core.composition.Composition objects.

        Args:
            endpoints: the list of compositions, as strings (e.g. ["Al2O3"])
        """

        for i in range(len(endpoints)):
            endpoints[i] = Composition(endpoints[i])
            endpoints[i] = endpoints[i].reduced_composition

        self.endpoints = endpoints

        # objective function lives here
        self.objective_function = self.infer_objective_function()

    def infer_objective_function(self):
        """
        Infers the objective function (energy per atom or phase diagram) based
        on the composition space.

        Returns either "epa" or "pd".
        """

        if len(self.endpoints) == 1:
            return "epa"
        else:
            for point in self.endpoints:
                for next_point in self.endpoints:
                    if not point.almost_equals(next_point, 0.0, 0.0):
                        return "pd"
        return "epa"

    def get_all_elements(self):
        """
        Returns a list of all the elements
        (as pymatgen.core.periodic_table.Element objects) that are in the
        composition space.
        """

        # get each element type from the composition_space object
        elements = []
        for point in self.endpoints:
            elements = elements + point.elements

        # remove duplicates
        elements = list(set(elements))
        return elements

    def get_all_pairs(self):
        """
        Returns all possible pairs of elements in the composition space, as
        list of strings, where each string contains the symbols of two
        elements, separated by a space. Does not include self-pairs
        (e.g., "Cu Cu").
        """

        elements = self.get_all_elements()
        if len(elements) == 1:
            return []
        pairs = []
        for i in range(0, len(elements) - 1):
            for j in range(i + 1, len(elements)):
                pairs.append(str(elements[i].symbol + " " +
                                 elements[j].symbol))
        return pairs

    def get_all_swappable_pairs(self):
        """
        Computes all pairs of elements in the the composition space that are
        allowed to be swapped based on their electronegativities. Only pairs
        whose electronegativities differ by 1.1 or less can be swapped.

        Returns a list of strings, where each string contains the symbols of
        two elements, separated by a space. Does not include self-pairs
        (e.g., "Cu Cu").
        """

        all_pairs = self.get_all_pairs()
        allowed_pairs = []
        for pair in all_pairs:
            element1 = Element(pair.split()[0])
            element2 = Element(pair.split()[1])
            diff = abs(element1.X - element2.X)
            if diff <= 1.1:
                allowed_pairs.append(pair)
        return allowed_pairs


class StoppingCriteria(object):
    """
    Defines when the search should be stopped.
    """

    def __init__(self, stopping_parameters, composition_space):
        """
        Makes a StoppingCriteria, and sets default parameter values if
        necessary.

        Args:
            stopping_parameters: a dictionary of parameters

            composition_space: the CompositionSpace of the search
        """

        # set the defaults
        if len(composition_space.endpoints) == 1:
            # for fixed composition searches
            self.default_num_energy_calcs = 800
        elif len(composition_space.endpoints) == 2:
            # for binary phase diagram searches
            self.default_num_energy_calcs = 1000
        elif len(composition_space.endpoints) == 3:
            # for ternary phase diagram searches
            self.default_num_energy_calcs = 3000
        elif len(composition_space.endpoints) == 4:
            # for quaternary phase diagram searches
            self.default_num_energy_calcs = 6000

        # the energy per atom at which to stop
        self.default_epa_achieved = None
        # the Cell at which to stop when found
        self.default_found_cell = None
        # whether or not the stopping criteria are satisfied
        self.are_satisfied = False
        # to keep track of how many energy calculations have been done
        self.calc_counter = 0

        # set defaults if stopping_parameters equals 'default' or None
        if stopping_parameters in (None, 'default'):
            self.num_energy_calcs = self.default_num_energy_calcs
            self.epa_achieved = self.default_epa_achieved
            self.found_cell = self.default_found_cell
        # check each keyword to see if it's been included
        else:
            # value achieved
            if 'epa_achieved' in stopping_parameters:
                if stopping_parameters['epa_achieved'] in (None, 'default'):
                    self.epa_achieved = self.default_epa_achieved
                elif composition_space.objective_function == 'epa':
                    self.epa_achieved = stopping_parameters['epa_achieved']
            else:
                self.epa_achieved = self.default_epa_achieved

            # found structure
            if 'found_structure' in stopping_parameters:
                if stopping_parameters['found_structure'] in (None, 'default'):
                    self.found_cell = self.default_found_cell
                else:
                    self.path_to_structure_file = stopping_parameters[
                        'found_structure']
                    self.found_cell = Cell.from_file(
                            stopping_parameters['found_structure'])
            else:
                self.found_cell = self.default_found_cell

            # number energy calculations
            if 'num_energy_calcs' in stopping_parameters:
                if stopping_parameters['num_energy_calcs'] in (None,
                                                               'default'):
                    if self.epa_achieved is None and \
                            self.found_cell is None:
                        self.num_energy_calcs = self.default_num_energy_calcs
                    else:
                        self.num_energy_calcs = None
                else:
                    self.num_energy_calcs = stopping_parameters[
                        'num_energy_calcs']
            elif self.epa_achieved is None and self.found_cell is None:
                self.num_energy_calcs = self.default_num_energy_calcs
            else:
                self.num_energy_calcs = None

    def update_calc_counter(self):
        """
        If num_energy_calcs stopping criteria is being used, increments calc
        counter and updates are_satisfied if necessary.
        """

        if self.num_energy_calcs is not None:
            self.calc_counter += 1
            if self.calc_counter >= self.num_energy_calcs:
                self.are_satisfied = True

    def check_organism(self, organism, redundancy_guard, geometry):
        """
        If value_achieved or found_structure stopping criteria are used, checks
        if the relaxed organism satisfies them, and updates are_satisfied.

        Args:
            organism: a relaxed Organism whose value has been computed

            redundancy_guard: the RedundancyGuard of the search

            geometry: the Geometry of the search
        """

        # check the objective function value if needed
        if self.epa_achieved is not None:
            if organism.epa <= self.epa_achieved:
                self.are_satisfied = True

        # check the structure if needed
        if self.found_cell is not None:
            if geometry.shape == 'cluster':
                mol1 = Molecule(organism.cell.species,
                                organism.cell.cart_coords)
                mol2 = Molecule(self.found_cell.species,
                                self.found_cell.cart_coords)
                self.are_satisfied = \
                    redundancy_guard.molecule_matcher.fit(mol1, mol2)
            else:
                self.are_satisfied = \
                    redundancy_guard.structure_matcher.fit(organism.cell,
                                                           self.found_cell)


class DataWriter(object):
    """
    For writing useful data to a file in the course of a search.
    """

    def __init__(self, file_path, composition_space):
        """
        Makes a DataWriter.

        Args:
            file_path: the path to the file where the data is to be written

            composition_space: the CompositionSpace of the search
        """

        self.file_path = file_path
        with open(self.file_path, 'a') as data_file:
            data_file.write('Composition space endpoints: ')
            for endpoint in composition_space.endpoints:
                data_file.write(' {}'.format(
                    endpoint.reduced_formula.replace(' ', '')))
            data_file.write('\n\n')
            data_file.write('id\t\t composition\t total energy\t\t epa\t\t\t '
                            'num calcs\t best value\n\n')

    def write_data(self, organism, num_calcs, progress):
        """
        Writes data to self.file_path in the format:

            id    composition    total energy    epa    num calcs    best value

        Args:
            organism: the Organism whose data to write

            num_calcs: the number of calculations finished so far,
                including failures

            progress: for fixed-composition searches, the best epa seen by the
                algorithm so far. For phase diagram searches, the area/volume
                of the convex hull, or None if no convex hull could be
                constructed.
        """

        format_string = '{0}\t\t {1}\t\t {2:.6f}\t\t {3:.6f}\t\t {4}\t\t'
        if progress is None:
            format_string = format_string + ' None\n'
        else:
            format_string = format_string + ' {5:.6f}\n'

        with open(self.file_path, 'a') as data_file:
            data_file.write(format_string.format(
                organism.id, organism.composition.formula.replace(' ', ''),
                organism.total_energy, organism.epa, num_calcs, progress))
