# coding: utf-8
# Copyright(c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Population module:

1. InitialPopulation: represents the initial population of organisms,
        produced from one or more of the organism creators

2. Pool: represents the population of organisms, after the initial population

"""

from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram

import os
import copy
import math
import numpy as np
from collections import deque
from scipy.spatial.qhull import ConvexHull


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
        # the parameters for the selection distribution and composition fitness
        # weight
        # note: these are just placeholders - they get set in
        # objects_maker.make_objects()
        self.selection = []
        self.comp_fitness_weight = []
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

        # check that the initial population contains at least three organisms
        if len(organisms_list) < 3:
            print('The initial population must contain at least three '
                  'organisms.')
            print('Quitting...')
            quit()

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

        # transform the pdentries and put them in a dictionary, with the
        # organism id's as the keys
        transformed_pdentries = {}
        for org_id in pdentries:
            transformed_pdentries[org_id] = compound_pd.transform_entries(
                    [pdentries[org_id]], composition_space.endpoints)[0][0]

        # put the values in a dictionary, with the organism id's as the keys
        values = {}
        for org_id in pdentries:
            values[org_id] = compound_pd.get_e_above_hull(
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
        if best_value == worst_value:  # in case they're all on the convex hull
            for organism in best_organisms:
                organism.fitness = 1.0
        else:
            for organism in best_organisms:
                organism.fitness = (organism.value - worst_value)/(
                    best_value - worst_value)

        # assign fitnesses of zero to organisms not in this subset
        for organism in self.to_list():
            if organism not in best_organisms:
                organism.fitness = 0.0

    def compute_relative_fitnesses(self, ref_organism, composition_space):
        """
        Calculates and assigns relative fitnesses to all the organisms in the
        pool.

        For fixed-composition searches, the relative fitness is equivalent to
        the regular fitness, except the relative fitness of ref_organism is
        set to zero.

        For phase diagram searches, the relative fitness is taken to be the
        average of the regular fitness and the composition fitness, where the
        composition fitness is defined as 1 minus the normalized distance
        between an organism and the ref_organism in composition space. An
        exception to this occurs when the ref_organism and the organism are
        both at the same endpoint of the composition space; in this case, the
        composition fitness of the organism is set to zero. The relative
        fitness of ref_organism is always set to zero.

        Args:
            ref_organism: relative fitnesses are computed w.r.t. this Organism

            composition_space: the CompositionSpace of the search
        """

        pool_list = self.to_list()

        if composition_space.objective_function == 'epa':
            for organism in pool_list:
                organism.relative_fitness = organism.fitness
            ref_organism.relative_fitness = 0.0
        elif composition_space.objective_function == 'pd':

            # compute the weight for the composition fitness based on the
            # distance between the ref organism's composition and the center of
            # the composition space (closer to center -> smaller weight)
            dist_from_center = self.get_composition_distance(
                ref_organism.composition_vector, composition_space.center)
            # normalize distance from center to lie between 0 and 1
            normalized_dist_from_center = (
                dist_from_center/composition_space.max_dist_from_center)
            # the weight for the composition fitness
            comp_fit_weight = self.comp_fitness_weight.max_weight*math.pow(
                normalized_dist_from_center, self.comp_fitness_weight.power)

            # compute the relative fitnesses from the weighted average of the
            # regular and composition fitnesses
            for organism in pool_list:
                if organism.fitness > 0.0:
                    comp_fitness = 1 - self.get_composition_distance(
                        organism.composition_vector,
                        ref_organism.composition_vector)
                    # in case both the organisms are at the same endpoint
                    if ref_organism.is_at_endpoint(composition_space) and \
                            ref_organism.composition.reduced_composition.almost_equals(
                                organism.composition.reduced_composition):
                        organism.relative_fitness = (1 - comp_fit_weight
                                                     )*organism.fitness
                    else:
                        organism.relative_fitness = \
                            comp_fit_weight*comp_fitness + (
                                1 - comp_fit_weight)*organism.fitness
                else:
                    organism.relative_fitness = 0.0
            # set the relative fitness of the reference organism to 0
            ref_organism.relative_fitness = 0.0

    def get_composition_distance(self, comp_vect1, comp_vect2):
        """
        Returns the normalized distance between two composition vectors, which
        is defined as the L1 norm of their difference, divided by 2 to
        normalize (so the maximum distance is 1 and the minimum distance is 0).

        Args:
            comp_vect1: the first composition vector, as a numpy array

            comp_vect2: the second composition vector, as a numpy array
        """

        # compute the different between the two composition vectors
        diff_vector = np.subtract(comp_vect1, comp_vect2)
        # compute the L1 norm of the difference vector
        return 0.5*np.linalg.norm(diff_vector, ord=1)

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

    def compute_relative_selection_probs(self):
        """
        Calculates and assigns relative selection probabilities to all the
        organisms in the pool.

        Precondition: the organisms in the pool all have up-to-date relative
            fitnesses
        """

        # get a list of the organisms in the pool, sorted by increasing
        # relative fitness
        organisms = self.to_list()
        organisms.sort(key=lambda x: x.relative_fitness, reverse=False)

        # get the denominator of the expression for selection probability
        denom = 0
        for organism in organisms:
            denom = denom + math.pow(organism.relative_fitness,
                                     self.selection.power)

        # compute the relative selection probability and interval location for
        # each organism in best_organisms. The interval location is defined
        # here as the starting point of the interval
        relative_selection_loc = 0
        for organism in organisms:
            organism.relative_selection_prob = math.pow(
                organism.relative_fitness, self.selection.power)/denom
            organism.relative_selection_loc = relative_selection_loc
            relative_selection_loc = relative_selection_loc + \
                organism.relative_selection_prob

    def select_organism(self, random, composition_space, excluded_org=None):
        """
        Randomly selects an organism from the pool based on selection
        probabilities (either standard or relative). Returns the selected
        organism.

        Args:
            random: a copy of Python's built in PRNG

            composition_space: the CompositionSpace of the search

            excluded_org: an Organism to exclude from being selected
        """

        pool_list = self.to_list()

        # if not excluding an organism, then select based on standard selection
        # probabilities
        if excluded_org is None:
            rand = random.random()
            for organism in pool_list:
                if rand >= organism.selection_loc and rand < (
                        organism.selection_loc + organism.selection_prob):
                    return organism
        # if excluding an organism, first compute relative selection
        # probabilities, then select based on those
        else:
            self.compute_relative_fitnesses(excluded_org, composition_space)
            self.compute_relative_selection_probs()
            rand = random.random()
            for organism in pool_list:
                if rand >= organism.relative_selection_loc and rand < (
                        organism.relative_selection_loc +
                        organism.relative_selection_prob):
                    return organism

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
