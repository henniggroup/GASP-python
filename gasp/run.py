# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Run module:

This module is run to do a genetic algorithm structure search.

"""

import general
import objects_maker
import parameters_printer

import copy
import threading
import random
import yaml
import sys
import os


# get dictionaries from the input file (in yaml format)
input_file = os.path.abspath(sys.argv[1])
try:
    with open(input_file, 'r') as f:
        parameters = yaml.load(f)
except:
    print('Error reading input file.')
    print('Quitting...')

# make the objects needed by the algorithm
objects_dict = objects_maker.make_objects(parameters)

# get the objects from the dictionary for convenience
run_dir_name = objects_dict['run_dir_name']
organism_creators = objects_dict['organism_creators']
num_calcs_at_once = objects_dict['num_calcs_at_once']
composition_space = objects_dict['composition_space']
constraints = objects_dict['constraints']
geometry = objects_dict['geometry']
developer = objects_dict['developer']
redundancy_guard = objects_dict['redundancy_guard']
stopping_criteria = objects_dict['stopping_criteria']
energy_calculator = objects_dict['energy_calculator']
pool = objects_dict['pool']
variations = objects_dict['variations']
id_generator = objects_dict['id_generator']

# make the run directory and move into it
garun_dir = str(os.getcwd()) + '/' + run_dir_name
os.mkdir(garun_dir)
os.chdir(garun_dir)

# make the temp subdirectory where the energy calculations will be done
os.mkdir(garun_dir + '/temp')

# print the search parameters to a file in the run directory
parameters_printer.print_parameters(objects_dict)

whole_pop = []
num_finished_calcs = 0
threads = []
initial_population = general.InitialPopulation(run_dir_name)

# To temporarily hold relaxed organisms. The key to each relaxed organism is
# the index of the Thread in the list threads that did the energy calculation.
relaxed_organisms = {}

# populate the initial population
for creator in organism_creators:
    print('Making {} organisms with {}'.format(creator.number, creator.name))
    while not creator.is_finished:

        # start initial batch of energy calculations
        if len(threads) < num_calcs_at_once:
            # make a new organism - keep trying until we get one
            new_organism = creator.create_organism(
                id_generator, composition_space, constraints, random)
            while new_organism is None and not creator.is_finished:
                new_organism = creator.create_organism(
                    id_generator, composition_space, constraints, random)
            if new_organism is not None:  # loop above could return None
                geometry.unpad(new_organism, constraints)
                if developer.develop(new_organism, composition_space,
                                     constraints, geometry, pool):
                    redundant_organism = redundancy_guard.check_redundancy(
                        new_organism, whole_pop)
                    if redundant_organism is None:  # no redundancy
                        # add a copy to whole_pop so the organisms in whole_pop
                        # don't change upon relaxation
                        whole_pop.append(copy.deepcopy(new_organism))
                        stopping_criteria.update_calc_counter()
                        geometry.pad(new_organism)
                        thread_index = len(threads)
                        thread = threading.Thread(
                            target=energy_calculator.do_energy_calculation,
                            args=[new_organism, relaxed_organisms,
                                  thread_index])
                        thread.start()
                        threads.append(thread)

        # process finished calculations and start new ones
        else:
            for thread in threads:
                if not thread.is_alive():
                    num_finished_calcs += 1
                    thread_index = threads.index(thread)
                    relaxed_organism = relaxed_organisms[thread_index]
                    relaxed_organisms[thread_index] = None

                    # take care of relaxed organism
                    if relaxed_organism is not None:
                        geometry.unpad(relaxed_organism, constraints)
                        if developer.develop(relaxed_organism,
                                             composition_space, constraints,
                                             geometry, pool):
                            redundant_organism = \
                                redundancy_guard.check_redundancy(
                                    relaxed_organism, whole_pop)
                            if redundant_organism is not None:  # redundancy
                                if redundant_organism.is_active and \
                                        redundant_organism.epa > \
                                        relaxed_organism.epa:
                                    initial_population.replace_organism(
                                        redundant_organism, relaxed_organism,
                                        composition_space)
                                    initial_population.print_progress(
                                        composition_space)
                                    print('Number of energy calculations so '
                                          'far: {} '.format(
                                              num_finished_calcs))
                            else:  # no redundancy
                                stopping_criteria.check_organism(
                                    relaxed_organism)
                                initial_population.add_organism(
                                    relaxed_organism, composition_space)
                                whole_pop.append(relaxed_organism)
                                initial_population.print_progress(
                                    composition_space)
                                print('Number of energy calculations so far: '
                                      '{} '.format(num_finished_calcs))
                                if creator.is_successes_based and \
                                        relaxed_organism.made_by == \
                                        creator.name:
                                    creator.update_status()

                    # make another organism for the initial population
                    started_new_calc = False
                    while not started_new_calc and not creator.is_finished:
                        new_organism = creator.create_organism(
                            id_generator, composition_space,
                            constraints, random)
                        while new_organism is None and not creator.is_finished:
                            new_organism = creator.create_organism(
                                id_generator, composition_space, constraints,
                                random)
                        if new_organism is not None:
                            geometry.unpad(new_organism, constraints)
                            if developer.develop(new_organism,
                                                 composition_space,
                                                 constraints, geometry, pool):
                                redundant_organism = \
                                    redundancy_guard.check_redundancy(
                                        new_organism, whole_pop)
                                if redundant_organism is None:  # no redundancy
                                    whole_pop.append(
                                        copy.deepcopy(new_organism))
                                    stopping_criteria.update_calc_counter()
                                    geometry.pad(new_organism)
                                    thread = threading.Thread(
                                        target=energy_calculator.do_energy_calculation,
                                        args=(new_organism, relaxed_organisms,
                                              thread_index))
                                    thread.start()
                                    threads[thread_index] = thread
                                    started_new_calc = True

# process all the calculations that were still running when the last creator
# finished
num_to_get = num_calcs_at_once  # how many threads we have left to handle
handled_indices = []  # the indices of the threads we've already handled
while num_to_get > 0:
    for thread in threads:
        thread_index = threads.index(thread)
        if not thread.is_alive() and thread_index not in handled_indices:
            num_finished_calcs += 1
            relaxed_organism = relaxed_organisms[thread_index]
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            relaxed_organisms[thread_index] = None

            # take care of relaxed organism
            if relaxed_organism is not None:
                geometry.unpad(relaxed_organism, constraints)
                if developer.develop(relaxed_organism, composition_space,
                                     constraints, geometry, pool):
                    redundant_organism = redundancy_guard.check_redundancy(
                        relaxed_organism, whole_pop)
                    if redundant_organism is not None:  # redundancy
                        if redundant_organism.is_active and \
                                redundant_organism.epa > relaxed_organism.epa:
                            initial_population.replace_organism(
                                redundant_organism, relaxed_organism,
                                composition_space)
                            initial_population.print_progress(
                                composition_space)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))
                    else:  # no redundancy
                        stopping_criteria.check_organism(relaxed_organism)
                        initial_population.add_organism(relaxed_organism,
                                                        composition_space)
                        whole_pop.append(relaxed_organism)
                        initial_population.print_progress(composition_space)
                        print('Number of energy calculations so far: '
                              '{} '.format(num_finished_calcs))

# populate the pool with the initial population
pool.add_initial_population(initial_population, composition_space)

# To temporarily hold relaxed organisms. The key to each relaxed organism is
# the index of the Thread in the list threads that did the energy calculation.
relaxed_organisms = {}

offspring_generator = general.OffspringGenerator()
threads = []

# create the initial batch of offspring organisms and submit them for energy
# calculations
for _ in range(num_calcs_at_once):
    unrelaxed_offspring = offspring_generator.make_offspring_organism(
        random, pool, variations, geometry, id_generator, whole_pop, developer,
        redundancy_guard, composition_space, constraints)
    whole_pop.append(copy.deepcopy(unrelaxed_offspring))
    geometry.pad(unrelaxed_offspring)
    stopping_criteria.update_calc_counter()
    thread_index = len(threads)
    thread = threading.Thread(
        target=energy_calculator.do_energy_calculation,
        args=[unrelaxed_offspring, relaxed_organisms, thread_index])
    thread.start()
    threads.append(thread)

# process finished calculations and start new ones
while not stopping_criteria.are_satisfied:
    for thread in threads:
        if not thread.is_alive():
            num_finished_calcs += 1
            thread_index = threads.index(thread)
            relaxed_offspring = relaxed_organisms[thread_index]
            relaxed_organisms[thread_index] = None

            # take care of relaxed offspring organism
            if relaxed_offspring is not None:
                geometry.unpad(relaxed_offspring, constraints)
                if developer.develop(relaxed_offspring, composition_space,
                                     constraints, geometry, pool):
                    # check for redundancy with the organisms in the pool first
                    redundant_organism = redundancy_guard.check_redundancy(
                        relaxed_offspring, pool.to_list())
                    if redundant_organism is not None:  # redundancy
                        if redundant_organism.epa > relaxed_offspring.epa:
                            pool.replace_organism(redundant_organism,
                                                  relaxed_offspring,
                                                  composition_space)
                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            pool.print_summary(composition_space)
                            pool.print_progress(composition_space)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))
                    # check for redundancy with all the organisms in whole_pop
                    else:
                        redundant_organism = redundancy_guard.check_redundancy(
                            relaxed_offspring, whole_pop)
                    if redundant_organism is None:  # no redundancy
                        stopping_criteria.check_organism(relaxed_offspring)
                        pool.add_organism(relaxed_offspring, composition_space)
                        whole_pop.append(relaxed_offspring)

                        # check if we've added enough new offspring organisms
                        # to the pool that we can remove the initial population
                        # organisms from the front (right end) of the queue.
                        if pool.num_adds == pool.size:
                            print('Removing the initial population from the '
                                  'pool ')
                            for _ in range(len(
                                    initial_population.initial_population)):
                                removed_org = pool.queue.pop()
                                removed_org.is_active = False
                                print('Removing organism {} from the '
                                      'pool '.format(removed_org.id))

                        # if the initial population organisms have already been
                        # removed from the pool's queue, then just need to pop
                        # one organism from the front (right end) of the queue.
                        elif pool.num_adds > pool.size:
                            removed_org = pool.queue.pop()
                            removed_org.is_active = False
                            print('Removing organism {} from the pool '.format(
                                removed_org.id))

                        pool.compute_fitnesses()
                        pool.compute_selection_probs()
                        pool.print_summary(composition_space)
                        pool.print_progress(composition_space)
                        print('Number of energy calculations so far: '
                              '{} '.format(num_finished_calcs))

            # make another offspring organism
            if not stopping_criteria.are_satisfied:
                unrelaxed_offspring = \
                    offspring_generator.make_offspring_organism(
                        random, pool, variations, geometry, id_generator,
                        whole_pop, developer, redundancy_guard,
                        composition_space, constraints)
                whole_pop.append(copy.deepcopy(unrelaxed_offspring))
                geometry.pad(unrelaxed_offspring)
                stopping_criteria.update_calc_counter()
                thread = threading.Thread(
                    target=energy_calculator.do_energy_calculation,
                    args=[unrelaxed_offspring, relaxed_organisms,
                          thread_index])
                thread.start()
                threads[thread_index] = thread

# process all the calculations that were still running when the last creator
# finished
num_to_get = num_calcs_at_once  # how many threads we have left to handle
handled_indices = []  # the indices of the threads we've already handled
while num_to_get > 0:
    for thread in threads:
        thread_index = threads.index(thread)
        if not thread.is_alive() and thread_index not in handled_indices:
            num_finished_calcs += 1
            relaxed_offspring = relaxed_organisms[thread_index]
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            relaxed_organisms[thread_index] = None

            # take care of relaxed offspring organism
            if relaxed_offspring is not None:
                geometry.unpad(relaxed_offspring, constraints)
                if developer.develop(relaxed_offspring, composition_space,
                                     constraints, geometry, pool):
                    # check for redundancy with the organisms in the pool first
                    redundant_organism = redundancy_guard.check_redundancy(
                        relaxed_offspring, pool.to_list())
                    if redundant_organism is not None:  # redundancy
                        if redundant_organism.epa > relaxed_offspring.epa:
                            pool.replace_organism(redundant_organism,
                                                  relaxed_offspring,
                                                  composition_space)
                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            pool.print_summary(composition_space)
                            pool.print_progress(composition_space)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))
                    # check for redundancy with all the organisms in whole_pop
                    else:
                        redundant_organism = redundancy_guard.check_redundancy(
                            relaxed_offspring, whole_pop)
                    if redundant_organism is None:  # no redundancy
                        pool.add_organism(relaxed_offspring, composition_space)
                        whole_pop.append(relaxed_offspring)
                        removed_org = pool.queue.pop()
                        removed_org.is_active = False
                        print('Removing organism {} from the pool '.format(
                            removed_org.id))
                        pool.compute_fitnesses()
                        pool.compute_selection_probs()
                        pool.print_summary(composition_space)
                        pool.print_progress(composition_space)
                        print('Number of energy calculations so far: '
                              '{} '.format(num_finished_calcs))
