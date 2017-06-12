# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Run module:

This module is run to do a genetic algorithm structure search.

Usage: python run.py /path/to/gasp/input/file

"""

from gasp import general
from gasp import objects_maker
from gasp import parameters_printer

import copy
import threading
import random
import yaml
import sys
import os


def main():
    # get dictionaries from the input file (in yaml format)
    if len(sys.argv) < 2:
        print('No input file given.')
        print('Quitting...')
        quit()
    else:
        input_file = os.path.abspath(sys.argv[1])

    try:
        with open(input_file, 'r') as f:
            parameters = yaml.load(f)
    except:
        print('Error reading input file.')
        print('Quitting...')
        quit()

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
    surrogate_model = objects_dict['surrogate_model']
    sm_output_interpreter = objects_dict['sm_output_interpreter']

    # make the run directory and move into it
    garun_dir = str(os.getcwd()) + '/' + run_dir_name
    os.mkdir(garun_dir)
    os.chdir(garun_dir)

    # make the temp subdirectory where the energy calculations will be done
    os.mkdir(garun_dir + '/temp')

    # print the search parameters to a file in the run directory
    parameters_printer.print_parameters(objects_dict)

    # make the data writer
    data_writer = general.DataWriter(garun_dir + '/run_data',
                                     composition_space)

    whole_pop = []
    num_finished_calcs = 0
    threads = []
    initial_population = general.InitialPopulation(run_dir_name)

    # To temporarily hold relaxed organisms. The key to each relaxed organism
    # is the index of the Thread in the list threads that did the energy
    # calculation.
    relaxed_organisms = {}

    # populate the initial population
    for creator in organism_creators:
        print('Making {} organisms with {}'.format(creator.number,
                                                   creator.name))
        while not creator.is_finished and not stopping_criteria.are_satisfied:

            # start initial batch of energy calculations
            if len(threads) < num_calcs_at_once:
                # make a new organism - keep trying until we get one
                new_organism = creator.create_organism(
                    id_generator, composition_space, constraints, random)
                while new_organism is None and not creator.is_finished:
                    new_organism = creator.create_organism(
                        id_generator, composition_space, constraints, random)
                if new_organism is not None:  # loop above could return None
                    geometry.unpad(new_organism.cell, constraints)
                    if developer.develop(new_organism, composition_space,
                                         constraints, geometry, pool):
                        redundant_organism = redundancy_guard.check_redundancy(
                            new_organism, whole_pop, geometry)
                        if redundant_organism is None:  # no redundancy
                            # add a copy to whole_pop so the organisms in
                            # whole_pop don't change upon relaxation
                            whole_pop.append(copy.deepcopy(new_organism))
                            geometry.pad(new_organism.cell)
                            # check the surrogate model's prediction
                            predicted_epa = surrogate_model.predict_epa(
                                new_organism)
                            if sm_output_interpreter.is_below_threshold(
                                    new_organism, predicted_epa,
                                    composition_space,
                                    initial_population.initial_population):
                                predicted_cell = \
                                    surrogate_model.predict_relaxed_cell(
                                        new_organism)
                                if predicted_cell is not None:
                                    new_organism.cell = predicted_cell
                                # do the energy calculation
                                stopping_criteria.update_calc_counter()
                                index = len(threads)
                                thread = threading.Thread(
                                    target=energy_calculator.do_energy_calculation,
                                    args=[new_organism, relaxed_organisms,
                                          index, composition_space])
                                thread.start()
                                threads.append(thread)

            # process finished calculations and start new ones
            else:
                for index, thread in enumerate(threads):
                    if not thread.is_alive():
                        num_finished_calcs += 1
                        relaxed_organism = relaxed_organisms[index]
                        relaxed_organisms[index] = None

                        # take care of relaxed organism
                        if relaxed_organism is not None:
                            geometry.unpad(relaxed_organism.cell, constraints)
                            if developer.develop(relaxed_organism,
                                                 composition_space,
                                                 constraints, geometry, pool):
                                redundant_organism = \
                                    redundancy_guard.check_redundancy(
                                        relaxed_organism, whole_pop, geometry)
                                if redundant_organism is not None:  # redundant
                                    if redundant_organism.is_active and \
                                            redundant_organism.epa > \
                                            relaxed_organism.epa:
                                        initial_population.replace_organism(
                                            redundant_organism,
                                            relaxed_organism,
                                            composition_space)
                                        progress = \
                                            initial_population.get_progress(
                                                composition_space)
                                        data_writer.write_data(
                                            relaxed_organism,
                                            num_finished_calcs, progress)
                                        print('Number of energy calculations '
                                              'so far: {} '.format(
                                                  num_finished_calcs))
                                else:  # not redundant
                                    stopping_criteria.check_organism(
                                        relaxed_organism, redundancy_guard,
                                        geometry)
                                    initial_population.add_organism(
                                        relaxed_organism, composition_space)
                                    whole_pop.append(relaxed_organism)
                                    progress = \
                                        initial_population.get_progress(
                                            composition_space)
                                    data_writer.write_data(
                                        relaxed_organism, num_finished_calcs,
                                        progress)
                                    print('Number of energy calculations so '
                                          'far: {} '.format(
                                              num_finished_calcs))
                                    # give a copy of the padded organism to the
                                    # surrogate model
                                    org_copy = copy.deepcopy(relaxed_organism)
                                    geometry.pad(org_copy.cell)
                                    surrogate_model.add_relaxed_organism(
                                        org_copy)
                                    # update creator status if needed
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
                            while new_organism is None and not \
                                    creator.is_finished:
                                new_organism = creator.create_organism(
                                    id_generator, composition_space,
                                    constraints, random)
                            if new_organism is not None:
                                geometry.unpad(new_organism.cell, constraints)
                                if developer.develop(new_organism,
                                                     composition_space,
                                                     constraints, geometry,
                                                     pool):
                                    redundant_organism = \
                                        redundancy_guard.check_redundancy(
                                            new_organism, whole_pop, geometry)
                                    if redundant_organism is None:
                                        # not redundant
                                        whole_pop.append(
                                            copy.deepcopy(new_organism))
                                        geometry.pad(new_organism.cell)
                                        # check the surrogate model prediction
                                        predicted_epa = surrogate_model.predict_epa(
                                            new_organism)
                                        if sm_output_interpreter.is_below_threshold(
                                                new_organism, predicted_epa,
                                                composition_space,
                                                initial_population.initial_population):
                                            predicted_cell = surrogate_model.predict_relaxed_cell(
                                                new_organism)
                                            if predicted_cell is not None:
                                                new_organism.cell = predicted_cell
                                            # do the energy calculation
                                            stopping_criteria.update_calc_counter()
                                            new_thread = threading.Thread(
                                                target=energy_calculator.do_energy_calculation,
                                                args=[new_organism,
                                                      relaxed_organisms,
                                                      index,
                                                      composition_space])
                                            new_thread.start()
                                            threads[index] = new_thread
                                            started_new_calc = True

    # depending on how the loop above exited, update bookkeeping
    if not stopping_criteria.are_satisfied:
        num_finished_calcs = num_finished_calcs - 1

    # process all the calculations that were still running when the last
    # creator finished
    num_to_get = num_calcs_at_once  # number of threads left to handle
    handled_indices = []  # the indices of the threads we've already handled
    while num_to_get > 0:
        for index, thread in enumerate(threads):
            if not thread.is_alive() and index not in handled_indices:
                num_finished_calcs += 1
                relaxed_organism = relaxed_organisms[index]
                num_to_get = num_to_get - 1
                handled_indices.append(index)
                relaxed_organisms[index] = None

                # take care of relaxed organism
                if relaxed_organism is not None:
                    geometry.unpad(relaxed_organism.cell, constraints)
                    if developer.develop(relaxed_organism, composition_space,
                                         constraints, geometry, pool):
                        redundant_organism = redundancy_guard.check_redundancy(
                            relaxed_organism, whole_pop, geometry)
                        if redundant_organism is not None:  # redundant
                            if redundant_organism.is_active and \
                                    redundant_organism.epa > \
                                    relaxed_organism.epa:
                                initial_population.replace_organism(
                                    redundant_organism, relaxed_organism,
                                    composition_space)
                                progress = initial_population.get_progress(
                                    composition_space)
                                data_writer.write_data(relaxed_organism,
                                                       num_finished_calcs,
                                                       progress)
                                print('Number of energy calculations so far: '
                                      '{} '.format(num_finished_calcs))
                        else:  # no redundancy
                            stopping_criteria.check_organism(
                                relaxed_organism, redundancy_guard, geometry)
                            initial_population.add_organism(relaxed_organism,
                                                            composition_space)
                            whole_pop.append(relaxed_organism)
                            progress = initial_population.get_progress(
                                    composition_space)
                            data_writer.write_data(relaxed_organism,
                                                   num_finished_calcs,
                                                   progress)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))
                            # give a copy of the padded organism to the
                            # surrogate model
                            org_copy = copy.deepcopy(relaxed_organism)
                            geometry.pad(org_copy.cell)
                            surrogate_model.add_relaxed_organism(org_copy)

    # check if the stopping criteria were already met when making the initial
    # population
    if stopping_criteria.are_satisfied:
        quit()

    # populate the pool with the initial population
    pool.add_initial_population(initial_population, composition_space)

    # To temporarily hold relaxed organisms. The key to each relaxed organism
    # is the index of the Thread in the list threads that did the energy
    # calculation.
    relaxed_organisms = {}

    offspring_generator = general.OffspringGenerator()
    threads = []

    # create the initial batch of offspring organisms and submit them for
    # energy calculations
    for _ in range(num_calcs_at_once):
        unrelaxed_offspring = offspring_generator.make_offspring_organism(
            random, pool, variations, geometry, id_generator, whole_pop,
            developer, redundancy_guard, composition_space, constraints)
        whole_pop.append(copy.deepcopy(unrelaxed_offspring))
        geometry.pad(unrelaxed_offspring.cell)
        # check the surrogate model's prediction
        predicted_epa = surrogate_model.predict_epa(new_organism)
        if sm_output_interpreter.is_below_threshold(
                new_organism, predicted_epa, composition_space,
                initial_population.initial_population):
            predicted_cell = surrogate_model.predict_relaxed_cell(new_organism)
            if predicted_cell is not None:
                new_organism.cell = predicted_cell
            # do the energy calculation
            stopping_criteria.update_calc_counter()
            index = len(threads)
            new_thread = threading.Thread(
                target=energy_calculator.do_energy_calculation,
                args=[unrelaxed_offspring, relaxed_organisms,
                      index, composition_space])
            new_thread.start()
            threads.append(new_thread)

    # process finished calculations and start new ones
    while not stopping_criteria.are_satisfied:
        for index, thread in enumerate(threads):
            if not thread.is_alive():
                num_finished_calcs += 1
                relaxed_offspring = relaxed_organisms[index]
                relaxed_organisms[index] = None

                # take care of relaxed offspring organism
                if relaxed_offspring is not None:
                    geometry.unpad(relaxed_offspring.cell, constraints)
                    if developer.develop(relaxed_offspring, composition_space,
                                         constraints, geometry, pool):
                        # check for redundancy with the the pool first
                        redundant_organism = redundancy_guard.check_redundancy(
                            relaxed_offspring, pool.to_list(), geometry)
                        if redundant_organism is not None:  # redundant
                            if redundant_organism.epa > relaxed_offspring.epa:
                                pool.replace_organism(redundant_organism,
                                                      relaxed_offspring,
                                                      composition_space)
                                pool.compute_fitnesses()
                                pool.compute_selection_probs()
                                pool.print_summary(composition_space)
                                progress = pool.get_progress(composition_space)
                                data_writer.write_data(relaxed_offspring,
                                                       num_finished_calcs,
                                                       progress)
                                print('Number of energy calculations so far: '
                                      '{} '.format(num_finished_calcs))
                        # check for redundancy with all the organisms
                        else:
                            redundant_organism = \
                                redundancy_guard.check_redundancy(
                                    relaxed_offspring, whole_pop, geometry)
                        if redundant_organism is None:  # not redundant
                            stopping_criteria.check_organism(
                                relaxed_offspring, redundancy_guard, geometry)
                            pool.add_organism(relaxed_offspring,
                                              composition_space)
                            whole_pop.append(relaxed_offspring)

                            # check if we've added enough new offspring
                            # organisms to the pool that we can remove the
                            # initial population organisms from the front
                            # (right end) of the queue.
                            if pool.num_adds == pool.size:
                                print('Removing the initial population from '
                                      'the pool ')
                                for _ in range(len(
                                        initial_population.initial_population)):
                                    removed_org = pool.queue.pop()
                                    removed_org.is_active = False
                                    print('Removing organism {} from the '
                                          'pool '.format(removed_org.id))

                            # if the initial population organisms have already
                            # been removed from the pool's queue, then just
                            # need to pop one organism from the front (right
                            # end) of the queue.
                            elif pool.num_adds > pool.size:
                                removed_org = pool.queue.pop()
                                removed_org.is_active = False
                                print('Removing organism {} from the '
                                      'pool '.format(removed_org.id))

                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            pool.print_summary(composition_space)
                            progress = pool.get_progress(composition_space)
                            data_writer.write_data(relaxed_offspring,
                                                   num_finished_calcs,
                                                   progress)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))
                            # give a copy of the padded organism to the
                            # surrogate model
                            org_copy = copy.deepcopy(relaxed_organism)
                            geometry.pad(org_copy.cell)
                            surrogate_model.add_relaxed_organism(org_copy)

                # make another offspring organism
                if not stopping_criteria.are_satisfied:
                    unrelaxed_offspring = \
                        offspring_generator.make_offspring_organism(
                            random, pool, variations, geometry, id_generator,
                            whole_pop, developer, redundancy_guard,
                            composition_space, constraints)
                    whole_pop.append(copy.deepcopy(unrelaxed_offspring))
                    geometry.pad(unrelaxed_offspring.cell)
                    # check the surrogate model's prediction
                    predicted_epa = surrogate_model.predict_epa(new_organism)
                    if sm_output_interpreter.is_below_threshold(
                            new_organism, predicted_epa, composition_space,
                            initial_population.initial_population):
                        predicted_cell = surrogate_model.predict_relaxed_cell(
                            new_organism)
                        if predicted_cell is not None:
                            new_organism.cell = predicted_cell
                        # do the energy calculation
                        stopping_criteria.update_calc_counter()
                        new_thread = threading.Thread(
                            target=energy_calculator.do_energy_calculation,
                            args=[unrelaxed_offspring, relaxed_organisms,
                                  index, composition_space])
                        new_thread.start()
                        threads[index] = new_thread

    # process all the calculations that were still running when the
    # stopping criteria were achieved
    num_to_get = num_calcs_at_once  # how many threads we have left to handle
    handled_indices = []  # the indices of the threads we've already handled
    while num_to_get > 0:
        for index, thread in enumerate(threads):
            if not thread.is_alive() and index not in handled_indices:
                num_finished_calcs += 1
                relaxed_offspring = relaxed_organisms[index]
                num_to_get -= 1
                handled_indices.append(index)
                relaxed_organisms[index] = None

                # take care of relaxed offspring organism
                if relaxed_offspring is not None:
                    geometry.unpad(relaxed_offspring.cell, constraints)
                    if developer.develop(relaxed_offspring, composition_space,
                                         constraints, geometry, pool):
                        # check for redundancy with the pool first
                        redundant_organism = redundancy_guard.check_redundancy(
                            relaxed_offspring, pool.to_list(), geometry)
                        if redundant_organism is not None:  # redundant
                            if redundant_organism.epa > relaxed_offspring.epa:
                                pool.replace_organism(redundant_organism,
                                                      relaxed_offspring,
                                                      composition_space)
                                pool.compute_fitnesses()
                                pool.compute_selection_probs()
                                pool.print_summary(composition_space)
                                progress = pool.get_progress(composition_space)
                                data_writer.write_data(relaxed_offspring,
                                                       num_finished_calcs,
                                                       progress)
                                print('Number of energy calculations so far: '
                                      '{} '.format(num_finished_calcs))
                        # check for redundancy with all the organisms
                        else:
                            redundant_organism = \
                                redundancy_guard.check_redundancy(
                                    relaxed_offspring, whole_pop, geometry)
                        if redundant_organism is None:  # not redundant
                            pool.add_organism(relaxed_offspring,
                                              composition_space)
                            whole_pop.append(relaxed_offspring)
                            removed_org = pool.queue.pop()
                            removed_org.is_active = False
                            print('Removing organism {} from the pool '.format(
                                removed_org.id))
                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            pool.print_summary(composition_space)
                            progress = pool.get_progress(composition_space)
                            data_writer.write_data(relaxed_offspring,
                                                   num_finished_calcs,
                                                   progress)
                            print('Number of energy calculations so far: '
                                  '{} '.format(num_finished_calcs))


if __name__ == "__main__":
    main()
