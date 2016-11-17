# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Main module:

This module is run to do a genetic algorithm structure search.

"""

import yaml
import sys
import os
import general
import objects_maker
import parameters_printer
import copy
import threading 
import random


# get the path to the input file (in yaml format)
input_file = os.path.abspath(sys.argv[1]) 

# parse the input file as nested dictionaries
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

# list to hold copies of all the valid organisms made by the algorithm
whole_pop = []

# the number of energy calculations completed. This gets incremented for both successful and unsuccessful energy calculations
num_finished_calcs = 0

# create the initial population
initial_population = general.InitialPopulation(run_dir_name)

# list of threads to do the energy calculations
threads = []  

# dictionary to temporarily hold the relaxed organisms. The key to each relaxed organism is the index of the Thread in the list threads that did the energy calculation
relaxed_organisms = {} 

for creator in organism_creators:
    print('Making {} organisms with {}'.format(creator.number, creator.name))
    while not creator.is_finished:
            
        # start by doing num_calcs_at_once energy calculations, each on its own thread
        if len(threads) < num_calcs_at_once:
        
            # make a new organism - keep trying until we get one
            new_organism = creator.create_organism(id_generator, composition_space, constraints, random) 
            while new_organism == None and not creator.is_finished:
                new_organism = creator.create_organism(id_generator, composition_space, constraints, random)
            
            # this check is needed because the loop above could be exited with a None organism for attempts-based creators
            if new_organism != None:
                # unpad the organism (does nothing for bulk search)
                geometry.unpad(new_organism, constraints)
                # develop the organism
                developed_org = developer.develop(new_organism, composition_space, constraints, geometry, pool)
                if developed_org != None: # successful development
                    # check for redundancy
                    redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)
                    if redundant_org == None: # no redundancy
                        whole_pop.append(copy.deepcopy(developed_org)) # we want to add copies to whole_pop so the organisms in whole_pop don't change upon relaxation, etc.
                        stopping_criteria.update_calc_counter()  # if num calcs is one of the stopping criteria, this updates it   
                        geometry.pad(developed_org) # for bulk search, this does nothing except rotate into principal directions
                        thread_index = len(threads) # the index this thread will have once it's appended to the list threads
                        thread = threading.Thread(target=energy_calculator.do_energy_calculation, args=[developed_org, relaxed_organisms, thread_index])
                        thread.start()
                        threads.append(thread)
                        
        else:
            # check for dead threads
            for thread in threads:
                if not thread.is_alive():
                    # increment the number of completed energy calculations
                    num_finished_calcs += 1
                    # get the relaxed structure from the dictionary
                    thread_index = threads.index(thread)
                    relaxed_org = relaxed_organisms[thread_index]
                    # remove the relaxed organism from the dictionary since we got it out
                    relaxed_organisms[thread_index] = None
                    # if the relaxed organism is not None, then do development and redundancy checking
                    if relaxed_org != None:
                        
                        geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                        developed_org = developer.develop(relaxed_org, composition_space, constraints, geometry, pool)
                        if developed_org != None:
                            redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)
                            if redundant_org != None:
                                if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                                    initial_population.replace_organism(redundant_org, developed_org, composition_space)
                                    # print out info on the best organism and how many calcs have been done so far
                                    initial_population.print_progress(composition_space)
                                    print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                    
                            else: 
                                stopping_criteria.check_organism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                                initial_population.add_organism(developed_org, composition_space)
                                whole_pop.append(developed_org)
                                # print out info on the best organism and how many calcs have been done so far
                                initial_population.print_progress(composition_space)
                                print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                # if this organism was made by the current creator and the creator is success-based, then update the status of the creator
                                if creator.is_successes_based and developed_org.made_by == creator.name: 
                                    creator.update_status()
                        
                    started_new_calc = False
                    while not started_new_calc and not creator.is_finished:
                        # make another organism for the initial population
                        new_organism = creator.create_organism(id_generator, composition_space, constraints, random)
                        while new_organism == None and not creator.is_finished:
                            new_organism = creator.create_organism(id_generator, composition_space, constraints, random)
            
                        # this check is needed because the loop above could be exited with a None organism for attempts-based creators
                        if new_organism != None:
                            # unpad the organism (does nothing for bulk search)
                            geometry.unpad(new_organism, constraints)
                            # develop the organism
                            developed_org = developer.develop(new_organism, composition_space, constraints, geometry, pool)
                            if developed_org != None:
                                # check for redundancy
                                redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)
                                if redundant_org == None: # no redundancy
                                    whole_pop.append(copy.deepcopy(developed_org))
                                    stopping_criteria.update_calc_counter() 
                                    geometry.pad(developed_org) # for bulk search, this does nothing 
                                    thread = threading.Thread(target=energy_calculator.do_energy_calculation, args=(developed_org, relaxed_organisms, thread_index))
                                    thread.start()
                                    # replace the dead thread in threads with the new one
                                    threads[thread_index] = thread
                                    # set the flag so we know we've started a new energy calculation
                                    started_new_calc = True
                                

# get the output of all the calculations that were still running when the last creator finished
num_to_get = num_calcs_at_once # how many threads we have left to handle
handled_indices = [] # the indices of the threads we've already handled
while num_to_get > 0:
    for thread in threads:
        thread_index = threads.index(thread)
        if not thread.is_alive() and thread_index not in handled_indices:
            # increment the number of completed energy calculations
            num_finished_calcs += 1
            # get the relaxed structure from the dictionary
            relaxed_org = relaxed_organisms[thread_index]
            # record that this thread has been handled
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = developer.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            initial_population.replace_organism(redundant_org, developed_org, composition_space)
                            # print out info on the best organism and how many calcs have been done so far
                            initial_population.print_progress(composition_space)
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                    else: 
                        stopping_criteria.check_organism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        initial_population.add_organism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        # print out info ont he best organism and how many calcs have been done so far
                        initial_population.print_progress(composition_space)
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                        

# add the initial population to the pool (this prints out a summary of the initial population)
pool.add_initial_population(initial_population, composition_space)

threads = []
relaxed_organisms = {} # dictionary to temporarily hold the relaxed organisms. The key to each relaxed organism is the index of the Thread in the list threads that did the energy calculation
offspring_generator = general.OffspringGenerator()

# create the initial batch of offspring organisms and submit them for energy calculations
for _ in range(num_calcs_at_once):
    unrelaxed_offspring = offspring_generator.make_offspring_organism(random, pool, variations, geometry, id_generator, whole_pop, developer, redundancy_guard, composition_space, constraints)
    whole_pop.append(copy.deepcopy(unrelaxed_offspring))
    geometry.pad(unrelaxed_offspring) 
    stopping_criteria.update_calc_counter()
    thread_index = len(threads) # the index this thread will have once it's appended to the list threads
    thread = threading.Thread(target=energy_calculator.do_energy_calculation, args=[unrelaxed_offspring, relaxed_organisms, thread_index])
    thread.start()
    threads.append(thread)
    
# continue the search until stopping criteria are met
while not stopping_criteria.are_satisfied:
    # check for dead threads
    for thread in threads:
        if not thread.is_alive():
            # increment the number of completed energy calculations
            num_finished_calcs += 1
            # get the relaxed structure from the dictionary
            thread_index = threads.index(thread)
            relaxed_org = relaxed_organisms[thread_index]
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = developer.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    # check for redundancy with the organisms in the pool first
                    redundant_org = redundancy_guard.check_redundancy(developed_org, pool.to_list())
                    if redundant_org != None:
                        if redundant_org.epa > developed_org.epa:
                            # replace the organism
                            pool.replace_organism(redundant_org, developed_org, composition_space)
                            # recompute fitnesses and selection probabilities
                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            # print out a summary of the pool
                            pool.print_summary(composition_space)
                            # print out the progress of the search - either the best value (for epa) or the volume of the convex hull (for pd)
                            pool.print_progress(composition_space)
                            # print out how many energy calculations have been done so far
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                    # if not redundant with an organism in the pool, check for redundancy with all the organisms in whole_pop
                    else:
                        redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)     
                    # if no redundancy in either the pool or whole_pop, then add the new organism to the pool
                    if redundant_org == None: 
                        stopping_criteria.check_organism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        pool.add_organism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        
                        # check if we've added enough new offspring organisms to the pool that we can remove the initial
                        # population organisms from the front (right end) of the queue. 
                        if pool.num_adds == pool.size:
                            print('Removing the initial population from the pool ')
                            for _ in range(len(initial_population.initial_population)):
                                removed_org = pool.queue.pop()
                                removed_org.is_active = False
                                print('Removing organism {} from the pool '.format(removed_org.id))
                        
                        # if the initial population organisms have already been removed from the pool's queue, then just need 
                        # to pop one organism from the front (right end) of the queue.
                        elif pool.num_adds > pool.size:
                            removed_org = pool.queue.pop()
                            removed_org.is_active = False
                            print('Removing organism {} from the pool '.format(removed_org.id))

                        # recompute fitnesses and selection probabilities
                        pool.compute_fitnesses()
                        pool.compute_selection_probs()
                            
                        # print out a summary of the pool
                        pool.print_summary(composition_space)
                        pool.print_progress(composition_space)
                        # print out how many energy calculations have been done so far
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                 
            # make another offspring organism if the stopping criteria aren't satisfied
            if not stopping_criteria.are_satisfied:
                unrelaxed_offspring = offspring_generator.make_offspring_organism(random, pool, variations, geometry, id_generator, whole_pop, developer, redundancy_guard, composition_space, constraints)
                whole_pop.append(copy.deepcopy(unrelaxed_offspring))
                geometry.pad(unrelaxed_offspring) 
                stopping_criteria.update_calc_counter()
                thread = threading.Thread(target=energy_calculator.do_energy_calculation, args=[unrelaxed_offspring, relaxed_organisms, thread_index])
                thread.start()
                # replace the dead thread in threads with the new one
                threads[thread_index] = thread


# get the output of all the calculations that were still running when the stopping criteria were met
num_to_get = num_calcs_at_once # how many threads we have left to handle
handled_indices = [] # the indices of the threads we've already handled
while num_to_get > 0:
    for thread in threads:
        thread_index = threads.index(thread)
        if not thread.is_alive() and thread_index not in handled_indices:
            # increment the number of completed energy calculations
            num_finished_calcs += 1
            # get the relaxed structure from the dictionary
            relaxed_org = relaxed_organisms[thread_index]
            # record that this thread has been handled
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = developer.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    # check for redundancy with the organisms in the pool first
                    redundant_org = redundancy_guard.check_redundancy(developed_org, pool.to_list())
                    if redundant_org != None:
                        if redundant_org.epa > developed_org.epa:
                            # replace the organism
                            pool.replace_organism(redundant_org, developed_org, composition_space)
                            # recompute fitnesses and selection probabilities
                            pool.compute_fitnesses()
                            pool.compute_selection_probs()
                            # print out a summary of the pool
                            pool.print_summary(composition_space)
                            # print out the progress of the search - either the best value (for epa) or the volume of the convex hull (for pd)
                            pool.print_progress(composition_space)
                            # print out how many energy calculations have been done so far
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                    # if not redundant with an organism in the pool, check for redundancy with all the organisms in whole_pop
                    else:
                        redundant_org = redundancy_guard.check_redundancy(developed_org, whole_pop)     
                    # if no redundancy in either the pool or whole_pop, then add the new organism to the pool
                    if redundant_org == None: 
                        pool.add_organism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        removed_org = pool.queue.pop()
                        removed_org.is_active = False
                        print('Removing organism {} from the pool '.format(removed_org.id))
                        # recompute fitnesses and selection probabilities
                        pool.compute_fitnesses()
                        pool.compute_selection_probs()
                        # print out a summary of the pool
                        pool.print_summary(composition_space)
                        pool.print_progress(composition_space)
                        # print out how many energy calculations have been done so far
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))