import yaml
import sys
import os
import classes
import objects_maker
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
    print('Quitting')

# just for testing. Normally we'll read the input file as an argument
#with open('/Users/benjaminrevard/GASPy/gaspy/src/gaspy_input.yaml', 'r') as f:
#    parameters = yaml.load(f)

# this line is just for testing. Normally the code will be executed in the folder where the search is to be done...
#os.chdir('/Users/benjaminrevard/testing/gaspy_testing') 

# make the objects needed by the algorithm
objects_dict = objects_maker.makeObjects(parameters)

# get the objects from the dictionary
run_dir_name = objects_dict['run_dir_name']
organism_creators = objects_dict['organism_creators']
num_calcs_at_once = objects_dict['num_calcs_at_once']
composition_space = objects_dict['composition_space']
constraints = objects_dict['constraints']
geometry = objects_dict['geometry']
development = objects_dict['development']
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

# if we're doing a phase diagram search, print out the composition space endpoints
if composition_space.objective_function == 'pd':
    composition_space.printEndpoints()

# print the search parameters to a file in the run directory
objects_maker.printParameters(objects_dict)

# list to hold copies of all the valid organisms made by the algorithm
whole_pop = []

# the number of energy calculations completed. This only gets incremented for successful energy calcs (ones that returned a structure and an energy)
num_finished_calcs = 0

# create the initial population
initial_population = classes.InitialPopulation(run_dir_name)

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
            new_organism = creator.createOrganism(id_generator, composition_space, constraints, random) 
            while new_organism == None and not creator.is_finished:
                new_organism = creator.createOrganism(id_generator, composition_space, constraints, random)
            
            # this check is needed because the loop above could be exited with a None organism for attempts-based creators
            if new_organism != None:
                # unpad the organism (does nothing for bulk search)
                geometry.unpad(new_organism, constraints)
                # develop the organism
                developed_org = development.develop(new_organism, composition_space, constraints, geometry, pool)
                if developed_org != None: # successful development
                    # check for redundancy
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org == None: # no redundancy
                        whole_pop.append(copy.deepcopy(developed_org)) # we want to add copies to whole_pop so the organisms in whole_pop don't change upon relaxation, etc.
                        stopping_criteria.updateCalcCounter()  # if num calcs is one of the stopping criteria, this updates it   
                        geometry.pad(developed_org) # for bulk search, this does nothing except rotate into principal directions
                        thread_index = len(threads) # the index this thread will have once it's appended to the list threads
                        thread = threading.Thread(target=energy_calculator.doEnergyCalculation, args=[developed_org, relaxed_organisms, thread_index])
                        thread.start()
                        threads.append(thread)
                        
        else:
            # check for dead threads
            for thread in threads:
                if not thread.is_alive():
                    # get the relaxed structure from the dictionary
                    thread_index = threads.index(thread)
                    relaxed_org = relaxed_organisms[thread_index]
                    # remove the relaxed organism from the dictionary since we got it out
                    relaxed_organisms[thread_index] = None
                    # if the relaxed organism is not None, then do development and redundancy checking
                    if relaxed_org != None:
                        num_finished_calcs += 1
                        geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                        developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, pool)
                        if developed_org != None:
                            redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                            if redundant_org != None:
                                if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                                    initial_population.replaceOrganism(redundant_org, developed_org, composition_space)
                                    # print out info on the best organism and how many calcs have been done so far
                                    initial_population.printProgress(composition_space)
                                    print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                    
                            else: 
                                stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                                initial_population.addOrganism(developed_org, composition_space)
                                whole_pop.append(developed_org)
                                # print out info on the best organism and how many calcs have been done so far
                                initial_population.printProgress(composition_space)
                                print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                # if this organism was made by the current creator and the creator is success-based, then update the status of the creator
                                if creator.is_successes_based and developed_org.made_by == creator.name: 
                                    creator.updateStatus()
                        
                    started_new_calc = False
                    while not started_new_calc and not creator.is_finished:
                        # make another organism for the initial population
                        new_organism = creator.createOrganism(id_generator, composition_space, constraints, random)
                        while new_organism == None and not creator.is_finished:
                            new_organism = creator.createOrganism(id_generator, composition_space, constraints, random)
            
                        # this check is needed because the loop above could be exited with a None organism for attempts-based creators
                        if new_organism != None:
                            # unpad the organism (does nothing for bulk search)
                            geometry.unpad(new_organism, constraints)
                            # develop the organism
                            developed_org = development.develop(new_organism, composition_space, constraints, geometry, pool)
                            if developed_org != None:
                                # check for redundancy
                                redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                                if redundant_org == None: # no redundancy
                                    whole_pop.append(copy.deepcopy(developed_org))
                                    stopping_criteria.updateCalcCounter() 
                                    geometry.pad(developed_org) # for bulk search, this does nothing 
                                    thread = threading.Thread(target=energy_calculator.doEnergyCalculation, args=(developed_org, relaxed_organisms, thread_index))
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
            # get the relaxed structure from the dictionary
            relaxed_org = relaxed_organisms[thread_index]
            # record that this thread has been handled
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                # update the number of energy calculations
                num_finished_calcs += 1
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            initial_population.replaceOrganism(redundant_org, developed_org, composition_space)
                            # print out info on the best organism and how many calcs have been done so far
                            initial_population.printProgress(composition_space)
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                    else: 
                        stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        initial_population.addOrganism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        # print out info ont he best organism and how many calcs have been done so far
                        initial_population.printProgress(composition_space)
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                        

# add the initial population to the pool (this prints out a summary of the initial population)
pool.addInitialPopulation(initial_population, composition_space)

threads = []
relaxed_organisms = {} # dictionary to temporarily hold the relaxed organisms. The key to each relaxed organism is the index of the Thread in the list threads that did the energy calculation
offspring_generator = classes.OffspringGenerator()

# create the initial batch of offspring organisms and submit them for energy calculations
for _ in range(num_calcs_at_once):
    unrelaxed_offspring = offspring_generator.makeOffspringOrganism(random, pool, variations, geometry, id_generator, whole_pop, development, redundancy_guard, composition_space, constraints)
    whole_pop.append(copy.deepcopy(unrelaxed_offspring))
    geometry.pad(unrelaxed_offspring) 
    stopping_criteria.updateCalcCounter()
    thread_index = len(threads) # the index this thread will have once it's appended to the list threads
    thread = threading.Thread(target=energy_calculator.doEnergyCalculation, args=[unrelaxed_offspring, relaxed_organisms, thread_index])
    thread.start()
    threads.append(thread)
    
# continue the search until stopping criteria are met
while not stopping_criteria.are_satisfied:
    # check for dead threads
    for thread in threads:
        if not thread.is_alive():
            # get the relaxed structure from the dictionary
            thread_index = threads.index(thread)
            relaxed_org = relaxed_organisms[thread_index]
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                num_finished_calcs += 1
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            pool.replaceOrganism(redundant_org, developed_org, composition_space)
                            # print out a summary of the pool
                            pool.printSummary(composition_space)
                            # print out the progress of the search - either the best value (for epa) or the volume of the convex hull (for pd)
                            pool.printProgress(composition_space)
                            # print out how many energy calculations have been done so far
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                    else: 
                        stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        pool.addOrganism(developed_org, composition_space)
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
                            
                        # print out a summary of the pool
                        pool.printSummary(composition_space)
                        pool.printProgress(composition_space)
                        # print out how many energy calculations have been done so far
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                                 
            # make another offspring organism if the stopping criteria aren't satisfied
            if not stopping_criteria.are_satisfied:
                unrelaxed_offspring = offspring_generator.makeOffspringOrganism(random, pool, variations, geometry, id_generator, whole_pop, development, redundancy_guard, composition_space, constraints)
                whole_pop.append(copy.deepcopy(unrelaxed_offspring))
                geometry.pad(unrelaxed_offspring) 
                stopping_criteria.updateCalcCounter()
                thread = threading.Thread(target=energy_calculator.doEnergyCalculation, args=[unrelaxed_offspring, relaxed_organisms, thread_index])
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
            # get the relaxed structure from the dictionary
            relaxed_org = relaxed_organisms[thread_index]
            # record that this thread has been handled
            num_to_get = num_to_get - 1
            handled_indices.append(thread_index)
            # remove the relaxed organism from the dictionary since we got it out
            relaxed_organisms[thread_index] = None
            # if the relaxed organism is not None, then do development and redundancy checking
            if relaxed_org != None:
                num_finished_calcs += 1
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, pool)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            pool.replaceOrganism(redundant_org, developed_org, composition_space)
                            pool.printSummary(composition_space)
                            pool.printProgress(composition_space)
                            print('Number of energy calculations so far: {} '.format(num_finished_calcs))
                            
                    else: 
                        pool.addOrganism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        removed_org = pool.queue.pop()
                        removed_org.is_active = False
                        print('Removing organism {} from the pool '.format(removed_org.id))
                        # print out a summary of the pool
                        pool.printSummary(composition_space)
                        pool.printProgress(composition_space)
                        # print out how many energy calculations have been done so far
                        print('Number of energy calculations so far: {} '.format(num_finished_calcs))
    
    
    
    








                        

# for testing
# print out the energies of all the structures in the initial population, and write out there structures to poscar files so we can look at them
#for organism in initial_population.initial_population:
#    print('Organism {} has energy per atom {}'.format(organism.id, organism.epa))
#    organism.structure.to('poscar', '/n/srv/brevard/testing/gaspy_testing/garun_test/{}.vasp'.format(organism.id))




# get an organism from the org creator (only one in this case)
#organism = None
#while organism == None:
    # generate an organism
#    organism = organism_creators[0].createOrganism(id_generator, composition_space, constraints, random)
    # develop the organism
#    if organism != None:
#        geometry.unpad(organism, constraints)
#        developed_organism = development.develop(organism, composition_space, constraints, geometry, None)
#        if developed_organism == None:
#            organism = None


# pad the organism
#geometry.pad(developed_organism)
#print(developed_organism.structure)
#print('')
# call the doEnergyCalculation method
#relaxed_organism = energy_calculator.doEnergyCalculation(developed_organism) 
#if relaxed_organism != None:
#    print(relaxed_organism.total_energy)
#    print(relaxed_organism.structure)
#    whole_pop.append(relaxed_organism)



# test the doVariation method. The first arg is just a placeholder for the pool
#offspring = variations[0].doVariation(None, random, geometry, id_generator)

# write out the offspring structure to a file so we can look at it
#offspring.structure.to('poscar', '/n/srv/brevard/structures/permutation.vasp')
#print('')
#print(offspring.structure)




 
    
    
              

# write out developed structure to cif or poscar so I can look at it
#developed_organism.structure.to('poscar', '{}_dev.vasp'.format(developed_organism.id))
# pad the structure
#geometry.pad(developed_organism)
# write out the padded structure to a file so I can look at it - check atom positions were shifted properly
#developed_organism.structure.to('poscar', '{}_padded.vasp'.format(developed_organism.id))
# unpad the structure
#geometry.unpad(developed_organism, constraints)
# write out unpadded structure to a file so I can look at it - check atom positions were shifted properly
#developed_organism.structure.to('poscar', '{}_unpadded.vasp'.format(developed_organism.id))



#print(stopping_criteria.num_energy_calcs)
#print(stopping_criteria.value_achieved)
#print(stopping_criteria.found_structure)


#print(developed_organism.id)
#print(developed_organism.structure)
#print('')
#print(constraints.per_species_mids)

#energy_calculator.doEnergyCalculation(developed_organism)

# TODO: now try to create an organism with one of the creators, and then pass it to the doEnergyCalculation method of GulpEnergyCalculator





#print(energy_calculator.potential)

#for i in range(len(energy_calculator.potential)):
#    print(energy_calculator.potential[i])



















# just for testing

#foc = organism_creators[0]
#roc = organism_creators[1]

#print(foc.files)
#print(roc.number)

#file_org = foc.createOrganism()
#for i in range(len(foc.files)):
#    if file_org == None:
#        print('pymatgen couldnt read the structure in')
#    else:
#        print(file_org.structure)
#        print("")











#random_org = organism_creators[0].createOrganism(composition_space, constraints)
#if random_org != None:
#    print(random_org.structure)
#else:
#    print("createOrg returned None - stupid volume scaling didn't work again")


#print(organism_creators[0].number)
#print(organism_creators[0].volume)






# make an organism 
#lattice1 = [[5, 0, 0], [0, 5, 0], [0, 0, 5]]
#species1 = ["C", "Si"]
#coordinates1 = [[0.3, 0.5, 0.5],[0.7, 0.5, 0.5]]
#structure1 = classes.Structure(lattice1, species1, coordinates1)
#org1 = classes.Organism(structure1)

# make another organism 
#lattice2 = [[10, 0, 0], [0, 5, 0], [0, 0, 5]]
#species2 = ["C", "Si", "C", "Si"]
#coordinates2 = [[0.15, 0.5, 0.5],[0.35, 0.5, 0.5], [0.65, 0.5, 0.5], [0.85, 0.5, 0.5]]
#structure2 = classes.Structure(lattice2, species2, coordinates2)
#org2 = classes.Organism(structure2)

#d1 = development.develop(org1, composition_space, constraints, geometry, None)
#d2 = development.develop(org2, composition_space, constraints, geometry, None)

#print(d1.structure)
#print("")
#print(d2.structure)

# make-shift whole_pop list
#whole_pop = []
#whole_pop.append(d1)

# check for redundancy
#print("")
#redundancy_guard.checkRedundancy(d2, whole_pop)

# see if the structures got changed
#print("")
#print("")
#print(d1.structure)
#print("")
#print(d2.structure)

# testing atomic masses and solid state densities
#carbon = Element('C')
#silicon = Element('Si')

#print("")
#print(carbon.atomic_mass)
#print("")
#print(silicon.atomic_mass)
#print("")
#print(carbon.density_of_solid)
#print("")
#rint(silicon.density_of_solid)


#undeveloped = Poscar(org1.structure)
#print(undeveloped.get_string())
#print("")

#developed = development.develop(org1, None)

#if developed != None:
#    print(Poscar(developed.structure).get_string())
#else:
#    print("org failed development")








# make the creator objects

#for i in parameters['InitialPopulation']:
#    if i == 'random':
#        print('make RandomOrgCreator')
#    elif i == 'fromPoscars':
#        print('make PoscarOrgCreator')
        
    #print parameters['InitialPopulation'][i]

# construct the creators and put them in a list
#creators = []


