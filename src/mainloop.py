# Outline of the algorithm
import yaml
import sys
import os
import classes
import copy
import threading 
import random
import time
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.util.convergence import id_generator
from _ast import Pass

# get the pash to the input file (in yaml format)
#input_file = os.path.abspath(sys.argv[1]) 

# parse the input file as nested dictionaries
# TODO: this should be in try-except clause, since opening the file is something that could fail (e.g., wrong path given, etc.)
#with open(input_file, 'r') as f:
#    parameters = yaml.load(f)

#with open('/Users/benjaminrevard/GASPy/gaspy/src/gaspy_input.yaml', 'r') as f:
#    parameters = yaml.load(f)


with open('/n/srv/brevard/python_GA/gaspy/src/gaspy_input.yaml', 'r') as f:
    parameters = yaml.load(f)

    
# make the composition space object
if 'CompositionSpace' in parameters:
    composition_space = classes.CompositionSpace(parameters['CompositionSpace'])
else:
    print('Input file must contain a "CompositionSpace" block.')
    print("Quitting...")
    quit()
    
# make the constraints object
if 'Constraints' in parameters:
    constraints = classes.Constraints(parameters['Constraints'], composition_space)
else:
    # if no Constraints block is given in the input file, then just use default values for everything
    constraints = classes.Constraints('default', composition_space)
    
# make the geometry object
if 'Geometry' in parameters:
    geometry = classes.Geometry(parameters['Geometry'])
else:
    # if no Geometry block is given in the input file, then just use default values for everything
    geometry = classes.Geometry('default')
    
# make the niggli object
if 'Niggli' in parameters:
    if parameters['Niggli'] == None or parameters['Niggli'] == 'default':
        niggli = True
    else:
        niggli = parameters['Niggli']
else:
    # if no Niggli block is given in the input file, set it to True
    niggli = True
    
# make the scale density object
if 'ScaleDesnsity' in parameters:
    if parameters['ScaleDensity'] == None or parameters['ScaleDensity'] == 'default':
        scale_density = True
    else:
        scale_density = parameters['ScaleDensity']
else:
    # if no ScaleDensity block is given in the input file, set it to True
    scale_density = True
   
# make the development object
development = classes.Development(niggli, scale_density)

# make the redundancy guard object
if 'RedundancyGuard' in parameters:
    redundancy_guard = classes.RedundancyGuard(parameters['RedundancyGuard'])
else:
    redundancy_guard = classes.RedundancyGuard('default')
    
# make the id generator
id_generator = classes.IDGenerator()
    
# make the organism creators
organism_creators = []
if 'InitialPopulation' in parameters:
    # the random organism creator
    if 'random' in parameters['InitialPopulation']:
        random_organism_creator = classes.RandomOrganismCreator(parameters['InitialPopulation']['random'], composition_space)
        organism_creators.append(random_organism_creator)
    # the from files organism creator
    if 'from_files' in parameters['InitialPopulation']:
        # if nothing is given after the from_files flag
        if parameters['InitialPopulation']['from_files'] == None:
            print('The path to the folder containing the files must be provided. Please use the "path_to_folder" flag.')
            print('Quitting...')
            quit()
        elif 'path_to_folder' in parameters['InitialPopulation']['from_files']:
            given_path = parameters['InitialPopulation']['from_files']['path_to_folder']
            # check if no path was given after path_to_folder flag
            if parameters['InitialPopulation']['from_files']['path_to_folder'] == None:
                print('The path to the folder containing the files for the initial population must be provided. Please give the path after the "path_to_folder" flag.')
                print('Quitting...')
                quit()
            elif not os.path.exists(given_path):
                print('The given folder containing structures for the initial population does not exist.')
                print('Quitting...')
                quit()
            # if the folder exists, check that it contains files
            elif len([f for f in os.listdir(given_path) if os.path.isfile(os.path.join(given_path, f))]) == 0:
                print('The given folder containing structures for the initial population does not contain any files.')
                print('Quitting...')
                quit()
            # otherwise, the directory exists and contains files. Instantiate the files organism creator
            else:
                files_organism_creator = classes.FileOrganismCreator(given_path)
                organism_creators.append(files_organism_creator)
        # in case an incorrect flag is given after the from_files flag 
        else:
            print('Incorrect flag given after "from_files" in the InitialPopulation block. Please use the "path_to_folder" flag.')
            print('Quitting...')
            quit()
        # TODO: if other organism creators are used, they should be instantiated here
# if no method specified for making the initial population, then make a random one by default
else:
    random_organism_creator = classes.RandomOrganismCreator('default', composition_space)
    organism_creators.append(random_organism_creator)

# If more than one organism creator, sort them so that the attempts-based ones are at the front and the successes-based ones are at the back
if len(organism_creators) > 1:
    organism_creators.sort(key=lambda x: x.is_successes_based)
    
# the number of energy calculations to run at a time
if 'NumCalcsAtOnce' in parameters:
    if parameters['NumCalcsAtOnce'] == None or parameters['NumCalcsAtOnce'] == 'default':
        num_calcs_at_once = 1
    else:
        num_calcs_at_once = parameters['NumCalcsAtOnce']
else:
    # if no NumCalcsAtOnce block is given in the input file, set it to 1
    num_calcs_at_once = 1
    
# get the run title, if specified
if 'RunTitle' in parameters:
    if parameters['RunTitle'] == None:
        run_title = ''
    else:
        run_title = parameters['RunTitle']
else:
    # if no RunTitle block is given in the input file, set it to an empty string
    run_title = ''
    
# using the run title, construct the name of the run directory
if run_title == '':
    run_dir_name = 'garun'
else:
    run_dir_name = 'garun' + '_' + str(run_title)

# the energy calculator
if 'EnergyCode' not in parameters:
    print('A method for calculating energy must be provided. Please use the "EnergyCode" flag.')
    print('Quitting...')
    quit()
elif parameters['EnergyCode'] == None:
    print('An energy code must be specified after the "EnergyCode" flag.')
    print('Quitting...')
    quit()
elif 'gulp' in parameters['EnergyCode']:
    if parameters['EnergyCode']['gulp'] == None:
        print('No gulp header or potential files given. Please use the "header_file" and "potential_file" flags.')
        print('Quitting...')
        quit()
    else:
        # get the header file
        if 'header_file' not in parameters['EnergyCode']['gulp']:
            print('A gulp header file must be provided. Please use the "header_file" flag.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['header_file'] == None:
            print('No gulp header file given after the "header_file" flag. Please provide one.')
            print('Quitting...')
            quit()
        else:
            # check that the given header file exists
            header_file_path = parameters['EnergyCode']['gulp']['header_file']
            if not os.path.exists(header_file_path):
                print('The given gulp header file does not exist.')
                print('Quitting...')
                quit() 
        # get the potential file 
        if 'potential_file' not in parameters['EnergyCode']['gulp']:
            print('A gulp potential file must be provided. Please use the "potential_file" flag.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['potential_file'] == None:
            print('No gulp potential file given after the "potential_file" flag. Please provide one.')
            print('Quitting...')
            quit()   
        else:
            # check that the given potential file exists
            potential_file_path = parameters['EnergyCode']['gulp']['potential_file']
            if not os.path.exists(potential_file_path):
                print('The given gulp potential file does not exist.')
                print('Quitting...')
                quit()
        # if we made it this far, then both the header and potential file exist, so make the gulp energy calculator
        energy_calculator = classes.GulpEnergyCalculator(header_file_path, potential_file_path, run_dir_name)
        
elif 'lammps' in parameters['EnergyCode']:
    # TODO: read in the stuff for setting up lammps calcs
    pass
elif 'vasp' in parameters['EnergyCode']:
    # TODO: read in the stuff for vasp calcs
    pass
# TODO: add other energy codes here
else:
    print('The given energy code name is invalid.')
    print('Quitting...')
    quit()
    
# the stopping criteria, if specified
if 'StoppingCriteria' in parameters:
    # if the entire block is empty or set to default, then set everything to default
    if parameters['StoppingCriteria'] == None or parameters['StoppingCriteria'] == 'default':
        stopping_criteria = classes.StoppingCriteria(None, composition_space)
    # need to check that if the found structure option is used and a file path is given, that it's valid
    elif 'found_structure' in parameters['StoppingCriteria']:
        if parameters['StoppingCriteria']['found_structure'] == None or parameters['StoppingCriteria']['found_structure'] == 'default':
            stopping_criteria = classes.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
        else:
            # check that the given path is valid
            given_path = parameters['StoppingCriteria']['found_structure']
            if not os.path.exists(given_path):
                print('The file containing the structure to find does not exist.')
                print('Quitting...')
                quit()
            # check that the file has the correct suffix or prefix
            elif given_path.endswith('.cif') or given_path.startswith('POSCAR'):
                # check that file can be read properly
                try:
                    found_struct = Structure.from_file(given_path)
                    stopping_criteria = classes.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
                # couldn't read the structure from the given file
                except ValueError:
                    print('Error reading the structure to find from the given file.')
                    print('Quitting...')
                    quit()
            else:
                print('File containing structure to find must be in POSCAR or cif format and begin with POSCAR or end with .cif, respectively.')
                print('Quitting...')
                quit()
    else:
        stopping_criteria = classes.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
# use defaults if no StoppingCriteria block has been used in the input file
else:
    stopping_criteria = classes.StoppingCriteria(None, composition_space)

# list to hold the variation objects
variations = []

# the default fractions for the variations
default_mating_fraction = 0.8
default_structure_mut_fraction = 0.1
default_num_stoichs_mut_fraction = 0.1
default_permutation_fraction = 0.0

# TODO: get rid of code duplication below
# create the variation objects
if 'Variations' not in parameters:
    # no variations specified, so make default ones
    mating = classes.Mating({'fraction': default_mating_fraction})
    structure_mut = classes.StructureMut({'fraction': default_structure_mut_fraction})
    num_stoichs_mut = classes.NumStoichsMut({'fraction': default_num_stoichs_mut_fraction})
    permutation = classes.Permutation({'fraction': default_permutation_fraction}, composition_space)
    
    # add them to the list
    variations.append(mating)
    variations.append(structure_mut)
    variations.append(num_stoichs_mut)
    variations.append(permutation)
    
elif parameters['Variations'] == None or parameters['Variations'] == 'default':
    # no variations specified, so make default ones
    mating = classes.Mating({'fraction': default_mating_fraction})
    structure_mut = classes.StructureMut({'fraction': default_structure_mut_fraction})
    num_stoichs_mut = classes.NumStoichsMut({'fraction': default_num_stoichs_mut_fraction})
    permutation = classes.Permutation({'fraction': default_permutation_fraction},  composition_space)
    
    # add them to the list
    variations.append(mating)
    variations.append(structure_mut)
    variations.append(num_stoichs_mut)
    variations.append(permutation)
    
else:
    # make each variation that's been specified
    # mating
    if 'Mating' not in parameters['Variations']:
        # Mating sub block was not used, so don't make a Mating variation
        pass
    elif parameters['Variations']['Mating'] == None:
        # Mating sub block was left blank
        print('If the "Mating" flag is used, its "fraction" flag must also be set.')
        print('Quittig...')
        quit()
    else:
        # Mating sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
        if parameters['Variations']['Mating']['fraction'] == None or parameters['Variations']['Mating']['fraction'] == 'default':
            print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the Mating variation.')
            print('Quitting...')
            quit()
        else:
            # make a Mating variation object with the supplied parameters
            mating = classes.Mating(parameters['Variations']['Mating'])
            variations.append(mating) 
               
    # structure mutation
    if 'StructureMut' not in parameters['Variations']:
        # StructureMut sub block was not used, so don't make a StructureMut variation
        pass
    elif parameters['Variations']['StructureMut'] == None:
        # StructureMut sub block was left blank
        print('If the "StructureMut" flag is used, its "fraction" flag must also be set.')
        print('Quittig...')
        quit()
    else:
        # StructureMut sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
        if parameters['Variations']['StructureMut']['fraction'] == None or parameters['Variations']['StructureMut']['fraction'] == 'default':
            print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the StructureMut variation.')
            print('Quitting...')
            quit()
        else:
            # make a StructureMut variation object with the supplied parameters
            structure_mut = classes.StructureMut(parameters['Variations']['StructureMut'])
            variations.append(structure_mut)
               
    # mutating the number of stoichiometries worth of atoms in the cell
    if 'NumStoichsMut' not in parameters['Variations']:
        # NumStoichsMut sub block was not used, so don't make a NumStoichsMut variation
        pass
    elif parameters['Variations']['NumStoichsMut'] == None:
        # NumStoichsMut block was left blank
        print('If the "NumStoichsMut" flag is used, its "fraction" flag must also be set.')
        print('Quittig...')
        quit()
    else:
        # NumStoichsMut sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
        if parameters['Variations']['NumStoichsMut']['fraction'] == None or parameters['Variations']['NumStoichsMut']['fraction'] == 'default':
            print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the NumStoichsMut variation.')
            print('Quitting...')
            quit()
        else:
            # make a NumStoichsMut variation object with the supplied parameters
            num_stoichs_mut = classes.NumStoichsMut(parameters['Variations']['NumStoichsMut']) 
            variations.append(num_stoichs_mut) 
             
    # permutation (swapping atoms)
    if 'Permutation' not in parameters['Variations']:
        # Permutation sub block was not used, so don't make a Permutation variation
        pass
    elif parameters['Variations']['Permutation'] == None:
        # Permutation sub block was left blank
        print('If the "Permutation" flag is used, its "fraction" flag must also be set.')
        print('Quitting...')
        quit()
    else:
        # Permutation sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
        if parameters['Variations']['Permutation']['fraction'] == None or parameters['Variations']['Permutation']['fraction'] == 'default':
            print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the Permutation variation.')
            print('Quitting...')
        else:
            # make a Permutation variation object with the supplied parameters
            permutation = classes.Permutation(parameters['Variations']['Permutation'], composition_space)
            variations.append(permutation)
        
# check that at least one variation has been used. This shouldn't happen...
if len(variations) == 0:
    print("At least one variation must be used. Either leave entire 'Variations' block blank to use default variations, or specify at least one variation within the 'Variations' block.")
    print('Quitting...')
    quit()

# check that the variations fraction variables sum to 1
frac_sum = 0
for variation in variations:
    frac_sum = frac_sum + variation.fraction
if frac_sum != 1.0:
    print("The Variations' fraction values must sum to 1.")
    print('Quitting...')
    quit()

# create the pool object
if 'Pool' not in parameters:
    pool = classes.Pool(None, composition_space, run_dir_name)
else:
    pool = classes.Pool(parameters['Pool'], composition_space, run_dir_name)
    
# create the selection object
if 'Selection' not in parameters:
    selection = classes.SelectionProbDist(None, pool.size)
else:
    selection = classes.SelectionProbDist(parameters['Selection'], pool.size)
    
# give the selection object to the pool
pool.selection = selection


################### actual algorithm is below here #######################


os.chdir('/n/srv/brevard/testing/gaspy_testing') # this line is just for testing. Normally the code will be executed in the folder where the search is to be done...
# make the run directory
os.mkdir(str(os.getcwd()) + '/' + run_dir_name)
# make the temp subdirectory where the energy calculations will be run
os.mkdir(str(os.getcwd()) + '/' + run_dir_name + '/temp')

# list to hold copies of all the valid organisms made by the algorithm
whole_pop = []

# create the initial population
initial_population = classes.InitialPopulation(run_dir_name)
threads = []  # list of threads to do the energy calculations. TODO: should we use a dictionary for this instead?
relaxed_organisms = {} # dictionary to temporarily hold the relaxed organisms. The key to each relaxed organism is the index of the Thread in the list threads that did the energy calculation
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
                developed_org = development.develop(new_organism, composition_space, constraints, geometry, None)
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
                        geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                        developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, None)
                        if developed_org != None:
                            redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                            if redundant_org != None:
                                if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                                    initial_population.replaceOrganism(redundant_org, developed_org)
                            else: 
                                stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                                initial_population.addOrganism(developed_org)
                                whole_pop.append(developed_org)
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
                            developed_org = development.develop(new_organism, composition_space, constraints, geometry, None)
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
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, None)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            initial_population.replaceOrganism(redundant_org, developed_org)
                    else: 
                        stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        initial_population.addOrganism(developed_org)
                        whole_pop.append(developed_org)
                        

# add the initial population to the pool
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
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, None)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            pool.replaceOrganism(redundant_org, developed_org, composition_space)
                    else: 
                        stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                        pool.addOrganism(developed_org, composition_space)
                        whole_pop.append(developed_org)
                        
                        # check if we've added enough new offspring organisms to the pool that we can trim off the initial
                        # population organisms from the queue. 
                        # TODO: check the logic here
                        if pool.num_adds == (pool.size - pool.num_promoted):
                            print('Trimming the initial population from the pool')
                            for i in range(0, len(pool.queue) - (pool.size - pool.num_promoted)):
                                removed_org = pool.queue.pop()
                                removed_org.is_active = False
                
                        # check if we've already added more than enough offspring to the pool to allow trimming off the 
                        # orgs from the initial population, in which case it has already happened, and we just need to
                        # remove one organism from the front (right end) of the pool's queue
                        elif (pool.num_adds > (pool.size - pool.num_promoted)):
                            removed_org = pool.queue.pop()
                            removed_org.is_active = False
                                 
            # make another offspring organism
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
                geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing
                developed_org = development.develop(relaxed_org, composition_space, constraints, geometry, None)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None:
                        if redundant_org.is_active and redundant_org.epa > developed_org.epa:  
                            pool.replaceOrganism(redundant_org, developed_org, composition_space)
                    else: 
                        pool.addOrganism(developed_org, composition_space)
                        whole_pop.append(developed_org)




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


'''


# get the parameters from the input file, and store them as nested dictionaries
parameters_parser = ParametersParser(path_to_input_file)
dict_of_parameters = parameters_parser.parse()

# create objects that will be used throughout the search by passing the appropriate dictionary 
# of input data to the constructor.

# make the organism creator objects (from poscars, random, etc.) and put them in a list
organism_creators = []
for name in dict_of_parameters['Initial Population']
    Name = str(name)
    organism_creators.append(Name(name))
# TODO: sort the list of creators so that the num_attempts ones are at the front and the num_successes ones are at the back

# create variation objects (Mutation, Mating, etc.) and put them in a list
# don't actually need these until after we've finished making the initial population
variations = []
for name in dict_of_parameters['Variations']:
    Name = str(name)
    variations.append(Name(name))

# create other objects we'll need
id_generator = IDGenerator()    # id generator - do we need this up here? Not sure...
structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher(input)
redundancy_guard = RedundancyGuard(structure_matcher, other_inputs)
development = Development(niggli_params, constraint_params)
offspring_generator = OffspringGenerator(variations, development, redundancy_guard, 10000) # don't need this until after we've finished making the initial population
energy_calculator = "energy code name" + "EnergyCalculator"(energy_code_params)
waiting_queue = deque()         # queue to hold the organisms that have finished their energy calculations and are waiting to be added to the pool
whole_pop = []                  # holds every organism that has been submitted for energy calculation, both relaxed and unrelaxed
others?

# create the initial population
initial_population = InitialPopulation()
threads = []  # list of threads to do the energy calculations
while not stopping_criteria.are_satisfied:
    for creator in organism_creators:
        while not creator.is_finished:
            # start by making N threads
            if len(threads) < N:
        
                # keep trying until we get one
                new_organism = creator.createOrganism(id_generator, composition_space, constraints) 
                while new_organism == None and not creator.is_finished:
                    new_organism = creator.create_organsim # this will cause problems since the different createOrganism methods take different numbers of arguments... 
            
                # unpad the organism (does nothing for bulk search)
                # does this make sense for random organism creator?
                geometry.unpad(new_organism, constraints)
                # develop the organism
                developed_org = development.develop(new_organism)
                if developed_org != None:
                    # check for redundancy
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org == None: # no redundancy
                        whole_pop.append(copy.deepcopy(developed_org)) # we want to add copies to whole_pop so the organisms in whole_pop don't change upon relaxation, etc.
                        stopping_criteria.updateCalcCounter()  # if num calcs is one of the stopping criteria, this updates it
                        geometry.pad(developed_org) # for bulk search, this does nothing except rotate into principal directions
                        thread = Thread(target=energy_calculator.doEnergyCalculation, args=(developed_org))
                        thread.start()
                        threads.append(thread)
            else:
                # check for dead threads
                for thread in threads:
                    if not thread.isAlive:
                        # TODO: need to figure out how to get the return value from a dead thread
                        relaxed_org = thread.return_value  
                        if relaxed_org != None:
                            geometry.unpad(relaxed_org, constraints) # for bulk search, this does nothing except rotate into principal directions
                            developed_org = development.develop(relaxed_org)
                            if developed_org != None:
                                redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                                if redundant_org != None:
                                    if redundant_org.isActive and redundant_org.value > developed_org.value:
                                        initial_population.replaceOrganism(redundant_org, developed_org)
                                else:
                                    whole_pop.append(copy.deepcopy(developed_org))
                                    stopping_criteria.checkOrganism(developed_org) # if value achieved or found structure stopping criteria are used, this checks if it's met
                                    initial_population.add(developed_org)
                                    if creator.is_successes_based: # update status of success-based creators
                                        creator.updateStatus()
                                
                        # remove the dead thread and make another one
                        threads.remove(thread)
                        
                        # keep trying until we get one
                        new_organism = creator.createOrganism(id_generator, composition_space, constraints)
                        while new_organism == None and not creator.is_finished:
                            new_organism = creator.createOrganism(id_generator, composition_space, constraints)
            
                        # develop the organism
                        developed_org = development.develop(new_organism)
                        if developed_org != None:
                            # check for redundancy
                            redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                            if redundant_org == None: # no redundancy
                                whole_pop.append(copy.deepcopy(developed_org))
                                stopping_criteria.updateCalcCounter() 
                                geometry.pad(developed_org) # for bulk search, this does nothing except rotate into principal directions
                                thread = Thread(target=energy_calculator.doEnergyCalculation, args=(developed_org))
                                thread.start()
                                threads.append(thread)   


# create the pool
pool = Pool(initial_population, other args)
threads = []  # list of threads to do the energy calculations (overwrite the old list used in initial population)

# create the initial batch of N offspring organisms and submit them for energy calculations
for i in range(0, N):
    unrelaxed_offspring = offspring_generator.makeOffspringOrganism(pool, whole_pop)
    geometry.pad(unrelaxed_offspring) 
    stopping_criteria.updateCalcCounter()
    thread = Thread(target=energy_calculator.doEnergyCalculation, args=(new_organism))
    thread.start()
    threads.append(thread)

# continue the search until stopping criteria are met
while not stopping_criteria.are_satisfied:
    # check for dead threads
    for thread in threads:
        if not thread.isAlive:
            # TODO: need to figure out how to get the return value from a dead thread
            relaxed_org = thread.return_value  
            if relaxed_org != None:
                geometry.unpad(relaxed_org, constraints) 
                developed_org = development.develop(relaxed_org)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None and redundant_org.isActive and redundant_org.value > developed_org.value:
                            pool.replaceOrganism(redundant_org, developed_org)
                    else:
                        pool.add(developed_org)
                        whole_pop.append(copy.deepcopy(developed_org))
                        stopping_criteria.checkOrganism(developed_org)
                        # check if we've added enough new offspring orgs to the pool that we can trim off the initial
                        # population orgs from the queue. poolSize and numPromoted are user inputs.
                        if (pool.numAdds == (poolSize - numPromoted)):
                            for i in range(0, len(pool.queue) - (poolSize - numPromoted)):
                                pool.queue.popleft()
                
                        # check if we've already added more than enough offspring to the pool to allow trimming off the 
                        # orgs from the initial population, in which case it has already happened, and we just need to
                        # remove one organism from the pool
                        elif (pool.numAdds > (poolSize - numPromoted)):
                            pool.queue.popleft()
                                                
            # remove the dead thread and make another one
            threads.remove(thread)
            pool.calculateFitnesses()
            pool.calculateSelectionProbs()
            offspring = offspring_generator.makeOffspringOrganism(pool, whole_pop) # this handles development and redundancy, and appending a copy to whole_pop
            stopping_criteria.updateCalcCounter()
            geometry.pad(offspring) 
            thread = Thread(target=energy_calculator.doEnergyCalculation, args=(offspring))
            thread.start()
            threads.append(thread)

'''
