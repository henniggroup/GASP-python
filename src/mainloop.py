# Outline of the algorithm
import yaml
import sys
import os
import classes
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.periodic_table import Element

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
    print("Input file must contain a CompositionSpace block.")
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
                print('The path to the folder containing the files must be provided. Please give the path after the "path_to_folder" flag.')
                print('Quitting...')
                quit()
            elif not os.path.exists(given_path):
                print('The given folder does not exist.')
                print('Quitting...')
                quit()
            # if the folder exists, check that it contains files
            elif len([f for f in os.listdir(given_path) if os.path.isfile(os.path.join(given_path, f))]) == 0:
                print('The given folder does not contain any files.')
                print('Quitting...')
                quit()
            # otherwise, the directory exists and contains files. Instantiate the files organism creator
            else:
                files_organism_creator = classes.FileOrganismCreator(given_path)
                organism_creators.append(files_organism_creator)
        # in case an incorrect flag is given after the from_files flag 
        else:
            print('Incorrect flag given after "from_files". Please use the "path_to_folder" flag.')
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







# just for testing

foc = organism_creators[0]
roc = organism_creators[1]

print(foc.files)
print(roc.number)

file_org = foc.createOrganism()
for i in range(len(foc.files)):
    if file_org == None:
        print('pymatgen couldnt read the structure in')
    else:
        print(file_org.structure)
        print("")











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
variations = []
for name in dict_of_parameters['Variations']:
    Name = str(name)
    variations.append(Name(name))

# create other objects we'll need
id_generator = IDGenerator()    # id generator - do we need this up here? Not sure...
structure_matcher = pymatgen.analysis.structure_matcher.StructureMatcher(input)
redundancy_guard = RedundancyGuard(structure_matcher, other_inputs)
development = Development(niggli_params, constraint_params)
offspring_generator = OffspringGenerator(variations, development, redundancy_guard, 10000)
energy_calculator = "energy code name" + "EnergyCalculator"(energy_code_params)
waiting_queue = deque()         # queue to hold the organisms that have finished their energy calculations and are waiting to be added to the pool
whole_pop = []                  # holds every organism that has been submitted for energy calculation, both relaxed and unrelaxed
others?

# create the initial population
initial_population = InitialPopulation(whole_pop)
threads = []  # list of threads to do the energy calculations
for creator in organism_creators:
    while not creator.is_finished:
        # start by making N threads
        if len(threads) < N:
        
            # keep trying until we get one
            new_organism = creator.create_organism() 
            while new_organism == None and not creator.is_finished:
                new_organism = creator.create_organsim 
            
            # develop the organism
            developed_org = development.develop(new_organism)
            if developed_org != None:
                # check for redundancy
                redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                if redundant_org == None: # no redundancy
                    whole_pop.append(developed_org)
                    developed_org.structure = geometry.pad(developed_org.structure) # for bulk search, this does nothing except rotate into principal directions
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
                        relaxed_org.structure = geometry.unpad(relaxed_org.structure) # for bulk search, this does nothing except rotate into principal directions
                        developed_org = development.develop(relaxed_org)
                        if developed_org != None:
                            redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                            if redundant_org != None:
                                if redundant_org.isActive and redundant_org.value > developed_org.value:
                                    initial_population.replaceOrganism(redundant_org, developed_org)
                            else:
                                whole_pop.append(developed_org)
                                initial_population.add(developed_org)
                                if creator.when_stop == "successes": # update status of success-based creators
                                    creator.updateStatus()
                                
                    # remove the dead thread and make another one
                    threads.remove(thread)
                        
                    # keep trying until we get one
                    new_organism = creator.create_organism() 
                    while new_organism == None and not creator.is_finished:
                        new_organism = creator.create_organsim 
            
                    # develop the organism
                    developed_org = development.develop(new_organism)
                    if developed_org != None:
                        # check for redundancy
                        redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                        if redundant_org == None: # no redundancy
                            whole_pop.append(developed_org)
                            developed_org.structure = geometry.pad(developed_org.structure) # for bulk search, this does nothing except rotate into principal directions
                            thread = Thread(target=energy_calculator.doEnergyCalculation, args=(developed_org))
                            thread.start()
                            threads.append(thread)   


# create the pool
pool = Pool(initial_population, other args)
threads = []  # list of threads to do the energy calculations (overwrite the old list used in initial population)

# create the initial batch of N offspring organisms and submit them for energy calculations
for i in range(0, N):
    unrelaxed_offspring = offspring_generator.makeOffspringOrganism(pool, whole_pop)
    unrelaxed_offspring.structure = geometry.pad(unrelaxed_offspring.structure) 
    thread = Thread(target=energy_calculator.doEnergyCalculation, args=(new_organism))
    thread.start()
    threads.append(thread)

# continue the search until stopping criteria are met
while the stopping criteria are not met:
    # check for dead threads
    for thread in threads:
        if not thread.isAlive:
            # TODO: need to figure out how to get the return value from a dead thread
            relaxed_org = thread.return_value  
            if relaxed_org != None:
                relaxed_org.structure = geometry.unpad(relaxed_org.structure) 
                developed_org = development.develop(relaxed_org)
                if developed_org != None:
                    redundant_org = redundancy_guard.checkRedundancy(developed_org, whole_pop)
                    if redundant_org != None and redundant_org.isActive and redundant_org.value > developed_org.value:
                            pool.replaceOrganism(redundant_org, developed_org)
                    else:
                        pool.add(developed_org)
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
            offspring = offspring_generator.makeOffspringOrganism(pool, whole_pop) # this handles development and redundancy
            offspring.structure = geometry.pad(offspring.structure) 
            thread = Thread(target=energy_calculator.doEnergyCalculation, args=(offspring))
            thread.start()
            threads.append(thread)

'''
