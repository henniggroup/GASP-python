'''
Outline of the algorithm


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

# create variation objects (Mutation, Mating, etc.) and put them in a list
variations = []
for name in dict_of_parameters['Variations']:
    Name = str(name)
    variations.append(Name(name))

# create other objects we'll need
id_generator = IDGenerator()    # id generator - do we need this up here? Not sure...
structure_matcher = StructureMatcher(input)
redundancy_guard = RedundancyGuard(structure_matcher, other_inputs)
development = Development(niggli_params, constraint_params)
offspring_generator = OffspringGenerator(variations, development, redundancy_guard, 10000)
energy_calculator = "energy code name" + "EnergyCalulator"(energy_code_params)
waiting_queue = deque()         # queue to hold the organisms that have finished their energy calculations and are
                                # waiting to be added to the pool
whole_pop = []                  # holds every org seen, both relaxed and unrelaxed, to all redundancy checking
others?

# populate the initial population
# TODO: unpack/expand this a little more. We want to move some of the detail of the create_organisms
#     method up to this level. Think about how to wait after submitting N jobs, how to do constraints 
#     and redundancy checks, and how to keep a wholePop list of all the structures, both relaxed and 
#     unrelaxed.
initial_population = []
for creator in organism_creators:
    initial_population.extend(creator.create_organisms())

# create the pool
pool = Pool(initial_population, other args)

# create the initial batch of N offspring organisms and submit them for energy calculations
for i in range(0, N):
    unrelaxed_offspring = offspring_generator.makeOffspringOrganism(pool, whole_pop)
    whole_pop.append(unrelaxed_offspring)
    energy_calculator.doEnergyCalculation(unrelaxed_offspring) # this handles threads, development, whole_pop and waiting q  

# continue the search until stopping criteria are met
while the stopping criteria are not met:

    # fetch next offspring in line
    next_org = waiting_queue.popleft()  # or something similar
    
    # check if redundant with the pool 
    redundant_pool_org = redundancy_guard.checkPool(next_org, pool)
    if (redundant_pool_org != None):
        pool.replaceOrganism(redundant_pool_org, next_org) # TODO: only swap if new on is better...
    else:
        # add it to the pool
        pool.addOrganism(offspring)
        
        # check if we've added enough new offspring orgs to the pool that we can trim off the initial
        # population orgs from the queue. poolSize and numPromoted are user inputs.
        if (pool.numAdds == (poolSize - numPromoted)):
            for i in range(0, len(pool.queue) - (poolSize - numPromoted)):
                pool.queue.popleft()
                
        # check if we've already added more than enough offspring to the pool to allow trimming off the 
        # orgs from the initial population, in which case it has already happened.
        elif (pool.numAdds > (poolSize - numPromoted)):
            pool.queue.popleft()
            
        # now create a new organism from the pool and submit it for an energy calculation
        pool.calculateFitnesses()
        pool.calculateSelectionProbs()
        offspring = offspring_generator.makeOffspringOorganism(pool, whole_pop) # this handles development and redundancy
        energy_calculator.doEnergyCalculation(offspring)                        # use threads

'''
