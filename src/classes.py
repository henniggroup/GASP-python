from __future__ import division, unicode_literals, print_function

from pymatgen.core.structure import Structure
# from abc import abstractmethod, ABCMeta
from _pyio import __metaclass__
# import collections.deque
from _collections import deque
import random  # TODO: need to make sure all random numbers from the same PRNG
# TODO: import other needed stuff from pymatgen

'''
This module contains all the classes used by the algorithm.

TODO: should probably break this into several modules eventually.
'''

class IDGenerator(object):
    '''
    Generates successive integer ID numbers, starting from 1.
    
    This class is a singleton.
    '''
    
    def __init__(self):
        '''
        Creates an id generator.
        '''
        self.id = 0
    
    def makeID(self):
        '''
        Returns the next id number.
        '''
        self.id = self.id + 1
        return self.id



class Organism(object):
    '''
    An organism
    '''

    def __init__(self, structure, value=None, fitness=None, select_prob=None):
        '''
        Creates an organism
        
        Args:
            structure: The structure of this organism, as a pymatgen.core.structure.Structure
            
            value: The objective function value of this organism, which is either the energy
                per atom (for fixed-composition search) or the distance from the current best
                convex hull (for phase diagram search).
                
            fitness: The fitness of this organism. Ranges from 0 to 1, including both endpoints.
            
            select_prob: The selection probability of this organism. Ranges from 0 to 1, including
                both endpoints.
        '''
        # initialize instance variables
        self.structure = structure
        self.value = value
        self.fitness = fitness
        self.select_prob = select_prob 
        self._id = IDGenerator.makeID(); # unique id number for this organism. Should not be changed.
    
    # this keeps the id (sort of) immutable by causing an exception to be raised if the user tries to 
    # the set the id with org.id = some_id.
    @property
    def id(self):
        return self._id
    
    # TODO: maybe make setter methods for fitness and select_prob that check that they're between 0 and 1. 
    #    The Structure class has checks at initialization, so we don't need to check it again here.
    


class Pool(object):
    '''
    Pool to hold all the organisms that are currently candidates for selection to become parents.
    
    Composed of two parts: a promotion set containing the best few organisms, and a queue containing
        rest of the organisms in the pool.
        
    This class is a singleton.
    '''

    def __init__(self, pool_params_dict, initial_population, selection_probability_dict):
        '''
        Creates a pool of organisms
        
        Args:
            pool_params_dict: a dictionary containing the pool parameters: poolSize, numPromoted
        
            initial_population: a list of organism.Organism's that comprise the initial population. They must
                already have had their energies calculated.
            
            selection_probability_dict: a dictionary containing the parameters needed to calculating selection 
                probabilities: numParents and selectionPower
        '''
        # TODO: implement me
        # 1. calculate the fitnesses of the organisms in the initial population, based on their
        #    objective function values
        # 2. calculate the selection probabilities of the organisms, based on their fitnesses
        # 3. put the best few of them in the promotion set (exact number determined from input)
        # 4. put the rest in the queue. The order doesn't matter, since they'll all get thrown away
        #    at the same time
        # I think I can just use a list for the promotion set 
        # Use a deque for the queue
        self.promotionSet = [] # TODO: should be set to a list of the best few orgs
        self.queue = deque() # TODO: should be a deque containing the rest of the orgs
        self.selectionDist = [] # TODO: the parameters for the selection distribution
        self.numAdds = 0 # the number of organisms added to the pool (excluding the initial population)
    
    
    def addOrganism(self, org):
        '''
        Adds a new organism to the pool. 
        
        If the new organism better than one of the orgs currently in the promotion set, then it is added 
        to the promotion set, and the worst org in the promotion set is moved to the back of the queue. 
        Otherwise, the new org is just appended to the back of the queue.
        
        Args:
            org: the organism.Organism to add.
        '''
        # TODO: implement me. Look at deque methods...
        # 1. If doing a pd search, will need to transform value from epa to distance from current best convex hull
        # 2. Once value is updated (in necessary), decide whether to place in promotion set or queue, based on org values
        
        self.numAdds = self.numAdds + 1
        
    
    def replaceOrganism(self, old_org, new_org):
        '''
        Replaces an organism in the pool with a new organism. The new organism has the same location in the pool
        as the old one.
        
        Precondition: the old_org is a member of the current pool.
        
        Args:
            old_org: the organism in the pool to replace
            new_org: the new organism to replace the old one
        '''
        # TODO: implement me
        # 1. determine if old_org is in promotion set or queue
        # 2. do the replacement
        # 3. if the new_org is either the best or worst in the pool, will need to update fitnesses and selection probs
    
    
    def calculateFitnesses(self):
        '''
        Calculates and assigns fitnesses to all organisms in the pool.
        
        Precondition: the organisms in the pool all have valid values
        '''
        # TODO: implement me. 
        # There might be some tricks to speed this up, like:
        #    always keeping track of the best and worst values in the pool, 
        #    so we don't have to search for them each time
        #
        #    might only have to update the fitness of the newest addition to the pool, 
        #    if the best and worst values didn't change when it was added.
        
    
    def calculateSelectionProbs(self):
        '''
        Calculates and assigns selection probabilities to all the organisms in the pool.
        
        Precondition: the organisms in the pool all have valid values.
        '''
        # TODO: implement me
        # some of the same tricks as in calculateFitnesses possible here too
        
        
    def toList(self):
        '''
        Returns a list containing all the organisms in the pool.
        '''
        # TODO: implement me
        


class ParametersParser(object):
    '''
    A parser for the input parameters.
    
    This class is a singleton.
    '''
    
    def __init__(self, path_to_input_file):
        '''
        Creates a parser for the input parameters
        
        Args:
            path_to_input_file: the path to the file containing the search parameters
        '''
        self._input_file = path_to_input_file
    
    def parse(self):
        '''
        Parses the input file
        
        Returns a dictionary of dictionaries containing all the input parameters. The dictionary 
            names correspond to the section headers (e.g. Constraints, Variations). Note that some 
            sections will have multiple subsections, for example the Variations sections will have 
            Mutation, Mating, etc. subsections. For cases like these, the keys of the Variations dict 
            will be the subsection names, and the values will themselves be dicts containing the data 
            for the given subsection.
        '''
        # TODO: implement me
        #
        # The exact details of how the parsing is done will depend on how we decide to format the input file,
        # but the basic idea is:
        #
        #    1. Read the input file into a string or list of strings
        #    2. Make a dictionary to hold everything
        #    3. For each section header:
        #        create a key in the master dictionary with that name
        #        create a dictionary containing the parameters and make it the value for the key
        #        do this recursively if a subsection header is encountered instead of a list of data
        #    4. Return the master dictionary
        
        

class Variation(object):
    '''
    A general variation object for creating offspring organisms from parent organism(s)
    
    Not meant to be instantiated, but rather subclassed by particular Variations, like Mating,
    Mutation, etc.
    '''
    
    def doVariation(self):
        '''
        Creates an offspring organism from parent organism(s).
        
        Returns an organism.
        '''
        raise NotImplementedError("Please implement this method.")
    
    
    def getSelectionProb(self, variation_parameters):
        '''
        Returns the selection probability of this variation, as a float between 0 and 1
        
        Args:
            variation_parameters: a dictionary containing the parameters for the variation
        '''
        return variation_parameters["selectProb"]
    
    
    def selectParents(self, n, pool):
        '''
        Selects n distinct organisms from the pool.
        
        Returns a list containing n organisms.
        
        Args:
            n: how many organisms to select from the pool
            pool: the current pool of organisms
            
        Precondition: all the organisms in the pool have been assigned selection probabilities.
        '''
        # TODO: implement me



class Mating(Variation):
    '''
    A mating operator.
    
    This class is a singleton.
    '''

    def __init__(self, mating_params):
        '''
        Creates a mating operator
        
        Args:
            
            mating_params: The parameters for doing the mating operation, as a dictionary.
                
        '''
        # TODO: initialize the instance variables using the values in the dict. Maybe just keeping a 
        # copy of the dict as an instance variable is enough...
    
    def doVariation(self):
        '''
        Performs the mating operation, as described in ref. TODO
        
        Returns the resulting offspring as an organism.Organism
        '''
        # TODO: implement me
        # 1. select two parent organisms from the pool - selectParents(2, pool)
        # 2. combine the two parents to make an offspring structure, and return a new offspring organism 
        #    with that structure
        
        
        
class Mutation(Variation):
    '''
    A mutation operator.
    
    This class is a singleton.
    '''
    
    def __init__(self, mutation_params):
        '''
        Creates a mutation operator
        
        Args:
        
            mutation_params: The parameters for doing the mutation operation, as a dict.
        '''   
        # TODO: initialize the instance variables using the values in the dict. Maybe just keeping a
        #    copy of the dict as an instance variable is enough...
        
    
    def doVariation(self):
        '''
        Performs the mutation operation, as described in ref. TODO
        
        Returns the resulting offspring as an organism.Organism
        ''' 
        # TODO: implement me
        # 1. select a parent organism from the pool - selectParents(1, pool)
        # 2. do the mutation to make an offspring structure, and return a new offspring organism with 
        #    that structure
   
        

class Permutation(Variation):
    '''
    A permutation operator.
    
    This class is a singleton.
    '''
    
    def __init__(self, permutation_params):
        '''
        Creates a permutation operator
        
        Args:
        
        permutation_params: The parameters for doing the permutation operation, as a dict.
        '''
        # TODO: initialize the instance variables using the values in the dict. Maybe just keeping a
        #    copy of the dict as an instance variable is enough...
    
    
    def doVariation(self):
        '''
        Performs the permutation operation, as described in ref. TODO
        
        Returns the resulting offspring as an organism.Organism
        '''
        # TODO: implement me
        # 1. select a parent organism from the pool - selectParents(1, pool)
        # 2. do the permutation to make an offspring structure, and return an offspring organism with 
        #    that structure
        
        

class NumStoichsMut(Variation):
    '''
    An operator that creates an offspring organism by mutating the number of stoichiometries' worth of atoms 
    in the parent organism.
    
    This class is a singleton.
    '''
    
    def __init__(self, numstoichsmut_params):
        '''
        Creates a NumStoichsMut operator
        
        Args:
        
        numstoichsmut_params: The parameters for doing the numstoichsmut operation, as a dict.
        '''
        # TODO: initialize the instance variables using the values in the dict. Maybe just keeping a
        #    copy of the dict as an instance variable is enough...
    
    def doVariation(self, parent):
        '''
        Performs the numstoichsmut operation, as described in ref. TODO
        
        Returns the resulting offspring as an organism.Organism
        '''
        # TODO: implement me
        # 1. select a parent organism from the pool - selectParents(1, pool)
        # 2. do the numstoichsmutation to make an offspring structure, and return an offspring organism 
        #    with that structure
        
        
        
class OrganismCreator(object):
    '''
    Creates organisms for the initial population
    
    Not meant to be instantiated, but rather subclassed by particular Creators, like RandomOrganismCreator
    or PoscarsOrganismCreator.
    '''
    
    def create_organisms(self):
        '''
        Creates the organisms for the initial population.
        
        Returns a list of organisms whose energies have been calculated.
        
        Args:
            TODO: think about what data this method needs (e.g. objective function data, etc.) and add it 
                to the argument list.
        '''
        raise NotImplementedError("Please implement this method.")
        #
        # The general outline of this method (to be implemented in subclasses) is:
        #     1. Get the information needed to create the organisms (path to poscars, random vol, 
        #        stoichiometry, etc.)
        #     2. Get the information needed to submit an energy calc (energy code, objective function, etc.)
        #     3. For each organism:
        #            - create the organism
        #            - check against constraints
        #            - if it passes, add it to wholePop list and submit for energy evaluation
        #            - if not, throw it away
        


class RandomOrganismCreator(OrganismCreator):
    '''
    Creates random organisms for the initial population
    '''
    def __init__(self, needed_parameters):
        '''
        Creates a RandomOrganismCreator.
        
        Args:
            needed_parameters: all the parameters needed for creating the random organisms.
            This includes how many to make and whether to scale their volumes (and if yes, to what value).
            It will also need stoichiometry information so it knows what types of atoms to use.
        '''
    
    def create_organisms(self):
        '''
        Creates random organisms for the initial population.
        
        Returns a list of organisms whose energies have been calculated.
        
        Args:
        
            TODO: think about what data this method needs (e.g. objective function data, etc.) and add it 
                to the argument list.
        '''
        # This method will need to create random offspring organisms using the methods specified in 
        # needed_parameters and submit them in batches of N at a time.
        #
        # Will need to have methods to create random structures and submit organisms for energy calculations.
        #
        # Some of the needed data will be N, the objective function data, and data concerning the external code
        #
        # Once an organism is created and it passes the constraints, it should be added to the wholePop list
        # before relaxation. After relaxation, it should be added to the wholePop list again
        


class PoscarOrganismCreator(OrganismCreator):
    '''
    Creates organisms from poscar files for the initial population.
    '''
    def __init__(self, needed_parameters):
        '''
        Creates a RandomOrganismCreator.
        
        Args:
            needed_parameters: all the parameters needed for creating the random organisms.
            This includes how many to make and whether to scale their volumes (and if yes, to what value).
            It will also need stoichiometry information so it knows what types of atoms to use.
        '''
    
    def create_organisms(self):
        '''
        Creates organisms for the initial population from poscar files.
        
        Returns a list of organisms whose energies have been calculated.
        
        Args:
        
            TODO: think about what data this method needs (e.g. objective function data, etc.) and add it 
                to the argument list.
        '''
        # This method will need to create organisms from poscar files and submit them in batches of N at a time.
        #
        # Will need to have methods to read in a structure from a poscar file and submit organisms for energy 
        # calculations.
        #
        # Some of the needed data will be N, the objective function data, and data concerning the external code
        #
        # Once an organism is created and it passes the constraints, it should be added to the wholePop list
        # before relaxation. After relaxation, it should be added to the wholePop list again



class RedundancyGuard(object):
    '''
    A redundancy guard.
    
    This is a singleton class.
    '''
    
    def __init__(self, structure_matcher, other_params):
        '''
        Creates a redundancy guard.
        
        Args:
            structure_matcher: a StructureMatcher object for comparing organisms' structures.
            
            other_params: other optional parameters for identifying redundant structures, like
                d-value, etc.
        '''
        # TODO: implement me. Maybe just keeping copies of the relevant dictionaries as instance variable
        #    is enough...
    
    
    def checkStructures(self, org, list_of_organisms):
        '''
        Checks if an organism's structure is redundant with that of any organism in a list of organisms,
        according to the tolerances in the StructureMatcher.
        
        Returns the organism that org is redundant with. If it's not redundant, returns None
        
        Args:
            org: organism to check for redundancy.
            
            list_of_organsisms: list of organisms against which to check for redundancy.
        ''' 
        # TODO: implement me. Use StructureMatcher
    
    
    def checkPool(self, org, pool):
        '''
        Checks if an organism is redundant with with any organism in the pool. This includes 
        both structure comparisons and also d-value comparisons, if specified.
        
        Returns the organism from the pool that org is redundant with. If it's not redundant, returns None.
        
        Args:
            org: organism to check for redundancy
            
            pool: the current pool
        '''   
        # TODO: implement me.
        #
        # start by checking just the structures - call checkStructures(org, pool.toList())
        # if a duplicate is found AND the new structure has a lower objective function value,
        #    do a pool.replaceOrganism() and exit the method
        # 2. if not a redundant structure, then check the d-value constraint, if it's specified. 
        #        if a duplicate is found with this criteria, do a pool.replaceOrganism() and exit
        #        the method
        


class Development(object):
    '''
    A development object is used to develop an organism before evaluating its energy or adding it
    to the pool. Doesn't do redundancy checking.
    
    This is a singleton class.
    '''
    
    def __init__(self, niggli_reduction_params, structure_constraints_params):
        '''
        Creates a Development object.
        
        Args:
            niggli_reduction_params: dictionary containing data for Niggli cell reduction
                        
            structure_constraints: dictionary of structure constraints parameters
        '''
        # TODO: implement me
        # I think I can just keep the two dictionaries as instance variables
        
    
    def develop(self, organism):
        '''
        Develops an organism.
        
        Returns the developed organism, or None if the organism failed development
        
        Args:
            organism: the organism to develop
        '''
        # TODO: implement me
        # 
        # 1. Do Niggli cell reduction, if specified
        # 2. Check the structural constraints. Print error message if it failed


class OffspringGenerator(object):
    '''
    Used to generate offspring structures
    '''
    
    def __init__(self, variations, development, redundancy_guard, num_tries_limit):
        '''
        Args:
            variations: a list of Variation objects 
            
            development: the Development object (for cell reduction and structure constraints)
            
            redundancy_guard: the redundancyGuard object 
            
            num_tries_limit: the max number of times to attempt creating an offspring organism from a given variation
                before giving up and trying a different variation.
        '''
        self.variation = variations
        self.development = development
        self.redundancy_guard = redundancy_guard
        self.num_tries_limit = num_tries_limit
        
        
    def makeOffspringOrganism(self, pool, whole_pop):
        '''
        Generates a valid offspring organism using the variations.
        
        Returns an unrelaxed offspring organism.
        
        Args:
            pool: the current Pool
            
            whole_pop: the list containing all the orgs seen so far (both relaxed an unrelaxed)
        '''
        tried_variations = []
        while(len(self.variations) > len(tried_variations)):
            variation = self.selectVariation(tried_variations)
            num_tries = 0
            while (num_tries < self.num_tries_limit):
                offspring = variation.doVariation()
                if (self.development.develop(offspring) != None) and (self.redundancy_guard.checkStructures(offspring, whole_pop) == None):
                    return offspring
                else:
                    num_tries = num_tries + 1
            tried_variations.append(variation)
        print("Could not make valid offspring organism with any Variation.") # This point should never be reached
        
    
    def selectVariation(self, tried_variations):
        '''
        Selects a variation that hasn't been tried yet based on their selection probabilities
        
        Args:
            tried_variations: list of Variations that have already been unsuccessfully tried
        '''
    # TODO: implement me
    #
    #    while(true):
    #        variation = random_selection(variations)    # choose one randomly based on their selection probs.
    #        if (variation not in tried_variations):
    #            return variation
    
        
    
class EnergyCalculator(object):
    '''
    Handles calculating the energy of organisms.
    
    Not meant to be instantiated, but rather subclassed by particular Calculators, like VaspEnergyCalculator
    or GulpEnergyCalculator.
    '''
    
    def doEnergyCalculation(self, org):
        '''
        Calculates the energy of an organism.
        
        Args:
            org: the organism whose energy we want to calculate
        '''
        raise NotImplementedError("Please implement this method.")
        # TODO: think about how to impelement this. There are two main parts: preparing for the calculation (writing
        #    input files), and actually submitting it, by calling an external script or something. Once the calculation
        #    is finished, we need to develop the relaxed organism, add it to the whole_pop list, and add it to the
        #    waiting_queue.
        #
        #    All this should be done on it's own thread, so that multiple energy calculations can be running at once.
        #    It might be best to handle the threads inside the the method (not sure)
        #
        #    The goal is be able to call EnergyCalculator.doEnergyCalcualtion(unrelaxed_organism) and then have the
        #    control flow immediately return (i.e. not having to wait for the method to finish)
        #
        #    Note: when the energy calculation finishes, this method will need to have access to the current versions
        #        of the whole_pop list and the waiting_queue. Not sure best way to do that...
        
        # set up the energy calc (prepare input files, etc.)
        # do the energy calc
        # if it finished correctly: 
        #    relaxed_org = development.develop(relaxed_org) it and append updated org to the waiting queue
        #    if (relaxed_org != None):
        #        waiting_queue.append(relaxed_org)
        #    else:
        #        print("failed constraints") # say which one
        # else:
        #    print("Energy calculation failed. Discarding org") 
        
        

class VaspEnergyCalculator(EnergyCalculator):
    '''
    Calculates the energy of an organism using VASP.
    '''
    
    def __init__(self, vasp_code_params):
        '''
        Args:
            vasp_code_params: the parameters needed for preparing a vasp calculation (INCAR, KPOINTS, POTCAR files)
        '''
        # TODO: implement me. Just keeping the paths to the input files should be enough
    
    
    def doEnergyCalculation(self, org):
         '''
        Calculates the energy of an organism using VASP.
        
        Args:
            org: the organism whose energy we want to calculate
        '''
        # TODO: implement me
        # 1. prepare input files for calculation
        # 2. submit calculation, by running external script
        # 3. when finished
        #        - develop relaxed organism
        #        - add relaxed organism to whole_pop list
        #        - add relaxed organism to waiting_queue
    

        
        
        
        
########## area for casual testing ########## 

# make a structure object
lattice = [[5,0,0],[0,5,0],[0,0,5]]
species = ["C", "Si"]
coordinates = [[0.25,0.25,0.25],[0.75,0.75,0.75]]
structure1 = Structure(lattice, species, coordinates)

# make an organism object
org1 = Organism(structure1)

print(org1.structure)
print(org1.fitness)
print(org1.id)

org1.id = 6
print(org1.id)

########## end of testing area ##########    


'''
One way to store the parameters is in several groups, one for each related group of parameters

Possible groups include

    variations, which has four subgroups: mutation, mating, permutation, numstoichsmut
    
    initial population settings
    
    redundancy guard settings
    
    termination criteria
    
    objective function info
    
    things needed by the algorithm at all times, like
        - number of calcs to run concurrently
        - pool size
        - the selection probability function
        - volume scaling
        - niggli cell reduction


I think where possible, the data should be stored inside the objects that need it, like the variation objects.

Idea: read in the input file with minimal processing (just a string or something), then have each object parse 
the data it needs out of that string when it is initialized. Some of it might just go into lists or something.

Maybe use dictionaries to store the data

Ok, start by making a list of all the input file options, and think about how to logically divide that into 
dictionaries. Let's worry about parsing from the input file later, including what format it should be in.

Optional flags [Current value]:
   --help : displays this message with parameters' default values
   --verbosity n : verbosity of output [4]
Genetic Algorithm
   --t runTitle : name the run
   --outDir : specify the output directory
   --keepTempFiles <true|false>
   --saveStateEachIter <true|false>
   --popSize <n> : use a non-initial population size of n
   --promotion <n> : promote n best structures (or the whole CH) to next gen
   --compositionSpace <numElements> <Sym>*  :     System composition for phase diagram search
                   or <numElements> <Sym>* <amount of element1> <amount of element2> etc.
   --optimizeDensity <weight of adaptation> <num orgs to avg. over>
   --useRedundancyGuard <wholePopulation|perGeneration|both> <atomic misfit> <lattice misfit> <angle misfit> <use PBCs?>
   --useSurrogateModel <gulp header file> <potentials_specification file> <refit frequency>
   --endgameNumGens <n>
   --useNiggliReducedCell <true|false>
   --use2DNiggliReducedCell <true|false>
   --useSubstrate <true|false> <path to substrate POSCAR file
   --writeHartkeFile <boolean>
   --colorOutput <boolean>
   
Initial Population
   --initialPopulation <num> random givenVol <volumeperatom>
   --initialPopulation <num> random randomVol
   --initialPopulation <num> poscars <directory>
   --initialPopulation <num> manual
   --initialPopulation <num> units <numMols> <numAtoms_1>...<numAtoms_n> (<symbol_i> <x_i> <y_i> <z_i>)+ <numUnits_1>...<numUnits_n> <targetDensity> <densityTol> <unitsOnly?>
   --initialPopulation <num> supercell a b c maxll minll maxla minla maxh maxna minna <randomsocreator args>
   
Objective Functions
   --objectiveFunction cluster <padding length> <other obj fcn args from below...>
   --objectiveFunction surface <padding length> <other obj fcn args from below...>
   --objectiveFunction substrate <padding length> <other obj fcn args from below...>
   --objectiveFunction <epa/pd> gulp <gulp header file> <gulp potential file> <cautious?> <species needing a shell>
   --objectiveFunction <epa/pd> vasp <cautious?> <kpoints> <incar> <element potcar>+ 
   --objectiveFunction <epa/pd> ohmms <header> <footer> <cautious?>
   --objectiveFunction <epa/pd> lammps <potlFile> <units> <relax box?>
   --objectiveFunction <epa/pd> castep <cautious?> <kpointSpacing> <pressure> <paramFile> <element potcar>+ 
   --objectiveFunction <epa/pd> avogadro <avog header file>
   --objectiveFunction <epa/pd> dlpoly <loc> <potl>
   --objectiveFunction <epa/pd> mopac <execpath>
   --objectiveFunction <epa/pd> dftpp <dftpp_inputs> <cautious?> <element ppFile.fhi>*
   --objectiveFunction <epa/pd> generic
   --parallelize <numCalcsInParallel> <minPopSize>
   
Variation Algorithms
   --variation <percentage> <percentage> slicer <thicknessMean> <thicknessSigma> <majorShiftFrac> <minorShiftFrac> <maxAmplitude> <maxFreq> <growParents?> <doublingProb>
   --variation <percentage> <percentage> structureMut <rate> <sigmaAtoms> <sigmaLattice>
   --variation <percentage> <percentage> permutation <meanSwaps> <sigmaSwaps> <pairsToSwap (e.g. Mn-O)>
   --variation <percentage> <percentage> numStoichsMut <meanNumAtoms> <sigmaNumAtoms>
   
Selection Algorithms
   --selection probDist <numParents> <selectionPower>
   
Convergence Criteria
   --convergenceCriterion maxFunctionEvals <n>
   --convergenceCriterion maxNumGens <n>
   --convergenceCriterion maxNumGensWOImpr <n> <dValue>
   --convergenceCriterion valueAchieved <maximum acceptable energy>
   --convergenceCriterion foundStructure <CIF filename> <rGuard misfits>
   
Hard Constraints
   --minInteratomicDistance d : minimum interatomic distance (Angstroms)
   --perSpeciesMID <symbol1> <symbol2> <distance>
   --maxLatticeLength d : maximum lattice vector length (Angstroms)
   --minLatticeLength d : minimum lattice vector length (Angstroms)
   --maxLatticeAngle d : maximum lattice angle (Degrees)
   --minLatticeAngle d : minimum lattice angle (Degrees)
   --maxCellHeight d : maximum height of cell in z-direction
   --maxNumAtoms n
   --minNumAtoms n
   --minNumSpecies n
   --doNonnegativityConstraint <boolean>
   --dValue x : discard organisms within a value of x of each other
   
   
   
   
   Variations:         # section
       Mating:         # subsection
           percentage
           thicknessMean
           thicknessSigma
           majorShiftFrac
           minorShiftFrac
           growParents
           doublingProb
        Mutation:
            percentage
            fracPerturbed
            sigmaAtoms
            sigmaLattice
        Permutation:
            percentage
            meanSwaps
            sigmaSwaps
            pairsToSwap (probably a list of pairs or something)
        NumStoichsMut:
            percentage
            meanNumAtoms
            sigmaNumAtoms
        
    Selection:
        numParents
        selectionPower
        
    Constraints:
        minInteratomicDistances (probably a dict of pairs and distances or something)
        maxLatticeLength
        minLatticeLength
        maxLatticeAngle
        minLatticeAngle
        maxLayerThickness
        maxNumAtoms
        minNumAtoms
        minNumSpecies
        doNonnegativityConstraint (maybe don't need this?)
        dValue
        
    Convergence Criteria:
        maxFunctionEvals
   
   Objective Function:
       objectiveFunction <epa|pd>
       geometry <bulk|sheet|cluster>
       padding (ignored if geometry is bulk)
       energyCode <vasp|gulp|more stuff>
       energyFiles (depends on energyCode. Not sure best way to handle this)
        
        
    
       

'''
