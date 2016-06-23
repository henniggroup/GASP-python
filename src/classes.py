from __future__ import division, unicode_literals, print_function

from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.transformations.standard_transformations import RotationTransformation
from pymatgen.core.operations import SymmOp
# from abc import abstractmethod, ABCMeta
from _pyio import __metaclass__
# import collections.deque
from _collections import deque
import random  # TODO: need to make sure all random numbers from the same PRNG
import threading
import numpy as np
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

    def __init__(self, structure, value=None, fitness=None, select_prob=None, isActive=False):
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
                
            isActive: Whether this organism is currently part of the pool or initial population
        '''
        # initialize instance variables
        self.structure = structure
        self.value = value
        self.fitness = fitness
        self.select_prob = select_prob 
        self.isActive = isActive
      #  self._id = IDGenerator.makeID(); # unique id number for this organism. Should not be changed.
    
    # this keeps the id (sort of) immutable by causing an exception to be raised if the user tries to 
    # the set the id with org.id = some_id.
    @property
    def id(self):
        return self._id
    
    # TODO: maybe make setter methods for fitness and select_prob that check that they're between 0 and 1. 
    #    The Structure class has checks at initialization, so we don't need to check it again here.
    
    def rotateToPrincipalDirections(self):
        '''
        Rotates the organism's structure into the principal directions. That is, a is parallel to the Cartesian x-axis 
        and b lies in the Cartesian x-y plane.
        
        '''
        # 1. rotate about the z-axis until a is vertically aligned with the x-axis
        # 2. rotate about y-axis until a is parallel to the x-axis
        # 3. rotate about x-axis until b lies in the x-y plane
        # 4. make sure all vectors are pointing in positive directions
        
        # rotate about the z-axis
        rotation = RotationTransformation([0, 0, 1], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[0][1], self.structure.lattice.matrix[0][0]))
        self.structure = rotation.apply_transformation(self.structure)
        # rotate about the y-axis
        rotation = RotationTransformation([0, 1, 0], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[0][2], self.structure.lattice.matrix[0][0]))
        self.structure = rotation.apply_transformation(self.structure)
        # rotate about the x-axis
        rotation = RotationTransformation([1, 0, 0], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[1][2], self.structure.lattice.matrix[1][1]))
        self.structure = rotation.apply_transformation(self.structure)
        # make sure they are all pointing in positive directions
        if self.structure.lattice.matrix[0][0] < 0:
            # rotate about y-axis to make a positive
            rotation = RotationTransformation([0, 1, 0], 180)
            self.structure = rotation.apply_transformation(self.structure)
        if self.structure.lattice.matrix[1][1] < 0:
            # rotate about x-axis to make b positive
            rotation = RotationTransformation([1, 0, 0], 180)
            self.structure = rotation.apply_transformation(self.structure)
        if self.structure.lattice.matrix[2][2] < 0:
            # mirror c across the x-y plane to make it positive
            reflection = SymmOp([1, 0, 0], [0, 0, 0])
            # TODO: figure out how to apply this SymmOp object to do the reflection
            
            
            #print(self.structure.lattice.matrix[2][2])
            #x_component = self.structure.lattice.matrix[2][0]
            #y_component = self.structure.lattice.matrix[2][1]
            #z_component = self.structure.lattice.matrix[2][2]
            
            #self.structure.lattice = [self.structure.lattice[0], self.structure.lattice[1], [x_component, y_component, -z_component]]
            
            #self.structure.lattice.matrix[2][2] = -1*self.structure.lattice.matrix[2][2]
            #print(self.structure.lattice.matrix[2][2])
            


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
        Adds a new organism to the pool, and also to whole_pop.
        
        If the new organism better than one of the orgs currently in the promotion set, then it is added 
        to the promotion set, and the worst org in the promotion set is moved to the back of the queue. 
        Otherwise, the new org is just appended to the back of the queue.
        
        Args:
            org: the organism.Organism to add.
        '''
        # TODO: implement me. Look at deque methods...
        # 1. If doing a pd search, will need to transform value from epa to distance from current best convex hull
        # 2. Once value is updated (if necessary), decide whether to place in promotion set or queue, based on org values
        # 3. Add organism to whole_pop (whole_pop.append(org))
        # 4. Set org.isActive = True
        
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
        # 3. set old_org.isActive = False and newOrg.isActive = True
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
        
        
        
class Geometry(object):
    '''
    Represents the geometry data.
    
    This class is a singleton
    '''
    
    def __init__(self, geometry_parameters):
        '''
        Creates a geometry object
        
        Args:
            geometry_parameters: a dictionary of parameters
        '''
        # default values
        self.default_shape = 'bulk'
        self.default_padding = 10 # this will only be used for non-bulk shapes
        
        # if entire "Geometry block was set to default or left blank"
        self.geometry_parameters = geometry_parameters
        self.geometry_parameters = geometry_parameters
        if self.geometry_parameters == None or self.geometry_parameters == 'default':
            self.shape = self.default_shape
            self.padding = None
        else:
            # check each one and see if it's been left blank or set to default
            self.shape = geometry_parameters['shape']
            if self.shape == None or self.shape == 'default':
                self.shape = self.default_shape
                self.padding = None
            elif self.shape != 'bulk':
                self.padding = geometry_parameters['padding']
                if self.padding == None or self.padding == 'default':
                    self.padding = self.default_padding
                
    def pad(self, structure):
        '''
        Makes an organism's structure conform to the required shape. For sheet, wire and cluster geometries, this 
        means adding vacuum padding to the cell. For bulk, the structure is unchanged. Used to pad a structure prior to 
        an energy calculation.
        
        Returns a structure that has been modified to conform to the shape (most likely padded with vacuum).
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # rotate structure into principal directions
        # Call other methods based on the value of self.shape
        if self.shape == 'sheet':
            return self.padSheet(structure)
        elif self.shape == 'wire':
            return self.padWire(structure)
        elif self.shape == 'cluster':
            return self.padCluster(structure)
        else:
            return structure
        
    def padSheet(self, structure):
        '''
        Adds vertical vacuum padding to a sheet, and makes the c-lattice vector normal to the plane of the sheet.
        
        Returns a sheet structure that has been padded with vacuum.
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
        # 1. rotate structure to the principal directions
        #         TODO: see if pymatgen has a method for this
        # 2. replace c with it's vertical component (make it normal to plane of the sheet) (make sure atomic positions are preserved)
        
        # 3. reduce c such that it's equal to the layer thickness (max vertical distance between atoms in the cell)
        # 4. add padding to c
        # 5. return padded structure
        
    def padWire(self, structure):
        '''
        Adds vacuum padding around a wire.
        
        Returns a wire structure that has been padded with vacuum.
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
    
    def padCluster(self, structure):
        '''
        Adds vacuum padding around a cluster.
        
        Returns a cluster structure that has been padded with vacuum.
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
    
    
    def unpad(self, structure):
        '''
        Removes vacuum padding to return an organism's structure to a form used by the variations.
        
        Returns a structure that has had the vacuum padding removed.
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # rotate structure into principal directions
        # Call other methods based on the value of self.shape
        if self.shape == 'sheet':
            return self.unpadSheet(structure)
        elif self.shape == 'wire':
            return self.unpadWire(structure)
        elif self.shape == 'cluster':
            return self.unpadCluster(structure)
        else:
            return structure
        
    
    def unpadSheet(self, structure):
        '''
        Removes vertical vacuum padding from a sheet.
        
        Returns a sheet structure with the vertical vacuum padding removed
        
        Precondition: the sheet structure is represented with the c-lattice vector perpendicular to the plane of the sheet
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
        
    def unpadWire(self, structure):
        '''
        Removes vacuum padding around a wire.
        
        Returns a wire structure with the vertical vacuum padding removed
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
        
    def unpadCluster(self, structure):
        '''
        Removes vacuum padding around a cluster.
        
        Returns a wire structure with the vertical vacuum padding removed
        
        Args:
            structure: the structure of an organism, as a pymatgen.core.structure.Structure object
        '''
        # TODO: implement me
        
       
class CompositionSpace(object):
    '''
    Represents the composition space to be searched by the algorithm.
    
    This class is a singleton.
    '''
    
    def __init__(self, end_points):
        '''
        Creates a CompositionSpace object, which is list of pymatgen.core.composition objects
        
        TODO: possible cast the elements of endpoints to pymatgen.core.composition.Composition objects
              (I tried before but I kept getting this error: "TypeError: 'unicode' object is not callable")
        
        Args:
            end_points: the list dictionaries, with each one representing a compositions
        '''
        self.end_points = end_points
        
            
       
        
class OrganismCreator(object):
    '''
    Creates organisms for the initial population
    
    Not meant to be instantiated, but rather subclassed by particular Creators, like RandomOrganismCreator
    or PoscarsOrganismCreator.
    '''
    
    def createOrganism(self):
        '''
        Creates an organism for the initial population. Handles development and redundancy checking, and adds valid
        organism to whole_pop.
        
        Returns a developed, non-redundant organism, or None if one could not be created
        
        Args:
            TODO: think about what data this method needs (e.g. objective function data, etc.) and add it 
                to the argument list.
        '''
        raise NotImplementedError("Please implement this method.")
        #
        # The general outline of this method (to be implemented in subclasses) is:
        #     1. create an organism
        #     2. develop the organism
        #     3. if the organism fails development, make another one and try again
        #     4. check for redundancy 
        #     5. if the organism fails redundancy, make another one and try again 
        #     6. add successful organism to whole_pop
        


class RandomOrganismCreator(OrganismCreator):
    '''
    Creates random organisms for the initial population
    '''
    def __init__(self, num_to_make, volume, composition_space, development, redundancy_guard, num_made=0, when_stop="successes", is_finished=False):
        '''
        Creates a RandomOrganismCreator.
        
        Args:
            needed_parameters: all the parameters needed for creating the random organisms.
            This includes how many to make and whether to scale their volumes (and if yes, to what value).
            It will also need stoichiometry information so it knows what types of atoms to use.
            
            development: the Development object (for cell reduction and structure constraints)
            
            redundancy_guard: the redundancyGuard object 
            
            num_to_make: the number of organisms to make with this creator. TODO: this will probably included in needed_parameters
            
            when_stop: the criteria for when this creator is finished ("successes" or "attempts") TODO: this will probably included in needed_parameters
            
            num_successes: the number of organisms successfully added to the initial population from this creator
            
            is_finished: whether the creator has made enough organisms
        '''
        self.development = development
        self.redundancy_guard = redundancy_guard
        self.num_to_make = num_to_make
        self.num_made = num_made
        self.when_stop = when_stop
        self.is_finished = is_finished
    
    def createOrganism(self):
        '''
        Creates a random organism for the initial population. Handles development and redundancy checking, and adds valid
        organism to whole_pop.
        
        Returns a developed, non-redundant organism
        
        Args:
        
            TODO: think about what data this method needs and add it to the argument list.
        '''
        # This method will need to:
        #     1. create a random organism
        #     2. develop the organism
        #     3. if the organism fails development, make another one and try again
        #     4. check for redundancy 
        #     5. if the organism fails redundancy, make another one and try again
        #     6. add successful organism to whole_pop
        #     7. return the successful organism
        
    
    def updateStatus(self):
        '''
        Increments num_made, and if necessary, updates is_finished
        '''
        self.num_made = self.num_made + 1
        if self.num_made == self.num_to_make:
            self.is_finished = True
        
        


class PoscarOrganismCreator(OrganismCreator):
    '''
    Creates organisms from poscar files for the initial population.
    '''
    def __init__(self, needed_parameters, development, redundancy_guard, num_to_make, num_made=0, when_stop="attempts", is_finished=False):
        '''
        Creates a PoscarOrganismCreator.
        
        Args:
            needed_parameters: all the parameters needed for creating the random organisms.
            This includes how many to make and whether to scale their volumes (and if yes, to what value).
            It will also need stoichiometry information so it knows what types of atoms to use.
            
            development: the Development object (for cell reduction and structure constraints)
            
            redundancy_guard: the redundancyGuard object 
            
            num_to_make: the number of organisms to make with this creator. TODO: this will probably included in needed_parameters
            
            when_stop: the criteria for when this creator is finished ("successes" or "attempts") TODO: this will probably included in needed_parameters
            
            num_successes: the number of organisms successfully added to the initial population from this creator
            
            is_finished: whether the creator has made enough organisms
        '''
        self.development = development
        self.redundancy_guard = redundancy_guard
        self.num_to_make = num_to_make
        self.num_made = num_made
        self.when_stop = when_stop
        self.is_finished = is_finished

    
    def createOrganism(self):
        '''
        Creates an organism for the initial population from a poscar file. Handles development and redundancy checking, 
        and adds valid organism to whole_pop.
        
        Returns a developed, non-redundant organism, or None if one could not be created
        
        Args:
        
            TODO: think about what data this method needs and add it to the argument list.
        '''
        # This method will need to:
        #     1. create an organism from a poscar file
        #     2. develop the organism
        #     3. if the organism fails development, make another one from the next poscar file
        #     4. check for redundancy
        #     5. if the organism fails redundancy, make another one from the next poscar file
        #     6. add successful organism to whole_pop
        #     7. increment num_made, and update is_finished if needed
        #          self.updateStatus()
        
        
    def updateStatus(self):
        '''
        Increments num_made, and if necessary, updates is_finished
        '''
        self.num_made = self.num_made + 1
        if self.num_made == self.num_to_make:
            self.is_finished = True
        



class RedundancyGuard(object):
    '''
    A redundancy guard.
    
    This is a singleton class.
    '''
    
    def __init__(self, structure_matcher, d_value):
        '''
        Creates a redundancy guard.
        
        Args:
            structure_matcher: a StructureMatcher object for comparing organisms' structures.
            
            other_params: other optional parameters for identifying redundant structures, like
                d-value, etc.
        '''
        self.structure_matcher = structure_matcher
        self.d_value = d_value

        
    def checkRedundancy(self, new_organism, whole_pop):
        '''
        Checks for redundancy, both structural and if specified, value (d-value)
        
        Returns the organism with which new_organism is redundant, or None if no redundancy
        
        Args:
            new_organism: the organism to check for redundancy
            
            whole_pop: the list containing all organisms to check against
        '''
        for organism in whole_pop:
            if self.structure_matcher.fit(new_organism.structure, organism.structure) == True:
                print("message that new_organism failed structural redundancy")
                return organism
            if new_organism.value != None and self.d_value != None and abs(new_organism.value - organism.value) < self.d_value:
                print("message that new_organism failed value redundancy")
                return organism    
        return None    # should only get here if no organisms are redundant with the new organism
        


class Constraints(object):
    '''
    Represents the general constraints imposed on structures considered by the algorithm. 
    '''
    
    def __init__(self, constraints_parameters, composition_space):
        '''
        Sets the general constraints imposed on structures. Assigns default values if needed.
        
        Args:
            constraints_parameters: a dictionary of parameters
            
            composition_space: a CompositionSpace object describing the composition space to be searched
        '''
        # default values
        self.default_min_num_atoms = 2
        self.default_max_num_atoms = 20
        self.default_min_lattice_length = 0.5
        self.default_max_lattice_length = 20
        self.default_min_lattice_angle = 40
        self.default_max_lattice_angle = 140
        
        # set defaults if entire "Constraints" block was set to default or left blank
        self.constraints_parameters = constraints_parameters
        if self.constraints_parameters == None or self.constraints_parameters == 'default':
            self.set_all_to_defaults()
        else:
            # check each one and see if it's been left blank or set to default
            self.min_num_atoms = constraints_parameters['min_num_atoms']
            if self.min_num_atoms == None or self.min_num_atoms == 'default':
                self.min_num_atoms = self.default_min_num_atoms
                
            self.max_num_atoms = constraints_parameters['max_num_atoms']
            if self.max_num_atoms == None or self.max_num_atoms == 'default':
                self.max_num_atoms = self.default_max_num_atoms
                
            self.min_lattice_length = constraints_parameters['min_lattice_length']
            if self.min_lattice_length == None or self.min_lattice_length == 'default':
                self.min_lattice_length = self.default_min_lattice_length
            
            self.max_lattice_length = constraints_parameters['max_lattice_length']
            if self.max_lattice_length == None or self.max_lattice_length == 'default':
                self.max_lattice_length = self.default_max_lattice_length
                
            self.min_lattice_angle = constraints_parameters['min_lattice_angle']
            if self.min_lattice_angle == None or self.min_lattice_angle == 'default':
                self.min_lattice_angle = self.default_min_lattice_angle
                
            self.max_lattice_angle = constraints_parameters['max_lattice_angle']
            if self.max_lattice_angle == None or self.max_lattice_angle == 'default':
                self.max_lattice_angle = self.default_max_lattice_angle
        
            self.per_species_mids = constraints_parameters['per_species_mids']
            if self.per_species_mids == None or self.per_species_mids == 'default':
                self.set_all_mids_to_defaults(composition_space)     
            else:
                # check each pair and set to default if needed
                for key in self.per_species_mids:
                    if self.per_species_mids[key] == None or self.per_species_mids[key] == 'default':
                        elements = key.split()
                        radius1 = Element(elements[0]).atomic_radius
                        radius2 = Element(elements[1]).atomic_radius
                        self.per_species_mids[key] = 0.8*(radius1 + radius2) 
                
                
    def set_all_to_defaults(self, composition_space):
        '''
        Sets all general constraints (those in Constraints block of input file) to default values
        
        Args:
            composition_space: the composition space object
        '''
        self.min_num_atoms = self.min_num_atoms
        self.max_num_atoms = self.max_num_atoms
        self.min_lattice_length = self.min_lattice_length
        self.max_lattice_length = self.max_lattice_length
        self.min_lattice_angle = self.default_min_lattice_angle
        self.max_lattice_angle = self.default_max_lattice_angle
        self.set_all_mids_to_defaults(composition_space) 
        
        
    def set_all_mids_to_defaults(self, composition_space):
        '''
        Sets all the per-species mids to default values based on atomic radii
        
        Args:
            composition_space: the composition space object
        '''
        # get each element type from the composition_space object
        elements = []
        for point in composition_space:
            for key in point:
                elements.append(key)
        # remove duplicates from the list of elements
        elements = list(set(elements))    
        # get the atomic radius for each element type and put them in a dictionary
        atomic_radii = {}
        for element in elements:
            atomic_radii[element] = Element(element).atomic_radius
            # compute the mid for each pair of elements from their atomic radii as 80% of the sum of the radii
        self.per_species_mids = {}
        for element in elements:
            for next_element in elements:      
                self.per_species_mids[element + " " + next_element] = 0.8*(atomic_radii[element] + atomic_radii[next_element])
        # TODO: maybe remove duplicates, or maybe just leave them - could be helpful to have all representations 
       
        


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
        # 1. Check composition to see if it's in the composition space (could be a little tricky for pd searches...)
        # 2. Do Niggli cell reduction, if specified
        # 3. Check the structural constraints. Do most likely to fail first
        #    - per-species MIDs
        #    - min and max num atoms
        #    - min and max lattice angles
        #    - min and max lattice lengths


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
        Generates a valid offspring organism using the variations and adds it to whole_pop.
        
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
                offspring = self.development.develop(offspring)
                if (offspring != None) and (self.redundancy_guard.checkRedundancy(offspring, whole_pop) == None):
                    whole_pop.append(offspring)
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
        Calculates the energy of an organism
        
        Returns an organism that has been parsed from the output files of the energy code, or None if the calculation 
        failed. Does not do development or redundancy checking.
        
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
        
        

class VaspEnergyCalculator(object):
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
        Calculates the energy of an organism using VASP
        
        Returns an organism that has been parsed from the output files of the energy code, or None if the calculation 
        failed. Does not do development or redundancy checking.
        
        Args:
            org: the organism whose energy we want to calculate
        '''
        # TODO: implement me
        # 1. prepare input files for calculation
        # 2. submit calculation (by running external script)
        # 3. when external script returns, parse organism from the energy code output files
        


class InitialPopulation():
    '''
    The initial population of organisms
    '''
    
    def __init__(self, whole_pop):
        '''
        Args:
            whole_pop: the list containing the organisms seen by the algorithm for redundancy checking
        '''
        self.initial_population = []
    
    
    def addOrganism(self, org, whole_pop):
         '''
        Adds a relaxed organism to the initial population and updates whole_pop.
        
        Args:
            org: the organism whose energy we want to calculate
            
            whole_pop: the list containing all the organisms that the algorithm has submitted for energy calculations
        '''
      #  initial_population.append(org)
      #  org.isActive = True
      #  whole_pop.append(org)
      
      
    def replaceOrganism(self, old_org, new_org):
        '''
        Replaces an organism in the initial population with a new organism.
        
        Precondition: the old_org is a current member of the initial population
        
        Args:
            old_org: the organism in the initial population to replace
            new_org: the new organism to replace the old one
        '''
        # TODO: implement me
        # 1. do the replacement
        # 2. set old_org.isActive = False and newOrg.isActive = True
        # 3. whole_pop.append(new_org)
    

        
        
        
        
########## area for casual testing ########## 

# make a structure object
lattice = [[1,0.5,0], [0.5,1,0], [0,0,-1]]
species = ["C", "Si"]
coordinates = [[0.25,0.25,0.25],[0.75,0.75,0.75]]
structure1 = Structure(lattice, species, coordinates)

# make an organism object
org1 = Organism(structure1)

print(org1.structure.lattice)

org1.rotateToPrincipalDirections()
print("")
print(org1.structure.lattice)

#print(org1.structure)
#print(org1.fitness)
#print(org1.id)

#org1.id = 6
#print(org1.id)

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
        maxClusterDiameter
        maxWireDiameter
        maxLayerThickness
        maxNumAtoms
        minNumAtoms
        minNumSpecies
        doNonnegativityConstraint (maybe don't need this?)
        dValue
        
    Convergence Criteria:
        maxFunctionEvals
        maybe others
   
   Objective Function:
       objectiveFunction <epa|pd>
       geometry <bulk|sheet|wire|cluster>
           # these should be ignored if geometry is bulk
           padding # how much vacuum to add (around cluster/wire/sheet)
           maxSize # how big it can be (diameter for cluster and wire, thickness for sheet. measured from atom-to-atom)
       energyCode <vasp|gulp|others>
       energyFiles (depends on energyCode. Not sure best way to handle this...)
        
'''
