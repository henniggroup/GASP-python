from __future__ import division, unicode_literals, print_function

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.sites import Site
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_matcher import ElementComparator
from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.transformations.standard_transformations import RotationTransformation
import pymatgen.command_line.gulp_caller as gulp_caller

from _pyio import __metaclass__
from _collections import deque
import os
from os import listdir, mkdir, getcwd
from os.path import isfile, join
import numpy as np
import copy
import math
import warnings
from subprocess import check_output, CalledProcessError, STDOUT
from numpy import inf, Inf

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
        self.id += 1
        return self.id



class Organism(object):
    '''
    An organism
    '''
    def __init__(self, structure, id_generator, maker):
        '''
        Creates an organism
        
        Args:
            structure: The structure of this organism, as a pymatgen.core.structure.Structure
            
            id_generator: the instance of IDGenerator used to assign id numbers to all organisms
            
            maker: the name of algorithm that made the organism, as a string. Either a creator or variation
        '''
        # initialize instance variables
        self.structure = structure
        self.composition = self.structure.composition
        self.total_energy = None
        self.epa = None
        self.value = None # the objective function value of this organism, which is either the energy per atom (for fixed-composition search) or the distance from the current best convex hull (for phase diagram search).
        self.fitness = None # the fitness of this organism. Ranges from 0 to 1, including both endpoints
        self.selection_prob = None # the selection probability of this organism. Ranges from 0 to 1, including both endpoints
        self.selection_loc = None # location on the number line (between 0 and 1) where the interval of length selection_prob starts
        self.is_active = False # whether this organism is currently part of the pool or initial population
        self._id = id_generator.makeID(); # unique id number for this organism. Should not be changed.
        self.made_by = maker   # the name of the algorithm that made this organism - either one of the creators or one of the variations. Should be set each time a new organism is made
    
    # this keeps the id (sort of) immutable by causing an exception to be raised if the user tries to 
    # the set the id with org.id = some_id.
    @property
    def id(self):
        return self._id
    
    
    def rotateToPrincipalDirections(self):
        '''
        Rotates the organism's structure into the principal directions. That is, a is parallel to the Cartesian x-axis, 
        b lies in the Cartesian x-y plane and the z-component of c is positive.
        
        Note: this method doesn't change the fractional coordinates of the sites. However, the Cartesian coordinates may be changed
        '''
        # rotate about the z-axis to align a vertically with the x-axis
        rotation = RotationTransformation([0, 0, 1], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[0][1], self.structure.lattice.matrix[0][0]))
        self.structure = rotation.apply_transformation(self.structure)
        # rotate about the y-axis to make a parallel to the x-axis
        rotation = RotationTransformation([0, 1, 0], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[0][2], self.structure.lattice.matrix[0][0]))
        self.structure = rotation.apply_transformation(self.structure)
        # rotate about the x-axis to make b lie in the x-y plane
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
            # mirror c across the x-y plane to make it positive - have to build a new lattice to do this
            # the components of a
            ax = self.structure.lattice.matrix[0][0]
            ay = self.structure.lattice.matrix[0][1]
            az = self.structure.lattice.matrix[0][2]
            # the components of b
            bx = self.structure.lattice.matrix[1][0]
            by = self.structure.lattice.matrix[1][1]
            bz = self.structure.lattice.matrix[1][2]
            # the components of c
            cx = self.structure.lattice.matrix[2][0]
            cy = self.structure.lattice.matrix[2][1]
            cz = -1*self.structure.lattice.matrix[2][2]
            
            self.structure.modify_lattice(Lattice([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]))
        
        
    def rotateCParallelToZ(self):
        '''
        Rotates the organism's structure such that the c lattice vector is parallel to the z axis.
        
        Note: this method doesn't change the fractional coordinates of the sites. However, the Cartesian coordinates may be changed
        '''
        # rotate about the z-axis until c lies in the x-z plane
        rotation = RotationTransformation([0, 0, 1], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[2][1], self.structure.lattice.matrix[2][0]))
        self.structure = rotation.apply_transformation(self.structure)
        # rotate about the y-axis to make c parallel to the z-axis
        rotation = RotationTransformation([0, 1, 0], 180 - (180/np.pi)*np.arctan2(self.structure.lattice.matrix[2][0], self.structure.lattice.matrix[2][2]))
        self.structure = rotation.apply_transformation(self.structure)
        # make sure c is pointing along the positive z-axis
        if self.structure.lattice.matrix[2][2] < 0:
            # rotate 180 degrees about the x-axis
            rotation = RotationTransformation([1, 0, 0], 180)
            self.structure = rotation.apply_transformation(self.structure)
        
     
    def translateAtomsIntoCell(self):
        '''
        Translates all the atoms into the cell, so that their fractional coordinates are between 0 and 1
        '''
        for i in range(0, len(self.structure.frac_coords)):
            translation_vector = []
            for j in range(0, len(self.structure.frac_coords[i])):
                # compute the needed shift in this dimension
                dim = self.structure.frac_coords[i][j]
                if dim > 1.0:
                    translation_vector.append(-1*int(dim))
                elif dim < 0.0:
                    translation_vector.append(-1*int(dim) + 1)
                else:
                    # no shift needed in this dimension
                    translation_vector.append(0)
            # translate the atom by the translation vector
            self.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        
        
    def reduceSheetCell(self):
        '''
        Applies Niggli cell reduction to a sheet structure. 
        
        The idea is to make c vertical and add lots of vertical vacuum so that the standard reduction algorithm only changes the a and b lattice vectors
        
        TODO: pymatgen's Niggli cell reduction algorithm sometimes moves the atoms' relative positions a little (I've seen up to 0.5 A...). 
        '''
        # rotate into principal directions
        self.rotateToPrincipalDirections()
        # get the species and their Cartesian coordinates
        species = self.structure.species
        cartesian_coords = self.structure.cart_coords
        # get the non-zero components of the a and b lattice vectors, and the vertical component of the c lattice vector
        ax = self.structure.lattice.matrix[0][0]
        bx = self.structure.lattice.matrix[1][0]
        by = self.structure.lattice.matrix[1][1]
        cz = self.structure.lattice.matrix[2][2]
        # make a new lattice with a ton of vertical vacuum (add 100 Angstroms)
        padded_lattice = Lattice([[ax, 0.0, 0.0], [bx, by, 0.0], [0.0, 0.0, cz + 100]])
        # make a new structure with the padded lattice and Cartesian coordinates
        padded_structure = Structure(padded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        # do cell reduction on the padded structure (the c lattice vector should still be parallel to z, and a and b should still lie in x-y plane)
        reduced_structure = padded_structure.get_reduced_structure()
        # unpad the reduced structure
        rspecies = reduced_structure.species
        rcartesian_coords = reduced_structure.cart_coords
        rax = reduced_structure.lattice.matrix[0][0]
        ray = reduced_structure.lattice.matrix[0][1]
        rbx = reduced_structure.lattice.matrix[1][0]
        rby = reduced_structure.lattice.matrix[1][1]
        unpadded_lattice = Lattice([[rax, ray, 0.0], [rbx, rby, 0.0], [0.0, 0.0, cz]])
        # set the organism's structure to the unpadded one, and make sure all the atoms are located inside the cell
        self.structure = Structure(unpadded_lattice, rspecies, rcartesian_coords, coords_are_cartesian=True)
        self.translateAtomsIntoCell()
        
        
    def getBoundingBox(self, cart_coords = True):
        '''
        Returns the smallest and largest coordinates in each dimension of all the sites in the organism's structure
        
        Args:
            frac_coords: whether to give the result in Cartesian or fractional coordinates
            
        TODO: maybe there's a cleaner way to do this...
        '''
        # get coordinates of the atoms in the structure
        if cart_coords:
            coords = self.structure.cart_coords
        else:
            coords = self.structure.frac_coords
        # find the largest and smallest coordinates in each dimenstion
        minx = Inf
        maxx = -Inf
        miny = Inf
        maxy = -Inf
        minz = Inf
        maxz = -Inf
        for coord in coords:
            if coord[0] < minx:
                minx = coord[0]
            if coord[0] > maxx:
                maxx = coord[0]
            if coord[1] < miny:
                miny = coord[1]
            if coord[1] > maxy:
                maxy = coord[1]
            if coord[2] < minz:
                minz = coord[2]
            if coord[2] > maxz:
                maxz = coord[2]
        bounds = [[minx, maxx], [miny, maxy], [minz, maxz]]
        return bounds
    
    

class InitialPopulation():
    '''
    The initial population of organisms
    '''
    def __init__(self, run_dir_name):
        '''
        Creates an initial population
        
        Args:
            run_dir_name: the name (not path) of the garun directory where the search will be done
        '''
        self.run_dir_name = run_dir_name
        self.initial_population = []
    
    
    def addOrganism(self, organism_to_add):
        '''
        Adds a relaxed organism to the initial population and updates whole_pop.
        
        Args:
            organism_to_add: the organism to add to the pool
        '''
        # write the organism we're adding to a poscar file
        organism_to_add.structure.to('poscar', os.getcwd() + '/' + self.run_dir_name + '/POSCAR.' + str(organism_to_add.id))
        
        print('Adding organism {} to the initial population'.format(organism_to_add.id))
        self.initial_population.append(organism_to_add)
        organism_to_add.is_active = True
      
      
    def replaceOrganism(self, old_org, new_org):
        '''
        Replaces an organism in the initial population with a new organism.
        
        Precondition: the old_org is a current member of the initial population
        
        TODO: instead of finding the old organism in the initial population by searching for the right id, it might better to eventually implement a __eq__(self, other)
              method in the Organism class to determine when two organisms are the same...
        
        Args:
            old_org: the organism in the initial population to replace
            new_org: the new organism to replace the old one
        '''
        # write the new organism to a poscar file
        new_org.structure.to('poscar', os.getcwd() + '/' + self.run_dir_name + '/POSCAR.' + str(new_org.id))
        
        print('Replacing organism {} with organism {} in the initial population'.format(old_org.id, new_org.id))
        # find the organism in self.initial_population with the same id as old_org
        for org in self.initial_population:
            if org.id == old_org.id:
                self.initial_population.remove(org)
                org.is_active = False
        # add the new organism to the initial population
        self.initial_population.append(new_org)
        new_org.is_active = True
        
    
    def printBestOrg(self, composition_space):
        '''
        Print out the value and number of the best organism in the initial population
        
        Args:
            composition_space: the CompositionSpace object
            
        TODO: right now, this only works for epa searches. For pd searches, we will want to print out the area/volume of the convex hull, and also which organisms are on the convex hull.
              Might be a little tricky because it won't be possible to construct a convex hull unless we have at least one organism from each of the composition endpoints. Will have to
              check for that somehow
        '''
        if composition_space.objective_function == 'epa':
            # sort the organisms in order of increasing epa
            sorted_list = copy.deepcopy(self.initial_population)
            sorted_list.sort(key = lambda x: x.epa, reverse = False)
            # print out the best one
            print('Organism {} is the best and has value {} eV/atom'.format(sorted_list[0].id, sorted_list[0].epa))
        elif composition_space.objective_function == 'pd':
            pass
            # TODO: implement me 
            
            

class Pool(object):
    '''
    Pool to hold all the organisms that are currently candidates for selection to become parents.
    
    Composed of two parts: a promotion set containing the best few organisms, and a queue containing
        the rest of the organisms in the pool.
    '''
    def __init__(self, pool_params, composition_space, run_dir_name):
        '''
        Creates a pool of organisms
        
        Args:
            pool_params: the parameters of the pool, as a dictionary 
            
            composition_space: the CompositionSpace object
            
            run_dir_name: the name (not path) of the garun directory
        '''
        # the default values
        if composition_space.objective_function == 'epa':
            self.default_size = 21
        elif composition_space.objective_function == 'pd':
            self.default_size = 30
        self.default_num_promoted = 3
        
        # if entire Pool block was set to default or left blank
        if pool_params == None or pool_params == 'default':
            self.size = self.default_size
            self.num_promoted = self.default_num_promoted
        else:     
            # check if size parameter was used
            if 'size' not in pool_params:
                self.size = self.default_size
            elif pool_params['size'] == None or pool_params['size'] == 'default':
                self.size = self.default_size
            else:
                self.size = pool_params['size']
            
            # check if num_promoted parameter was used
            if 'num_promoted' not in pool_params:
                self.num_promoted = self.default_num_promoted
            elif pool_params['num_promoted'] == None or pool_params['num_promoted'] == 'default':
                self.num_promoted = self.default_num_promoted
            else:
                self.num_promoted = pool_params['num_promoted']
        
        self.promotion_set = [] # list to hold the best few organisms in the pool
        self.queue = deque() # deque containing the rest of the organisms in the pool
        self.selection = [] # list to hold the parameters for the selection distribution
        self.num_adds = 0 # the number of organisms added to the pool (excluding the initial population)
        self.run_dir_name = run_dir_name # the name (not  path) of the garun directory
    
    
    def addInitialPopulation(self, initial_population, composition_space):
        '''
        Adds the organisms of the initial population to the pool.
        
        Args:
            initial_population: an InitialPopulation object containing the relaxed organisms of the initial population
            
            composition_space: the CompositionSpace object
            
        Precondition: the pool is empty (no organisms have been previously added to it)
        '''
        print('Populating the pool with the initial population')
        
        # get the list of organisms in the initial population
        organisms_list = initial_population.initial_population
        
        # compute the values of the organisms in the initial population
        self.computeInitialPopValues(organisms_list, composition_space)
        
        # sort the list of initial population organisms in order of increasing value
        organisms_list.sort(key = lambda x: x.value, reverse = False)
        
        # populate the promotion set and queue
        if composition_space.objective_function == 'epa':
            for i in range(len(organisms_list)):
                if i < self.num_promoted:
                    self.promotion_set.append(organisms_list[i])
                else:
                    self.queue.appendleft(organisms_list[i])
        elif composition_space.objective_function == 'pd':
            for organism in organisms_list:
                if organism.value < 0.0001: # TODO: this should be zero. Change it after testing
                    self.promotion_set.append(organism)
                else:
                    self.queue.appendleft(organism)
        
        # compute the fitnesses of the organisms in the initial population
        self.computeFitnesses()
        
        # compute the selection probabilities of all the organisms in the initial population
        self.computeSelectionProbs()
        
        # print out some info on the organisms in the initial population
        print('Summary of the initial population:')
        for organism in organisms_list:
            print('Organism {} has value {} and fitness {}'.format(organism.id, organism.value, organism.fitness))
            
        # for epa seaches, print out info on the best organism in the initial population
        if composition_space.objective_function == 'epa':
            print('Organism {} is the best and has value {} eV/atom'.format(organisms_list[0].id, organisms_list[0].value))
            
        elif composition_space.objective_function == 'pd':
            # print out the organisms that are on the convex hull (in the promotion set)
            for organism in self.promotion_set:
                print('Organsims {} is on the convex hull and has value {}'.format(organism.value))
                
            # TODO: compute and print out the area/volume of the convex hull here
        
    
    def computeInitialPopValues(self, organisms_list, composition_space):
        '''
        Computes and assigns the values of all the organisms in the initial population
        
        Args:
            organisms_list: a list of the relaxed organisms of the initial population
            
            composition_space: the CompositionSpace object
        '''
        if composition_space.objective_function == 'epa':
            for organism in organisms_list:
                organism.value = organism.epa
        elif composition_space.objective_function == 'pd':
            # create a PDEntry object for each organism in the list of organisms 
            pdentries = {}
            for organism in organisms_list:
                pdentries[organism.id] = PDEntry(organism.total_energy, organism.composition) 
       
            # create a compound phase diagram object from the pdentries
            pdentries_list = []
            for organism_id in pdentries:
                pdentries_list.append(pdentries[organism_id])
            compound_pd = CompoundPhaseDiagram(pdentries_list, composition_space.endpoints)
        
            # create a phase diagram analyzer object from the compound phase diagram
            pd_analyzer = PDAnalyzer(compound_pd)
        
            # compute the value of every organism in the list
            for organism in organisms_list:
                organism.value = pd_analyzer.get_e_above_hull(PDEntry(organism.total_energy, organism.composition))
        
    
    def addOrganism(self, organism_to_add, composition_space):
        '''
        Adds a new organism to the pool
        
        For epa, if the new organism better than one of the orgs currently in the promotion set, then it is added 
        to the promotion set, and the worst org in the promotion set is moved to the back of the queue. 
        Otherwise, the new org is just appended to the back of the queue.
        
        For pd, if the new organism is on the convex hull, then it is added to the promotion set and any organisms 
        in the promotion set that are no longer on the convex hull are moved the queue. Otherwise, the new organism
        is added to the back of the queue.
        
        Args:
            organism: the Organism to add to the pool
            
            composition_space: the CompositionSpace object
        '''
        print('Adding organism {} to the pool'.format(organism_to_add.id))
        
        # increment the number of adds to the pool
        self.num_adds = self.num_adds + 1
        
        # write the organism we're adding to a poscar file
        organism_to_add.structure.to('poscar', os.getcwd() + '/' + self.run_dir_name + '/POSCAR.' + str(organism_to_add.id))
        
        # the organism that we're adding to the pool is now active
        organism_to_add.is_active = True
        
        # compute the value of the new organism (this updates the values of the other organisms in the pool if needed)
        self.computeValue(organism_to_add, composition_space)
        
        if composition_space.objective_function == 'epa':
            # check if the new organism's value is better than the worst one in the promotion set
            worst_organism = self.getWorstInPromotionSet()
            if organism_to_add.value < worst_organism.value:
                # add the new organism to the promotion set
                self.promotion_set.append(organism_to_add)
                # move the worst organism from the promotion set to the back (left end) of the queue
                self.promotion_set.remove(worst_organism)
                self.queue.appendleft(worst_organism)
            else:
                # append the new organism to the back (left end) of the queue
                self.queue.appendleft(organism_to_add)
                  
        elif composition_space.objective_function == 'pd':
            # check if any of the organisms in the promotion set are no longer on the convex hull
            threshhold = 0.0001 # just in case the get_e_above_hull method doesn't work quite right. Can probably remove this once we confirm
            organisms_to_move = []
            # get a list of organisms in the promotion set that are no longer on the convex hull
            for organism in self.promotion_set:
                if organism.value >= threshhold:
                    organisms_to_move.append(organism)
            # move all the organism that are no longer on the convex hull from the promotion set to the queue       
            for organism_to_move in organisms_to_move:    
                self.promotion_set.remove(organism_to_move)
                self.queue.appendleft(organism_to_move)
                    
            # determine whether to add the new organism to the promotion set or the queue
            if organism_to_add.value < threshhold:
                self.promotion_set.append(organism_to_add)
            else:
                self.queue.appendleft(organism_to_add)
                    
        # compute fitnesses
        self.computeFitnesses()
        
        # compute selection probabilities
        self.computeSelectionProbs()
        
    
    def getWorstInPromotionSet(self):
        '''
        Returns the organism in the promotion set with the worst (largest) objective function value
        '''
        worst_organism = self.promotion_set[0]
        for organism in self.promotion_set:
            if organism.value > worst_organism.value:
                worst_organism = organism
        return worst_organism
        
    
    def replaceOrganism(self, old_org, new_org, composition_space):
        '''
        Replaces an organism in the pool with a new organism. The new organism has the same location in the pool
        as the old one.
        
        Precondition: the old_org is a member of the current pool.
        
        Args:
            old_org: the organism in the pool to replace
            
            new_org: the new organism to replace the old one
        '''
        print('Replacing organism {} with organism {} in the pool'.format(old_org.id, new_org.id))
        
        # write the new organism to a poscar file
        new_org.structure.to('poscar', os.getcwd() + '/' + self.run_dir_name + '/POSCAR.' + str(new_org.id))
        
        # compute the value of the new organism
        self.computeValue(new_org, composition_space)
        
        if old_org in self.promotion_set:
            self.promotion_set.remove(old_org)
            self.promotion_set.append(new_org)
        elif old_org in self.queue:
            # cast the queue to a list to do the replacement
            queue_list = list(self.queue)
            # insert the new organism next to the old organism in the list
            queue_list.insert(queue_list.index(old_org), new_org)
            # remove the old organism from the list
            queue_list.remove(old_org)
            # cast it back to a deque
            self.queue = deque(queue_list)
            
        # update is_active flags on both organisms
        old_org.is_active = False
        new_org.is_active = True
        
        # recompute fitnesses and selection probabilities of all organisms in the pool
        # TODO: don't actually have to update all their fitnesses every time...
        self.computeFitnesses()
        self.computeSelectionProbs()
    
    
    def computeValue(self, organism, composition_space):
        '''
        Computes the objective function value of an organism. This doesn't change the locations of any of the organisms in the pool, and it
        doesn't add or remove an organism from the pool.
        
        Args:
            organism: the Organism whose value we need to compute. Not yet added to the pool
            
            composition_space: the CompositionSpace object
        '''
        # if epa search, then value is just epa
        if composition_space.objective_function == 'epa':
            organism.value = organism.epa
        # if pd search, then need to get the convex hull to compute value
        elif composition_space.objective_function == 'pd':
            self.computePDValue(organism, composition_space)
        
    
    def computePDValue(self, organism, composition_space):
        '''
        Computes an organism's value, as the distance from the current best convex hull. If the organism lies on the convex hull, then recomputes of the values
        of all the organisms in the pool relative to the new best convex hull.
        
        Args:
            organism: the Organism whose value we need to compute
            
            composition_space: the CompositionSpace object
        '''
        # create a PDEntry object for the new organism and for each organism in the pool and place them in a dictionary, where the keys are the organism id's
        pdentries = {}
        pdentries[organism.id] = PDEntry(organism.total_energy, organism.composition)
        pool_list = self.toList()
        for pool_org in pool_list:
            pdentries[pool_org.id] = PDEntry(pool_org.total_energy, pool_org.composition) 
       
        # create a compound phase diagram object from the pdentries
        pdentries_list = []
        for organism_id in pdentries:
            pdentries_list.append(pdentries[organism_id])
        compound_pd = CompoundPhaseDiagram(pdentries_list, composition_space.endpoints)
        
        # create a phase diagram analyzer object from the compound phase diagram
        pd_analyzer = PDAnalyzer(compound_pd)
        
        # compute the value of the new organism, and also of every organism in the pool
        organism.value = pd_analyzer.get_e_above_hull(PDEntry(organism.total_energy, organism.composition))
        for pool_org in pool_list:
            pool_org.value = pd_analyzer.get_e_above_hull(PDEntry(pool_org.total_energy, pool_org.composition))
        
    
    def computeFitnesses(self):
        '''
        Calculates and assigns fitnesses to all organisms in the pool.
        
        Precondition: the organisms in the pool all have up-to-date values
        
        TODO: ways to speed this up:
                keep track of the best and worst values in the pool so we don't have to search for them each time
                the fitnesses of all the organisms only change if the best or worst value changed 
        '''
        # get the best and worst values
        best_value = self.getBestInPool().value
        worst_value = self.getWorstInPool().value
        
        # update all the fitnesses
        for organism in self.toList():
            organism.fitness = (organism.value - worst_value)/(best_value - worst_value)


    def getBestInPool(self):
        '''
        Returns the organism in the pool with the best (smallest) objective function value
        '''
        # only need to look in the promotion set because that's where the best ones live
        best_organism = self.promotion_set[0]
        for organism in self.promotion_set:
            if organism.value < best_organism.value:
                best_organism = organism
        return best_organism
    
    
    def getWorstInPool(self):
        '''
        Returns the organism in the pool with the worst (largest) objective function value
        '''
        # only need to look in the queue because that's where the worst ones live
        worst_organism = self.queue[0]
        for organism in self.queue:
            if organism.value > worst_organism.value:
                worst_organism = organism
        return worst_organism
        
    
    def computeSelectionProbs(self):
        '''
        Calculates and assigns selection probabilities to all the organisms in the pool.
        
        Precondition: the organisms in the pool all have up-to-date fitnesses
        
        TODO: there might be ways to speed this up, similar to computing fitnesses...
        '''
        # get a list of the best organisms in the pool, sorted by increasing fitness
        best_organisms = self.getNBestOrganisms(self.selection.num_parents)
        
        # get the denominator of the expression for fitness
        denom = 0
        for organism in best_organisms:
            denom = denom + math.pow(organism.fitness, self.selection.power)
            
        # compute the selection probability and interval location for each organism in best_organisms
        selection_loc = 0
        for organism in best_organisms:
            organism.selection_prob = math.pow(organism.fitness, self.selection.power)/denom
            organism.selection_loc = selection_loc
            # increment the selection location by this organism's selection probability
            selection_loc = selection_loc + organism.selection_prob
            
        # for any organisms not in best_organisms, set the selection probability and location to zero
        for organism in self.toList():
            if organism not in best_organisms:
                organism.selection_prob = 0
                organism.selection_loc = 0
        
       
    def getNBestOrganisms(self, N):
        '''
        Returns a list containing the N best organisms in the pool, sorted in order of increasing fitness
        
        Precondition: all the organisms in the pool have up-to-date fitnesses
        
        Args:
            N: the number of best organisms to get
        '''
        # get a list of all the organisms in the pool
        pool_list = self.toList()
        
        # sort it in order of increasing fitness
        pool_list.sort(key = lambda x: x.fitness, reverse = False)
        
        # get a sublist consisting of the last N elements of pool_list
        return pool_list[len(pool_list) - N:]
        
        
    def selectOrganisms(self, n, random):
        '''
        Randomly selects n distinct organisms from the pool based on their selection probabilities.
        
        Returns a list containing n organisms.
        
        Args:
            n: how many organisms to select from the pool
            
            random: Python's built in PRNG
            
        Precondition: all the organisms in the pool have been assigned selection probabilities.
        '''
        # list to hold the selected organisms
        selected_orgs = []
        pool_list = self.toList()
        # keep going until we've selected enough
        while True:
            rand = random.random()
            # find the organism in the pool with the corresponding selection probability
            for organism in pool_list:
                if rand >= organism.selection_loc and rand < (organism.selection_loc + organism.selection_prob):
                    # check that we haven't already selected this one
                    if organism not in selected_orgs:
                        selected_orgs.append(organism)
                        # check if we've selected enough
                        if len(selected_orgs) == n:
                            return selected_orgs
       
       
    def printSummary(self, composition_space):
        '''
        Prints out a summary of the organisms in the pool
        
        Args:
            composition_space: the CompositionSpace object
        '''
        # get a list of all the organisms in the pool
        pool_list = self.toList()
        
        # sort it in order of decreasing fitness
        pool_list.sort(key = lambda x: x.fitness, reverse = True)
        
        # print out the value and fitness of each organism
        print('Summary of the pool:')
        for organism in pool_list:
            print('Organism {} has value {} and fitness {}'.format(organism.id, organism.value, organism.fitness))
        
        # for epa searches, print out info on the best organism in the pool, and also the average value of organisms in the pool
        if composition_space.objective_function == 'epa':
            print('Organism {} is the best and has value {} eV/atom'.format(pool_list[0].id, pool_list[0].value))
            
        elif composition_space.objective_function == 'pd':
            # print out the organisms that are on the convex hull (in the promotion set)
            for organism in self.promotion_set:
                print('Organsims {} is on the convex hull and has value {}'.format(organism.value))
                
            # TODO: compute and print out the area/volume of the convex hull here
                
                
    def toList(self):
        '''
        Returns a list containing all the organisms in the pool.
        '''
        return self.promotion_set + list(self.queue)
        
        
  
class OffspringGenerator(object):
    '''
    This class handles generating offspring organisms from the pool and the variations. 
    Basically a place for the makeOffspringOrganism method to live
    '''
    def __init__(self):
        pass
  
     
    def makeOffspringOrganism(self, random, pool, variations, geometry, id_generator, whole_pop, development, redundancy_guard, composition_space, constraints):
        '''
        Returns a developed, non-redundant, unrelaxed offspring organism generated with one of the variations
        
        Args:
            TODO: seems like a ridiculous number of arguments...
            
        TODO: outline of how this method works
        '''
        tried_variations = [] # list to hold the variations that have already been tried
        max_num_tries = 10000 # the maximum number of times to try making a valid offspring with a given variation before giving up and trying a different variation
        while(len(variations) > len(tried_variations)):
            variation = self.selectVariation(random, tried_variations, variations)
            num_tries = 0
            while num_tries < max_num_tries:
                offspring = variation.doVariation(pool, random, geometry, id_generator)
                offspring = development.develop(offspring, composition_space, constraints, geometry, pool)
                # don't need to worry about replacement, since we're only dealing with unrelaxed organisms here
                if (offspring != None) and (redundancy_guard.checkRedundancy(offspring, whole_pop) == None):
                    return offspring
                else:
                    num_tries = num_tries + 1
            tried_variations.append(variation)
        # This point should never be reached
        print("Could not a make valid offspring organism with any Variation") 
        print('Quitting...')
        quit()
        
    
    def selectVariation(self, random, tried_variations, variations):
        '''
        Selects a variation that hasn't been tried yet based on their selection probabilities
        
        Args:
            random:
        
            tried_variations: list of Variations that have already been unsuccessfully tried
            
            variations: list of all the Variations used in this search
            
        TODO: description
        '''
        while True:
            rand = random.random()
            lower_bound = 0
            for i in range(len(variations)):
                if rand > lower_bound and rand <= (lower_bound + variations[i].fraction) and variations[i] not in tried_variations:
                    return variations[i]
                else:
                    lower_bound = lower_bound + variations[i].fraction
        


class SelectionProbDist(object):
    '''
    Defines the probability distribution over the fitnesses of the organisms in the pool that determines their selection probabilities
    
    Comprised of two numbers: the number of organisms to select from, and the selection power
    '''
    def __init__(self, selection_params, pool_size):
        '''
        Creates a selection probability distribution
        
        Args:
            selection_params: the parameters defining the distribution, as a dictionary
            
            pool_size: the size of the Pool (how many organisms it contains)
        '''
        # default values
        self.default_num_parents = pool_size
        self.default_power = 1
        
        # if entire Selection block was left blank or set to default
        if selection_params == None or selection_params == 'default':
            self.num_parents = self.default_num_parents
            self.power = self.default_power
        else:
            # check the num_parents parameter
            if 'num_parents' not in selection_params:
                self.num_parents = self.default_num_parents
            elif selection_params['num_parents'] == None or selection_params['num_parents'] == 'default':
                self.num_parents = self.default_num_parents
            # check that the given number isn't larger than the pool size
            elif selection_params['num_parents'] > pool_size:
                self.num_parents = self.default_num_parents
            else:
                self.num_parents = selection_params['num_parents']
                
            # check the selection_power parameter
            if 'power' not in selection_params:
                self.power = self.default_power
            elif selection_params['power'] == None or selection_params['power'] == 'default':
                self.power = self.default_power
            else:
                self.power = selection_params['power']
        
        

class Variation(object):
    '''
    A general variation object for creating offspring organisms from parent organism(s)
    
    Not meant to be instantiated, but rather subclassed by particular Variations, like Mating,
    Mutation, etc.
    
    TODO: is this superclass even needed?
    '''
    def doVariation(self):
        '''
        Creates an offspring organism from parent organism(s).
        
        Returns an organism.
        '''
        raise NotImplementedError("Please implement this method.")



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
               
        Precondition: the 'fraction' parameter in mating_params is not optional, and it is assumed that mating_params contains this parameter 
        '''
        # the name of this variation
        self.name = 'mating'
        
        # parse the fraction from the parameters. This argument is not optional, so it doesn't have a default value here
        self.fraction = mating_params['fraction']
        
        # default values
        # TODO: are these good default values?
        self.default_mu_cut_loc = 0.5  # the average (fractional) location along the randomly chosen lattice vector to make the cut
        self.default_sigma_cut_loc = 0.5 # the standard deviation of the (fractional) location along the randomly chosen lattice vector to make the cut
        self.default_shift_prob = 1.0 # the probability of randomly shifting the atoms along the lattice vector of the cut
        self.default_doubling_prob = 0.1 # the probability that one of the parents will be doubled before doing the variation
        self.default_grow_parents = True # whether or not to grow the smaller parent (by taking a supercell) to the approximate size the larger parent before doing the variation
        
        # TODO: maybe parse these parameters with a loop instead...
        # the mean of the cut location 
        if 'mu_cut_loc' not in mating_params:
            # use the default value if the flag hasn't been used
            self.mu_cut_loc = self.default_mu_cut_loc
        elif mating_params['mu_cut_loc'] == None or mating_params['mu_cut_loc'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.mu_cut_loc = self.default_mu_cut_loc
        else:
            # otherwise, parse the value from the parameters
            self.mu_cut_loc = mating_params['mu_cut_loc']
            
        # the standard deviation of the cut location 
        if 'sigma_cut_loc' not in mating_params:
            # use the default value if the flag hasn't been used
            self.sigma_cut_loc = self.default_sigma_cut_loc
        elif mating_params['sigma_cut_loc'] == None or mating_params['sigma_cut_loc'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.sigma_cut_loc = self.default_sigma_cut_loc
        else:
            # otherwise, parse the value from the parameters
            self.sigma_cut_loc = mating_params['sigma_cut_loc']
            
        # the probability of shifting the atoms
        if 'shift_prob' not in mating_params:
            # use the default value if the flag hasn't been used
            self.shift_prob = self.default_shift_prob
        elif mating_params['shift_prob'] == None or mating_params['shift_prob'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.shift_prob = self.default_shift_prob
        else:
            # otherwise, parse the value from the parameters
            self.shift_prob = mating_params['shift_prob']
            
        # the probability of doubling one of the parents before doing the variation
        if 'doubling_prob' not in mating_params:
            # use the default value if the flag hasn't been used
            self.doubling_prob = self.default_doubling_prob
        elif mating_params['doubling_prob'] == None or mating_params['doubling_prob'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.doubling_prob = self.default_doubling_prob
        else:
            # otherwise, parse the value from the parameters
            self.doubling_prob = mating_params['doubling_prob']
            
        # whether to grow the smaller of the parents before doing the variation
        if 'grow_parents' not in mating_params:
            # use the default value if the flag hasn't been used
            self.grow_parents = self.default_grow_parents
        elif mating_params['grow_parents'] == None or mating_params['grow_parents'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.grow_parents = self.default_grow_parents
        else:
            # otherwise, parse the value from the parameters
            self.grow_parents = mating_params['grow_parents']
            
    
    def doVariation(self, pool, random, geometry, id_generator):
        '''
        Performs the mating operation
        
        Returns the resulting offspring as an Organism
        
        Args:
            pool: the Pool of Organisms
            
            random: Python's built in PRNG
            
            geometry: the Geometry object
            
            id_generator: the IDGenerator
            
        Description:
        
            Creates an offspring organism by combining chunks cut from two parent organisms.
            
                1. Selects two organisms from the pool to act as parents, and makes copies of them.
            
                2. Optionally doubles one of the parents. This occurs with probability self.doubling_prob, and if it happens, the parent 
                   with the smallest cell volume is doubled.
                   
                3. Optionally grows one the parents to the approximate size of the other parent, if self.grow_parents is True. The number 
                   of times to double the smaller parent is determined by the ratio of the cell volumes of the parents. 
                   
                4. Randomly selects one of the three lattice vectors to slice.
                
                5. Determines the cut location (in fractional coordinates) by drawing from a Gaussian with mean self.mu_cut_loc and standard
                   deviation self.sigma_cut_loc.
                   
                6. In each parent, optionally shift the atoms (in fractional space, with probability self.shift_prob) by an amount drawn from 
                   a uniform distribution along the direction of the lattice vector to cut. For non-bulk geometries, shift only occurs if the 
                   lattice vector to cut is along a periodic direction. 
                   
                7. Copy the sites from the first parent organism with fractional coordinate less than the randomly chosen cut location along
                   the randomly chosen lattice vector to the offspring organism, and do the same for the second parent organism, except copy
                   the sites with fractional coordinate greater than the cut location.
        '''
        # for testing
        #structure_1 = Structure.from_file('/n/srv/brevard/structures/POSCAR.Cu_b')
        #structure_2 = Structure.from_file('/n/srv/brevard/structures/POSCAR.Ni_a')
        #parent_orgs = []
        #parent_orgs.append(Organism(structure_1, id_generator))
        #parent_orgs.append(Organism(structure_2, id_generator))
        
        # select two parent organisms from the pool
        parent_orgs = pool.selectOrganisms(2, random)
        
        # make deep copies of the parent organisms
        parent_1 = copy.deepcopy(parent_orgs[0])
        parent_2 = copy.deepcopy(parent_orgs[1])
        
        # optionally double one of the parents
        if random.random() < self.doubling_prob:
            # pick the smallest parent (based on cell volume)
            vol_1 = parent_1.structure.lattice.volume
            vol_2 = parent_2.structure.lattice.volume
            if vol_1 < vol_2:
                self.doubleParent(parent_1, geometry, random)
            else:
                self.doubleParent(parent_2, geometry, random)
            
        # grow the smaller parent if specified
        if self.grow_parents:
            # pick the smallest parent (based on cell volume)
            vol_1 = parent_1.structure.lattice.volume
            vol_2 = parent_2.structure.lattice.volume
            if vol_1 < vol_2:
                volume_ratio = vol_2/vol_1
                parent_to_grow = parent_1 
            else:
                volume_ratio = vol_1/vol_2
                parent_to_grow = parent_2
            # compute how many times to double the smaller parent based on ratio of cell volumes
            num_doubles = self.getNumDoubles(volume_ratio)
            # double the smaller parent the computed number of times
            for _ in range(num_doubles):
                self.doubleParent(parent_to_grow, geometry, random)
        
        # lists to hold the species of the sites contributed from each parent
        species_from_parent_1 = []
        species_from_parent_2 = []
        
        # loop needed here because sometimes no atoms get contributed from one of the parents, so have to try again
        while len(species_from_parent_1) == 0 or len(species_from_parent_2) == 0:
            # lists to hold the species and fractional coordinates of the sites contributed from each parent. Have to set them to empty each time the loop is entered
            species_from_parent_1 = []
            frac_coords_from_parent_1 = []
            species_from_parent_2 = []
            frac_coords_from_parent_2 = []
            
            # randomly select the lattice vector to cut
            cut_vector_index = random.randint(0, 2)
            
            # draw the random cut location from a Gaussian, and make sure it's between 0 and 1
            cut_location = random.gauss(self.mu_cut_loc, self.sigma_cut_loc)
            while cut_location > 1 or cut_location < 0:
                cut_location = random.gauss(self.mu_cut_loc, self.sigma_cut_loc)
            
            # possibly shift the atoms in the first parent along the cut vector
            if random.random() < self.shift_prob:
                shifted_parent_1 = self.doRandomShift(parent_1, cut_vector_index, geometry, random)
        
            # possibly shift the atoms in the second parent along the cut vector
            if random.random() < self.shift_prob:
                shifted_parent_2 = self.doRandomShift(parent_2, cut_vector_index, geometry, random)
            
            # get the species and fractional coordinates of each site in parent 1 with fractional coordinate along the cut vector less than the cut location
            for site in shifted_parent_1.structure.sites:
                if site.frac_coords[cut_vector_index] < cut_location:
                    species_from_parent_1.append(site.species_and_occu)
                    frac_coords_from_parent_1.append(site.frac_coords)
                
            # get the species and fractional coordinates of each site in parent 2 with fractional coordinate along the cut vector greater than the cut location
            for site in shifted_parent_2.structure.sites:
                if site.frac_coords[cut_vector_index] > cut_location:
                    species_from_parent_2.append(site.species_and_occu)
                    frac_coords_from_parent_2.append(site.frac_coords)
                
        # combine the information for the sites contributed by each parent
        offspring_species = species_from_parent_1 + species_from_parent_2
        offspring_frac_coords = frac_coords_from_parent_1 + frac_coords_from_parent_2
                
        # compute the lattice vectors of the offspring by taking the average of the lattice vectors of the parents
        offspring_lengths = 0.5*(np.array(shifted_parent_1.structure.lattice.abc) + np.array(shifted_parent_2.structure.lattice.abc)) 
        offspring_angles = 0.5*(np.array(shifted_parent_1.structure.lattice.angles) + np.array(shifted_parent_2.structure.lattice.angles)) 
        offspring_lattice = Lattice.from_lengths_and_angles(offspring_lengths, offspring_angles)    
        
        # make the offspring structure from the offspring lattice, species and fractional coordinates
        offspring_structure = Structure(offspring_lattice, offspring_species, offspring_frac_coords) 
        
        # make the offspring organism from the offspring structure, and return it
        offspring = Organism(offspring_structure, id_generator, self.name)
        
        # print out a message
        print("Creating offspring organism {} from parent organisms {} and {} with the mating variation".format(offspring.id, parent_orgs[0].id, parent_orgs[1].id))
        
        return offspring
    
    
    def getNumDoubles(self, volume_ratio):
        '''
        Returns the number of times to double a cell based on the given volume ratio. Essentially maps the volume ratio to a step function that 
        approximates the base 2 logarithm.
        
        Args:
            volume_ratio: the ratio of the volumes of the two parent organisms. Assumed to be less than 48.
        '''
        if volume_ratio < 1.5:
            return 0
        elif volume_ratio >= 1.5 and volume_ratio < 3:
            return 1
        elif volume_ratio >=3 and volume_ratio < 6:
            return 2
        elif volume_ratio >= 6 and volume_ratio < 12:
            return 3
        elif volume_ratio >= 12 and volume_ratio < 24:
            return 4
        elif volume_ratio >= 24 and volume_ratio < 48:
            return 5 
    
    
    def doubleParent(self, organism, geometry, random):
        '''
        Takes a supercell of the organism. For bulk geometries, the supercell is taken in the direction of the organism's shortest lattice vector. 
        For non-bulk geometries, the supercell is taken in the direction of a randomly chosen lattice vector. Modifies the structure of the organism.
        
        Args:
            organism: the Organism to take the supercell of
            
            geometry: a Geometry object
            
        TODO: does this method belong in the Organism class instead?
        '''
        if geometry.shape == 'bulk':
            # get the index of the smallest lattice vector of the organism
            lattice_lengths = organism.structure.lattice.abc
            smallest_vector = min(lattice_lengths)
            doubling_index = lattice_lengths.index(smallest_vector)
        else:
            # if geometry is not bulk, then  pick a random lattice vector to double
            doubling_index = random.choice([0, 1, 2])
        # take a supercell of the organism along the smallest lattice vector
        scaling_factors = [1, 1, 1]
        scaling_factors[doubling_index] = 2
        organism.structure.make_supercell(scaling_factors)
        
    
    def doRandomShift(self, organism, lattice_vector_index, geometry, random):
        '''
        Makes a copy of the organism, and shifts all the atoms in the copy by random amount (drawn from uniform distribution) along the 
        lattice vector specified by the given index. After shifting the atoms by the random amount along the specified lattice vector, calls 
        Organism.translateAtomsIntoCell() to make sure all the atoms are located inside the cell.
        
        Note: this method checks the geometry, and will not do the shift if the specified lattice vector is not in a periodic direction of 
        the geometry of the structure because that would destroy some of the local structure of the non-bulk structure.
        
        Args:
            organism: the Organism whose Structure to change by shifting the atoms
            
            lattice_vector_index: the index (0, 1 or 2) of the lattice vector along which to shift the atoms
            
            geometry: the Geometry object
            
            random: Python's built in PRNG
        '''
       
        # if shape is cluster, then no shift
        if geometry.shape == 'cluster':
            pass
        # if shape is wire, then don't shift if not c lattice vector
        elif geometry.shape == 'wire' and lattice_vector_index != 2:
            pass
        # if shape is sheet, then don't shift if c lattice vector
        elif geometry.shape == 'sheet' and lattice_vector_index == 2:
            pass
        # otherwise, shift all the atoms along the specified lattice vector by a random amount
        else:
            # make a copy of the organism
            shifted_org = copy.deepcopy(organism)
            # randomly select the fractional amount to shift the atoms along the lattice vector
            shift_vector = random.random()*shifted_org.structure.lattice.matrix[lattice_vector_index]
            site_indices = []
            for i in range(0, len(shifted_org.structure.sites)):
                site_indices.append(i)
            shifted_org.structure.translate_sites(site_indices, shift_vector, frac_coords=False, to_unit_cell=False)
            # translate the atoms back into the cell in case the shifts moved some of them outside of it
            shifted_org.translateAtomsIntoCell()
            return shifted_org
        
        
        
class StructureMut(Variation):
    '''
    A operator to create an offspring organism by mutating the structure of a parent organism.
    '''
    def __init__(self, structure_mut_params):
        '''
        Creates a mutation operator
        
        Args:
            structure_mut_params: The parameters for doing the mutation operation, as a dictionary
            
        Precondition: the 'fraction' parameter in structure_mut_params is not optional, and it is assumed that structure_mut_params contains this parameter
        '''    
        # the name of this variation
        self.name = 'structure mutation'
        
        # parse the fraction from the parameters. This argument is not optional, so it doesn't have a default value here
        self.fraction = structure_mut_params['fraction']
        
        # the default values
        # TODO: are these good default values?
        self.default_frac_atoms_perturbed = 1.0            # the fraction of atoms to perturb (on average)
        self.default_sigma_atomic_coord_perturbation = 0.5 # the standard deviation of the perturbation of an atomic coordinate (in Angstroms)
        self.default_max_atomic_coord_perturbation = 5.0   # the maximum allowed perturbation of an atomic coordinate (in Angstroms)
        self.default_sigma_strain_matrix_element = 0.2     # the standard deviation of the non-identity components of the elements of the strain matrix
        # note: the non-identity components of the strain matrix elements are constrained to be between -1 and 1
        
        # TODO: maybe parse these parameters with a loop instead...
        # the fraction of atoms to perturb, on average
        if 'frac_atoms_perturbed' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        elif structure_mut_params['frac_atoms_perturbed'] == None or structure_mut_params['frac_atoms_perturbed'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.frac_atoms_perturbed = self.default_frac_atoms_perturbed
        else:
            # otherwise, parse the value from the parameters
            self.frac_atoms_perturbed = structure_mut_params['frac_atoms_perturbed']
            
        # the standard deviation of the perturbation of each atomic coordinate
        if 'sigma_atomic_coord_perturbation' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_atomic_coord_perturbation = self.default_sigma_atomic_coord_perturbation
        elif structure_mut_params['sigma_atomic_coord_perturbation'] == None or structure_mut_params['sigma_atomic_coord_perturbation'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.sigma_atomic_coord_perturbation = self.default_sigma_atomic_coord_perturbation
        else:
            # otherwise, parse the value from the parameters
            self.sigma_atomic_coord_perturbation = structure_mut_params['sigma_atomic_coord_perturbation']
            
        # the maximum allowed magnitude of perturbation of each atomic coordinate
        if 'max_atomic_coord_perturbation' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.max_atomic_coord_perturbation = self.default_max_atomic_coord_perturbation
        elif structure_mut_params['max_atomic_coord_perturbation'] == None or structure_mut_params['max_atomic_coord_perturbation'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.max_atomic_coord_perturbation = self.default_max_atomic_coord_perturbation
        else:
            # otherwise, parse the value from the parameters
            self.max_atomic_coord_perturbation = structure_mut_params['max_atomic_coord_perturbation']
            
        # the standard deviation of the magnitude of the non-identity components of the elements in the strain matrix
        if 'sigma_strain_matrix_element' not in structure_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_strain_matrix_element = self.default_sigma_strain_matrix_element
        elif structure_mut_params['sigma_strain_matrix_element'] == None or structure_mut_params['sigma_strain_matrix_element'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.sigma_strain_matrix_element = self.default_sigma_strain_matrix_element
        else:
            # otherwise, parse the value from the parameters
            self.sigma_strain_matrix_element = structure_mut_params['sigma_strain_matrix_element']
        
    
    def doVariation(self, pool, random, geometry, id_generator):
        '''
        Performs the structural mutation operation
        
        Returns the resulting offspring as an Organism
        
         Args:
            pool: the Pool of Organisms
            
            random: Python's built in PRNG
            
            geometry: the Geometry object
            
            id_generator: the IDGenerator
            
        Description:
        
            Creates an offspring organism by perturbing the atomic positions and lattice vectors of the parent structure
            
                1. Selects a parent organism from the pool and makes a copy of it
                
                2. Perturbs the atomic coordinates of each site with probability self.frac_atoms_perturbed. The perturbation of each atomic
                   coordinate is drawn from a Gaussian with mean zero and standard deviation self.sigma_atomic_coord_perturbation. The 
                   magnitude of each atomic coordinate perturbation is constrained to not exceed self.max_atomic_coord_perturbation
                   
                3. The lattice vectors are perturbed by taking the product of each lattice vector with a strain matrix. The strain matrix is 
                   defined as
            
                        I + E
                
                   where I is the 3x3 identity matrix and E is a 3x3 perturbation matrix whose elements are distinct and drawn from a Gaussian 
                   with mean zero and standard deviation self.sigma_strain_matrix_element and are constrained to lie between -1 and 1 
        ''' 
        # just for testing, read in a structure from a file
        #structure = Structure.from_file('/n/srv/brevard/structures/POSCAR.NaCl')
        
        # select a parent organism from the pool
        parent_org = pool.selectOrganisms(1, random)
        
        # make a deep copy of the structure of the parent organism, so that the subsequent mutation doesn't affect the structure of the parent
        structure = copy.deepcopy(parent_org[0].structure)
        
        # for each site in the structure, determine whether to perturb it or not
        for site in structure.sites:
            if random.random() < self.frac_atoms_perturbed:
                # generate three random perturbation magnitudes (in Cartesian coordinates), one for each atomic coordinate, and make sure they aren't larger than the max allowed perturbation
                # the first one
                nudge_x = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_x) > self.max_atomic_coord_perturbation:
                    nudge_x = random.gauss(0, self.sigma_atomic_coord_perturbation)
                # the second one
                nudge_y = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_y) > self.max_atomic_coord_perturbation:
                    nudge_y = random.gauss(0, self.sigma_atomic_coord_perturbation)
                # the third one
                nudge_z = random.gauss(0, self.sigma_atomic_coord_perturbation)
                while np.abs(nudge_z) > self.max_atomic_coord_perturbation:
                    nudge_z = random.gauss(0, self.sigma_atomic_coord_perturbation)
                # make the perturbation vector
                perturbation_vector = [nudge_x, nudge_y, nudge_z]
                # translate the site be the computed perturbation vector
                structure.translate_sites(structure.sites.index(site), perturbation_vector, frac_coords=False, to_unit_cell=False)
            
        # compute the random non-identity components of the nine elements of the strain matrix
        epsilons = []
        for _ in range(9):
            epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            # make sure it's in [-1, 1]
            while np.abs(epsilon) > 1:
                epsilon = random.gauss(0, self.sigma_strain_matrix_element)
            epsilons.append(epsilon)
            
        # construct the strain matrix (I + epsilon_ij), and randomly assign positive or negative directions to the non-identity components
        row_1 = [1 + epsilons[0], epsilons[1], epsilons[2]]
        row_2 = [epsilons[3], 1 + epsilons[4], epsilons[5]]
        row_3 = [epsilons[6], epsilons[7], 1 + epsilons[8]]
        strain_matrix = np.array([row_1, row_2, row_3])
        
        # compute new lattice vectors by applying the strain matrix to the lattice vectors
        new_a = strain_matrix.dot(structure.lattice.matrix[0])
        new_b = strain_matrix.dot(structure.lattice.matrix[1])
        new_c = strain_matrix.dot(structure.lattice.matrix[2])
        new_lattice = Lattice([new_a, new_b, new_c])
        
        # assign the new lattice to the structure (this doesn't change the sites' fractional coordinates, but it does change their Cartesian coordinates)
        structure.modify_lattice(new_lattice)
        
        # make a new offspring organism
        offspring = Organism(structure, id_generator, self.name)
        
        # make sure all the site are within the cell (some of them could have been pushed outside by the atomic coordinate perturbations)
        offspring.translateAtomsIntoCell()
        
        # print out a message
        print("Creating offspring organism {} from parent organism {} with the structure mutation variation".format(offspring.id, parent_org[0].id))
        
        # return the offspring
        return offspring 
   
        
        
class NumStoichsMut(Variation):
    '''
    An operator that creates an offspring organism by mutating the number of stoichiometries' worth of atoms 
    in the parent organism.
    '''
    def __init__(self, num_stoichs_mut_params):
        '''
        Creates a NumStoichsMut operator
        
        Args:    
            num_stoichs_mut_params: The parameters for doing the NumStoichsMut operation, as a dictionary
            
        Precondition: the 'fraction' parameter in num_stoichs_mut_params is not optional, and it is assumed that num_stoichs_mut_params contains this parameter
        '''
        # the name of this variation
        self.name = 'number of stoichiometries mutation'
        
        # parse the fraction from the parameters. This argument is not optional, so it doesn't have a default value here
        self.fraction = num_stoichs_mut_params['fraction']
        
        # the default values
        # TODO: are these good default values?
        self.default_mu_num_adds = 0     # the average number of stoichimetries to add
        self.default_sigma_num_adds = 1  # the standard deviation of the number of stoichiometries to add
        self.default_scale_volume = True # whether to scale the volume of the offspring to equal that of the parent
        
        # the average number of stoichiometries to add
        if 'mu_num_adds' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.mu_num_adds = self.default_mu_num_adds
        elif num_stoichs_mut_params['mu_num_adds'] == None or num_stoichs_mut_params['mu_num_adds'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.mu_num_adds = self.default_mu_num_adds
        else:
            # otherwise, parse the value from the parameters
            self.mu_num_adds = num_stoichs_mut_params['mu_num_adds']
        
        # the standard deviation of the number of stoichiometries
        if 'sigma_num_adds' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.sigma_num_adds = self.default_sigma_num_adds
        elif num_stoichs_mut_params['sigma_num_adds'] == None or num_stoichs_mut_params['sigma_num_adds'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.sigma_num_adds = self.default_sigma_num_adds
        else:
            # otherwise, parse the value from the parameters
            self.sigma_num_adds = num_stoichs_mut_params['sigma_num_adds']
            
        # whether to scale the volume of the offspring
        if 'scale_volume' not in num_stoichs_mut_params:
            # use the default value if the flag hasn't been used
            self.scale_volume = self.default_scale_volume
        elif num_stoichs_mut_params['scale_volume'] == None or num_stoichs_mut_params['scale_volume'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.scale_volume = self.default_scale_volume
        else:
            # otherwise, parse the value from the parameters
            self.scale_volume = num_stoichs_mut_params['scale_volume']
        
    
    def doVariation(self, pool, random, geometry, id_generator):
        '''
        Performs the number of stoichiometries mutation operation
        
        Returns the resulting offspring as an Organism
        
        Args:
            pool: the Pool of Organisms 
            
            random: Python's built in PRNG
            
            geometry: the Geometry object
            
            id_generator: the IDGenerator
            
        Description:
        
            Creates an offspring organism by adding or removing a random number of stoichiometries' worth of atoms to or from the parent structure.
            
                1. Selects a parent organism from the pool and makes a copy of it
                
                2. Computes the number of stoichiometries to add or remove by drawing from a Gaussian with mean self.mu_num_adds and standard deviation
                   self.sigma_num_adds and rounding the result to the nearest integer
                   
                3. Computes the number of atoms of each type to add or remove, and does the additions or removals
                
                4. If self.scale_volume is True, scales the new structure to have the same volume per atom as the parent 
        '''
        # just for testing, read in a structure from a file
        #structure = Structure.from_file('/n/srv/brevard/structures/POSCAR.NaCl')
        #parent_org = Organism(structure, id_generator)
        
        # select a parent organism from the pool
        parent_org = pool.selectOrganisms(1, random)
        
        # make a deep copy of the structure of the parent organism
        structure = copy.deepcopy(parent_org[0].structure)
        # get the total number of atoms in the parent
        parent_num_atoms = len(structure.sites)
        # get the reduced composition of the parent
        reduced_composition = structure.composition.reduced_composition
        # get the volume per atom of the parent
        vol_per_atom = structure.lattice.volume/len(structure.sites)
        
        # loop needed here in case no atoms got added or removed, in which case need to try again. 
        # this happens if the randomly chosen number of atoms to remove exceeds the number of atoms in the cell
        while len(structure.sites) == parent_num_atoms:
        
            # compute the number of stoichiometries to add or remove, and make sure it's not zero
            num_add = int(round(random.gauss(self.mu_num_adds, self.sigma_num_adds)))
            while num_add == 0:
                num_add = int(round(random.gauss(self.mu_num_adds, self.sigma_num_adds)))
            
            # compute the number of each type of atom to add (or remove), and also the total number of atoms to add or remove
            amounts_to_add = {}
            total_add = 0
            for key in reduced_composition:
                amounts_to_add[key] = int(num_add*reduced_composition[key])
                total_add = total_add + int(num_add*reduced_composition[key])
        
            # if num_add is positive, put the new atoms in the cell at random locations
            if num_add > 0:
                for key in amounts_to_add:
                    for _ in range(amounts_to_add[key]):
                        frac_coords = [random.random(), random.random(), random.random()]
                        structure.append(Specie(key, 0), frac_coords) 
                # to remove the oxidation state of 0 we had to specify above in the structure.append() method
                structure.remove_oxidation_states()
                # to make sure all the atoms of the same element are listed consecutively, so they get printed consecutively in the poscar file
                structure.sort()
        
            # if num_adds is negative and the structure contains enough atoms, randomly remove the right number of atoms of each element 
            elif num_add < 0 and -1*total_add < len(structure.sites):
                site_indices_to_remove = [] # a list to hold the indices of the sites to remove
                for key in amounts_to_add:
                    for _ in range(0, -1*amounts_to_add[key]):
                        # pick a random site in the structure that has the right element and that hasn't already been picked
                        random_site = random.choice(structure.sites)
                        while str(random_site.specie.symbol) != str(key) or structure.sites.index(random_site) in site_indices_to_remove:
                            random_site = random.choice(structure.sites)
                        # record the index to of the site that has been designated for removal
                        site_indices_to_remove.append(structure.sites.index(random_site))
                # remove the chosen sites
                structure.remove_sites(site_indices_to_remove)
        
        # optionally scale the volume after the atoms have been added or removed
        if self.scale_volume:
            structure.scale_lattice(vol_per_atom*len(structure.sites))
        
        # create a new organism from the structure and return it
        offspring = Organism(structure, id_generator, self.name)
        
        # print out a message
        print("Creating offspring organism {} from parent organism {} with the number of stoichiometries mutation variation".format(offspring.id, parent_org[0].id))
        
        return offspring
        
        
        
class Permutation(Variation):
    '''
    A permutation operator.
    
    This class is a singleton.
    '''
    def __init__(self, permutation_params, composition_space):
        '''
        Creates a permutation operator
        
        Args:
            permutation_params: The parameters for doing the permutation operation, as a dictionary
            
            composition_space: the CompositionSpace object
            
        Precondition: the 'fraction' parameter in permutation_params is not optional, and it is assumed that permutation_params contains this parameter
        '''
        # the name of this variation
        self.name = 'permutation'
        
        # parse the fraction from the parameters. This argument is not optional, so it doesn't have a default value here
        self.fraction = permutation_params['fraction']
        
        # the max number of times to try getting a parent organism from the pool that can undergo at least one swap
        # this could be a problem if, e.g., the composition space is a single pure element but the Permutation variation has non-zero fraction
        self.max_num_selects = 1000
        
        # the default values
        # TODO: are these good default values?
        self.default_mu_num_swaps = 2                                    # the average number of pairs to swap
        self.default_sigma_num_swaps = 1                                 # the standard deviation of pairs to swap
        self.default_pairs_to_swap = composition_space.get_all_pairs()   # which atomic pairs to swap
        
        # the average number of swaps
        if 'mu_num_swaps' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.mu_num_swaps = self.default_mu_num_swaps
        elif permutation_params['mu_num_swaps'] == None or permutation_params['mu_num_swaps'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.mu_num_swaps = self.default_mu_num_swaps
        else:
            # otherwise, parse the value from the parameters
            self.mu_num_swaps = permutation_params['mu_num_swaps']
            
        # the standard deviation of the number of swaps
        if 'sigma_num_swaps' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.sigma_num_swaps = self.default_sigma_num_swaps
        elif permutation_params['sigma_num_swaps'] == None or permutation_params['sigma_num_swaps'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.sigma_num_swaps = self.default_sigma_num_swaps
        else:
            # otherwise, parse the value from the parameters
            self.sigma_num_swaps = permutation_params['sigma_num_swaps']
            
        # which pairs of atoms to swap
        if 'pairs_to_swap' not in permutation_params:
            # use the default value if the flag hasn't been used
            self.pairs_to_swap = self.default_pairs_to_swap
        elif permutation_params['pairs_to_swap'] == None or permutation_params['pairs_to_swap'] == 'default':
            # use the default value if the flag was left blank or set to 'default'
            self.pairs_to_swap = self.default_pairs_to_swap
        else:
            # otherwise, parse the value from the parameters
            self.pairs_to_swap = permutation_params['pairs_to_swap']
    
    
    def doVariation(self, pool, random, geometry, id_generator):
        '''
        Performs the permutation operation
        
        Returns the resulting offspring as an Organism, or None if no offspring could be created
        
        Args:
            pool: the Pool of Organisms
        
            random: Python's built in PRNG
            
            geometry: the Geometry object
            
            id_generator: the IDGenerator
            
        Description:
            
            Creates and offspring organism by swapping the elements of some of the sites in the parent structure.
            
                1. Selects a parent organism from the pool that is able to have at least one of the allowed swaps done on it
                
                2. Computes the number of swaps to try to do by drawing from a Gaussian with mean self.mu_num_swaps and standard
                   deviation self.sigma_num_swaps and rounding to the nearest integer
                   
                3. Tries to do the computed number of allowed swaps by randomly selecting an allowed pair to swap and then randomly 
                   selecting sites in the structure with elements of the allowed pair. This is repeated until either the computed 
                   number of swaps have been done or no more swaps are possible with the parent structure
                    
        '''
        # just for testing, read in a structure from a file
        #structure = Structure.from_file('/n/srv/brevard/structures/POSCAR.AlFeCo2')
        #parent_org = Organism(structure, id_generator)
        
        # select a parent organism from the pool
        parent_org = pool.selectOrganisms(1, random)
        # make a deep copy of the structure of the parent organism
        structure = copy.deepcopy(parent_org[0].structure)
        # keep trying until we get a parent that has at least one possible swap
        possible_swaps = self.getPossibleSwaps(structure)
        num_selects = 0
        while len(possible_swaps) == 0 and num_selects < self.max_num_selects:
            # select a parent organism from the pool
            parent_org = pool.selectOrganisms(1, random)
            # make a deep copy of the structure of the parent organism
            structure = copy.deepcopy(parent_org.structure)
            possible_swaps = self.getPossibleSwaps(structure)
            num_selects = num_selects + 1
            
        # if the maximum number of selections have been made, then this isn't working and it's time to stop
        if num_selects >= self.max_num_selects:
            return None
            
        # compute a positive random number of swaps to do
        num_swaps = int(round(random.gauss(self.mu_num_swaps, self.sigma_num_swaps)))
        while num_swaps <= 0:
            num_swaps = int(round(random.gauss(self.mu_num_swaps, self.sigma_num_swaps)))
            
        # try to select the computed number of swaps
        num_swaps_selected = 0 # how many swaps have been selected to do 
        pair_indices = [] # list hold the indices of the pair of atoms in each swap
        structure_to_check = copy.deepcopy(structure) # copy that we can remove atoms from to recompute possible swaps
        # keep getting more swaps until either we've got enough or no more swaps are possible
        while num_swaps_selected < num_swaps and len(possible_swaps) > 0:
            # pick a random pair to swap that we know is possible
            swap = random.choice(possible_swaps)
            symbols = swap.split()
            # keep trying until we find a site with the first element in the pair
            site_1 = random.choice(structure_to_check.sites)
            while str(site_1.specie.symbol) != symbols[0]:
                site_1 = random.choice(structure_to_check.sites)
            # keep trying until we find a site with the second element in the pair
            site_2 = random.choice(structure_to_check.sites)
            while str(site_2.specie.symbol) != symbols[1]:
                site_2 = random.choice(structure_to_check.sites)
            # record the indices (w.r.t. to the unchanged structure) for this pair of sites
            pair_index = [structure.sites.index(site_1), structure.sites.index(site_2)]
            pair_indices.append(pair_index)
            num_swaps_selected = num_swaps_selected + 1
            # remove these two sites from the structure to check 
            structure_to_check.remove_sites([structure_to_check.sites.index(site_1), structure_to_check.sites.index(site_2)])
            # update the possible swaps
            possible_swaps = self.getPossibleSwaps(structure_to_check)
        
        # do the swaps with the selected pairs
        for pair_index in pair_indices:
            species_1 = structure.sites[pair_index[0]].specie
            species_2 = structure.sites[pair_index[1]].specie
            structure.replace(pair_index[0], species_2)
            structure.replace(pair_index[1], species_1)
        
        # make a new organism from the structure and return it
        offspring = Organism(structure, id_generator, self.name)
        
        # print out a message
        print("Creating offspring organism {} from parent organism {} with the permutation variation".format(offspring.id, parent_org[0].id))
        
        return offspring
        
        
    def getPossibleSwaps(self, structure):
        '''
        Returns a list of swaps that are possible to do, based on what atoms are in the cell and which pairs are list in self.pairs_to_swap.
        The returned list is a sublist of self.pairs_to_swap. Does not change the structure.
        
        Args:
            structure: the Structure object to check.
        '''
        possible_pairs = [] # list to hold the possible pairs that could be swapped in the structure
        for pair in self.pairs_to_swap:
            symbols = pair.split()
            # check if the structure contains the elements in the chosen pair
            has_element_1 = False
            has_element_2 = False
            for site in structure.sites:
                if str(site.specie.symbol) == symbols[0]:
                    has_element_1 = True
                if str(site.specie.symbol) == symbols[1]:
                    has_element_2 = True
            # append the pair if both elements were found in the structure
            if has_element_1 and has_element_2:
                possible_pairs.append(pair)
        return possible_pairs
        
        
        
class Geometry(object):
    '''
    Represents the geometry data, including any geometry-specific constraints (max_size, etc.)
    ''' 
    def __init__(self, geometry_parameters):
        '''
        Creates a geometry object
        
        Args:
            geometry_parameters: a dictionary of parameters
        '''
        # default values
        self.default_shape = 'bulk'
        self.default_max_size = inf
        self.default_min_size = -inf
        self.default_padding = 10 # this will only be used for non-bulk shapes
        
        # if entire Geometry block was set to default or left blank
        if geometry_parameters == None or geometry_parameters == 'default':
            self.shape = self.default_shape
            self.max_size = self.default_max_size
            self.min_size = self.default_min_size
            self.padding = None
        else:     
            # check each one and see if it's been left blank or set to default, or not included at all 
            if 'shape' in geometry_parameters:
                # if no shape was given, assume bulk and set everything to default
                if geometry_parameters['shape'] == None or geometry_parameters['shape'] == 'default':
                    self.shape = self.default_shape
                    self.max_size = self.default_max_size
                    self.padding = None
                else:
                    self.shape = geometry_parameters['shape']
                    # set max size, and check if was left blank or set to None or 'default'
                    if 'max_size' in geometry_parameters:
                        if geometry_parameters['max_size'] == None or geometry_parameters['max_size'] == 'default':
                            self.max_size = self.default_max_size
                        else:
                            self.max_size = geometry_parameters['max_size']
                    else:
                        self.max_size = self.default_max_size
                    # set min size, and check if was left blank or set to None or 'default'
                    if 'min_size' in geometry_parameters:
                        if geometry_parameters['min_size'] == None or geometry_parameters['min_size'] == 'default':
                            self.min_size = self.default_min_size
                        else:
                            self.min_size = geometry_parameters['min_size']
                    else:
                        self.min_size = self.default_min_size
                    # set padding, and check if was left blank or set to None or 'default'
                    if 'padding' in geometry_parameters:
                        if geometry_parameters['padding'] == None or geometry_parameters['padding'] == 'default':
                            self.padding = self.default_padding
                        else:
                            self.padding = geometry_parameters['padding']
                    else:
                        self.padding = self.default_padding
            # if shape field was missing, assume bulk and set default values
            else:
                self.shape = self.default_shape
                self.max_size = self.default_max_size
                self.min_size = self.default_min_size
                self.padding = None          
     
                
    def pad(self, organism):
        '''
        Makes an organism's structure conform to the required shape. For sheet, wire and cluster geometries, this 
        means adding vacuum padding to the cell. For bulk, the structure is unchanged. Used to pad a structure prior to 
        an energy calculation. Changes the structure of an organism.
        
        Args:
            organism: the Organism to pad
        '''
        # call other methods based on the value of self.shape
        if self.shape == 'sheet':
            self.padSheet(organism)
        elif self.shape == 'wire':
            self.padWire(organism)
        elif self.shape == 'cluster':
            self.padCluster(organism)
     
        
    def padSheet(self, organism):
        '''
        Adds vertical vacuum padding to a sheet, and makes the c-lattice vector normal to the plane of the sheet. 
        The atoms are shifted up to the center of the padded sheet. Changes an organism's structure
        
        Args:
            organism: an Organism object
        ''' 
        # rotate into principal directions
        organism.rotateToPrincipalDirections()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the layer thickness of the sheet
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        minz = cart_bounds[2][0]
        maxz = cart_bounds[2][1]
        layer_thickness = maxz - minz
        # get the non-zero components of the a and b lattice vectors
        ax = organism.structure.lattice.matrix[0][0]
        bx = organism.structure.lattice.matrix[1][0]
        by = organism.structure.lattice.matrix[1][1]
        # make a new lattice with c vertical and with length layer_thickness + padding
        padded_lattice = Lattice([[ax, 0.0, 0.0], [bx, by, 0.0], [0.0, 0.0, layer_thickness + self.padding]])
        # make a new structure with the padded lattice, species, and Cartesian coordinates
        padded_structure = Structure(padded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        organism.structure = padded_structure
        # shift the atoms vertically so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
        
        
    def padWire(self, organism):
        '''
        Makes c lattice vector parallel to z-axis, and adds vacuum padding around a wire in the x and y directions by replacing a and b lattice vectors with padded vectors along 
        the x and y axes, respectively. The atoms are shifted to the center of the padded cell. Changes the structure of an organism.
        
        Args:
            organism: the organism to pad 
        '''
        # rotate c parallel to z-axis
        organism.rotateCParallelToZ()
        # get the species and Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the wire in the x and y directions
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        # get the non-zero component of the c lattice vector
        cz = organism.structure.lattice.matrix[2][2]
        # make a new lattice with padding in the x and y directions
        padded_lattice = Lattice([[x_extent + self.padding, 0, 0], [0, y_extent + self.padding, 0], [0, 0, cz]])
        # make a new structure with the padded lattice, species, and Cartesian coordinates
        padded_structure = Structure(padded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        organism.structure = padded_structure
        # shift the atoms horizontally so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
    
    
    def padCluster(self, organism):
        '''
        Adds vacuum padding around a cluster. Changes the structure of an organism.
        
        Args:
            organism: the organism to pad 
        '''
        # get the species and Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the cluster in each direction
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        z_min = cart_bounds[2][0]
        z_max = cart_bounds[2][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min
        # make a new orthorhombic lattice
        padded_lattice = Lattice([[x_extent + self.padding, 0, 0], [0, y_extent + self.padding, 0], [0, 0, z_extent + self.padding]])
        # make a new structure with the padded lattice, species, and Cartesian coordinates
        padded_structure = Structure(padded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        organism.structure = padded_structure
        # shift all the atoms so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
    
    
    def unpad(self, organism, constraints):
        '''
        Removes vacuum padding to return an organism's structure to a form used by the variations. Does nothing if shape is bulk. 
        Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad
            
            constraints: a Constraints object
        '''
        # call other methods based on the value of self.shape
        if self.shape == 'sheet':
            self.unpadSheet(organism, constraints)
        elif self.shape == 'wire':
            self.unpadWire(organism, constraints)
        elif self.shape == 'cluster':
            self.unpadCluster(organism, constraints)
        
    
    def unpadSheet(self, organism, constraints):
        '''
        Removes vertical vacuum padding from a sheet, leaving only enough to satisfy the per-species MID constraints, and makes the c-lattice vector normal to 
        the plane of the sheet (if it isn't already). Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad 
            
            constraints: a Constraints object
        '''
        # rotate into principal directions
        organism.rotateToPrincipalDirections()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the layer thickness of the sheet
        layer_thickness = self.getLayerThickness(organism)
        # get the largest per-species MID (add just a little to prevent corner cases...)
        max_mid = constraints.get_max_mid() + 0.01
        # get the non-zero components of the a and b lattice vectors
        ax = organism.structure.lattice.matrix[0][0]
        bx = organism.structure.lattice.matrix[1][0]
        by = organism.structure.lattice.matrix[1][1]
        # make a new lattice with c vertical and with length layer_thickness + padding
        unpadded_lattice = Lattice([[ax, 0.0, 0.0], [bx, by, 0.0], [0.0, 0.0, layer_thickness + max_mid]])
        # make a new structure with the unpadded lattice, species, and Cartesian coordinates
        unpadded_structure = Structure(unpadded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        # set the organism's structure to the unpadded one
        organism.structure = unpadded_structure
        # shift the atoms vertically a little so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
        
          
    def unpadWire(self, organism, constraints):
        '''
        Removes vacuum padding around a wire. Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad
            
            constraints: a Constraints object
        '''
        # rotate c parallel to z-axis
        organism.rotateCParallelToZ()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the wire in the x and y directions
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        # get the non-zero component of the c lattice vector
        cz = organism.structure.lattice.matrix[2][2]
        # get the largest per-species MID (add just a little to prevent corner cases...)
        max_mid = constraints.get_max_mid() + 0.01
        # make a new orthorhombic lattice with the vacuum removed
        unpadded_lattice = Lattice([[x_extent + max_mid, 0.0, 0.0], [0, y_extent + max_mid, 0.0], [0.0, 0.0, cz]])
        # make a new structure with the unpadded lattice, species, and Cartesian coordinates
        unpadded_structure = Structure(unpadded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        # set the organism's structure to the unpadded one
        organism.structure = unpadded_structure
        # shift the atoms horizontally so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
 
        
    def unpadCluster(self, organism, constraints):
        '''
        Removes vacuum padding around a cluster, leaving only enough to satisfy the per-species MID constraints. 
        Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad
            
            constraints: a Constraints object
        '''
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the cluster in each direction
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        x_min = cart_bounds[0][0]
        x_max = cart_bounds[0][1]
        y_min = cart_bounds[1][0]
        y_max = cart_bounds[1][1]
        z_min = cart_bounds[2][0]
        z_max = cart_bounds[2][1]
        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min
        # get the largest per-species MID (add just a little to prevent corner cases...)
        max_mid = constraints.get_max_mid() + 0.01
        # make a new orthorhombic lattice with the vacuum removed
        unpadded_lattice = Lattice([[x_extent + max_mid, 0.0, 0.0], [0, y_extent + max_mid, 0.0], [0.0, 0.0, z_extent + max_mid]])
        # make a new structure with the padded lattice, species, and Cartesian coordinates
        unpadded_structure = Structure(unpadded_lattice, species, cartesian_coords, coords_are_cartesian=True)
        # set the organism's structure to the unpadded one
        organism.structure = unpadded_structure
        # shift the atoms a little so they're in the center
        frac_bounds = organism.getBoundingBox(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translateAtomsIntoCell()
        
        
    def getSize(self, organism):
        '''
        Returns the size of an organism with a non-bulk shape. Returns 0 if the shape is bulk.
        
        Args:
            organism: the Organism whose size to get
        '''
        # call other methods based on the value of self.shape
        if self.shape == 'bulk':
            return 0
        elif self.shape == 'sheet':
            return self.getLayerThickness(organism)
        elif self.shape == 'wire':
            return self.getWireDiameter(organism)
        elif self.shape == 'cluster':
            return self.getClusterDiameter(organism)
       
        
    def getLayerThickness(self, organism):
        '''
        Returns the layer thickness of a sheet structure, which is the maximum vertical distance between atoms in the cell.
        
        Assumes that the organism has already been rotated into the principal directions, and that plane of the sheet is parallel to the a-b facet.
        '''
        cart_bounds = organism.getBoundingBox(cart_coords = True)
        layer_thickness = cart_bounds[2][1] - cart_bounds[2][0]
        return layer_thickness
    
    
    def getWireDiameter(self, organism):
        '''
        Returns the diameter of a wire structure, defined as the maximum distance between atoms projected to the x-y plane.
        
        Assumes that the organism has already been put into wire format (c lattice vector is parallel to z-axis), and that all sites are located inside the cell (i.e., have 
        fractional coordinates between 0 and 1). Generally called after Geometry.unpadWire has already been called on the organism.
        ''' 
        max_distance = 0
        for site_i in organism.structure.sites:
            # have to make Site versions of each PeriodicSite so that the computed distance won't include periodic images, and to project each site to the x-y plane
            non_periodic_site_i = Site(site_i.species_and_occu, [site_i.coords[0], site_i.coords[1], 0.0])
            for site_j in organism.structure.sites:
                non_periodic_site_j = Site(site_j.species_and_occu, [site_j.coords[0], site_j.coords[1], 0.0])
                distance = non_periodic_site_i.distance(non_periodic_site_j)
                if distance > max_distance:
                    max_distance = distance
        return max_distance 
        
    
    def getClusterDiameter(self, organism):
        '''
        Returns the diameter of a cluster structure, defined as the maximum distance between atoms in the cell
        
        Assumes that all sites are located inside the cell (i.e., have fractional coordinates between 0 and 1). Generally called after
        Geometry.unpadCluster has already been called on the organism.
        '''
        max_distance = 0
        for site_i in organism.structure.sites:
            # have to make Site versions of each PeriodicSite so that the computed distance won't include periodic images
            non_periodic_site_i = Site(site_i.species_and_occu, site_i.coords)
            for site_j in organism.structure.sites:
                non_periodic_site_j = Site(site_j.species_and_occu, site_j.coords)
                distance = non_periodic_site_i.distance(non_periodic_site_j)
                if distance > max_distance:
                    max_distance = distance
        return max_distance 
        
       
       
class CompositionSpace(object):
    '''
    Represents the composition space to be searched by the algorithm.
    '''
    def __init__(self, endpoints):
        '''
        Creates a CompositionSpace object, which is list of pymatgen.core.composition.Composition objects
        
        Args:
            endpoints: the list dictionaries mapping element symbols to amounts, with each one representing a composition
        '''
        for i in range(len(endpoints)):
            endpoints[i] = Composition(endpoints[i])
            endpoints[i] = endpoints[i].reduced_composition()
            
        self.endpoints = endpoints
        
        # for now, let's have the objective live here
        self.objective_function = self.inferObjectiveFunction()
    
    
    def inferObjectiveFunction(self):
        '''
        Infers the objective function (energy per atom or phase diagram) based on the composition space
        
        Returns either "epa" or "pd"
        '''
        # if only one composition, then it must be an epa search
        if len(self.endpoints) == 1:
            return "epa"
        # otherwise, compare all the compositions and see if any of them are different
        else:
            for point in self.endpoints:
                for next_point in self.endpoints:
                    if not point.almost_equals(next_point, 0.0, 0.0):
                        return "pd"
        # should only get here if there are multiple identical compositions in end_points (which would be weird)
        return "epa" 
    
    
    def get_all_elements(self):
        '''
        Returns a list of all the elements (as pymatgen.core.periodic_table.Element objects) that are in the composition space
        '''
        # get each element type from the composition_space object
        elements = []
        for point in self.endpoints:
            for key in point:
                elements.append(key)
        # remove duplicates from the list of elements
        elements = list(set(elements)) 
        return elements  
            
            
    def get_all_pairs(self):
        '''
        Returns all possible pairs of elements in the composition space, as list of strings, where each string contains the symbols of two elements, separated by a space. 
        
        Does not include self-pairs (e.g., "Cu Cu")
        '''
        # get all the Element objects
        elements = self.get_all_elements()
        # if only one type of element, then no pairs, so return an empty list
        if len(elements) == 1:
            return []
        # list to hold the pairs of symbols
        pairs = []
        # get all the possible distinct pairs
        for i in range(0, len(elements) - 1):
            for j in range(i + 1, len(elements)):
                pairs.append(str(elements[i].symbol + " " + elements[j].symbol))
        return pairs
                
        
        
class OrganismCreator(object):
    '''
    Creates organisms for the initial population
    
    Not meant to be instantiated, but rather subclassed by particular Creators, like RandomOrganismCreator
    or PoscarsOrganismCreator.
    
    TODO: is this even necessary? All it specifies is that a creator should have a createOrganism method that returns an organism or None. 
    '''
    
    def createOrganism(self):
        '''
        Creates an organism for the initial population.
        
        Returns an organism, or None if one could not be created
        '''
        raise NotImplementedError("Please implement this method.")
        


class RandomOrganismCreator(OrganismCreator):
    '''
    Creates random organisms for the initial population
    '''
    def __init__(self, random_org_parameters, composition_space):
        '''
        Creates a RandomOrganismCreator.
        
        Args:
            random_org_parameters: the parameters for generating random organisms
            
            composition_space: a CompositionSpace object   
            
            random: Python's PRNG
        '''
        # the name of this creator
        self.name = 'random organism creator'
        # the default number of random organisms to make
        # TODO: are these good values?
        if composition_space.objective_function == 'epa':
            self.default_number = 28
        elif composition_space.objective_function == 'pd':
            self.default_number = 40
        # the default volume scaling behavior
        self.default_volume = 'from_elemental_densities'
        
        # if entire random_org_parameters is None or 'default', then set to defaults
        if random_org_parameters == None or random_org_parameters == 'default':
            self.number = self.default_number
            self.volume = self.default_volume
        
        # otherwise, parse the parameters and set to defaults if necessary
        else:
            if 'number' in random_org_parameters:
                if random_org_parameters['number'] == None or random_org_parameters['number'] == 'default':
                    self.number = self.default_number
                else:
                    self.number = random_org_parameters['number']
            else:
                # if no 'number' tag, then just use the default
                self.number = self.default_number
            
            # get the volume to scale them to 
            if 'volume' in random_org_parameters:
                if random_org_parameters['volume'] == None or random_org_parameters['volume'] == 'default':
                    self.volume = self.default_volume
                else:
                    self.volume = random_org_parameters['volume']
            else:
                # if no 'volume' tag given, then just do the default
                self.volume = self.default_volume
                
        # variables to keep track of how many have been made, when to stop, and if this creator is finished   
        # for a random organism creator, num_made is defined as the number of organisms made that have been added to the initial population 
        self.num_made = 0
        self.is_successes_based = True
        self.is_finished = False
        
        # For elements that are normally gases, pymatgen doesn't have densities. These are the densities (in kg/m^3) for the ones that are 
        # missing from pymatgen's database. All values are taken from Materials Project, except for At and Rn, which were taken from Wikipedia
        self.missing_densities = {'H': 150, 'He': 421, 'N': 927, 'O': 1974, 'F': 1972, 'Ne': 1680, 'Cl': 829, 'Ar': 1478, 'Br': 3535, 'Kr': 2030, 'Xe': 2611, 'At': 6350, 'Rn': 4400}
    
    
    def createOrganism(self, id_generator, composition_space, constraints, random):
        '''
        Creates a random organism for the initial population.
        
        Returns a random organism, or None if an error was encountered during volume scaling
        
        Args:
            id_generator: an IDGenerator object
            
            composition_space: a CompositionSpace object
            
            constraints: a Constraints object 
            
            random: copy of Python's built in PRNG
        '''
        # make three random lattice vectors that satisfy the length constraints
        a = constraints.min_lattice_length + random.random()*(constraints.max_lattice_length - constraints.min_lattice_length)
        b = constraints.min_lattice_length + random.random()*(constraints.max_lattice_length - constraints.min_lattice_length)
        c = constraints.min_lattice_length + random.random()*(constraints.max_lattice_length - constraints.min_lattice_length)
        
        # make three random lattice angles that satisfy the angle constraints
        alpha = constraints.min_lattice_angle + random.random()*(constraints.max_lattice_angle - constraints.min_lattice_angle)
        beta = constraints.min_lattice_angle + random.random()*(constraints.max_lattice_angle - constraints.min_lattice_angle)
        gamma = constraints.min_lattice_angle + random.random()*(constraints.max_lattice_angle - constraints.min_lattice_angle)
        
        # build the random lattice
        random_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        
        # get a list of elements for the random organism
        if composition_space.objective_function == 'epa':
            reduced_formula = composition_space.endpoints[0].reduced_composition
            num_atoms_in_formula = reduced_formula.num_atoms
            max_num_formulas = int(constraints.max_num_atoms/num_atoms_in_formula)
            min_num_formulas = int(constraints.min_num_atoms/num_atoms_in_formula)
            if min_num_formulas == 0:
                min_num_formulas = 1
            # get a random number of formula units, and the resulting random number of atoms
            random_num_formulas = random.randint(min_num_formulas, max_num_formulas)
            num_atoms = int(random_num_formulas*num_atoms_in_formula)
            # add the right number of each element
            elements = []
            for element in reduced_formula:
                # for some reason, reduced_formula[element] is a float for elemental structures, so have to cast it to an int below
                for _ in range(random_num_formulas*int(reduced_formula[element])):
                    elements.append(element) 
        elif composition_space.objective_function == 'pd':
            # TODO: this doesn't ensure the organism will be in the composition space. If it's not, it will just fail development, but there might be a better way...
            num_atoms = random.randint(constraints.min_num_atoms, constraints.max_num_atoms)
            allowed_elements = constraints.get_all_elements(composition_space)
            elements = []
            for _ in range(num_atoms):
                elements.append(random.choice(allowed_elements))
        
        # for each element, generate a set of random fractional coordinates
        # TODO: this doesn't ensure the structure will satisfy the per-species mids, and in fact most won't. It's ok because they'll just fail development, but there might be a better way...
        # TODO: also, this doesn't ensure the structure will satisfy the max size constraint in Geometry. There could be a way to do this by applying more stringent constraints on how
        #       the random lattice vectors and angles are chosen when we have a non-bulk geometry....
        random_coordinates = []
        for _ in range(num_atoms):
            random_coordinates.append([random.random(), random.random(), random.random()])
        
        # make a random structure from the random lattice, random species, and random coordinates
        random_structure = Structure(random_lattice, elements, random_coordinates)
        
        # optionally scale the volume to the weighted average of the volumes per atom of the elemental structures
        if self.volume == 'from_elemental_densities':
            # compute volumes per atom (in Angstrom^3) of each element in the random organism
            reduced_composition = random_structure.composition.reduced_composition
            volumes_per_atom = {}
            for element in reduced_composition:
                # physical properties and conversion factors
                atomic_mass = float(element.atomic_mass) # in amu
                mass_conversion_factor = 1.660539040e-27 # converts amu to kg
                density = element.density_of_solid # in kg/m^3
                # if the density is not in pymatgen's database, then get it from the missing densities
                if density == None:
                    density = self.missing_densities[element.symbol]
                density = float(density)
                length_conversion_factor = 1.0e10 # converts meters to Angstroms 
            
                # compute the volume (in Angstrom^3) per atom of this element
                # take the log of the product to prevent numerical issues
                log_volume_per_atom = np.log(mass_conversion_factor) + np.log(atomic_mass) - np.log(density) + 3.0*np.log(length_conversion_factor)
                volume_per_atom = np.exp(log_volume_per_atom)                
                volumes_per_atom[element] = volume_per_atom
        
            # compute the desired volume per atom by taking the weighted average of the volumes per atom of the constituent elements
            weighted_sum = 0
            for element in reduced_composition:
                # the number of this type of element times it's calculated volume per atom
                weighted_sum = weighted_sum + reduced_composition[element]*volumes_per_atom[element]
        
            # normalize by the total number of atoms to get the weighted average
            mean_vpa = weighted_sum/reduced_composition.num_atoms
        
            # scale the volume of the random organism to satisfy the computed mean volume per atom
            # The structure.scale_lattice method sometimes throws divide-by-zero runtime warnings. To keep these from getting printed to the output, we're temporarily
            # suppressing warnings here. When the divide-by-zero error happens, the lattice either gets scaled to a huge volume, or else scaling fails and the lattice 
            # vectors are set to NaN. After doing the scaling, the first lattice vector is checked to see if either of these error occured.
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                # check if the volume scaling worked
                random_structure.scale_lattice(mean_vpa*len(random_structure.sites))
                if str(random_structure.lattice.a) == 'nan' or random_structure.lattice.a > 100:
                    return None          
        
        elif self.volume == 'random':
            # no volume scaling
            pass
        
        else:
            # scale to the given volume per atom
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                random_structure.scale_lattice(self.volume*len(random_structure.sites))
                # check if the volume scaling worked 
                if str(random_structure.lattice.a) == 'nan' or random_structure.lattice.a > 100:
                    return None 
        
        # return a random organism with the scaled random structure
        random_org = Organism(random_structure, id_generator, self.name)
        print('Random organism creator making organism {}'.format(random_org.id))
        return random_org
    
    
    def updateStatus(self):
        '''
        Increments num_made, and if necessary, updates is_finished
        '''
        self.num_made = self.num_made + 1
        print('Organisms left for {}: {}'.format(self.name, self.number - self.num_made))
        if self.num_made == self.number:
            self.is_finished = True
                


class FileOrganismCreator(OrganismCreator):
    '''
    Creates organisms from files (poscar or cif) for the initial population.
    '''
    def __init__(self, path_to_folder):
        '''
        Creates a FileOrganismCreator.
        
        Args:
            path_to_folder: the path to the folder containing the files from which to make organisms
                            Precondition: the folder exists and contains files
        '''
        # the name of this creator
        self.name = 'file organism creator'
        # all the files in the given folder
        self.path_to_folder = path_to_folder
        self.files = [f for f in listdir(self.path_to_folder) if isfile(join(self.path_to_folder, f))]
        self.number = len(self.files)
      
        # variables to keep track of how many have been made, when to stop, and if this creator is finished   
        # for a file organism creator, num_made is defined as the number of attempts to make organisms from files (usually the number of files provided)
        self.num_made = 0
        self.is_successes_based = False
        self.is_finished = False
    
    
    def createOrganism(self, id_generator, composition_space, constraints, random):
        '''
        Creates an organism for the initial population from a poscar or cif file. 
        
        Returns an organism, or None if one could not be created
        
        Args:
            id_generator: an IDGenerator object
             
            composition_space: a CompositionSpace object
            
            constraints: a Constraints object
            
            random: Python's built in PRNG
            
        TODO: the last three arguments are never actually used in this method, but I included them so the method has the same arguments as RandomOrganismCreator.creatorOrganism()
              to allow the createOrganism method to be called on both RandomOrganismCreator and FileOrganismCreator without having to know in advance which one it is.
              Maybe there's a better way to deal with this...
        '''
        # TODO: This is kind of annoying. Maybe contact pymatgen and ask if they can add support for files ending in .POSCAR instead of only files starting with POSCAR 
        if self.files[self.num_made - 1].endswith('.cif') or self.files[self.num_made - 1].startswith('POSCAR'):
            try:
                new_struct = Structure.from_file(str(self.path_to_folder) + "/" + str(self.files[self.num_made - 1]))
                new_org = Organism(new_struct, id_generator, self.name)
                print('Making organism {} from file: {}'.format(new_org.id, self.files[self.num_made - 1]))
                # update status each time the method is called, since this is an attempts-based creator
                self.updateStatus()
                return new_org
            # return None if a structure couldn't be read from a file
            except:
                print('Error reading structure from file: {}'.format(self.files[self.num_made - 1])) 
                self.updateStatus()
                return None
        else:
            print('File {} has invalid extension - file must end with .cif or begin with POSCAR'.format(self.files[self.num_made - 1]))
            self.updateStatus()
            return None
        
        
    def updateStatus(self):
        '''
        Increments num_made, and if necessary, updates is_finished
        '''
        self.num_made = self.num_made + 1
        print('Organisms left for {}: {}'.format(self.name, self.number - self.num_made))
        if self.num_made == len(self.files):
            self.is_finished = True
        


class RedundancyGuard(object):
    '''
    A redundancy guard.
    
    This is a singleton class.
    '''
    
    def __init__(self, redundancy_parameters):
        '''
        Creates a redundancy guard.
        
        TODO: I think the pymatgen structure matcher assumes periodic boundary conditions, but this might not make sense for all geometries...
        
        Args:
            redundancy parameters: a dictionary of parameters
        '''
        # TODO: are these sensible defaults?
        # default lattice length tolerance, in fractional coordinates (pymatgen uses 0.2 as default...)
        self.default_lattice_length_tol = 0.05
        # default lattice angle tolerance, in degrees (pymatgen uses 5 as default...)
        self.default_lattice_angle_tol = 2
        # default site tolerance, in fraction of average free length per atom (pymatgen uses 0.3 as default...)
        self.default_site_tol = 0.1
        # whether to transform to primitive cells before comparing
        self.default_use_primitive_cell = True
        # whether to check if structures are equivalent to supercells of each other
        self.default_attempt_supercell = True
        # the epa difference interval
        self.default_epa_diff = 0.0
        
        # parse the parameters, and set to defaults if necessary
        if redundancy_parameters == None or redundancy_parameters == 'default':
            self.set_all_to_defaults()
        else:
            # check each flag to see if it's been included, and if so, whether it has been set to default or left blank
            # lattice length tolerance
            if 'lattice_length_tol' in redundancy_parameters:
                if redundancy_parameters['lattice_length_tol'] == None or redundancy_parameters['lattice_length_tol'] == 'default':
                    self.lattice_length_tol = self.default_lattice_length_tol
                else:
                    self.lattice_length_tol = redundancy_parameters['lattice_length_tol']
            else:
                self.lattice_length_tol = self.default_lattice_length_tol
                
            # lattice angle tolerance
            if 'lattice_angle_tol' in redundancy_parameters:
                if redundancy_parameters['lattice_angle_tol'] == None or redundancy_parameters['lattice_angle_tol'] == 'default':
                    self.lattice_angle_to = self.default_lattice_angle_tol
                else:
                    self.lattice_angle_to = redundancy_parameters['lattice_angle_tol']
            else:
                self.lattice_angle_to = self.default_lattice_angle_tol
                
            # site tolerance
            if 'site_tol' in redundancy_parameters:
                if redundancy_parameters['site_tol'] == None or redundancy_parameters['site_tol'] == 'default':
                    self.site_tol = self.default_site_tol
                else:
                    self.site_tol = redundancy_parameters['site_tol']
            else:
                self.site_tol = self.default_site_tol
            
            # whether to use primitive cells
            if 'use_primitive_cell' in redundancy_parameters:
                if redundancy_parameters['use_primitive_cell'] == None or redundancy_parameters['use_primitive_cell'] == 'default':
                    self.use_primitive_cell = self.default_use_primitive_cell
                else:
                    self.use_primitive_cell = redundancy_parameters['use_primitive_cell']
            else:
                self.use_primitive_cell = self.default_use_primitive_cell
            
            # whether to try matching supercells
            if 'attempt_supercell' in redundancy_parameters:
                if redundancy_parameters['attempt_supercell'] == None or redundancy_parameters['attempt_supercell'] == 'default':
                    self.attempt_supercell = self.default_attempt_supercell
                else:
                    self.attempt_supercell = redundancy_parameters['attempt_supercell']
            else:
                self.attempt_supercell = self.default_attempt_supercell
                
            # d-value
            if 'epa_diff' in redundancy_parameters:
                if redundancy_parameters['epa_diff'] == None or redundancy_parameters['epa_diff'] == 'default':
                    self.epa_diff = self.default_epa_diff
                else:
                    self.epa_diff = redundancy_parameters['epa_diff']
            else:
                self.epa_diff = self.default_epa_diff
        
        # make the StructureMatcher object
        # The first False is to prevent the matcher from scaling the volumes, and the second False is to prevent subset matching
        self.structure_matcher = StructureMatcher(self.lattice_length_tol, self.site_tol, self.lattice_angle_to, self.use_primitive_cell, False, self.attempt_supercell, False, ElementComparator())
    
        
    def set_all_to_defaults(self):
        '''
        Sets all the redundancy parameters to default values
        '''
        self.lattice_length_tol = self.default_lattice_length_tol
        self.lattice_angle_to = self.default_lattice_angle_tol
        self.site_tol = self.default_site_tol
        self.use_primitive_cell = self.default_use_primitive_cell
        self.attempt_supercell = self.default_attempt_supercell
        self.epa_diff = self.default_epa_diff
    
        
    def checkRedundancy(self, new_organism, whole_pop):
        '''
        Checks for redundancy, both structural and if specified, epa (d-value)
        
        Returns the organism with which new_organism is redundant, or None if no redundancy
        
        Args:
            new_organism: the organism to check for redundancy
            
            whole_pop: the list containing all organisms to check against
        '''
        # if the new organism hasn't been relaxed, then check its structure against that of every organism in whole pop
        if new_organism.epa == None:
            for organism in whole_pop:
                # need to check id's because copies of both relaxed and unrelaxed organisms are added to whole pop
                if new_organism.id != organism.id:
                # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure, organism.structure):
                        print("Organism {} failed structural redundancy - looks like organism {}".format(new_organism.id, organism.id))
                        return organism
                
        # if the new organism has been relaxed, then only check its structure against those of organisms in whole pop that have also been relaxed
        else:
            for organism in whole_pop:
                # need to check id's because copies of both relaxed and unrelaxed organisms are added to whole pop
                if new_organism.id != organism.id and organism.epa != None:
                    # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure, organism.structure):
                        print("Organism {} failed structural redundancy - looks like organism {}".format(new_organism.id, organism.id))
                        return organism
            
                    # if specified, check if their epa's match within the epa difference interval
                    if self.epa_diff > 0:
                        if abs(new_organism.epa - organism.epa) < self.epa_diff:
                            print("Organism {} failed energy per atom redundancy - looks like organism {}".format(new_organism.id, organism.id))
                            return organism    
        # should only get here if no organisms are redundant with the new organism
        return None    
        


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
        # TODO: are these reasonable?
        self.default_min_num_atoms = 2
        self.default_max_num_atoms = 50
        self.default_min_lattice_length = 0.5
        self.default_max_lattice_length = 20
        self.default_min_lattice_angle = 40
        self.default_max_lattice_angle = 140
        self.default_allow_endpoints = True
        
        # set defaults if constraints_parameters equals 'default' or None
        if constraints_parameters == None or constraints_parameters == 'default':
            self.set_all_to_defaults(composition_space)
        else:
            # check each flag to see if it's been included, and if so, whether it has been set to default or left blank
            # min number of atoms
            if 'min_num_atoms' in constraints_parameters:
                if constraints_parameters['min_num_atoms'] == None or constraints_parameters['min_num_atoms'] == 'default':
                    self.min_num_atoms = self.default_min_num_atoms
                else:
                    self.min_num_atoms = constraints_parameters['min_num_atoms']
            else:
                self.min_num_atoms = self.default_min_num_atoms    
                
            # max number of atoms   
            if 'max_num_atoms' in constraints_parameters:
                if constraints_parameters['max_num_atoms'] == None or constraints_parameters['max_num_atoms'] == 'default':
                    self.max_num_atoms = self.default_max_num_atoms
                else:
                    self.max_num_atoms = constraints_parameters['max_num_atoms']
            else:
                self.max_num_atoms = self.default_max_num_atoms    
                
            # min lattice length    
            if 'min_lattice_length' in constraints_parameters:
                if constraints_parameters['min_lattice_length'] == None or constraints_parameters['min_lattice_length'] == 'default':
                    self.min_lattice_length = self.default_min_lattice_length
                else:
                    self.min_lattice_length = constraints_parameters['min_lattice_length']
            else:
                self.min_lattice_length = self.default_min_lattice_length     
                 
            # max lattice length    
            if 'max_lattice_length' in constraints_parameters:
                if constraints_parameters['max_lattice_length'] == None or constraints_parameters['max_lattice_length'] == 'default':
                    self.max_lattice_length = self.default_max_lattice_length
                else:
                    self.max_lattice_length = constraints_parameters['max_lattice_length']
            else:
                self.max_lattice_length = self.default_max_lattice_length 
             
            # min lattice angle    
            if 'min_lattice_angle' in constraints_parameters:
                if constraints_parameters['min_lattice_angle'] == None or constraints_parameters['min_lattice_angle'] == 'default':
                    self.min_lattice_angle = self.default_min_lattice_angle
                else:
                    self.min_lattice_angle = constraints_parameters['min_lattice_angle']
            else:
                self.min_lattice_angle = self.default_min_lattice_angle
             
            # max lattice angle    
            if 'max_lattice_angle' in constraints_parameters:
                if constraints_parameters['max_lattice_angle'] == None or constraints_parameters['max_lattice_angle'] == 'default':
                    self.max_lattice_angle = self.default_max_lattice_angle
                else:
                    self.max_lattice_angle = constraints_parameters['max_lattice_angle']
            else:
                self.max_lattice_angle = self.default_max_lattice_angle    
             
            # allowing endpoints   
            if 'allow_endpoints' in constraints_parameters:
                if constraints_parameters['allow_endpoints'] == None or constraints_parameters['allow_endpoints'] == 'default':
                    self.allow_endpoints = self.default_allow_endpoints
                else:
                    self.allow_endpoints = constraints_parameters['allow_endpoints']
            else:
                self.allow_endpoints = self.default_allow_endpoints
                  
            # the per-species min interatomic distances
            if 'per_species_mids' in constraints_parameters:
                if constraints_parameters['per_species_mids'] != None and constraints_parameters['per_species_mids'] != 'default':
                    self.per_species_mids = constraints_parameters['per_species_mids'] 
                    # check each pair that's been specified to see if it needs a default mid
                    for key in self.per_species_mids:
                        if self.per_species_mids[key] == None or self.per_species_mids[key] == 'default':
                            elements = key.split()
                            radius1 = Element(elements[0]).atomic_radius
                            radius2 = Element(elements[1]).atomic_radius
                            self.per_species_mids[key] = 0.75*(radius1 + radius2)
                    # check to see if any pairs are missing, and if so, add them and set to default values
                    self.set_some_mids_to_defaults(composition_space)
                # if the per_species_mids block has been left blank or set to default, then set all the pairs to defaults
                else:
                    self.set_all_mids_to_defaults(composition_space)
            # if the per_species_mids block wasn't set in the input file, then set all the pairs to defaults
            else:
                self.set_all_mids_to_defaults(composition_space)
            
                
    def set_all_to_defaults(self, composition_space):
        '''
        Sets all general constraints (those in Constraints block of input file) to default values
        
        Args:
            composition_space: the composition space object
        '''
        self.min_num_atoms = self.default_min_num_atoms
        self.max_num_atoms = self.default_max_num_atoms
        self.min_lattice_length = self.default_min_lattice_length
        self.max_lattice_length = self.default_max_lattice_length
        self.min_lattice_angle = self.default_min_lattice_angle
        self.max_lattice_angle = self.default_max_lattice_angle
        self.allow_endpoints = self.default_allow_endpoints
        self.set_all_mids_to_defaults(composition_space) 
        
        
    def set_all_mids_to_defaults(self, composition_space):
        '''
        Sets all the per-species mids to default values based on atomic radii
        
        Args:
            composition_space: the composition space object
        '''
        # get each element type from the composition_space object
        elements = composition_space.get_all_elements()   
        # compute the per_species_mids based on atomic radii
        self.per_species_mids = {}
        for i in range(0, len(elements)):
            for j in range(i, len(elements)):
                self.per_species_mids[str(elements[i].symbol + " " + elements[j].symbol)] = 0.75*(elements[i].atomic_radius + elements[j].atomic_radius)
        
    
    def set_some_mids_to_defaults(self, composition_space):
        '''
        Compares all the possible pairs of elements to what is contained in self.per_species_mids. If any pairs are missing, adds them with default values.
        
        Args:
            composition_space: a CompositionSpace object
        ''' 
        # get each element type from the composition_space object
        elements = composition_space.get_all_elements()
        # list to hold lists of missing pairs
        missing_pairs = []
        # scan through every possible pair to check if it's already included in self.per_species_mids
        for i in range(0, len(elements)):
            for j in range(i, len(elements)):
                # check both orders
                test_key1 = elements[i].symbol + " " + elements[j].symbol
                test_key2 = elements[j].symbol + " " + elements[i].symbol
                if test_key1 not in self.per_species_mids and test_key2 not in self.per_species_mids:
                    missing_pairs.append(test_key1)
                        
        # calculate the per species mids for all the missing pairs and add them to self.per_species_mids
        for pair in missing_pairs:
            p = pair.split()
            self.per_species_mids[str(pair)] = 0.75*(Element(p[0]).atomic_radius + Element(p[1]).atomic_radius)
    
    
    def get_max_mid(self):
        '''
        Returns largest per-species minimum interatomic distance constraint
        '''
        max_mid = 0
        for key in self.per_species_mids:
            if self.per_species_mids[key] > max_mid:
                max_mid = self.per_species_mids[key]
        return max_mid     
        


class Development(object):
    '''
    A development object is used to develop an organism before evaluating its energy or adding it
    to the pool. Doesn't do redundancy checking.
    
    This is a singleton class.
    '''
    def __init__(self, niggli, scale_density):
        '''
        Creates a Development object.
        
        Args:
            niggli: a boolean indicating whether or not to do Niggli cell reduction
            
            scale_density: a boolean indicating whether or not to scale the density
        '''
        self.niggli = niggli
        self.scale_density = scale_density
        
    
    def develop(self, organism, composition_space, constraints, geometry, pool):
        '''
        Develops an organism.
        
        Returns the developed organism, or None if the organism failed development
        
        TODO: it might make more sense to return a flag indicating whether the organism survived development, since this method modifies the organism...
        
        Args:
            organism: the organism to develop
            
            composition_space: a CompositionSpace object
            
            constraints: a Constraints object
            
            geometry: a Geometry object
            
            pool: the current pool. If this method is called before a pool exists (e.g., while making the initial population)
                  then pass None as the argument instead.
        '''
        # check max num atoms constraint
        if len(organism.structure.sites) > constraints.max_num_atoms:
            print("Organism {} failed max number of atoms constraint".format(organism.id))
            return None
            
        # check min num atoms constraint
        if len(organism.structure.sites) < constraints.min_num_atoms:
            print("Organism {} failed min number of atoms constraint".format(organism.id))
            return None
        
        # check if the organism has the right composition for fixed-composition searches
        if composition_space.objective_function == "epa":
            # compare the reduced compositions to ensure a valid comparison
            reduced_composition = composition_space.endpoints[0].reduced_composition
            org_reduced_composition = organism.composition.reduced_composition
            if not reduced_composition.almost_equals(org_reduced_composition):
                print("Organism {} has incorrect composition".format(organism.id))
                return None
        
        # check if the organism is in the composition space for phase-diagram searches
        # This is kind of hacky, but the idea is to use the CompoundPhaseDiagram.transform_entries method to do 
        # the heavy lifting of determining whether a composition lies in the composition space 
        elif composition_space.objective_function == "pd":
            # cast all the endpoints to PDEntries, and just make up some energies
            pdentries = []
            for endpoint in composition_space.endpoints:
                pdentries.append(PDEntry(endpoint, -10))
            # also cast the organism we want to check to a PDEntry
            pdentries.append(PDEntry(organism.composition, -10))
            # construct the CompoundPhaseDiagram object that we'll use to check if the organism is in the composition space
            composition_checker = CompoundPhaseDiagram(pdentries, composition_space.endpoints)
            # use the CompoundPhaseDiagram to check if the organism is in the composition space by seeing how many entries it returns
            if len(composition_checker.transform_entries(pdentries, composition_space.endpoints)[0]) == len(composition_space.endpoints):
                print("Organism {} is outside the composition space".format(organism.id))
                return None
            else:
                # check the endpoints if specified and if we're not making the initial population
                if constraints.allow_endpoints == False and pool != None:
                    for endpoint in composition_space.endpoints:
                        if endpoint.almost_equals(organism.composition.reduced_composition):
                            print("Organism {} is at a composition endpoint".format(organism.id))
                            return None
                        
        # optionally do Niggli cell reduction 
        # sometimes pymatgen's reduction routine fails, so we check for that. 
        # TODO: what should we do when it fails? For now we're just rejecting it...
        if self.niggli:
            if geometry.shape == "bulk":
                # do normal Niggli cell reduction
                try:
                    organism.structure = organism.structure.get_reduced_structure()
                    organism.rotateToPrincipalDirections()
                except:
                    return None
            elif geometry.shape == "sheet":
                # do the sheet Niggli cell reduction
                try:
                    organism.reduceSheetCell()
                    organism.rotateToPrincipalDirections()
                except:
                    return None
            # TODO: call special cell reduction for other geometries here if needed (doesn't makes sense for wires or clusters)
                     
        # optionally scale the density to the average of the densities of the organisms in the promotion set 
        if self.scale_density and composition_space.objective_function == "epa" and pool != None and organism.epa == None:
            # get average volume per atom of the organisms in the promotion set
            vpa_sum = 0
            for org in pool.promotion_set:
                vpa_sum = vpa_sum + org.structure.volume/len(org.structure.sites)
            vpa_mean = vpa_sum/len(pool.promotion_set)
            # compute the new volume per atom
            num_atoms = len(organism.structure.sites)
            new_vol = vpa_mean*num_atoms
            # scale to the new volume
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                organism.structure.scale_lattice(new_vol)
                # check if the volume scaling worked 
                if str(organism.structure.lattice.a) == 'nan' or organism.structure.lattice.a > 100:
                    return None 
            
        # check the max and min lattice length constraints
        lengths = organism.structure.lattice.abc
        for length in lengths:
            if length > constraints.max_lattice_length:
                print("Organism {} failed max lattice length constraint".format(organism.id))
                return None
            elif length < constraints.min_lattice_length:
                print("Organism {} failed min lattice length constraint".format(organism.id))
                return None
            
        # check the max and min lattice angle constraints
        angles = organism.structure.lattice.angles
        for angle in angles:
            if angle > constraints.max_lattice_angle:
                print("Organism {} failed max lattice angle constraint".format(organism.id))
                return None
            elif angle < constraints.min_lattice_angle:
                print("Organism {} failed min lattice angle constraint".format(organism.id))
                return None
            
        # check the per-species minimum interatomic distance constraints
        species_symbols = organism.structure.symbol_set
        for site in organism.structure.sites:
            for species_symbol in species_symbols:
                # get the mid for this particular pair. We don't know the ordering in per_species_mids, so try both
                test_key1 = species_symbol + " " + site.specie.symbol
                test_key2 = site.specie.symbol + " " + species_symbol
                if test_key1 in constraints.per_species_mids:
                    mid = constraints.per_species_mids[test_key1]
                elif test_key2 in constraints.per_species_mids:
                    mid = constraints.per_species_mids[test_key2]  
                # get all the sites within a sphere of radius mid centered on the current site
                neighbors = organism.structure.get_neighbors(site, mid)
                # check each neighbor in the sphere to see if it has the forbidden type
                for neighbor in neighbors:
                    if neighbor[0].specie.symbol == species_symbol:
                        print("Organism {} failed per-species minimum interatomic distance constraint".format(organism.id))
                        return None
            
        # check the max size constraint (can only fail for non-bulk geometries)
        if geometry.getSize(organism) > geometry.max_size:
            print("Organism {} failed max size constraint".format(organism.id))
            return None
        
        # check the min size constraint (can only fail for non-bulk geometries)
        if geometry.getSize(organism) < geometry.min_size:
            print("Organism {} failed min size constraint".format(organism.id))
            return None
        
        # return the organism if it survived
        return organism
    
        
    
class EnergyCalculator(object):
    '''
    Handles calculating the energy of organisms.
    
    Not meant to be instantiated, but rather subclassed by particular Calculators, like VaspEnergyCalculator
    or GulpEnergyCalculator.
    
    TODO: is the superclass even necessary?
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
        
        

class VaspEnergyCalculator(object):
    '''
    Calculates the energy of an organism using VASP.
    '''
    def __init__(self, vasp_code_params, run_dir_name):
        '''
        Args:
            vasp_code_params: the parameters needed for preparing a vasp calculation (INCAR, KPOINTS, POTCAR files)
            
            run_dir_name: the name of the garun directory containing all the output of the search (just the name, not the path)
        '''
        # TODO: implement me. Just keeping the paths to the input files should be enough
    
    
    def doEnergyCalculation(self, org):
        '''
        Calculates the energy of an organism using VASP
        
        Returns an organism that has been parsed from the output files of the energy code, or None if the calculation 
        failed. Does not do development or redundancy checking.
        
        Args:
            organism: the organism whose energy we want to calculate
            
        TODO: use custodian package for error handling
        '''
        # TODO: implement me
        # 1. prepare input files for calculation
        # 2. submit calculation (by running external script)
        # 3. when external script returns, parse organism from the energy code output files
    


class GulpEnergyCalculator(object):
    '''
    Calculates the energy of an organism using GULP.
    '''
    def __init__(self, header_file, potential_file, run_dir_name):
        '''
        Args:
            header_file: the path to the gulp header file
            
            potential_file: the path to the gulp potential file
            
            run_dir_name: the name of the garun directory containing all the output of the search (just the name, not the path)
            
        Precondition: the header and potential files exist and are valid
        '''
        # read the gulp header file and store it as a string
        with open (header_file, "r") as gulp_header_file:
            self.header = gulp_header_file.read()
            
        # read the gulp potential file and store it as a string
        with open (potential_file, "r") as gulp_potential_file:
            self.potential = gulp_potential_file.read()
            
        # the name of the garun directory
        self.run_dir_name = run_dir_name
        
        # for processing gulp input and output
        self.gulp_io = gulp_caller.GulpIO()
        
        # whether the anions and cations are polarizable in the gulp potential
        self.anions_shell, self.cations_shell = self.getShells()
        
    
    def getShells(self):
        '''
        Determines whether the anions and cations have shells by looking at the potential file
        
        Returns two booleans indicating whether the anions and cations have shells, respectively
        '''
        # get the symbols of the elements with shells
        shells = []
        for line in self.potential.split('\n'):
            if 'shel' in line:
                line_parts = line.split()
                # the element symbol is always right before the 'shel' keyword
                shells.append(str(line_parts[line_parts.index('shel') - 1]))
        
        # remove duplicates
        shells = list(set(shells))
        
        # determine whether the elements with shells are anions and/or cations
        anions_shell = False
        cations_shell = False
        for symbol in shells:
            element = Element(symbol)
            if element in gulp_caller._anions:
                anions_shell = True
            elif element in gulp_caller._cations:
                cations_shell = True
            
        return anions_shell, cations_shell
    
    
    def doEnergyCalculation(self, organism, dictionary, key):
        '''
        Calculates the energy of an organism using GULP, and stores the relaxed organism in the provided dictionary at the provided key. 
        If the calculation fails, stores None in the dictionary instead.
        
        Args:
            organism: the organism whose energy we want to calculate
            
            dictionary: a dictionary in which to store the relaxed organism
            
            key: the key specifying where to store the relaxed organism in the dictionary
            
        Precondition: the garun directory and temp subdirectory exist
        
        TODO: think about ways to make this more robust - look at Will's code
        TODO: might be better to eventually use the custodian package for error handling...
        '''
        # make the job directory
        job_dir_path = str(getcwd()) + '/' + str(self.run_dir_name) + '/temp/' + str(organism.id)
        mkdir(job_dir_path)
        
        # get the structure in gulp input format
        # TODO: might be better to make my own implementation of this so abritrary elements can have shells (not just all anions and/or all cations)
        structure_lines = self.gulp_io.structure_lines(organism.structure, anion_shell_flg = self.anions_shell, cation_shell_flg = self.cations_shell, symm_flg = False)
        
        # make the gulp input from the structrure, header and potential
        gulp_input = self.header + structure_lines + self.potential
        
        # print gulp input to a file for user's reference 
        gin_path = job_dir_path + '/' + str(organism.id) + '.gin'
        gin_file = open(gin_path, 'w')
        gin_file.write(gulp_input)
        gin_file.close()
        
        print('Starting GULP calculation on organism {}'.format(organism.id))
        
        # run the gulp calculation by running a 'callgulp' script as a subprocess. If errors are thrown, print them to the gulp output file
        try:
            gulp_output = check_output(['callgulp', gin_path], stderr = STDOUT)
        except CalledProcessError as e:
            # print the output of the bad gulp calc to a file for the user's reference
            gout_file = open(job_dir_path + '/' + str(organism.id) + '.gout', 'w')
            gout_file.write(e.output)
            gout_file.close()
            # error message and set dictionary element to None before returning
            print('Error running GULP on organism {}'.format(organism.id))
            dictionary[key] = None
            return
        
        # print gulp output to a file for user's reference 
        gout_file = open(job_dir_path + '/' + str(organism.id) + '.gout', 'w')
        gout_file.write(gulp_output)
        gout_file.close()
       
        # check if not converged (part of this is copied from pymatgen)
        conv_err_string = "Conditions for a minimum have not been satisfied"
        gradient_norm = self.getGradNorm(gulp_output)
        if conv_err_string in gulp_output and gradient_norm > 0.1:
            print('The GULP calculation on organism {} did not converge'.format(organism.id))
            dictionary[key] = None
            return
       
        # parse the relaxed structure from the gulp output
        try:
            # TODO: will have to change this line if pymatgen fixes the gulp parser
            relaxed_structure = self.get_relaxed_structure(gulp_output)
        except:
            print('Error reading structure of organism {} from GULP output'.format(organism.id))
            dictionary[key] = None
            return
        
        # parse the total energy from the gulp output
        try:
            total_energy = self.getEnergy(gulp_output)
        except:
            print('Error reading energy of organism {} from GULP output'.format(organism.id))
            dictionary[key] = None
            return
        
        # parse the number of atoms used by gulp from the gulp output (sometimes gulp takes a supercell)
        num_atoms = self.getNumAtoms(gulp_output)
            
        # for conveniance/testing...
        #relaxed_structure.to('poscar', job_dir_path + '/' + str(organism.id) + '.vasp')
            
        # assign the relaxed structure and energy to the organism, and compute the epa
        organism.structure = relaxed_structure
        organism.epa = total_energy/num_atoms
        organism.total_energy = organism.epa*organism.structure.num_sites
        
        # print out the enerngy per atom of the organism
        print('Setting energy of organism {} to {} eV/atom'.format(organism.id, organism.epa))
        
        # store the relaxed organism in the specified dictionary with the specified key 
        dictionary[key] = organism
        
       
    def getGradNorm(self, gout):
        '''
        Parses the final gradient norm from the Gulp output
        
        Args:
            gout: the Gulp output, as a string
        ''' 
        output_lines = gout.split("\n")
        for line in output_lines:
            if "Final Gnorm," in line:
                # the gradient norm is the fourth element in the line
                line_parts = line.split()
                return float(line_parts[3])  
       
       
    def getEnergy(self, gout):
        '''
        Parses the final energy from the Gulp output
        
        Args:
            gout: the Gulp output, as a string
        '''
        output_lines = gout.split("\n")
        for line in output_lines:
            if "Final energy" in line:
                # the energy is the fourth element in the line
                line_parts = line.split()
                return float(line_parts[3])
       
       
    def getNumAtoms(self, gout):
        '''
        Parses the number of atoms used by Gulp in the calculation
        
        Args:
            gout: the Gulp output, as a string
        '''
        output_lines = gout.split("\n")
        for line in output_lines:
            if "Total number atoms" in line:
                # the number is the last element in the list
                line_parts = line.split()
                return int(line_parts[-1])
       
        
    # this method is copied from GulpIO.get_relaxed_structure
    # I modified it slightly to get it to work
    # TODO: if pymatgen fixes this method, then I can delete this. Alternatively, I could submit a pull request with my fix    
    def get_relaxed_structure(self, gout):
        #Find the structure lines
        structure_lines = []
        cell_param_lines = []
        output_lines = gout.split("\n")
        no_lines = len(output_lines)
        i = 0
        # Compute the input lattice parameters
        while i < no_lines:
            line = output_lines[i]
            if "Full cell parameters" in line:
                i += 2
                line = output_lines[i]
                a = float(line.split()[8])
                alpha = float(line.split()[11])
                line = output_lines[i + 1]
                b = float(line.split()[8])
                beta = float(line.split()[11])
                line = output_lines[i + 2]
                c = float(line.split()[8])
                gamma = float(line.split()[11])
                i += 3
                break
            elif "Cell parameters" in line:
                i += 2
                line = output_lines[i]
                a = float(line.split()[2])
                alpha = float(line.split()[5])
                line = output_lines[i + 1]
                b = float(line.split()[2])
                beta = float(line.split()[5])
                line = output_lines[i + 2]
                c = float(line.split()[2])
                gamma = float(line.split()[5])
                i += 3
                break
            else:
                i += 1

        while i < no_lines:
            line = output_lines[i]
            if "Final fractional coordinates of atoms" in line or "Final asymmetric unit coordinates" in line:
                # read the site coordinates in the following lines
                i += 6
                line = output_lines[i]
                while line[0:2] != '--':
                    structure_lines.append(line)
                    i += 1
                    line = output_lines[i]
                    # read the cell parameters
                i += 9
                line = output_lines[i]
                if "Final cell parameters" in line:
                    i += 3
                    for del_i in range(6):
                        line = output_lines[i + del_i]
                        cell_param_lines.append(line)
                break
            else:
                i += 1

        #Process the structure lines
        if structure_lines:
            sp = []
            coords = []
            for line in structure_lines:
                fields = line.split()
                if fields[2] == 'c':
                    sp.append(fields[1])
                    coords.append(list(float(x) for x in fields[3:6]))
        else:
            raise IOError("No structure found")

        if cell_param_lines:
            a = float(cell_param_lines[0].split()[1])
            b = float(cell_param_lines[1].split()[1])
            c = float(cell_param_lines[2].split()[1])
            alpha = float(cell_param_lines[3].split()[1])
            beta = float(cell_param_lines[4].split()[1])
            gamma = float(cell_param_lines[5].split()[1])
        latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        return Structure(latt, sp, coords)
            
        
        
class StoppingCriteria(object):
    '''
    Defines when the search should be stopped.
    '''
    def __init__(self, stopping_parameters, composition_space):
        '''
        Args:
            stopping_parameters: a dictionary of parameters
            
            composition_space: a CompositionSpace object
        '''
        # set the defaults
        if composition_space.objective_function == 'epa':
            self.default_num_energy_calcs = 500
        elif composition_space.objective_function == 'pd':
            self.default_num_energy_calcs = 1000
        self.default_value_achieved = None
        self.default_found_structure = None
        
        # whether or not the stopping criteria are satisfied
        self.are_satisfied = False
        # to keep track of how many energy calculations have been done
        self.calc_counter = 0
        
        # set defaults if stopping_parameters equals 'default' or None
        if stopping_parameters == None or stopping_parameters == 'default':
            self.num_energy_calcs = self.default_num_energy_calcs
            self.value_achieved = self.default_value_achieved
            self.found_structure = self.default_found_structure
        else:
            # check each flag to see if it's been included, and if so, whether it has been set to default or left blank  
            # value achieved
            if 'value_achieved' in stopping_parameters:
                if stopping_parameters['value_achieved'] == None or stopping_parameters['value_achieved'] == 'default':
                    self.value_achieved = self.default_value_achieved
                else:
                    self.value_achieved = stopping_parameters['value_achieved']
            else:
                self.value_achieved = self.default_value_achieved
                
            # found structure
            if 'found_structure' in stopping_parameters:
                if stopping_parameters['found_structure'] == None or stopping_parameters['found_structure'] == 'default':
                    self.found_structure = self.default_found_structure
                else:
                    # read the structure from the file
                    self.found_structure = Structure.from_file(stopping_parameters['found_structure'])     
            else:
                self.found_structure = self.default_found_structure
                
            # num energy calcs
            if 'num_energy_calcs' in stopping_parameters:
                if stopping_parameters['num_energy_calcs'] == None or stopping_parameters['num_energy_calcs'] == 'default':
                    # only use default if other two methods haven't been specified
                    if self.value_achieved == None and self.found_structure == None:
                        self.num_energy_calcs = self.default_num_energy_calcs
                    else:
                        self.num_energy_calcs = None
                else:
                    self.num_energy_calcs = stopping_parameters['num_energy_calcs']
            # only use default if other two methods haven't been specified
            elif self.value_achieved == None and self.found_structure == None:    
                self.num_energy_calcs = self.default_num_energy_calcs 
            else:
                self.num_energy_calcs = None
          
            
    def updateCalcCounter(self):
        '''
        If num_energy_calcs stopping criteria is being used, increments calc counter and updates are_satisfied if necessary.
        '''
        if self.num_energy_calcs != None:
            self.calc_counter = self.calc_counter + 1
            if self.calc_counter >= self.num_energy_calcs:
                self.are_satisfied = True
                
    
    def checkOrganism(self, organism):
        '''
        If value_achieved or found_structure stopping criteria are used, checks if the relaxed organism satisfies them, and if so, updates are_satisfied.
        
        Args:
            organism: a relaxed organism whose value has been computed
        '''
        if self.value_achieved != None:
            if organism.value <= self.value_achieved:
                self.are_satisfied = True
        if self.found_structure != None:
            # TODO: check the tolerances for the structure matching algorithm - pymatgen's defualts might be too loose...
            if self.found_structure.matches(organism.structure):
                self.are_satisfied = True
        
                
           
        
    

        
        
        
        
########## area for casual testing ########## 

# make a structure object
#lattice = [[1,0.5,0], [0.5,1,0], [0,0,-1]]
#species = ["C", "Si"]
#coordinates = [[0.25,0.25,0.25],[0.75,0.75,0.75]]
#structure1 = Structure(lattice, species, coordinates)

# make an organism object
#org1 = Organism(structure1)

#print(org1.structure.lattice)

#org1.rotateToPrincipalDirections()
#print("")
#print(org1.structure.lattice)

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
