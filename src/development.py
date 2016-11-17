# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Development module:

This module contains the classes used when developing an organism 
before and after it is submitted for an energy calculation.

1. Constraints: contains the constraints placed on structures
2. Developer: develops organisms
3. RedundancyGuard: checks if an organism is redundant 
4. Geometry: contains geometry-specific data for non-bulk searches, including any additional constraints

"""

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element, DummySpecie
from pymatgen.core.sites import Site
from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_matcher import ElementComparator

import warnings
import numpy as np


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
        self.default_max_num_atoms = 50
        self.default_min_lattice_length = 0.5
        self.default_max_lattice_length = 20
        self.default_min_lattice_angle = 40
        self.default_max_lattice_angle = 140
        self.default_mid_factor = 0.7
        
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
                  
            # the per-species minimum interatomic distances
            if 'per_species_mids' in constraints_parameters:
                if constraints_parameters['per_species_mids'] != None and constraints_parameters['per_species_mids'] != 'default':
                    self.per_species_mids = constraints_parameters['per_species_mids'] 
                    # check each pair that's been specified to see if it needs a default mid
                    for key in self.per_species_mids:
                        if self.per_species_mids[key] == None or self.per_species_mids[key] == 'default':
                            elements = key.split()
                            radius1 = Element(elements[0]).atomic_radius
                            radius2 = Element(elements[1]).atomic_radius
                            self.per_species_mids[key] = self.default_mid_factor*(radius1 + radius2)
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
                self.per_species_mids[str(elements[i].symbol + " " + elements[j].symbol)] = self.default_mid_factor*(elements[i].atomic_radius + elements[j].atomic_radius)
        
    
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
            self.per_species_mids[str(pair)] = self.default_mid_factor*(Element(p[0]).atomic_radius + Element(p[1]).atomic_radius)
    
    
    def get_max_mid(self):
        '''
        Returns largest per-species minimum interatomic distance constraint
        '''
        max_mid = 0
        for key in self.per_species_mids:
            if self.per_species_mids[key] > max_mid:
                max_mid = self.per_species_mids[key]
        return max_mid     
        


class Developer(object):
    '''
    A development object is used to develop an organism before evaluating its energy or adding it
    to the pool. Doesn't do redundancy checking.
    '''
    def __init__(self, development_parameters, geometry): 
        '''
        Creates a Development object.
        
        Args:
            niggli: a boolean indicating whether or not to do Niggli cell reduction
            
            scale_density: a boolean indicating whether or not to scale the density
            
            geometry: a Geometry object
        '''
        # defaults
        self.default_niggli = True
        
        if geometry.shape == 'bulk':
            self.default_scale_density = True
        else:
            self.default_scale_density = False
        
        # set defaults if development_parameters equals 'default' or None
        if development_parameters == None or development_parameters == 'default':
            self.niggli = self.default_niggli
            self.scale_density = self.default_scale_density
        
        # otherwise, parse the parameters and set to defaults if necessary
        else:
            if 'niggli' in development_parameters:
                if development_parameters['niggli'] == None or development_parameters['niggli'] == 'default':
                    self.niggli = self.default_niggli
                else:
                    self.niggli = development_parameters['niggli']
            else:
                # if no 'niggli' flag, then just use the default
                self.niggli = self.default_niggli
                
            if 'scale_density' in development_parameters:
                if development_parameters['scale_density'] == None or development_parameters['scale_density'] == 'default':
                    self.scale_density = self.default_scale_density
                else:
                    self.scale_density = development_parameters['scale_density']
            else:
                # if no 'scale_density' flag, then just use the default
                self.scale_density = self.default_scale_density
        
    
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
            
            pool: the Pool object
        '''
        # check max num atoms constraint
        if len(organism.structure.sites) > constraints.max_num_atoms:
            print("Organism {} failed max number of atoms constraint ".format(organism.id))
            return None
            
        # check min num atoms constraint
        if len(organism.structure.sites) < constraints.min_num_atoms:
            print("Organism {} failed min number of atoms constraint ".format(organism.id))
            return None
        
        # check if the organism has the right composition for fixed-composition searches
        if composition_space.objective_function == "epa":
            # compare the reduced compositions to ensure a valid comparison
            reduced_composition = composition_space.endpoints[0].reduced_composition
            org_reduced_composition = organism.composition.reduced_composition
            if not reduced_composition.almost_equals(org_reduced_composition):
                print("Organism {} has incorrect composition ".format(organism.id))
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
                print("Organism {} is outside the composition space ".format(organism.id))
                return None
            # check the endpoints if we're not making the initial population
            elif len(pool.to_list()) > 0:
                for endpoint in composition_space.endpoints:
                    if endpoint.almost_equals(organism.composition.reduced_composition):
                        print("Organism {} is at a composition endpoint ".format(organism.id))
                        return None
                        
        # optionally do Niggli cell reduction 
        # sometimes pymatgen's reduction routine fails, so we check for that. 
        if self.niggli:
            if geometry.shape == 'bulk':
                # do normal Niggli cell reduction
                try:
                    organism.structure = organism.structure.get_reduced_structure()
                    organism.rotate_to_principal_directions()
                except:
                    print('Niggli cell reduction failed on organism {} during development '.format(organism.id))
                    return None
            elif geometry.shape == 'sheet':
                # do the sheet Niggli cell reduction
                try:
                    organism.reduce_sheet_cell()
                    organism.rotate_to_principal_directions()
                except:
                    print('2D Niggli cell reduction failed on organism {} during development '.format(organism.id))
                    return None
            # TODO: call special cell reduction for other geometries here if needed (doesn't makes sense for wires or clusters)
                     
        # optionally scale the volume per atom of unrelaxed structures
        if self.scale_density and len(pool.promotion_set) > 0 and organism.epa == None:
            # scale to the average of the volumes per atom of the organisms in the promotion set, and increase by 10%
            if composition_space.objective_function == 'epa':
                # get average volume per atom of the organisms in the promotion set
                vpa_sum = 0
                for org in pool.promotion_set:
                    vpa_sum += org.structure.volume/len(org.structure.sites)
                # take the mean, and increase by 10%
                vpa_mean = 1.1*(vpa_sum/len(pool.promotion_set))
                # compute the new volume
                num_atoms = len(organism.structure.sites)
                new_vol = vpa_mean*num_atoms
                # scale to the new volume
                with warnings.catch_warnings(): # this is to suppress the warnings produced if the scale_lattice method fails
                    warnings.simplefilter('ignore')
                    organism.structure.scale_lattice(new_vol)
                    # check if the volume scaling worked 
                    if str(organism.structure.lattice.a) == 'nan' or organism.structure.lattice.a > 100:
                        print('Volume scaling failed on organism {} during development '.format(organism.id))
                        return None 
                    
            # scale to the weighted average of the volumes per atom of the structures on the convex hull that this one would decompose to      
            elif composition_space.objective_function == 'pd':
                # make a compound phase diagram from the organisms in the pool
                pdentries = []
                for org in pool.promotion_set:
                    pdentries.append(PDEntry(org.composition, org.total_energy))
                compound_pd = CompoundPhaseDiagram(pdentries, composition_space.endpoints)
                
                # make a pdanalyzer object 
                pdanalyzer = PDAnalyzer(compound_pd)
                # transform the organism's composition (just use a dummy value for the energy)
                transformed_entry = compound_pd.transform_entries([PDEntry(organism.composition, 10)], composition_space.endpoints)
                # get a list containing the transformed species and amounts
                transformed_list = str(transformed_entry[0][0]).split()
                # remove unneeded items from the list
                del transformed_list[0]
                popped = ''
                while popped != 'with':
                    popped = transformed_list.pop()    
                # separate the dummy species symbols from the amounts
                symbols = []
                amounts = []
                for entry in transformed_list:
                    split_entry = entry.split('0+')
                    symbols.append(split_entry[0])
                    amounts.append(float(split_entry[1]))
                # make a dictionary mapping dummy species to amounts
                dummy_species_amounts = {}
                for i in range(len(symbols)):
                    dummy_species_amounts[DummySpecie(symbol=symbols[i])] = amounts[i]   
                # create a composition object with dummy species
                dummy_comp = Composition(dummy_species_amounts)
                # use the pdanalyzer to decompose of the organism's (dummy) composition
                decomp = pdanalyzer.get_decomposition(dummy_comp)
                # parse the original compositions and amounts from the decomposition
                fractions = []
                comps = []
                for item in decomp:
                    fractions.append(decomp[item])
                    first_split = str(item).split(',')
                    second_split = first_split[0].split()
                    while second_split[0] != 'composition':
                        del second_split[0]
                    del second_split[0]
                    # build the composition string
                    comp_string = ''
                    for symbol in second_split:
                        comp_string += str(symbol)
                    comps.append(Composition(comp_string))
                # get weighted average volume per atom of the organisms in the decomposition
                vpa_mean = 0
                for i in range(len(comps)):
                    for org in pool.promotion_set:
                        if (comps[i].reduced_composition).almost_equals(org.composition.reduced_composition):
                            vpa_mean += (org.structure.volume/len(org.structure.sites))*fractions[i]
                # increase the weighted average volume per atom by 10% 
                vpa_mean = 1.1*vpa_mean
                # compute the new volume
                num_atoms = len(organism.structure.sites)
                new_vol = vpa_mean*num_atoms
                # scale to the new volume
                with warnings.catch_warnings(): # this is to suppress the warnings produced if the scale_lattice method fails
                    warnings.simplefilter('ignore')
                    organism.structure.scale_lattice(new_vol)
                    # check if the volume scaling worked 
                    if str(organism.structure.lattice.a) == 'nan' or organism.structure.lattice.a > 100:
                        print('Volume scaling failed on organism {} during development '.format(organism.id))
                        return None 
                        
        # check the max and min lattice length constraints
        lengths = organism.structure.lattice.abc
        for length in lengths:
            if length > constraints.max_lattice_length:
                print("Organism {} failed max lattice length constraint ".format(organism.id))
                return None
            elif length < constraints.min_lattice_length:
                print("Organism {} failed min lattice length constraint ".format(organism.id))
                return None
            
        # check the max and min lattice angle constraints
        angles = organism.structure.lattice.angles
        for angle in angles:
            if angle > constraints.max_lattice_angle:
                print("Organism {} failed max lattice angle constraint ".format(organism.id))
                return None
            elif angle < constraints.min_lattice_angle:
                print("Organism {} failed min lattice angle constraint ".format(organism.id))
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
                        print("Organism {} failed per-species minimum interatomic distance constraint ".format(organism.id))
                        return None
            
        # check the max size constraint (can only fail for non-bulk geometries)
        if geometry.get_size(organism) > geometry.max_size:
            print("Organism {} failed max size constraint ".format(organism.id))
            return None
        
        # check the min size constraint (can only fail for non-bulk geometries)
        if geometry.get_size(organism) < geometry.min_size:
            print("Organism {} failed min size constraint ".format(organism.id))
            return None
        
        # return the organism if it survived
        return organism

    

class RedundancyGuard(object):
    '''
    A redundancy guard.
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
        self.structure_matcher = StructureMatcher(self.lattice_length_tol, self.site_tol, self.lattice_angle_tol, self.use_primitive_cell, False, self.attempt_supercell, False, ElementComparator())
    
        
    def set_all_to_defaults(self):
        '''
        Sets all the redundancy parameters to default values
        '''
        self.lattice_length_tol = self.default_lattice_length_tol
        self.lattice_angle_tol = self.default_lattice_angle_tol
        self.site_tol = self.default_site_tol
        self.use_primitive_cell = self.default_use_primitive_cell
        self.attempt_supercell = self.default_attempt_supercell
        self.epa_diff = self.default_epa_diff
    
        
    def check_redundancy(self, new_organism, orgs_list):
        '''
        Checks for redundancy, both structural and if specified, epa (d-value)
        
        Returns the organism with which new_organism is redundant, or None if no redundancy
        
        Args:
            new_organism: the organism to check for redundancy
            
            orgs_list: the list containing all organisms to check against
        '''
        # if the new organism hasn't been relaxed, then check its structure against that of every organism in the list
        if new_organism.epa == None:
            for organism in orgs_list:
                # need to check id's because copies of both relaxed and unrelaxed organisms are added to the whole pop list of organisms
                if new_organism.id != organism.id:
                # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure, organism.structure):
                        print("Organism {} failed structural redundancy - looks like organism {} ".format(new_organism.id, organism.id))
                        return organism
                
        # if the new organism has been relaxed, then only check its structure against those of organisms in the list that have also been relaxed
        else:
            for organism in orgs_list:
                # need to check id's because copies of both relaxed and unrelaxed organisms are added to whole pop list of organisms
                if new_organism.id != organism.id and organism.epa != None:
                    # check if their structures match
                    if self.structure_matcher.fit(new_organism.structure, organism.structure):
                        print("Organism {} failed structural redundancy - looks like organism {} ".format(new_organism.id, organism.id))
                        return organism
            
                    # if specified, check if their epa's match within the epa difference interval
                    if self.epa_diff > 0:
                        if abs(new_organism.epa - organism.epa) < self.epa_diff:
                            print("Organism {} failed energy per atom redundancy - looks like organism {} ".format(new_organism.id, organism.id))
                            return organism    
        # should only get here if no organisms are redundant with the new organism
        return None 
    


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
        self.default_max_size = np.inf
        self.default_min_size = -np.inf
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
            self.pad_sheet(organism)
        elif self.shape == 'wire':
            self.pad_wire(organism)
        elif self.shape == 'cluster':
            self.pad_cluster(organism)
     
        
    def pad_sheet(self, organism):
        '''
        Adds vertical vacuum padding to a sheet, and makes the c-lattice vector normal to the plane of the sheet. 
        The atoms are shifted up to the center of the padded sheet. Changes an organism's structure
        
        Args:
            organism: an Organism object
        ''' 
        # rotate into principal directions
        organism.rotate_to_principal_directions()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the layer thickness of the sheet
        cart_bounds = organism.get_bounding_box(cart_coords = True)
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
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
        
        
    def pad_wire(self, organism):
        '''
        Makes c lattice vector parallel to z-axis, and adds vacuum padding around a wire in the x and y directions by replacing a and b lattice vectors with padded vectors along 
        the x and y axes, respectively. The atoms are shifted to the center of the padded cell. Changes the structure of an organism.
        
        Args:
            organism: the organism to pad 
        '''
        # rotate c parallel to z-axis
        organism.rotate_c_parallel_to_z()
        # get the species and Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the wire in the x and y directions
        cart_bounds = organism.get_bounding_box(cart_coords = True)
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
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
    
    
    def pad_cluster(self, organism):
        '''
        Adds vacuum padding around a cluster. Changes the structure of an organism.
        
        Args:
            organism: the organism to pad 
        '''
        # get the species and Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the cluster in each direction
        cart_bounds = organism.get_bounding_box(cart_coords = True)
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
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        for i in range(0, len(organism.structure.sites)):
            organism.structure.translate_sites(i, translation_vector, frac_coords = True, to_unit_cell = False)
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
    
    
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
            self.unpad_sheet(organism, constraints)
        elif self.shape == 'wire':
            self.unpad_wire(organism, constraints)
        elif self.shape == 'cluster':
            self.unpad_cluster(organism, constraints)
        
    
    def unpad_sheet(self, organism, constraints):
        '''
        Removes vertical vacuum padding from a sheet, leaving only enough to satisfy the per-species MID constraints, and makes the c-lattice vector normal to 
        the plane of the sheet (if it isn't already). Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad 
            
            constraints: a Constraints object
        '''
        # rotate into principal directions
        organism.rotate_to_principal_directions()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the layer thickness of the sheet
        layer_thickness = self.get_layer_thickness(organism)
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
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
        # slightly shift the atoms vertically so they lie in the (vertical) center of the cell
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0, 0, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector, frac_coords = True, to_unit_cell = False)
        
          
    def unpad_wire(self, organism, constraints):
        '''
        Removes vacuum padding around a wire. Changes the structure of an organism.
        
        Args:
            organism: the Organism to unpad
            
            constraints: a Constraints object
        '''
        # rotate c parallel to z-axis
        organism.rotate_c_parallel_to_z()
        # get the species and their Cartesian coordinates
        species = organism.structure.species
        cartesian_coords = organism.structure.cart_coords
        # get the size of the wire in the x and y directions
        cart_bounds = organism.get_bounding_box(cart_coords = True)
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
        # translate all the atoms so they're in the cell (needed in case new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
        # shift the atoms horizontally so they're in the center
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.0]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector, frac_coords = True, to_unit_cell = False)
        
        
    def unpad_cluster(self, organism, constraints):
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
        cart_bounds = organism.get_bounding_box(cart_coords = True)
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
        # translate all the atoms so they're in the cell (needed in case the new lattice caused some of them to lie outside the cell)
        organism.translate_atoms_into_cell()
        # shift the atoms a little so they're in the center
        frac_bounds = organism.get_bounding_box(cart_coords = False)
        x_center = frac_bounds[0][0] + (frac_bounds[0][1] - frac_bounds[0][0])/2
        y_center = frac_bounds[1][0] + (frac_bounds[1][1] - frac_bounds[1][0])/2
        z_center = frac_bounds[2][0] + (frac_bounds[2][1] - frac_bounds[2][0])/2
        translation_vector = [0.5 - x_center, 0.5 - y_center, 0.5 - z_center]
        site_indices = [i for i in range(len(organism.structure.sites))]
        organism.structure.translate_sites(site_indices, translation_vector, frac_coords = True, to_unit_cell = False)
        
        
    def get_size(self, organism):
        '''
        Returns the size of an organism with a non-bulk shape. Returns 0 if the shape is bulk.
        
        Args:
            organism: the Organism whose size to get
        '''
        # call other methods based on the value of self.shape
        if self.shape == 'bulk':
            return 0
        elif self.shape == 'sheet':
            return self.get_layer_thickness(organism)
        elif self.shape == 'wire':
            return self.get_wire_diameter(organism)
        elif self.shape == 'cluster':
            return self.get_cluster_diameter(organism)
       
        
    def get_layer_thickness(self, organism):
        '''
        Returns the layer thickness of a sheet structure, which is the maximum vertical distance between atoms in the cell.
        
        Assumes that the organism has already been rotated into the principal directions, and that plane of the sheet is parallel to the a-b facet.
        '''
        cart_bounds = organism.get_bounding_box(cart_coords = True)
        layer_thickness = cart_bounds[2][1] - cart_bounds[2][0]
        return layer_thickness
    
    
    def get_wire_diameter(self, organism):
        '''
        Returns the diameter of a wire structure, defined as the maximum distance between atoms projected to the x-y plane.
        
        Assumes that the organism has already been put into wire format (c lattice vector is parallel to z-axis), and that all sites are located inside the cell (i.e., have 
        fractional coordinates between 0 and 1). Generally called after Geometry.unpad_wire has already been called on the organism.
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
        
    
    def get_cluster_diameter(self, organism):
        '''
        Returns the diameter of a cluster structure, defined as the maximum distance between atoms in the cell
        
        Assumes that all sites are located inside the cell (i.e., have fractional coordinates between 0 and 1). Generally called after
        Geometry.unpad_cluster has already been called on the organism.
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
