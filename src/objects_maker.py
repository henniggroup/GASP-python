# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Objects Maker module:

This module contains a function for creating the (singleton) objects 
used by the genetic algorithm during the search.

"""

import os
import general
import variations
import energy_code_interfaces
import organism_creators
import development
from pymatgen.core.structure import Structure


def make_objects(parameters):
    '''
    Constructs the needed objects for the genetic algorithm search
    
    Returns a dictionary containing the objects
    
    Args:
        parameters: the dictionary produced by calling yaml.load() on the input file
    '''
    # initial the dictionary to hold everything
    objects_dict = {}
 
    # make the composition space object
    if 'CompositionSpace' in parameters:
        composition_space = general.CompositionSpace(parameters['CompositionSpace'])
    else:
        print('Input file must contain a "CompositionSpace" block.')
        print("Quitting...")
        quit()
        
    # put the composition space in the dictionary
    objects_dict['composition_space'] = composition_space
    
    # make the constraints object
    if 'Constraints' in parameters:
        constraints = development.Constraints(parameters['Constraints'], composition_space)
    else:
        # if no Constraints block is given in the input file, then just use default values for everything
        constraints = development.Constraints('default', composition_space)
        
    # put the constraints in the dictionary
    objects_dict['constraints'] = constraints
    
    # make the geometry object
    if 'Geometry' in parameters:
        geometry = development.Geometry(parameters['Geometry'])
    else:
        # if no Geometry block is given in the input file, then just use default values for everything
        geometry = development.Geometry('default')
        
    # put the geometry in the dictionary
    objects_dict['geometry'] = geometry
    
    # make the development object
    if 'Development' in parameters:
        developer = development.Developer(parameters['Development'], geometry)
    else:
        # if no Development block is given in the input file, then just use default values for everything
        developer = development.Developer('default', geometry)
        
    # put the development in the dictionary
    objects_dict['developer'] = developer

    # make the redundancy guard object
    if 'RedundancyGuard' in parameters:
        redundancy_guard = development.RedundancyGuard(parameters['RedundancyGuard'])
    else:
        redundancy_guard = development.RedundancyGuard('default')
        
    # put the redundancy guard in the dictionary
    objects_dict['redundancy_guard'] = redundancy_guard
    
    # make the id generator
    id_generator = general.IDGenerator()
    
    # put the id generator in the dictionary
    objects_dict['id_generator'] = id_generator
    
    # make the organism creators
    initial_organism_creators = []
    # if no method is specified for making the initial population, throw an error if we're doing a phase diagram search
    if 'InitialPopulation' not in parameters: 
        if composition_space.objective_function == 'pd':
            print('For phase diagram searches, reference structures at each endpoint of the composition space must be provided in the initial population.')
            print('Please use the "from_files" flag in the InitialPopulation block to provide the reference structures.')
            print('Quitting...')
            quit()
        # otherwise, make a random one by default
        else:
            random_organism_creator = organism_creators.RandomOrganismCreator('default', composition_space)
            initial_organism_creators.append(random_organism_creator)
            
    # if the InitialPopulation block is given but left blank or default, throw an error if we're doing a phase diagram search
    elif parameters['InitialPopulation'] == None or parameters['InitialPopulation'] == 'default':
        if composition_space.objective_function == 'pd':
            print('For phase diagram searches, reference structures at each endpoint of the composition space must be provided in the initial population.')
            print('Please use the "from_files" flag in the InitialPopulation block to provide the reference structures.')
            print('Quitting...')
            quit()
        # otherwise, make a random one by default
        else:
            random_organism_creator = organism_creators.RandomOrganismCreator('default', composition_space)
            initial_organism_creators.append(random_organism_creator)
            
    else:
        # the random organism creator
        if 'random' in parameters['InitialPopulation']:
            random_organism_creator = organism_creators.RandomOrganismCreator(parameters['InitialPopulation']['random'], composition_space)
            initial_organism_creators.append(random_organism_creator)
        
        # the from files organism creator
        if 'from_files' not in parameters['InitialPopulation']:
            if composition_space.objective_function == 'pd':
                print('For phase diagram searches, reference structures at each endpoint of the composition space must be provided.')
                print('Please use the "from_files" flag in the InitialPopulation block to provide the reference structures.')
                print('Quitting...')
                quit()  
        # if nothing is given after the from_files flag
        elif parameters['InitialPopulation']['from_files'] == None:
            print('The path to the folder containing the files must be provided. Please use the "path_to_folder" flag.')
            print('Quitting...')
            quit()
        # if path_to_folder flag is not given
        elif 'path_to_folder' not in parameters['InitialPopulation']['from_files']:
            print('Incorrect flag given after "from_files" in the InitialPopulation block. Please use the "path_to_folder" flag.')
            print('Quitting...')
            quit()
        else:
            given_path = parameters['InitialPopulation']['from_files']['path_to_folder']
            # check if no path was given after path_to_folder flag
            if given_path == None:
                print('The path to the folder containing the files for the initial population must be provided. Please give the path after the "path_to_folder" flag.')
                print('Quitting...')
                quit()
            # check that the given path exists
            elif not os.path.exists(given_path):
                print('The given folder containing structures for the initial population does not exist.')
                print('Quitting...')
                quit()
            # if the folder exists, check that it contains files
            elif len([f for f in os.listdir(given_path) if os.path.isfile(os.path.join(given_path, f))]) == 0:
                print('The given folder containing structures for the initial population does not contain any files.')
                print('Quitting...')
                quit()
            else:   
                files_organism_creator = organism_creators.FileOrganismCreator(given_path)
                # check that the provided files cover all the composition space endpoints
                if composition_space.objective_function == 'pd':
                    # all the structures created from the provided files
                    structures = files_organism_creator.get_structures()
                    # list to hold the composition space endpoints for which structures have been provided
                    provided_endpoints = []
                    for endpoint in composition_space.endpoints:
                        for structure in structures:
                            if structure.composition.reduced_composition.almost_equals(endpoint.reduced_composition) and endpoint not in provided_endpoints:
                                provided_endpoints.append(endpoint)
                    # check if we got them all
                    for endpoint in composition_space.endpoints:
                        if endpoint not in provided_endpoints:
                            print('Error: valid structure files not provided to the initial population for all endpoints of the composition space.')
                            print('Quitting...')
                            quit()
                initial_organism_creators.append(files_organism_creator)
        # TODO: if other organism creators are used, they should be instantiated here

    # If more than one organism creator, sort them so that the attempts-based ones are at the front and the successes-based ones are at the back
    if len(initial_organism_creators) > 1:
        initial_organism_creators.sort(key=lambda x: x.is_successes_based)
        
    # put the list of organism creators in the dictionary
    objects_dict['organism_creators'] = initial_organism_creators
    
    # the number of energy calculations to run at a time
    if 'NumCalcsAtOnce' in parameters:
        if parameters['NumCalcsAtOnce'] == None or parameters['NumCalcsAtOnce'] == 'default':
            num_calcs_at_once = 1
        else:
            num_calcs_at_once = parameters['NumCalcsAtOnce']
    else:
        # if no NumCalcsAtOnce block is given in the input file, set it to 1
        num_calcs_at_once = 1
    objects_dict['num_calcs_at_once'] = num_calcs_at_once
    
    # get the run title, if specified
    if 'RunTitle' in parameters:
        if parameters['RunTitle'] == None:
            run_title = ''
        else:
            run_title = parameters['RunTitle']
    else:
        # if no RunTitle block is given in the input file, set it to an empty string
        run_title = ''
    
    # using the run title, construct the name (not path) of the run directory
    if run_title == '':
        run_dir_name = 'garun'
    else:
        run_dir_name = 'garun' + '_' + str(run_title)
    
    # put the name (not path) of the run directory in the dictionary
    objects_dict['run_dir_name'] = run_dir_name

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
            print('No GULP header or potential files given. Please use the "header_file" and "potential_file" flags.')
            print('Quitting...')
            quit()
        else:
            # get the header file
            if 'header_file' not in parameters['EnergyCode']['gulp']:
                print('A GULP header file must be provided. Please use the "header_file" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['gulp']['header_file'] == None:
                print('No GULP header file given after the "header_file" flag. Please provide one.')
                print('Quitting...')
                quit()
            else:
                # get the path to the header file
                header_file_path = parameters['EnergyCode']['gulp']['header_file']
                # check that the given header file exists
                if not os.path.exists(header_file_path):
                    print('The given GULP header file does not exist.')
                    print('Quitting...')
                    quit() 
            # get the potential file 
            if 'potential_file' not in parameters['EnergyCode']['gulp']:
                print('A GULP potential file must be provided. Please use the "potential_file" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['gulp']['potential_file'] == None:
                print('No GULP potential file given after the "potential_file" flag. Please provide one.')
                print('Quitting...')
                quit()   
            else:
                # get the path to the potential file
                potential_file_path = parameters['EnergyCode']['gulp']['potential_file']
                # check that the given potential file exists
                if not os.path.exists(potential_file_path):
                    print('The given GULP potential file does not exist.')
                    print('Quitting...')
                    quit()
            # if we made it this far, then both the header and potential file exist, so make the gulp energy calculator
            energy_calculator = energy_code_interfaces.GulpEnergyCalculator(header_file_path, potential_file_path, geometry)
        
    elif 'lammps' in parameters['EnergyCode']:
        if parameters['EnergyCode']['lammps'] == None:
            print('No LAMMPS input script given. Please use the "input_script" flag.')
            print('Quitting...')
            quit()
        else:
            # get the input script
            if 'input_script' not in parameters['EnergyCode']['lammps']:
                print('A LAMMPS input script must be provided. Please use the "header_file" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['lammps']['input_script'] == None:
                print('No LAMMPS input script given after the "input_script" flag. Please provide one.')
                print('Quitting...')
                quit()
            else:
                # get the path to the input script
                input_script_path = parameters['EnergyCode']['lammps']['input_script']
                # check that the given input script exists
                if not os.path.exists(input_script_path):
                    print('The given LAMMPS input script does not exist.')
                    print('Quitting...')
                    quit()
            # if we made it this far, then the input script exists, so make the lammps energy calculator
            energy_calculator = energy_code_interfaces.LammpsEnergyCalculator(input_script_path, geometry)
            
    elif 'vasp' in parameters['EnergyCode']:
        if parameters['EnergyCode']['vasp'] == None:
            print('No VASP input files given.')
            print('Quitting...')
            quit()
        else:
            # the INCAR file
            if 'incar' not in parameters['EnergyCode']['vasp']:
                print('An INCAR file must be provided. Please use the "incar" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['vasp']['incar'] == None:
                print('No INCAR file was given after the "incar" flag. Please provide one.')
                print('Quitting...')
                quit()
            else:
                # get the path to the INCAR file
                incar_path = parameters['EnergyCode']['vasp']['incar']
                # check that the given INCAR file exists
                if not os.path.exists(incar_path):
                    print('The given INCAR file does not exist.')
                    print('Quitting...')
                    quit()
            # the KPOINTS file
            if 'kpoints' not in parameters['EnergyCode']['vasp']:
                print('A KPOINTS file must be provided. Please use the "kpoints" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['vasp']['kpoints'] == None:
                print('No KPOINTS file was given after the "kpoints" flag. Please provide one.')
                print('Quitting...')
                quit()
            else:
                # get the path to the KPOINTS file
                kpoints_path = parameters['EnergyCode']['vasp']['kpoints']
                # check that the given KPOINTS file exists
                if not os.path.exists(kpoints_path):
                    print('The given KPOINTS file does not exist.')
                    print('Quitting...')
                    quit()
            # the POTCAR files
            if 'potcars' not in parameters['EnergyCode']['vasp']:
                print('POTCAR file(s) must be provided. Please use the "potcars" flag.')
                print('Quitting...')
                quit()
            elif parameters['EnergyCode']['vasp']['potcars'] == None:
                print('No POTCAR files were given after the "potcars" flag. Please provide them.')
                print('Quitting...')
                quit()
            else:
                # get the dictionary containing the paths to the POTCAR files of each element
                potcar_paths = parameters['EnergyCode']['vasp']['potcars']
                # check that enough POTCAR files have been provided
                elements_list = composition_space.get_all_elements()
                if len(potcar_paths) < len(elements_list):
                    print('Not enough POTCAR files provided - one must be given for each element in the composition space. Please provide them.')
                    print('Quitting...')
                    quit()
                # check that each element has been specified below the 'potcars' flag
                for element in elements_list:
                    if element.symbol not in potcar_paths:
                        print('No POTCAR file given for {}. Please provide one.'.format(element.symbol))
                        print('Quitting...')
                        quit()
                # for each element, check that a POTCAR file has been given and that it exists
                for key in potcar_paths:
                    if potcar_paths[key] == None:
                        print('No POTCAR file given for {}. Please provide one.'.format(key))
                        print('Quitting...')
                        quit()
                    elif not os.path.exists(potcar_paths[key]):
                        print('The POTCAR file given for {} does not exist.'.format(key))
                        print('Quitting...')
                        quit()
            # if we made it this far, then the all the files were provided properly and exist, so make the vasp energy calculator
            energy_calculator = energy_code_interfaces.VaspEnergyCalculator(incar_path, kpoints_path, potcar_paths, geometry)
    
    # TODO: add other energy codes here
    else:
        print('The given energy code name is invalid.')
        print('Quitting...')
        quit()
        
    # put the energy calculator in the dictionary
    objects_dict['energy_calculator'] = energy_calculator
    
    # the stopping criteria, if specified
    if 'StoppingCriteria' in parameters:
        # if the entire block is empty or set to default, then set everything to default
        if parameters['StoppingCriteria'] == None or parameters['StoppingCriteria'] == 'default':
            stopping_criteria = general.StoppingCriteria(None, composition_space)
        # need to check that if the found structure option is used and a file path is given, that it's valid
        elif 'found_structure' in parameters['StoppingCriteria']:
            if parameters['StoppingCriteria']['found_structure'] == None or parameters['StoppingCriteria']['found_structure'] == 'default':
                stopping_criteria = general.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
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
                        stopping_criteria = general.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
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
            stopping_criteria = general.StoppingCriteria(parameters['StoppingCriteria'], composition_space)
    # use defaults if no StoppingCriteria block has been used in the input file
    else:
        stopping_criteria = general.StoppingCriteria(None, composition_space)
    
    # put the stopping criteria in the dictionary
    objects_dict['stopping_criteria'] = stopping_criteria

    # list to hold the variation objects
    variations_list = []

    # the default fractions for the variations
    default_structure_mut_fraction = 0.1
    default_num_stoichs_mut_fraction = 0.1
    if len(composition_space.get_all_swappable_pairs()) > 0:
        default_mating_fraction = 0.7
        default_permutation_fraction = 0.1
    else:
        default_mating_fraction = 0.8
        default_permutation_fraction = 0.0

    # TODO: get rid of code duplication below
    # create the variation objects
    if 'Variations' not in parameters:
        # no variations specified, so make default ones
        mating =  variations.Mating({'fraction': default_mating_fraction})
        structure_mut = variations.StructureMut({'fraction': default_structure_mut_fraction})
        num_stoichs_mut = variations.NumStoichsMut({'fraction': default_num_stoichs_mut_fraction})
        permutation = variations.Permutation({'fraction': default_permutation_fraction}, composition_space)
    
        # add them to the list
        variations_list.append(mating)
        variations_list.append(structure_mut)
        variations_list.append(num_stoichs_mut)
        variations_list.append(permutation)
    
    elif parameters['Variations'] == None or parameters['Variations'] == 'default':
        # no variations specified, so make default ones
        mating = variations.Mating({'fraction': default_mating_fraction})
        structure_mut = variations.StructureMut({'fraction': default_structure_mut_fraction})
        num_stoichs_mut = variations.NumStoichsMut({'fraction': default_num_stoichs_mut_fraction})
        permutation = variations.Permutation({'fraction': default_permutation_fraction},  composition_space)
    
        # add them to the list
        variations_list.append(mating)
        variations_list.append(structure_mut)
        variations_list.append(num_stoichs_mut)
        variations_list.append(permutation)
    
    else:
        # make each variation that's been specified
        # mating
        if 'Mating' not in parameters['Variations']:
            # Mating sub block was not used, so don't make a Mating variation
            pass
        elif parameters['Variations']['Mating'] == None:
            # Mating sub block was left blank
            print('If the "Mating" flag is used, its "fraction" flag must also be set.')
            print('Quitting...')
            quit()
        else:
            # Mating sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
            if parameters['Variations']['Mating']['fraction'] == None or parameters['Variations']['Mating']['fraction'] == 'default':
                print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the Mating variation.')
                print('Quitting...')
                quit()
            else:
                # make a Mating variation object with the supplied parameters
                mating = variations.Mating(parameters['Variations']['Mating'])
                variations_list.append(mating) 
               
        # structure mutation
        if 'StructureMut' not in parameters['Variations']:
            # StructureMut sub block was not used, so don't make a StructureMut variation
            pass
        elif parameters['Variations']['StructureMut'] == None:
            # StructureMut sub block was left blank
            print('If the "StructureMut" flag is used, its "fraction" flag must also be set.')
            print('Quitting...')
            quit()
        else:
            # StructureMut sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
            if parameters['Variations']['StructureMut']['fraction'] == None or parameters['Variations']['StructureMut']['fraction'] == 'default':
                print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the StructureMut variation.')
                print('Quitting...')
                quit()
            else:
                # make a StructureMut variation object with the supplied parameters
                structure_mut = variations.StructureMut(parameters['Variations']['StructureMut'])
                variations_list.append(structure_mut)
               
        # mutating the number of stoichiometries worth of atoms in the cell
        if 'NumStoichsMut' not in parameters['Variations']:
            # NumStoichsMut sub block was not used, so don't make a NumStoichsMut variation
            pass
        elif parameters['Variations']['NumStoichsMut'] == None:
            # NumStoichsMut block was left blank
            print('If the "NumStoichsMut" flag is used, its "fraction" flag must also be set.')
            print('Quitting...')
            quit()
        else:
            # NumStoichsMut sub block was used is not empty. Check that the non-optional 'fraction' flag specifies a valid value
            if parameters['Variations']['NumStoichsMut']['fraction'] == None or parameters['Variations']['NumStoichsMut']['fraction'] == 'default':
                print('The "fraction" flag is not optional and must contain a valid entry (between 0 and 1) for the NumStoichsMut variation.')
                print('Quitting...')
                quit()
            else:
                # make a NumStoichsMut variation object with the supplied parameters
                num_stoichs_mut = variations.NumStoichsMut(parameters['Variations']['NumStoichsMut']) 
                variations_list.append(num_stoichs_mut) 
             
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
                permutation = variations.Permutation(parameters['Variations']['Permutation'], composition_space)
                variations_list.append(permutation)
        
    # check that at least one variation has been used. This shouldn't happen...
    if len(variations_list) == 0:
        print("At least one variation must be used. Either leave entire 'Variations' block blank to use default variations, or specify at least one variation within the 'Variations' block.")
        print('Quitting...')
        quit()

    # check that the variations fraction variables sum to 1
    frac_sum = 0.0
    for variation in variations_list:
        frac_sum = frac_sum + variation.fraction
    if frac_sum < 0.999 or frac_sum > 1.001:
        print("The Variations' fraction values must sum to 1.")
        print('Quitting...')
        quit()
        
    # put the list of variations in the dictionary
    objects_dict['variations'] = variations_list

    # create the pool object
    if 'Pool' not in parameters:
        pool = general.Pool(None, composition_space, run_dir_name)
    else:
        pool = general.Pool(parameters['Pool'], composition_space, run_dir_name)
    
    # create the selection object
    if 'Selection' not in parameters:
        selection = general.SelectionProbDist(None, pool.size)
    else:
        selection = general.SelectionProbDist(parameters['Selection'], pool.size)
    
    # give the selection object to the pool
    pool.selection = selection
    
    # put the pool in the dictionary
    objects_dict['pool'] = pool
    
    # return the dictionary containing all the objects
    return objects_dict    