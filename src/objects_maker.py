import os
import classes
from pymatgen.core.structure import Structure


def makeObjects(parameters):
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
        composition_space = classes.CompositionSpace(parameters['CompositionSpace'])
    else:
        print('Input file must contain a "CompositionSpace" block.')
        print("Quitting...")
        quit()
        
    # put the composition space in the dictionary
    objects_dict['composition_space'] = composition_space
    
    # make the constraints object
    if 'Constraints' in parameters:
        constraints = classes.Constraints(parameters['Constraints'], composition_space)
    else:
        # if no Constraints block is given in the input file, then just use default values for everything
        constraints = classes.Constraints('default', composition_space)
        
    # put the constraints in the dictionary
    objects_dict['constraints'] = constraints
    
    # make the geometry object
    if 'Geometry' in parameters:
        geometry = classes.Geometry(parameters['Geometry'])
    else:
        # if no Geometry block is given in the input file, then just use default values for everything
        geometry = classes.Geometry('default')
        
    # put the geometry in the dictionary
    objects_dict['geometry'] = geometry
    
    # make the development object
    if 'Development' in parameters:
        development = classes.Development(parameters['Development'])
    else:
        # if no Development block is given in the input file, then just use default values for everything
        development = classes.Development('default')
        
    # put the development in the dictionary
    objects_dict['development'] = development

    # make the redundancy guard object
    if 'RedundancyGuard' in parameters:
        redundancy_guard = classes.RedundancyGuard(parameters['RedundancyGuard'])
    else:
        redundancy_guard = classes.RedundancyGuard('default')
        
    # put the redundancy guard in the dictionary
    objects_dict['redundancy_guard'] = redundancy_guard
    
    # make the id generator
    id_generator = classes.IDGenerator()
    
    # put the id generator in the dictionary
    objects_dict['id_generator'] = id_generator
    
    # make the organism creators
    organism_creators = []
    # if no method is specified for making the initial population, throw an error if we're doing a phase diagram search
    if 'InitialPopulation' not in parameters: 
        if composition_space.objective_function == 'pd':
            print('For phase diagram searches, reference structures at each endpoint of the composition space must be provided in the initial population.')
            print('Please use the "from_files" flag in the InitialPopulation block to provide the reference structures.')
            print('Quitting...')
            quit()
        # otherwise, make a random one by default
        else:
            random_organism_creator = classes.RandomOrganismCreator('default', composition_space)
            organism_creators.append(random_organism_creator)
            
    # if the InitialPopulation block is given but left blank or default, throw an error if we're doing a phase diagram search
    elif parameters['InitialPopulation'] == None or parameters['InitialPopulation'] == 'default':
        if composition_space.objective_function == 'pd':
            print('For phase diagram searches, reference structures at each endpoint of the composition space must be provided in the initial population.')
            print('Please use the "from_files" flag in the InitialPopulation block to provide the reference structures.')
            print('Quitting...')
            quit()
        # otherwise, make a random one by default
        else:
            random_organism_creator = classes.RandomOrganismCreator('default', composition_space)
            organism_creators.append(random_organism_creator)
            
    else:
        # the random organism creator
        if 'random' in parameters['InitialPopulation']:
            random_organism_creator = classes.RandomOrganismCreator(parameters['InitialPopulation']['random'], composition_space)
            organism_creators.append(random_organism_creator)
        
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
                files_organism_creator = classes.FileOrganismCreator(given_path)
                # check that the provided files cover all the composition space endpoints
                if composition_space.objective_function == 'pd':
                    # all the structures created from the provided files
                    structures = files_organism_creator.getStructures()
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
                organism_creators.append(files_organism_creator)
        # TODO: if other organism creators are used, they should be instantiated here

    # If more than one organism creator, sort them so that the attempts-based ones are at the front and the successes-based ones are at the back
    if len(organism_creators) > 1:
        organism_creators.sort(key=lambda x: x.is_successes_based)
        
    # put the list of organism creators in the dictionary
    objects_dict['organism_creators'] = organism_creators
    
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
            energy_calculator = classes.GulpEnergyCalculator(header_file_path, potential_file_path)
        
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
            energy_calculator = classes.LammpsEnergyCalculator(input_script_path)
       
            
    elif 'vasp' in parameters['EnergyCode']:
        # TODO: read in the stuff for vasp calcs
        pass
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
    
    # put the stopping criteria in the dictionary
    objects_dict['stopping_criteria'] = stopping_criteria

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
                mating = classes.Mating(parameters['Variations']['Mating'])
                variations.append(mating) 
               
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
                structure_mut = classes.StructureMut(parameters['Variations']['StructureMut'])
                variations.append(structure_mut)
               
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
    frac_sum = 0.0
    for variation in variations:
        frac_sum = frac_sum + variation.fraction
    if frac_sum < 0.999 or frac_sum > 1.001:
        print("The Variations' fraction values must sum to 1.")
        print('Quitting...')
        quit()
        
    # put the list of variations in the dictionary
    objects_dict['variations'] = variations

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
    
    # put the pool in the dictionary
    objects_dict['pool'] = pool
    
    # return the dictionary containing all the objects
    return objects_dict


def printParameters(objects_dict):
    '''
    Prints out the parameters for the search to a file called 'search_parameters' inside the run directory
    
    Args:
        objects_dict: a dictionary of objects used by the algorithm, as returned by the make_objects method
    '''
    # get all the objects from the dictionary
    run_dir_name = objects_dict['run_dir_name']
    organism_creators = objects_dict['organism_creators']
    num_calcs_at_once = objects_dict['num_calcs_at_once']
    composition_space = objects_dict['composition_space']
    development = objects_dict['development']
    constraints = objects_dict['constraints']
    geometry = objects_dict['geometry']
    redundancy_guard = objects_dict['redundancy_guard']
    stopping_criteria = objects_dict['stopping_criteria']
    energy_calculator = objects_dict['energy_calculator']
    pool = objects_dict['pool']
    variations = objects_dict['variations']
    
    # make the file where the parameters will be printed
    parameters_file = open(os.getcwd() + '/search_parameters', 'w')
    
    # write the directory where the search will be done
    parameters_file.write('search output directory: ' + run_dir_name + '\n')
    parameters_file.write('\n')
    
    # write the endpoints of the composition space
    parameters_file.write('CompositionSpace: \n')
    for endpoint in composition_space.endpoints:
        parameters_file.write('    - ' + endpoint.reduced_formula + '\n')
    parameters_file.write('\n')
        
    # write the objective function
    parameters_file.write('ObjectiveFunction: ' + composition_space.objective_function + '\n')
    parameters_file.write('\n')
    
    # write the name of the energy code being used
    parameters_file.write('EnergyCode: ' + energy_calculator.name + '\n')
    parameters_file.write('\n')
    
    # write the number of energy calculations to run at once
    parameters_file.write('NumCalcsAtOnce: ' + str(num_calcs_at_once) + '\n')
    parameters_file.write('\n')
    
    # write the methods used to create the initial population
    parameters_file.write('InitialPopulation: \n')
    for creator in organism_creators:
        if creator.name == 'random organism creator':
            parameters_file.write('    random: \n')
            parameters_file.write('        number: ' + str(creator.number) + '\n')
            parameters_file.write('        volume: ' + str(creator.volume) + '\n')
        elif creator.name == 'file organism creator':
            parameters_file.write('    from_files: \n')
            parameters_file.write('        number: ' + str(creator.number) + '\n')
            parameters_file.write('        path_to_folder: ' + str(creator.path_to_folder) + '\n')
    parameters_file.write('\n')
    
    # write the selection probability distribution
    parameters_file.write('Selection: \n')
    parameters_file.write('    num_parents: ' + str(pool.selection.num_parents) + '\n')
    parameters_file.write('    power: ' + str(pool.selection.power) + '\n')
    parameters_file.write('\n')
    
    # write the pool info
    parameters_file.write('Pool: \n')
    parameters_file.write('    size: ' + str(pool.size) + '\n')
    parameters_file.write('    num_promoted: ' + str(pool.num_promoted) + '\n')
    parameters_file.write('\n')
    
    # write the variations info
    parameters_file.write('Variations: \n')
    for variation in variations:
        # don't print out info for variations with zero fraction since they aren't being used
        if variation.fraction == 0:
            pass
        else:
            if variation.name == 'mating':
                parameters_file.write('    Mating: \n')
                parameters_file.write('        fraction: ' + str(variation.fraction) + '\n')
                parameters_file.write('        mu_cut_loc: ' + str(variation.mu_cut_loc) + '\n')
                parameters_file.write('        sigma_cut_loc: ' + str(variation.sigma_cut_loc) + '\n')
                parameters_file.write('        shift_prob: ' + str(variation.shift_prob) + '\n')
                parameters_file.write('        doubling_prob: ' + str(variation.doubling_prob) + '\n')
                parameters_file.write('        grow_parents: ' + str(variation.grow_parents) + '\n')
                parameters_file.write('        merge_sites: ' + str(variation.merge_sites) + '\n')
            
            elif variation.name == 'structure mutation':
                parameters_file.write('    StructureMut: \n')
                parameters_file.write('        fraction: ' + str(variation.fraction) + '\n')
                parameters_file.write('        frac_atoms_perturbed: ' + str(variation.frac_atoms_perturbed) + '\n')
                parameters_file.write('        sigma_atomic_coord_perturbation: ' + str(variation.sigma_atomic_coord_perturbation) + '\n')
                parameters_file.write('        max_atomic_coord_perturbation: ' + str(variation.max_atomic_coord_perturbation) + '\n')
                parameters_file.write('        sigma_strain_matrix_element: ' + str(variation.sigma_strain_matrix_element) + '\n')
            
            elif variation.name == 'number of stoichiometries mutation':
                parameters_file.write('    NumStoichsMut: \n')
                parameters_file.write('        fraction: ' + str(variation.fraction) + '\n')
                parameters_file.write('        mu_num_adds: ' + str(variation.mu_num_adds) + '\n')
                parameters_file.write('        sigma_num_adds: ' + str(variation.sigma_num_adds) + '\n')
                parameters_file.write('        scale_volume: ' + str(variation.scale_volume) + '\n')
            
            elif variation.name == 'permutation':
                parameters_file.write('    Permutation: \n')
                parameters_file.write('        fraction: ' + str(variation.fraction) + '\n')
                parameters_file.write('        mu_num_swaps: ' + str(variation.mu_num_swaps) + '\n')
                parameters_file.write('        sigma_num_swaps = ' + str(variation.sigma_num_swaps) + '\n')
                parameters_file.write('        pairs_to_swap: \n')
                for pair in variation.pairs_to_swap:
                    parameters_file.write('        ' + pair + '\n')
    parameters_file.write('\n')
    
    # write the development 
    parameters_file.write('Development: \n')
    parameters_file.write('    niggli: ' + str(development.niggli) + '\n')
    parameters_file.write('    scale_density: ' + str(development.scale_density) + '\n')
    parameters_file.write('\n')
    
    # write the constraints
    parameters_file.write('Constraints: \n')
    parameters_file.write('    min_num_atoms: ' + str(constraints.min_num_atoms) + '\n')
    parameters_file.write('    max_num_atoms: ' + str(constraints.max_num_atoms) + '\n')
    parameters_file.write('    min_lattice_length: ' + str(constraints.min_lattice_length) + '\n')
    parameters_file.write('    max_lattice_length: ' + str(constraints.max_lattice_length) + '\n')
    parameters_file.write('    min_lattice_angle: ' + str(constraints.min_lattice_angle) + '\n')
    parameters_file.write('    max_lattice_angle: ' + str(constraints.max_lattice_angle) + '\n')
    parameters_file.write('    per_species_mids: \n')
    for pair in constraints.per_species_mids:
        parameters_file.write('        ' + pair + ': ' + str(float(constraints.per_species_mids[pair])) + '\n')
    parameters_file.write('\n')
    
    # write the redundancy guard
    parameters_file.write('RedundancyGuard: \n')
    parameters_file.write('    lattice_length_tol: ' + str(redundancy_guard.lattice_length_tol) + '\n')
    parameters_file.write('    lattice_angle_tol: ' + str(redundancy_guard.lattice_angle_tol) + '\n')
    parameters_file.write('    site_tol: ' + str(redundancy_guard.site_tol) + '\n')
    parameters_file.write('    use_primitive_cell: ' + str(redundancy_guard.use_primitive_cell) + '\n')
    parameters_file.write('    attempt_supercell: ' + str(redundancy_guard.attempt_supercell) + '\n')
    parameters_file.write('    epa_diff: ' + str(redundancy_guard.epa_diff) + '\n')
    parameters_file.write('\n')
    
    # write the geometry
    parameters_file.write('Geometry: \n')
    parameters_file.write('    shape: ' + geometry.shape + '\n')
    parameters_file.write('    max_size: ' + str(geometry.max_size) + '\n')
    parameters_file.write('    min_size: ' + str(geometry.min_size) + '\n')
    parameters_file.write('    padding: ' + str(geometry.padding) + '\n')
    parameters_file.write('\n')
    
    # write the stopping criteria
    parameters_file.write('StoppingCriteria: \n')
    if stopping_criteria.num_energy_calcs != None:
        parameters_file.write('    num_energy_calcs: ' + str(stopping_criteria.num_energy_calcs) + '\n')
    if stopping_criteria.value_achieved != None:
        parameters_file.write('    value_achieved: ' + str(stopping_criteria.value_achieved) + '\n')
    if stopping_criteria.found_structure != None:
        parameters_file.write('    found_structure: ' + stopping_criteria.path_to_structure_file + '\n')
    parameters_file.write('\n')
    