import os
import classes
from pymatgen.core.structure import Structure


def make_objects(parameters):
    '''
    Constructs the needed objects for the genetic algorithm search
    
    Returns a dictionary containing the objects
    
    Args:
        parameters: the dictionary produced by calling yaml.load on the input file
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
    
    # make the niggli object
    if 'Niggli' in parameters:
        if parameters['Niggli'] == None or parameters['Niggli'] == 'default':
            niggli = True
        else:
            niggli = parameters['Niggli']
    else:
        # if no Niggli block is given in the input file, set it to True
        niggli = True
        
    # put the niggli reduction in the dictionary
    objects_dict['niggli'] = niggli
    
    # make the scale density object
    if 'ScaleDesnsity' in parameters:
        if parameters['ScaleDensity'] == None or parameters['ScaleDensity'] == 'default':
            scale_density = True
        else:
            scale_density = parameters['ScaleDensity']
    else:
        # if no ScaleDensity block is given in the input file, set it to True
        scale_density = True
        
    # put the scale density in the dictionary
    objects_dict['scale_density'] = scale_density
   
    # make the development object
    development = classes.Development(niggli, scale_density)
    
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
    
    # using the run title, construct the name of the run directory
    if run_title == '':
        run_dir_name = 'garun'
    else:
        run_dir_name = 'garun' + '_' + str(run_title)
    
    # put the name of the run directory in the dictionary
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