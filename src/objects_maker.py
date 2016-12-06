# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Objects Maker module:

This module contains functions for creating the (singleton) objects
used by the genetic algorithm during the search.

"""

import general
import variations
import energy_code_interfaces
import organism_creators
import development

from pymatgen.core.structure import Structure

import os


def make_objects(parameters):
    """
    Constructs the needed objects for the genetic algorithm search.

    Returns a dictionary containing the objects.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file
    """

    # to hold all the objects
    objects_dict = {}

    # make the composition space object
    if 'CompositionSpace' in parameters:
        composition_space = general.CompositionSpace(
            parameters['CompositionSpace'])
    else:
        print('Input file must contain a "CompositionSpace" block.')
        print("Quitting...")
        quit()

    objects_dict['composition_space'] = composition_space

    # make the constraints object
    if 'Constraints' in parameters:
        constraints = development.Constraints(parameters['Constraints'],
                                              composition_space)
    else:
        constraints = development.Constraints('default', composition_space)

    objects_dict['constraints'] = constraints

    # make the geometry object
    if 'Geometry' in parameters:
        geometry = development.Geometry(parameters['Geometry'])
    else:
        geometry = development.Geometry('default')

    objects_dict['geometry'] = geometry

    # make the development object
    if 'Development' in parameters:
        developer = development.Developer(parameters['Development'], geometry)
    else:
        developer = development.Developer('default', geometry)

    objects_dict['developer'] = developer

    # make the redundancy guard object
    if 'RedundancyGuard' in parameters:
        redundancy_guard = development.RedundancyGuard(
            parameters['RedundancyGuard'])
    else:
        redundancy_guard = development.RedundancyGuard('default')

    objects_dict['redundancy_guard'] = redundancy_guard

    # make the id generator
    id_generator = general.IDGenerator()
    objects_dict['id_generator'] = id_generator

    # make the organism creators
    initial_organism_creators = make_organism_creators(parameters,
                                                       composition_space)

    # if more than one organism creator, sort them so that the attempts-based
    # ones are at the front and the successes-based ones are at the back
    if len(initial_organism_creators) > 1:
        initial_organism_creators.sort(key=lambda x: x.is_successes_based)

    objects_dict['organism_creators'] = initial_organism_creators

    # the number of energy calculations to run at a time
    if 'NumCalcsAtOnce' not in parameters:
        num_calcs_at_once = 1
    elif parameters['NumCalcsAtOnce'] in (None, 'default'):
        num_calcs_at_once = 1
    else:
        num_calcs_at_once = parameters['NumCalcsAtOnce']

    objects_dict['num_calcs_at_once'] = num_calcs_at_once

    # get the run title
    if 'RunTitle' not in parameters:
        run_dir_name = 'garun'
    elif parameters['RunTitle'] in (None, 'default'):
        run_dir_name = 'garun'
    else:
        run_dir_name = 'garun_' + str(parameters['RunTitle'])

    objects_dict['run_dir_name'] = run_dir_name

    # make the energy calculator
    energy_calculator = make_energy_calculator(parameters, geometry,
                                               composition_space)
    objects_dict['energy_calculator'] = energy_calculator

    # make the stopping criteria
    stopping_criteria = make_stopping_criteria(parameters, composition_space)
    objects_dict['stopping_criteria'] = stopping_criteria

    # default fractions for the variations
    default_variation_fractions = {}
    default_variation_fractions['structure_mut'] = 0.1
    default_variation_fractions['num_stoichs_mut'] = 0.1
    if len(composition_space.get_all_swappable_pairs()) > 0:
        default_variation_fractions['mating'] = 0.7
        default_variation_fractions['permutation'] = 0.1
    else:
        default_variation_fractions['mating'] = 0.8
        default_variation_fractions['permutation'] = 0.0

    # make the variations
    variations_list = make_variations(parameters, default_variation_fractions,
                                      composition_space)

    # check that at least one variation has been used
    if len(variations_list) == 0:
        print('At least one variation must be used. Either leave entire '
              '"Variations" block blank to use default variations, or specify '
              'at least one variation within the "Variations" block.')
        print('Quitting...')
        quit()

    # check that the variations' fraction variables sum to 1
    frac_sum = 0.0
    for variation in variations_list:
        frac_sum = frac_sum + variation.fraction
    if frac_sum < 0.999 or frac_sum > 1.001:
        print("The Variations' fraction values must sum to 1.")
        print('Quitting...')
        quit()

    objects_dict['variations'] = variations_list

    # make the pool and selection
    if 'Pool' not in parameters:
        pool = general.Pool(None, composition_space, run_dir_name)
    else:
        pool = general.Pool(parameters['Pool'], composition_space,
                            run_dir_name)

    if 'Selection' not in parameters:
        selection = general.SelectionProbDist(None, pool.size)
    else:
        selection = general.SelectionProbDist(parameters['Selection'],
                                              pool.size)

    pool.selection = selection
    objects_dict['pool'] = pool

    return objects_dict


def make_organism_creators(parameters, composition_space):
    """
    Returns a list containing organism creator objects.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        composition_space: the CompositionSpace of the search
    """

    if 'InitialPopulation' not in parameters:
        return make_default_organism_creator(composition_space)
    elif parameters['InitialPopulation'] in (None, 'default'):
        return make_default_organism_creator(composition_space)
    # make the specified creators
    else:
        initial_organism_creators = []

        # the random organism creator
        if 'random' in parameters['InitialPopulation']:
            random_organism_creator = organism_creators.RandomOrganismCreator(
                parameters['InitialPopulation']['random'], composition_space)
            initial_organism_creators.append(random_organism_creator)

        # the from files organism creator
        if 'from_files' not in parameters['InitialPopulation']:
            if composition_space.objective_function == 'pd':
                print('For phase diagram searches, reference structures at '
                      'each endpoint of the composition space must be '
                      'provided.')
                print('Please use the "from_files" keyword in the '
                      'InitialPopulation block to provide the reference '
                      'structures.')
                print('Quitting...')
                quit()
        # if nothing is given after the from_files keyword
        elif parameters['InitialPopulation']['from_files'] is None:
            print('The path to the folder containing the files must be '
                  'provided. Please use the "path_to_folder" keyword.')
            print('Quitting...')
            quit()
        # if path_to_folder keyword is not given
        elif 'path_to_folder' not in parameters['InitialPopulation'][
                'from_files']:
            print('Incorrect keyword given after "from_files" in the '
                  'InitialPopulation block. Please use the "path_to_folder" '
                  'keyword.')
            print('Quitting...')
            quit()
        else:
            given_path = parameters['InitialPopulation']['from_files'][
                    'path_to_folder']
            # if no path was given after path_to_folder keyword
            if given_path is None:
                print('The path to the folder containing the files for the '
                      'initial population must be provided. Please give the '
                      'path after the "path_to_folder" keyword.')
                print('Quitting...')
                quit()
            # if the given path does not exist
            elif not os.path.exists(given_path):
                print('The given folder containing structures for the initial '
                      'population does not exist.')
                print('Quitting...')
                quit()
            # if the folder exists, check that it contains files
            elif len([f for f in os.listdir(given_path) if
                      os.path.isfile(os.path.join(given_path, f))]) == 0:
                print('The given folder containing structures for the initial '
                      'population does not contain any files.')
                print('Quitting...')
                quit()
            else:
                files_organism_creator = organism_creators.FileOrganismCreator(
                    given_path)
                # check that the files cover all composition space endpoints
                if composition_space.objective_function == 'pd':
                    structures = files_organism_creator.get_structures()
                    provided_endpoints = []
                    for endpoint in composition_space.endpoints:
                        for structure in structures:
                            if structure.composition.reduced_composition.almost_equals(
                                    endpoint.reduced_composition) and \
                                    endpoint not in provided_endpoints:
                                provided_endpoints.append(endpoint)
                    # check if we got them all
                    for endpoint in composition_space.endpoints:
                        if endpoint not in provided_endpoints:
                            print('Error: valid structure files not provided '
                                  'to the initial population for all '
                                  'endpoints of the composition space.')
                            print('Quitting...')
                            quit()
                initial_organism_creators.append(files_organism_creator)

        # TODO: if other organism creators are used, they should be
        # instantiated here

        return initial_organism_creators


def make_default_organism_creator(composition_space):
    """
    Returns a list containing a RandomOrganismCreator, or quits.

    Args:
        composition_space: the CompositionSpace of the search
    """

    if composition_space.objective_function == 'pd':
            print('For phase diagram searches, reference structures at each '
                  'endpoint of the composition space must be provided in the '
                  'initial population.')
            print('Please use the "from_files" keyword in the '
                  'InitialPopulation block to provide the reference '
                  'structures.')
            print('Quitting...')
            quit()
    else:
        random_organism_creator = organism_creators.RandomOrganismCreator(
            'default', composition_space)
        return [random_organism_creator]


def make_energy_calculator(parameters, geometry, composition_space):
    """
    Returns an EnergyCode object corresponding to which energy code was
    specified in the input file. Quits if an energy code object cannot be made.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        geometry: the Geometry for the search

        composition_space: the CompositionSpace of the search
    """

    if 'EnergyCode' not in parameters:
        print('A method for calculating energy must be provided. Please use '
              'the "EnergyCode" keyword.')
        print('Quitting...')
        quit()
    elif parameters['EnergyCode'] is None:
        print('An energy code must be specified after the "EnergyCode" '
              'keyword.')
        print('Quitting...')
        quit()
    # for GULP
    elif 'gulp' in parameters['EnergyCode']:
        return make_gulp_energy_calculator(parameters, geometry)
    # for LAMMPS
    elif 'lammps' in parameters['EnergyCode']:
        return make_lammps_energy_calculator(parameters, geometry)
    # for VASP
    elif 'vasp' in parameters['EnergyCode']:
        return make_vasp_energy_calculator(parameters, composition_space,
                                           geometry)
    else:
        print('The given energy code name is invalid.')
        print('Quitting...')
        quit()


def make_gulp_energy_calculator(parameters, geometry):
    """
    Returns a GulpEnergyCalculator object, or quits if one cannot be made.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        geometry: the Geometry for the search
    """

    if parameters['EnergyCode']['gulp'] is None:
        print('No GULP header or potential files given. Please use the '
              '"header_file" and "potential_file" keywords.')
        print('Quitting...')
        quit()
    else:
        # get the header file
        if 'header_file' not in parameters['EnergyCode']['gulp']:
            print('A GULP header file must be provided. Please use the '
                  '"header_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['header_file'] is None:
            print('No GULP header file given after the "header_file" '
                  'keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            # get the path to the header file
            header_file_path = parameters['EnergyCode']['gulp'][
                'header_file']
            # check that the header file exists
            if not os.path.exists(header_file_path):
                print('The given GULP header file does not exist.')
                print('Quitting...')
                quit()
        # get the potential file
        if 'potential_file' not in parameters['EnergyCode']['gulp']:
            print('A GULP potential file must be provided. Please use the '
                  '"potential_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['gulp']['potential_file'] is None:
            print('No GULP potential file given after the '
                  '"potential_file" keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            # get the path to the potential file
            potential_file_path = parameters['EnergyCode']['gulp'][
                'potential_file']
            # check that the potential file exists
            if not os.path.exists(potential_file_path):
                print('The given GULP potential file does not exist.')
                print('Quitting...')
                quit()

        return energy_code_interfaces.GulpEnergyCalculator(
            header_file_path, potential_file_path, geometry)


def make_lammps_energy_calculator(parameters, geometry):
    """
    Returns a LammpsEnergyCalculator object, or quits if one cannot be made.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        geometry: the Geometry for the search
    """

    if parameters['EnergyCode']['lammps'] is None:
        print('No LAMMPS input script given. Please use the "input_script" '
              'keyword.')
        print('Quitting...')
        quit()
    else:
        # get the input script
        if 'input_script' not in parameters['EnergyCode']['lammps']:
            print('A LAMMPS input script must be provided. Please use the '
                  '"header_file" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['lammps']['input_script'] is None:
            print('No LAMMPS input script given after the "input_script" '
                  'keyword. Please provide one.')
            print('Quitting...')
            quit()
        else:
            # get the path to the input script
            input_script_path = parameters['EnergyCode']['lammps'][
                'input_script']
            # check that the input script exists
            if not os.path.exists(input_script_path):
                print('The given LAMMPS input script does not exist.')
                print('Quitting...')
                quit()

        return energy_code_interfaces.LammpsEnergyCalculator(
                input_script_path, geometry)


def make_vasp_energy_calculator(parameters, composition_space, geometry):
    """
    Returns a VaspEnergyCalculator object, or quits if one cannot be made.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        composition_space: the CompositionSpace of the search

        geometry: the Geometry for the search
    """

    if parameters['EnergyCode']['vasp'] is None:
        print('No VASP input files given.')
        print('Quitting...')
        quit()
    else:
        # the INCAR file
        if 'incar' not in parameters['EnergyCode']['vasp']:
            print('An INCAR file must be provided. Please use the "incar" '
                  'keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['incar'] is None:
            print('No INCAR file was given after the "incar" keyword. Please '
                  'provide one.')
            print('Quitting...')
            quit()
        else:
            # get the path to the INCAR file
            incar_path = parameters['EnergyCode']['vasp']['incar']
            # check that the INCAR file exists
            if not os.path.exists(incar_path):
                print('The given INCAR file does not exist.')
                print('Quitting...')
                quit()
        # the KPOINTS file
        if 'kpoints' not in parameters['EnergyCode']['vasp']:
            print('A KPOINTS file must be provided. Please use the '
                  '"kpoints" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['kpoints'] is None:
            print('No KPOINTS file was given after the "kpoints" keyword. '
                  'Please provide one.')
            print('Quitting...')
            quit()
        else:
            # get the path to the KPOINTS file
            kpoints_path = parameters['EnergyCode']['vasp']['kpoints']
            # check that the KPOINTS file exists
            if not os.path.exists(kpoints_path):
                print('The given KPOINTS file does not exist.')
                print('Quitting...')
                quit()
        # the POTCAR files
        if 'potcars' not in parameters['EnergyCode']['vasp']:
            print('POTCAR file(s) must be provided. Please use the '
                  '"potcars" keyword.')
            print('Quitting...')
            quit()
        elif parameters['EnergyCode']['vasp']['potcars'] is None:
            print('No POTCAR files were given after the "potcars" keyword. '
                  'Please provide them.')
            print('Quitting...')
            quit()
        else:
            # get the the paths to the POTCAR files of each element
            potcar_paths = parameters['EnergyCode']['vasp']['potcars']
            # check that enough POTCAR files have been provided
            elements_list = composition_space.get_all_elements()
            if len(potcar_paths) < len(elements_list):
                print('Not enough POTCAR files provided - one must be '
                      'given for each element in the composition space. '
                      'Please provide them.')
                print('Quitting...')
                quit()
            # check that each element has been specified below the
            # 'potcars' keyword
            for element in elements_list:
                if element.symbol not in potcar_paths:
                    print('No POTCAR file given for {}. Please provide '
                          'one.'.format(element.symbol))
                    print('Quitting...')
                    quit()
            # for each element, check that a POTCAR file has been given and
            # that it exists
            for key in potcar_paths:
                if potcar_paths[key] is None:
                    print('No POTCAR file given for {}. Please provide '
                          'one.'.format(key))
                    print('Quitting...')
                    quit()
                elif not os.path.exists(potcar_paths[key]):
                    print('The POTCAR file given for {} does not '
                          'exist.'.format(key))
                    print('Quitting...')
                    quit()

        return energy_code_interfaces.VaspEnergyCalculator(
                incar_path, kpoints_path, potcar_paths, geometry)


def make_stopping_criteria(parameters, composition_space):
    """
    Returns a StoppingCriteria object.

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        composition_space: the CompositionSpace of the search
    """

    if 'StoppingCriteria' not in parameters:
        return general.StoppingCriteria(None, composition_space)
    elif parameters['StoppingCriteria'] in (None, 'default'):
        return general.StoppingCriteria(None, composition_space)
    elif 'found_structure' in parameters['StoppingCriteria']:
        if parameters['StoppingCriteria']['found_structure'] in (None,
                                                                 'default'):
            return general.StoppingCriteria(parameters['StoppingCriteria'],
                                            composition_space)
        else:
            # check that the file exists
            given_path = parameters['StoppingCriteria']['found_structure']
            if not os.path.exists(given_path):
                print('The file containing the structure to find does not '
                      'exist.')
                print('Quitting...')
                quit()
            # check that the file has the correct suffix or prefix
            elif not (given_path.endswith('.cif') or given_path.startswith(
                    'POSCAR')):
                print('File containing structure to find must be in POSCAR or '
                      'cif format and begin with POSCAR. or end with .cif, '
                      'respectively.')
                print('Quitting...')
                quit()
            # check that file can be read properly
            else:
                try:
                    Structure.from_file(given_path)
                    return general.StoppingCriteria(
                        parameters['StoppingCriteria'], composition_space)
                except ValueError:
                    print('Error reading the structure to find from the given '
                          'file.')
                    print('Quitting...')
                    quit()
    else:
        return general.StoppingCriteria(parameters['StoppingCriteria'],
                                        composition_space)


def make_variations(parameters, default_fractions, composition_space):
    """
    Creates the variations, using default parameter values if needed.

    Returns a list containing the variation objects (Mating, StructureMut,
    NumStoichsMut and Permutation).

    Args:
        parameters: the dictionary produced by calling yaml.load() on the input
            file

        default_fractions: a dictionary containing the default fractions to use
             for each variation

        composition_space: the CompositionSpace of the search
    """

    if 'Variations' not in parameters:
        return make_default_variations(default_fractions, composition_space)
    elif parameters['Variations'] in (None, 'default'):
        return make_default_variations(default_fractions, composition_space)
    else:
        variations_list = []
        # mating
        if 'Mating' not in parameters['Variations']:
            pass
        elif parameters['Variations']['Mating'] is None:
            print('If the "Mating" keyword is used, its "fraction" keyword '
                  'must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['Mating']['fraction'] in (None,
                                                                  'default'):
                print('The "fraction" kwyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the Mating '
                      'variation.')
                print('Quitting...')
                quit()
            else:
                mating = variations.Mating(parameters['Variations']['Mating'])
                variations_list.append(mating)

        # structure mutation
        if 'StructureMut' not in parameters['Variations']:
            pass
        elif parameters['Variations']['StructureMut'] is None:
            print('If the "StructureMut" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['StructureMut']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'StructureMut variation.')
                print('Quitting...')
                quit()
            else:
                structure_mut = variations.StructureMut(
                    parameters['Variations']['StructureMut'])
                variations_list.append(structure_mut)

        # mutating the number of stoichiometries worth of atoms in the cell
        if 'NumStoichsMut' not in parameters['Variations']:
            pass
        elif parameters['Variations']['NumStoichsMut'] is None:
            print('If the "NumStoichsMut" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['NumStoichsMut']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'NumStoichsMut variation.')
                print('Quitting...')
                quit()
            else:
                num_stoichs_mut = variations.NumStoichsMut(
                    parameters['Variations']['NumStoichsMut'])
                variations_list.append(num_stoichs_mut)

        # permutation (swapping atoms)
        if 'Permutation' not in parameters['Variations']:
            pass
        elif parameters['Variations']['Permutation'] is None:
            print('If the "Permutation" keyword is used, its "fraction" '
                  'keyword must also be set.')
            print('Quitting...')
            quit()
        else:
            if parameters['Variations']['Permutation']['fraction'] in (
                    None, 'default'):
                print('The "fraction" keyword is not optional and must '
                      'contain a valid entry (between 0 and 1) for the '
                      'Permutation variation.')
                print('Quitting...')
            else:
                permutation = variations.Permutation(
                    parameters['Variations']['Permutation'], composition_space)
                variations_list.append(permutation)

        return variations_list


def make_default_variations(default_fractions, composition_space):
    """
    Creates the variations with default parameter values and the provided
    default fractions.

    Returns a list containing the variation objects (Mating, StructureMut,
    NumStoichsMut and Permutation).

    Args:
        default_fractions: a dictionary containing the default fractions to use
             for each variation

        composition_space: the CompositionSpace of the search
    """

    variations_list = []
    mating = variations.Mating({'fraction': default_fractions['mating']})
    structure_mut = variations.StructureMut(
                {'fraction': default_fractions['structure_mut']})
    num_stoichs_mut = variations.NumStoichsMut(
                {'fraction': default_fractions['num_stoichs_mut']})
    permutation = variations.Permutation(
                {'fraction': default_fractions['permutation']},
                composition_space)
    variations_list.append(mating)
    variations_list.append(structure_mut)
    variations_list.append(num_stoichs_mut)
    variations_list.append(permutation)
    return variations_list
