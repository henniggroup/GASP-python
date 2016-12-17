# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Parameters Printer module:

This module contains a function for printing the parameters
of a structure search to a file for the user's reference.

"""

import os


def print_parameters(objects_dict):
    """
    Prints out the parameters for the search to a file called 'ga_parameters'
    inside the garun directory.

    Args:
        objects_dict: a dictionary of objects used by the algorithm, as
            returned by the make_objects method
    """

    # get all the objects from the dictionary
    run_dir_name = objects_dict['run_dir_name']
    organism_creators = objects_dict['organism_creators']
    num_calcs_at_once = objects_dict['num_calcs_at_once']
    composition_space = objects_dict['composition_space']
    developer = objects_dict['developer']
    constraints = objects_dict['constraints']
    geometry = objects_dict['geometry']
    redundancy_guard = objects_dict['redundancy_guard']
    stopping_criteria = objects_dict['stopping_criteria']
    energy_calculator = objects_dict['energy_calculator']
    pool = objects_dict['pool']
    variations = objects_dict['variations']

    # make the file where the parameters will be printed
    with open(os.getcwd() + '/ga_parameters', 'w') as parameters_file:

        # write the directory where the search will be done
        parameters_file.write('search output directory: ' + run_dir_name +
                              '\n')
        parameters_file.write('\n')

        # write the endpoints of the composition space
        parameters_file.write('CompositionSpace: \n')
        for endpoint in composition_space.endpoints:
            parameters_file.write('    - ' +
                                  endpoint.reduced_formula.replace(' ', '') +
                                  '\n')
        parameters_file.write('\n')

        # write the name of the energy code being used
        parameters_file.write('EnergyCode: \n')
        parameters_file.write('    ' + energy_calculator.name + ': \n')
        if energy_calculator.name == 'gulp':
            parameters_file.write('        header_file: ' +
                                  energy_calculator.header_path + '\n')
            parameters_file.write('        potential_file: ' +
                                  energy_calculator.potential_path + '\n')
        elif energy_calculator.name == 'lammps':
            parameters_file.write('        input_script: ' +
                                  energy_calculator.input_script + '\n')
        elif energy_calculator.name == 'vasp':
            parameters_file.write('        incar: ' +
                                  energy_calculator.incar_file + '\n')
            parameters_file.write('        kpoints: ' +
                                  energy_calculator.kpoints_file + '\n')
            parameters_file.write('        potcars: \n')
            for key in energy_calculator.potcar_files:
                parameters_file.write('            ' + key + ': ' +
                                      energy_calculator.potcar_files[key] +
                                      '\n')
        parameters_file.write('\n')

        # write the number of energy calculations to run at once
        parameters_file.write('NumCalcsAtOnce: ' + str(num_calcs_at_once) +
                              '\n')
        parameters_file.write('\n')

        # write the methods used to create the initial population
        parameters_file.write('InitialPopulation: \n')
        for creator in organism_creators:
            if creator.name == 'random organism creator':
                parameters_file.write('    random: \n')
                parameters_file.write('        number: ' +
                                      str(creator.number) + '\n')
                parameters_file.write('        volume: ' +
                                      str(creator.volume) + '\n')
            elif creator.name == 'file organism creator':
                parameters_file.write('    from_files: \n')
                parameters_file.write('        number: ' +
                                      str(creator.number) + '\n')
                parameters_file.write('        path_to_folder: ' +
                                      str(creator.path_to_folder) + '\n')
        parameters_file.write('\n')

        # write the pool info
        parameters_file.write('Pool: \n')
        parameters_file.write('    size: ' + str(pool.size) + '\n')
        parameters_file.write('    num_promoted: ' + str(pool.num_promoted) +
                              '\n')
        parameters_file.write('\n')

        # write the selection probability distribution
        parameters_file.write('Selection: \n')
        parameters_file.write('    num_parents: ' +
                              str(pool.selection.num_parents) + '\n')
        parameters_file.write('    power: ' + str(pool.selection.power) + '\n')
        parameters_file.write('\n')

        # write the variations info
        parameters_file.write('Variations: \n')
        for variation in variations:
            if variation.fraction == 0:
                pass
            else:
                if variation.name == 'mating':
                    parameters_file.write('    Mating: \n')
                    parameters_file.write('        fraction: ' +
                                          str(variation.fraction) + '\n')
                    parameters_file.write('        mu_cut_loc: ' +
                                          str(variation.mu_cut_loc) + '\n')
                    parameters_file.write('        sigma_cut_loc: ' +
                                          str(variation.sigma_cut_loc) + '\n')
                    parameters_file.write('        shift_prob: ' +
                                          str(variation.shift_prob) + '\n')
                    parameters_file.write('        doubling_prob: ' +
                                          str(variation.doubling_prob) + '\n')
                    parameters_file.write('        grow_parents: ' +
                                          str(variation.grow_parents) + '\n')
                    parameters_file.write('        merge_cutoff: ' +
                                          str(variation.merge_cutoff) + '\n')

                elif variation.name == 'structure mutation':
                    parameters_file.write('    StructureMut: \n')
                    parameters_file.write('        fraction: ' +
                                          str(variation.fraction) + '\n')
                    parameters_file.write('        frac_atoms_perturbed: ' +
                                          str(variation.frac_atoms_perturbed) +
                                          '\n')
                    parameters_file.write(
                        '        sigma_atomic_coord_perturbation: ' +
                        str(variation.sigma_atomic_coord_perturbation) + '\n')
                    parameters_file.write(
                        '        max_atomic_coord_perturbation: ' +
                        str(variation.max_atomic_coord_perturbation) + '\n')
                    parameters_file.write(
                        '        sigma_strain_matrix_element: ' +
                        str(variation.sigma_strain_matrix_element) + '\n')

                elif variation.name == 'number of stoichiometries mutation':
                    parameters_file.write('    NumStoichsMut: \n')
                    parameters_file.write('        fraction: ' +
                                          str(variation.fraction) + '\n')
                    parameters_file.write('        mu_num_adds: ' +
                                          str(variation.mu_num_adds) + '\n')
                    parameters_file.write('        sigma_num_adds: ' +
                                          str(variation.sigma_num_adds) + '\n')
                    parameters_file.write('        scale_volume: ' +
                                          str(variation.scale_volume) + '\n')

                elif variation.name == 'permutation':
                    parameters_file.write('    Permutation: \n')
                    parameters_file.write('        fraction: ' +
                                          str(variation.fraction) + '\n')
                    parameters_file.write('        mu_num_swaps: ' +
                                          str(variation.mu_num_swaps) + '\n')
                    parameters_file.write('        sigma_num_swaps: ' +
                                          str(variation.sigma_num_swaps) +
                                          '\n')
                    parameters_file.write('        pairs_to_swap: \n')
                    for pair in variation.pairs_to_swap:
                        parameters_file.write('            - ' + pair + '\n')
        parameters_file.write('\n')

        # write the development info
        parameters_file.write('Development: \n')
        parameters_file.write('    niggli: ' + str(developer.niggli) + '\n')
        parameters_file.write('    scale_density: ' +
                              str(developer.scale_density) + '\n')
        parameters_file.write('\n')

        # write the constraints info
        parameters_file.write('Constraints: \n')
        parameters_file.write('    min_num_atoms: ' +
                              str(constraints.min_num_atoms) + '\n')
        parameters_file.write('    max_num_atoms: ' +
                              str(constraints.max_num_atoms) + '\n')
        parameters_file.write('    min_lattice_length: ' +
                              str(constraints.min_lattice_length) + '\n')
        parameters_file.write('    max_lattice_length: ' +
                              str(constraints.max_lattice_length) + '\n')
        parameters_file.write('    min_lattice_angle: ' +
                              str(constraints.min_lattice_angle) + '\n')
        parameters_file.write('    max_lattice_angle: ' +
                              str(constraints.max_lattice_angle) + '\n')
        parameters_file.write('    per_species_mids: \n')
        for pair in constraints.per_species_mids:
            parameters_file.write('        ' + pair + ': ' +
                                  str(float(
                                      constraints.per_species_mids[pair])) +
                                  '\n')
        parameters_file.write('\n')

        # write the redundancy guard info
        parameters_file.write('RedundancyGuard: \n')
        parameters_file.write('    lattice_length_tol: ' +
                              str(redundancy_guard.lattice_length_tol) + '\n')
        parameters_file.write('    lattice_angle_tol: ' +
                              str(redundancy_guard.lattice_angle_tol) + '\n')
        parameters_file.write('    site_tol: ' +
                              str(redundancy_guard.site_tol) + '\n')
        parameters_file.write('    use_primitive_cell: ' +
                              str(redundancy_guard.use_primitive_cell) + '\n')
        parameters_file.write('    attempt_supercell: ' +
                              str(redundancy_guard.attempt_supercell) + '\n')
        parameters_file.write('    epa_diff: ' +
                              str(redundancy_guard.epa_diff) + '\n')
        parameters_file.write('\n')

        # write the geometry info
        parameters_file.write('Geometry: \n')
        parameters_file.write('    shape: ' + geometry.shape + '\n')
        parameters_file.write('    max_size: ' + str(geometry.max_size) + '\n')
        parameters_file.write('    min_size: ' + str(geometry.min_size) + '\n')
        parameters_file.write('    padding: ' + str(geometry.padding) + '\n')
        parameters_file.write('\n')

        # write the stopping criteria
        parameters_file.write('StoppingCriteria: \n')
        if stopping_criteria.num_energy_calcs is not None:
            parameters_file.write('    num_energy_calcs: ' +
                                  str(stopping_criteria.num_energy_calcs) +
                                  '\n')
        if stopping_criteria.value_achieved is not None:
            parameters_file.write('    value_achieved: ' +
                                  str(stopping_criteria.value_achieved) + '\n')
        if stopping_criteria.found_structure is not None:
            parameters_file.write('    found_structure: ' +
                                  stopping_criteria.path_to_structure_file +
                                  '\n')
        parameters_file.write('\n')
