# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Test module:

This module contains the unit tests for GASP.

"""

from gasp import general, development, variations, population, \
    energy_calculators, organism_creators, objects_maker

from gasp import geometry as geo

from pymatgen.core.periodic_table import Element

import unittest
import random
import copy
import numpy as np
import os
import shutil
import tempfile


class TestIDGenerator(unittest.TestCase):
    def setUp(self):
        self.id_generator = general.IDGenerator()

    def test_initial_id_value(self):
        self.assertEqual(self.id_generator.id, 0)

    def test_incrementing_value(self):
        first_id = self.id_generator.make_id()
        second_id = self.id_generator.make_id()
        self.assertEqual(first_id + 1, second_id)


class TestMating(unittest.TestCase):
    def setUp(self):
        mating_params = {'fraction': 0.1}
        self.mating = variations.Mating(mating_params)
        # get the current working directory and make a copy of it
        self.starting_cwd = copy.deepcopy(os.getcwd())
        # make a temporary directory and move into it
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        # move back to the original working directory
        os.chdir(self.starting_cwd)
        # delete the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_double_parent(self):
        # make a cell to test
        species = [Element('Cu')]
        coords = [[0.5, 0.5, 0.5]]
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)

        # test with bulk geometry
        geometry = geo.Bulk()
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(2*old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a different lattice
        lattice = [[3, 0, 0], [0, 1, 0], [0, 0, 2]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(2*old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a different lattice
        lattice = [[3, 0, 0], [0, 2, 0], [0, 0, 1]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(2*old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a random lattice
        lattice = [
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)]]
        cell = general.Cell(lattice, species, coords)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test with sheet geometry
        geometry = geo.Sheet({'shape': 'sheet'})
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = cell.lattice.volume
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(2*old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a different lattice
        lattice = [[3, 0, 0], [0, 1, 0], [0, 0, 2]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = cell.lattice.volume
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(2*old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a different lattice
        lattice = [[3, 0, 0], [0, 2, 0], [0, 0, 1]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = cell.lattice.volume
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(2*old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a random lattice
        lattice = [
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)]]
        cell = general.Cell(lattice, species, coords)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test with wire geometry
        geometry = geo.Wire({'shape': 'wire'})
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = cell.lattice.volume
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(2*old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a different lattice
        lattice = [[3, 0, 0], [0, 2, 0], [0, 0, 1]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = cell.lattice.volume
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(2*old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test a cell with a random lattice
        lattice = [
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)]]
        cell = general.Cell(lattice, species, coords)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.volume, 5), round(2*old_volume, 5))
        self.assertEqual(cell.num_sites, 2*old_num_atoms)

        # test with cluster geometry and a random lattice
        geometry = geo.Cluster({'shape': 'cluster'})
        lattice = [
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)],
            [random.uniform(-5, 5), random.uniform(-5, 5),
             random.uniform(-5, 5)]]
        cell = general.Cell(lattice, species, coords)
        old_lattice_lengths = copy.deepcopy(cell.lattice.abc)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        self.mating.double_parent(cell, geometry)
        self.assertEqual(round(cell.lattice.a, 5),
                         round(old_lattice_lengths[0], 5))
        self.assertEqual(round(cell.lattice.b, 5),
                         round(old_lattice_lengths[1], 5))
        self.assertEqual(round(cell.lattice.c, 5),
                         round(old_lattice_lengths[2], 5))
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)

    def test_grow_parent_cell(self):
        # non-cluster geometry
        shape = random.choice(['bulk', 'sheet', 'wire'])
        if shape == 'bulk':
            geometry = geo.Bulk()
        elif shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})

        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        old_num_atoms_1 = copy.deepcopy(cell_1.num_sites)
        old_volume_1 = copy.deepcopy(cell_1.lattice.volume)

        # second parent same size as first
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent twice the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, 2*old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(2*old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent four times the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 6, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, 4*old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(4*old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent eight times the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 6, 0], [0, 0, 6]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, 8*old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(8*old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # cluster geometry
        geometry = geo.Cluster({'shape': 'cluster'})
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        old_num_atoms_1 = copy.deepcopy(cell_1.num_sites)
        old_volume_1 = copy.deepcopy(cell_1.lattice.volume)

        # second parent same size as first
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent twice the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent four times the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 6, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

        # second parent eight times the size of first parent
        species_1 = [Element('Cu')]
        coords_1 = [[0.5, 0.5, 0.5]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Al')]
        coords_2 = [[0.4, 0.4, 0.4]]
        lattice_2 = [[6, 0, 0], [0, 6, 0], [0, 0, 6]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        old_num_atoms_2 = copy.deepcopy(cell_2.num_sites)
        old_volume_2 = copy.deepcopy(cell_2.lattice.volume)
        self.mating.grow_parent_cell(cell_1, cell_2, geometry, random)
        self.assertEqual(cell_1.num_sites, old_num_atoms_1)
        self.assertEqual(cell_2.num_sites, old_num_atoms_2)
        self.assertEqual(round(cell_1.lattice.volume, 5),
                         round(old_volume_1, 5))
        self.assertEqual(round(cell_2.lattice.volume, 5),
                         round(old_volume_2, 5))

    def test_get_num_doubles(self):
        self.assertEqual(0, self.mating.get_num_doubles(-10))
        self.assertEqual(0, self.mating.get_num_doubles(-1))
        self.assertEqual(0, self.mating.get_num_doubles(0))
        self.assertEqual(0, self.mating.get_num_doubles(1))
        self.assertEqual(1, self.mating.get_num_doubles(1.5))
        self.assertEqual(1, self.mating.get_num_doubles(2))
        self.assertEqual(2, self.mating.get_num_doubles(3))
        self.assertEqual(2, self.mating.get_num_doubles(5))
        self.assertEqual(3, self.mating.get_num_doubles(6))
        self.assertEqual(3, self.mating.get_num_doubles(10))
        self.assertEqual(4, self.mating.get_num_doubles(12))
        self.assertEqual(4, self.mating.get_num_doubles(20))
        self.assertEqual(5, self.mating.get_num_doubles(24))
        self.assertEqual(5, self.mating.get_num_doubles(35))
        self.assertEqual(6, self.mating.get_num_doubles(48))
        self.assertEqual(6, self.mating.get_num_doubles(80))
        self.assertEqual(7, self.mating.get_num_doubles(96))
        self.assertEqual(7, self.mating.get_num_doubles(1000))

    def test_make_offspring_cell(self):
        # make two cells to test
        species_1 = [Element('Al'), Element('Al')]
        coords_1 = [[0.3, 0.3, 0.3], [0.7, 0.7, 0.7]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 4]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        species_2 = [Element('Cu'), Element('Cu')]
        coords_2 = [[0.7, 0.7, 0.3], [0.7, 0.3, 0.3]]
        lattice_2 = [[4, 0, 0], [0, 4, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)

        # make the geometry and constraints
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'bulk':
            geometry = geo.Bulk()
        elif shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)

        # make the offspring cell
        offspring_cell = self.mating.make_offspring_cell(
            cell_1, cell_2, geometry, constraints, random)

        # offspring cell should contain at least one atom from each parent
        # note: this test only works if parents don't contain same species
        has_specie_from_parent_1 = False
        has_specie_from_parent_2 = False
        for specie in offspring_cell.species:
            if specie in cell_1.species:
                has_specie_from_parent_1 = True
            elif specie in cell_2.species:
                has_specie_from_parent_2 = True
        self.assertTrue(has_specie_from_parent_1)
        self.assertTrue(has_specie_from_parent_2)

    def test_do_random_shift(self):
        # make a cell to test
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_volume = copy.deepcopy(cell.lattice.volume)
        old_num_atoms = copy.deepcopy(cell.num_sites)

        # test with bulk geometry and random lattice vector index
        geometry = geo.Bulk()
        vector_index = random.choice([0, 1, 2])
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        for site_coords in cell.frac_coords:
            for coord in site_coords:
                self.assertTrue(coord >= 0.0)
                self.assertTrue(coord < 1.0)

        # test with sheet geometry and lattice vector index for shift
        geometry = geo.Sheet({'shape': 'sheet'})
        vector_index = random.choice([0, 1])
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        for site_coords in cell.frac_coords:
            for coord in site_coords:
                self.assertTrue(coord >= 0.0)
                self.assertTrue(coord < 1.0)

        # test with sheet geometry and lattice vector index for no shift
        vector_index = 2
        old_frac_coords = cell.frac_coords
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(np.array_equal(cell.frac_coords, old_frac_coords))

        # test with wire geometry and lattice vector index for shift
        geometry = geo.Wire({'shape': 'wire'})
        vector_index = 2
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        for site_coords in cell.frac_coords:
            for coord in site_coords:
                self.assertTrue(coord >= 0.0)
                self.assertTrue(coord < 1.0)

        # test with wire geometry and lattice vector index for no shift
        vector_index = random.choice([0, 1])
        old_frac_coords = cell.frac_coords
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(np.array_equal(cell.frac_coords, old_frac_coords))

        # test with cluster geometry and random lattice vector index (no shift)
        geometry = geo.Cluster({'shape': 'cluster'})
        vector_index = random.choice([0, 1, 2])
        old_frac_coords = cell.frac_coords
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.mating.do_random_shift(cell, vector_index, geometry, random)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(np.array_equal(cell.frac_coords, old_frac_coords))

    def test_do_random_rotation(self):
        # make a cell to test
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_frac_coords = cell.frac_coords
        composition_space = general.CompositionSpace(['Cu'])
        constraints = development.Constraints(None, composition_space)

        # test with bulk geometry
        geometry = geo.Bulk()
        self.mating.do_random_rotation(cell, geometry, constraints, random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(np.array_equal(cell.frac_coords, old_frac_coords))

        # test with sheet geometry
        geometry = geo.Sheet({'shape': 'sheet'})
        old_frac_coords = cell.frac_coords
        self.mating.do_random_rotation(cell, geometry, constraints, random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(np.array_equal(cell.frac_coords, old_frac_coords))

        # test with wire geometry
        geometry = geo.Wire({'shape': 'wire'})
        self.mating.do_random_rotation(cell, geometry, constraints, random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        for site_coords in cell.frac_coords:
            for coord in site_coords:
                self.assertTrue(coord >= 0.0)
                self.assertTrue(coord < 1.0)

        # test with cluster geometry
        geometry = geo.Cluster({'shape': 'cluster'})
        self.mating.do_random_rotation(cell, geometry, constraints, random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        for site_coords in cell.frac_coords:
            for coord in site_coords:
                self.assertTrue(coord >= 0.0)
                self.assertTrue(coord < 1.0)

    def test_rotate_about_axis(self):
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        random_axis = random.choice([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        random_angle = random.uniform(0, 360)
        self.mating.rotate_about_axis(cell, random_axis, random_angle)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertEqual(round(cell.lattice.volume, 5), round(old_volume, 5))

    def test_merge_sites(self):
        # only one site
        species = [Element('Cu')]
        coords = [[0.5, 0.5, 0.5]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertEqual(round(new_cell.lattice.volume, 5),
                         round(old_volume, 5))

        # two sites close enough to merge but different species
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.49], [0.5, 0.5, 0.51]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(new_cell.num_sites, old_num_atoms)
        self.assertEqual(round(new_cell.lattice.volume, 5),
                         round(old_volume, 5))

        # two sites close enough to merge
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.49], [0.5, 0.5, 0.51]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(new_cell.num_sites, old_num_atoms - 1)

        # three sites close enough to merge
        species = [Element('Cu'), Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.49], [0.5, 0.5, 0.5], [0.5, 0.5, 0.51]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(new_cell.num_sites, old_num_atoms - 2)

        # three sites close enough to merge, but one is differenc species
        species = [Element('Cu'), Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.49], [0.5, 0.5, 0.5], [0.5, 0.5, 0.51]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(new_cell.num_sites, old_num_atoms - 1)

        # two duplicate sites
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        old_volume = copy.deepcopy(cell.lattice.volume)
        new_cell = self.mating.merge_sites(cell, geometry, constraints)
        self.assertEqual(new_cell.num_sites, old_num_atoms - 1)

    def test_do_variation(self):
        # test with fixed-composition search
        composition_space = general.CompositionSpace(['AlCu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.mating.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 1)

        # test with variable composition search
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make the two endpoint organisms
        species_1 = [Element('Al'), Element('Al')]
        coords_1 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        org_1 = general.Organism(cell_1, id_generator, 'maker',
                                 composition_space)
        org_1.epa = random.random() - 2
        org_1.total_energy = org_1.epa*org_1.cell.num_sites
        species_2 = [Element('Cu'), Element('Cu')]
        coords_2 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        org_2 = general.Organism(cell_2, id_generator, 'maker',
                                 composition_space)
        org_2.epa = random.random() - 2
        org_2.total_energy = org_2.epa*org_2.cell.num_sites
        initial_population.add_organism(org_1, composition_space)
        initial_population.add_organism(org_2, composition_space)
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                random_org.total_energy = \
                    random_org.epa*random_org.cell.num_sites
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.mating.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 1)


class TestStructureMut(unittest.TestCase):
    def setUp(self):
        structure_mut_params = {'fraction': 0.1}
        self.structure_mut = variations.StructureMut(structure_mut_params)
        # get the current working directory and make a copy of it
        self.starting_cwd = copy.deepcopy(os.getcwd())
        # make a temporary directory and move into it
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        # move back to the original working directory
        os.chdir(self.starting_cwd)
        # delete the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_perturb_atomic_coords(self):
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        shape = random.choice(['bulk', 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        constraints = development.Constraints(None, composition_space)
        geometry.unpad(cell, constraints)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        self.structure_mut.perturb_atomic_coords(cell, geometry, constraints,
                                                 random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(old_cell.composition.almost_equals(cell.composition))

    def test_perturb_lattice_vectors(self):
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        old_num_atoms = copy.deepcopy(cell.num_sites)
        self.structure_mut.perturb_lattice_vectors(cell, random)
        self.assertEqual(cell.num_sites, old_num_atoms)
        self.assertTrue(old_cell.composition.almost_equals(cell.composition))

    def test_do_variation(self):
        # test with fixed-composition search
        composition_space = general.CompositionSpace(['AlCu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.structure_mut.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 0)
        # this is only true if fixed-composition
        self.assertEqual(offspring.cell.composition.reduced_composition,
                         composition_space.endpoints[0])

        # test with variable composition search
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make the two endpoint organisms
        species_1 = [Element('Al'), Element('Al')]
        coords_1 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        org_1 = general.Organism(cell_1, id_generator, 'maker',
                                 composition_space)
        org_1.epa = random.random() - 2
        org_1.total_energy = org_1.epa*org_1.cell.num_sites
        species_2 = [Element('Cu'), Element('Cu')]
        coords_2 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        org_2 = general.Organism(cell_2, id_generator, 'maker',
                                 composition_space)
        org_2.epa = random.random() - 2
        org_2.total_energy = org_2.epa*org_2.cell.num_sites
        initial_population.add_organism(org_1, composition_space)
        initial_population.add_organism(org_2, composition_space)
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                random_org.total_energy = \
                    random_org.epa*random_org.cell.num_sites
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.structure_mut.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 0)


class TestNumAtomsMut(unittest.TestCase):
    def setUp(self):
        num_atoms_mut_params = {'fraction': 0.1}
        self.num_atoms_mut = variations.NumAtomsMut(num_atoms_mut_params)
        # get the current working directory and make a copy of it
        self.starting_cwd = copy.deepcopy(os.getcwd())
        # make a temporary directory and move into it
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        # move back to the original working directory
        os.chdir(self.starting_cwd)
        # delete the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_compute_num_adds(self):
        # test with fixed-composition and no removes possible
        composition_space = general.CompositionSpace(['AlCu'])
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = self.num_atoms_mut.compute_num_adds(cell, composition_space,
                                                       random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertTrue(num_adds > 0)

        # test with fixed-composition and 1 remove possible
        composition_space = general.CompositionSpace(['AlCu'])
        species = [Element('Cu'), Element('Cu'), Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = self.num_atoms_mut.compute_num_adds(cell, composition_space,
                                                       random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertTrue(num_adds > 0 or num_adds == -1)

        # test with variable composition and no removes possible
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        species = [Element('Cu')]
        coords = [[0.5, 0.5, 0.25]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = self.num_atoms_mut.compute_num_adds(cell, composition_space,
                                                       random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertTrue(num_adds > 0)

        # test with variable composition and no removes possible
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = self.num_atoms_mut.compute_num_adds(cell, composition_space,
                                                       random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertTrue(num_adds > 0 or num_adds == -1)

    def test_add_atoms_epa(self):
        # test with adding 1 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 1
        self.num_atoms_mut.add_atoms_epa(cell, num_adds, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites +
            num_adds*old_cell.composition.reduced_composition.num_atoms)

        # test with adding 2 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 2
        self.num_atoms_mut.add_atoms_epa(cell, num_adds, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites +
            num_adds*old_cell.composition.reduced_composition.num_atoms)

        # test with adding 10 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 10
        self.num_atoms_mut.add_atoms_epa(cell, num_adds, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites +
            num_adds*old_cell.composition.reduced_composition.num_atoms)

    def test_remove_atoms_epa(self):
        # test with removing 1 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al'),
                   Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 1
        self.num_atoms_mut.remove_atoms_epa(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites -
            num_removes*old_cell.composition.reduced_composition.num_atoms)

        # test with removing 2 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al'),
                   Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 2
        self.num_atoms_mut.remove_atoms_epa(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites -
            num_removes*old_cell.composition.reduced_composition.num_atoms)

        # test with removing 5 stoichiometry's worth of atoms
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Cu'),
                   Element('Cu'), Element('Cu')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 5
        self.num_atoms_mut.remove_atoms_epa(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.composition.reduced_composition,
                         old_cell.composition.reduced_composition)
        self.assertEqual(
            cell.num_sites, old_cell.num_sites -
            num_removes*old_cell.composition.reduced_composition.num_atoms)

    def test_add_atoms_pd(self):
        # test with adding 1 atom
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 1
        self.num_atoms_mut.add_atoms_pd(cell, num_adds, composition_space,
                                        random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites + num_adds)

        # test with adding 2 atoms
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 2
        self.num_atoms_mut.add_atoms_pd(cell, num_adds, composition_space,
                                        random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites + num_adds)

        # test with adding 10 atoms
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_adds = 10
        self.num_atoms_mut.add_atoms_pd(cell, num_adds, composition_space,
                                        random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites + num_adds)

    def test_remove_atoms_pd(self):
        # test with removing 1 atom
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al'),
                   Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 1
        self.num_atoms_mut.remove_atoms_pd(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites - num_removes)

        # test with removing 2 atoms
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al'),
                   Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 2
        self.num_atoms_mut.remove_atoms_pd(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites - num_removes)

        # test with removing 5 atoms
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al'),
                   Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.25, 0.25], [0.5, 0.75, 0.25],
                  [0.5, 0.25, 0.75], [0.25, 0.75, 0.75], [0.75, 0.75, 0.75]]
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        num_removes = 5
        self.num_atoms_mut.remove_atoms_pd(cell, num_removes, random)
        self.assertEqual(cell.lattice.volume, old_cell.lattice.volume)
        self.assertEqual(cell.num_sites, old_cell.num_sites - num_removes)

    def test_do_variation(self):
        # test with fixed-composition search
        composition_space = general.CompositionSpace(['AlCu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.num_atoms_mut.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 0)
        # this is only true if fixed-composition
        self.assertEqual(offspring.cell.composition.reduced_composition,
                         composition_space.endpoints[0])

        # test with variable composition search
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make the two endpoint organisms
        species_1 = [Element('Al'), Element('Al')]
        coords_1 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        org_1 = general.Organism(cell_1, id_generator, 'maker',
                                 composition_space)
        org_1.epa = random.random() - 2
        org_1.total_energy = org_1.epa*org_1.cell.num_sites
        species_2 = [Element('Cu'), Element('Cu')]
        coords_2 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        org_2 = general.Organism(cell_2, id_generator, 'maker',
                                 composition_space)
        org_2.epa = random.random() - 2
        org_2.total_energy = org_2.epa*org_2.cell.num_sites
        initial_population.add_organism(org_1, composition_space)
        initial_population.add_organism(org_2, composition_space)
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                random_org.total_energy = \
                    random_org.epa*random_org.cell.num_sites
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.num_atoms_mut.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        self.assertTrue(offspring is not None)
        self.assertTrue(offspring.cell.num_sites > 0)


class TestPermutation(unittest.TestCase):
    def setUp(self):
        permutation_params = {'fraction': 0.1}
        composition_space = general.CompositionSpace(['Cu', 'Al'])
        self.permutation = variations.Permutation(permutation_params,
                                                  composition_space)
        # get the current working directory and make a copy of it
        self.starting_cwd = copy.deepcopy(os.getcwd())
        # make a temporary directory and move into it
        self.temp_dir = tempfile.mkdtemp()
        os.chdir(self.temp_dir)

    def tearDown(self):
        # move back to the original working directory
        os.chdir(self.starting_cwd)
        # delete the temporary directory
        shutil.rmtree(self.temp_dir)

    def test_select_valid_parent(self):
        # need to make a Pool object to pass as an argument to the method
        composition_space = general.CompositionSpace(['AlCu'])
        constraints = development.Constraints(None, composition_space)
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        pool.add_initial_population(initial_population, composition_space)

        # make sure the selected parent has at least 1 possible swap
        valid_parent = self.permutation.select_valid_parent(
            pool, composition_space, random)
        if valid_parent is not None:
            possible_swaps = self.permutation.get_possible_swaps(
                valid_parent.cell)
            self.assertTrue(len(possible_swaps) > 0)

    def test_get_possible_swap(self):
        # test with no possible swaps
        species = [Element('Cu'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        swaps = self.permutation.get_possible_swaps(cell)
        self.assertEqual(len(swaps), 0)

        # test with one possible swap
        species = [Element('Al'), Element('Cu')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        swaps = self.permutation.get_possible_swaps(cell)
        self.assertEqual(len(swaps), 1)

    def test_get_indices_to_swap(self):
        # can get 1 swap, and asked for 1
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        num_swaps = 1
        indices_to_swap = self.permutation.get_indices_to_swap(cell, num_swaps,
                                                               random)
        self.assertTrue(len(indices_to_swap), num_swaps)

        # can get 2 swaps, and asked for 2
        species = [Element('Cu'), Element('Cu'), Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        num_swaps = 2
        indices_to_swap = self.permutation.get_indices_to_swap(cell, num_swaps,
                                                               random)
        self.assertTrue(len(indices_to_swap), num_swaps)

        # can only get 1 swap, but asked for 2
        species = [Element('Cu'), Element('Al')]
        coords = [[0.5, 0.5, 0.25], [0.5, 0.5, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        num_swaps = 2
        indices_to_swap = self.permutation.get_indices_to_swap(cell, num_swaps,
                                                               random)
        self.assertTrue(len(indices_to_swap), 1)

        # can get 1 swaps, and asked for 3
        species = [Element('Cu'), Element('Cu'), Element('Cu'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        num_swaps = 3
        indices_to_swap = self.permutation.get_indices_to_swap(cell, num_swaps,
                                                               random)
        self.assertTrue(len(indices_to_swap), 1)

    def test_swap_pairs(self):
        # swapping one pair
        species = [Element('Cu'), Element('Al'), Element('Cu'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        indices_to_swap = [[0, 1]]
        self.permutation.swap_pairs(cell, indices_to_swap)
        self.assertEqual(cell.composition, old_cell.composition)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(round(cell.lattice.volume),
                         round(old_cell.lattice.volume))
        self.assertEqual(cell.sites[0].specie, old_cell.sites[1].specie)
        self.assertEqual(cell.sites[1].specie, old_cell.sites[0].specie)

        # swapping two pairs
        species = [Element('Cu'), Element('Cu'), Element('Al'), Element('Al')]
        coords = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75],
                  [0.25, 0.75, 0.75]]
        lattice = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell = general.Cell(lattice, species, coords)
        old_cell = copy.deepcopy(cell)
        indices_to_swap = [[0, 1], [2, 3]]
        self.permutation.swap_pairs(cell, indices_to_swap)
        self.assertEqual(cell.composition, old_cell.composition)
        self.assertEqual(cell.num_sites, old_cell.num_sites)
        self.assertEqual(round(cell.lattice.volume),
                         round(old_cell.lattice.volume))
        self.assertEqual(cell.sites[0].specie, old_cell.sites[1].specie)
        self.assertEqual(cell.sites[1].specie, old_cell.sites[0].specie)
        self.assertEqual(cell.sites[2].specie, old_cell.sites[3].specie)
        self.assertEqual(cell.sites[3].specie, old_cell.sites[2].specie)

    def test_do_variation(self):
        # test with fixed-composition search
        composition_space = general.CompositionSpace(['AlCu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.permutation.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        if offspring is not None:
            self.assertTrue(offspring.cell.num_sites > 0)
            # this is only true if fixed-composition
            self.assertEqual(offspring.cell.composition.reduced_composition,
                             composition_space.endpoints[0])

        # test with variable composition search
        composition_space = general.CompositionSpace(['Al', 'Cu'])
        constraints = development.Constraints(None, composition_space)
        shape = random.choice(['bulk,' 'sheet', 'wire', 'cluster'])
        if shape == 'sheet':
            geometry = geo.Sheet({'shape': shape})
        elif shape == 'wire':
            geometry = geo.Wire({'shape': shape})
        elif shape == 'cluster':
            geometry = geo.Cluster({'shape': shape})
        else:
            geometry = geo.Bulk()
        initial_population = population.InitialPopulation('doesnt_matter')
        random_org_creator = organism_creators.RandomOrganismCreator(
            None, composition_space, constraints)
        id_generator = general.IDGenerator()
        # make the two endpoint organisms
        species_1 = [Element('Al'), Element('Al')]
        coords_1 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_1 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_1 = general.Cell(lattice_1, species_1, coords_1)
        org_1 = general.Organism(cell_1, id_generator, 'maker',
                                 composition_space)
        org_1.epa = random.random() - 2
        org_1.total_energy = org_1.epa*org_1.cell.num_sites
        species_2 = [Element('Cu'), Element('Cu')]
        coords_2 = [[0.25, 0.25, 0.25], [0.75, 0.75, 0.75]]
        lattice_2 = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        cell_2 = general.Cell(lattice_2, species_2, coords_2)
        org_2 = general.Organism(cell_2, id_generator, 'maker',
                                 composition_space)
        org_2.epa = random.random() - 2
        org_2.total_energy = org_2.epa*org_2.cell.num_sites
        initial_population.add_organism(org_1, composition_space)
        initial_population.add_organism(org_2, composition_space)
        # make a bunch of random organisms for the initial population
        for _ in range(30):
            random_org = random_org_creator.create_organism(
                id_generator, composition_space, constraints, random)
            if random_org is not None:
                random_org.epa = random.random() - 2  # just make up an epa
                random_org.total_energy = \
                    random_org.epa*random_org.cell.num_sites
                initial_population.add_organism(random_org, composition_space)
        pool = population.Pool(None, composition_space, 'doesnt_matter')
        selection = general.SelectionProbDist(None, 20)
        pool.selection = selection
        composition_fitness_weight = general.CompositionFitnessWeight(None)
        pool.comp_fitness_weight = composition_fitness_weight
        pool.add_initial_population(initial_population, composition_space)
        offspring = self.permutation.do_variation(
            pool, random, geometry, constraints, id_generator,
            composition_space)

        # check the offspring
        if offspring is not None:
            self.assertTrue(offspring.cell.num_sites > 0)

if __name__ == '__main__':
    unittest.main()
