# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Plotter module:

This module contains the Plotter class, which is used to plot various data
from the genetic algorithm structure search.

"""

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter

import matplotlib.pyplot as plt
import os


class Plotter(object):
    """
    Used to to plot various data from a structure search.
    """

    def __init__(self, data_file_path):
        """
        Makes a Plotter.

        Args:
            data_file_path: the path to file (called run_data) containing the
                data for the search
        """

        # get the input file contents
        input_file = os.path.abspath(data_file_path)
        try:
            with open(input_file) as input_data:
                self.lines = input_data.readlines()
        except:
            print('Error reading data file.')
            print('Quitting...')
            quit()

    def get_progress_plot(self):
        """
        Returns a plot of the best value versus the number of energy
        calculations, as a matplotlib plot object.
        """

        # set the font to Times, rendered with Latex
        plt.rc('font', **{'family': 'serif', 'serif': ['Times']})
        plt.rc('text', usetex=True)

        # parse the number of composition space endpoints
        endpoints_line = self.lines[0].split()
        endpoints = []
        for word in endpoints_line[::-1]:
            if word == 'endpoints:':
                break
            else:
                endpoints.append(word)
        num_endpoints = len(endpoints)

        if num_endpoints == 1:
            y_label = r'Best value (eV/atom)'
        elif num_endpoints == 2:
            y_label = r'Area of convex hull'
        else:
            y_label = r'Volume of convex hull'

        # parse the best values and numbers of energy calculations
        best_values = []
        num_calcs = []
        for i in range(4, len(self.lines)):
            line = self.lines[i].split()
            num_calcs.append(int(line[4]))
            best_values.append(line[5])

        # check for None best values
        none_indices = []
        for value in best_values:
            if value == 'None':
                none_indices.append(best_values.index(value))

        for index in none_indices:
            del best_values[index]
            del num_calcs[index]

        # make the plot
        plt.plot(num_calcs, best_values, color='blue', linewidth=2)
        plt.xlabel(r'Number of energy calculations', fontsize=22)
        plt.ylabel(y_label, fontsize=22)
        plt.tick_params(which='both', width=1, labelsize=18)
        plt.tick_params(which='major', length=8)
        plt.tick_params(which='minor', length=4)
        plt.xlim(xmin=0)
        return plt

    def plot_progress(self):
        """
        Plots the best value versus the number of energy calculations.
        """

        self.get_progress_plot().show()

    def get_system_size_plot(self):
        """
        Returns a plot of the system size versus the number of energy
        calculations, as a matplotlib plot object.
        """

        # set the font to Times, rendered with Latex
        plt.rc('font', **{'family': 'serif', 'serif': ['Times']})
        plt.rc('text', usetex=True)

        # parse the compositions and numbers of energy calculations
        compositions = []
        num_calcs = []
        for i in range(4, len(self.lines)):
            line = self.lines[i].split()
            compositions.append(line[1])
            num_calcs.append(int(line[4]))

        # get the numbers of atoms from the compositions
        nums_atoms = []
        for composition in compositions:
            comp = Composition(composition)
            nums_atoms.append(comp.num_atoms)

        # make the plot
        plt.plot(num_calcs, nums_atoms, 'D', markersize=5,
                 markeredgecolor='blue', markerfacecolor='blue')
        plt.xlabel(r'Number of energy calculations', fontsize=22)
        plt.ylabel(r'Number of atoms in the cell', fontsize=22)
        plt.tick_params(which='both', width=1, labelsize=18)
        plt.tick_params(which='major', length=8)
        plt.tick_params(which='minor', length=4)
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
        return plt

    def plot_system_size(self):
        """
        Plots the system size versus the number of energy calculations.
        """

        self.get_system_size_plot().show()

    def get_phase_diagram_plot(self):
        """
        Returns a phase diagram plot, as a matplotlib plot object.
        """

        # set the font to Times, rendered with Latex
        plt.rc('font', **{'family': 'serif', 'serif': ['Times']})
        plt.rc('text', usetex=True)

        # parse the composition space endpoints
        endpoints_line = self.lines[0].split()
        endpoints = []
        for word in endpoints_line[::-1]:
            if word == 'endpoints:':
                break
            else:
                endpoints.append(Composition(word))

        if len(endpoints) < 2:
            print('There must be at least 2 endpoint compositions to make a '
                  'phase diagram.')
            quit()

        # parse the compositions and total energies of all the structures
        compositions = []
        total_energies = []
        for i in range(4, len(self.lines)):
            line = self.lines[i].split()
            compositions.append(Composition(line[1]))
            total_energies.append(float(line[2]))

        # make a list of PDEntries
        pdentries = []
        for i in range(len(compositions)):
            pdentries.append(PDEntry(compositions[i], total_energies[i]))

        # make a CompoundPhaseDiagram
        compound_pd = CompoundPhaseDiagram(pdentries, endpoints)

        # make a PhaseDiagramPlotter
        pd_plotter = PDPlotter(compound_pd, show_unstable=100)
        return pd_plotter.get_plot(label_unstable=False)

    def plot_phase_diagram(self):
        """
        Plots the phase diagram.
        """

        self.get_phase_diagram_plot().show()
