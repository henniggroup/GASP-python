# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Surrogate Models module:

This module contains classes that correspond to surrogate models, which
estimate the energy per atom and/or relaxed structure of an organism before
its energy is actually computed with an external energy code.

With the exception of OutputInterpreter, all classes in this module must
implement the following three methods:

add_relaxed_organism()
predict_epa()
predict_relaxed_cell()

1. ShreyasModel: predicts energies using Shreyas' machine learning model

2. OutputInterpreter: interprets the output of the surrogate model

"""

from pymatgen.phasediagram.maker import CompoundPhaseDiagram
from pymatgen.phasediagram.entries import PDEntry
from pymatgen.phasediagram.analyzer import PDAnalyzer

# TODO: other imports here


class ShreyasModel(object):
    """
    Shreyas' machine learning surrogate model.
    """

    def __init__(self, model_params, temp_dir_path):
        """
        Makes a ShreyasModel, and sets default parameter values if necessary.

        Args:
            model_params: the parameters for the surrogate model, as a
                dictionary

            temp_dir_path: the path to temp directory containing the all the
                job directories where the energy calculations are done
        """

        # TODO: implement me
        pass

    def add_relaxed_organism(self, relaxed_organism):
        """
        Adds a new relaxed organism to the model's data set, and optionally
        re-trains the model with the additional data point.

        Args:
            relaxed_organism: an Organism that has been relaxed and whose
                energy has been computed.
        """

        # TODO: implement me
        pass

    def predict_epa(self, unrelaxed_organism):
        """
        Uses the model to predict the energy per atom an unrelaxed Organism
        would have after being relaxed.

        Returns the predicted energy per atom (in eV/atom), or None if no
        energy could be predicted.

        Args:
            unrelaxed_organism: the unrelaxed Organism whose energy to predict
        """

        # TODO: implement me
        pass

    def predict_relaxed_cell(self, unrelaxed_organism):
        """
        Uses the model to predict the Cell an unrelaxed Organism would have
        after being relaxed.

        Returns the predicted relaxed structure, as a Cell object, or None if
        no structure could be predicted.

        Args:
            unrelaxed_organism: the unrelaxed Organism whose relaxed Cell to
                predict
        """

        return None


class OutputInterpreter(object):
    """
    To interpret the output of surrogate models.
    """

    def __init__(self, interpreter_params):
        """
        Makes an OutputInterpreter, and sets default parameter values if
        necessary.

        Args:
            interpreter_params: the parameters for the output interpreter, as a
                dictionary
        """

        # default acceptance threshold, in eV/atom
        self.default_threshold = 0.3

        # the acceptance threshold
        if 'threshold' not in interpreter_params:
            self.threshold = self.default_threshold
        elif interpreter_params['threshold'] in (None, 'default'):
            self.threshold = self.default_threshold
        else:
            self.threshold = interpreter_params['threshold']

    def is_below_threshold(self, organism, predicted_epa, composition_space,
                           current_population):
        """
        Returns a boolean indicating whether an organism's predicted epa means
        that it lies below the threshold.

        Args:
            organism: the unrelaxed Organism to check

            predicted_epa: the energy per atom (eV/atom) that a surrogate model
                has predicted for the unrelaxed organism

            composition_space: the CompositionSpace of the search

            current_population: list of Organisms comprising the current
                population
        """

        # first check if the predicted epa is None
        if predicted_epa is None:
            return True

        # for fixed-composition searches
        if composition_space.objective_function == 'epa':
            best_epa = current_population.sort(key=lambda x: x.epa,
                                               reverse=False)[0].epa
            return (predicted_epa - best_epa) < self.threshold

        # for phase diagram searches
        if composition_space.objective_function == 'pd':
            # check if a phase diagram can be made with the current population
            if self.has_endpoints(composition_space, current_population) and \
                    self.has_non_endpoint(composition_space,
                                          current_population):

                # make a phase diagram
                pdentries = []
                for organism in current_population:
                    pdentries.append(PDEntry(
                        organism.composition, organism.total_energy))
                    compound_pd = CompoundPhaseDiagram(
                        pdentries, composition_space.endpoints)

                # create a pd analyzer object from the compound phase diagram
                pd_analyzer = PDAnalyzer(compound_pd)

                # make a pd entry for the candidate organism
                predicted_total_energy = organism.cell.num_sites*predicted_epa
                pdentry = PDEntry(organism.composition, predicted_total_energy)
                transformed_pdentry = compound_pd.transform_entries(
                    [pdentry], composition_space.endpoints)[0][0]

                # get the predicted distance above the convex hull
                predicted_e_above_hull = pd_analyzer.get_e_above_hull(
                    transformed_pdentry)
                return predicted_e_above_hull < self.threshold

            # if no phase diagram can be made, just return True
            else:
                return True

    def has_endpoints(self, composition_space, current_population):
        """
        Checks if there are organisms in the current population at all of
        endpoints of the composition space.

        Returns a boolean indicating whether or not there are organisms at all
        the endpoints.

        Args:
            composition_space: the CompositionSpace of the search

            current_population: the list of Organisms comprising the current
                population
        """

        for endpoint in composition_space.endpoints:
            has_endpoint = False
            for organism in current_population:
                if endpoint.almost_equals(
                        organism.composition.reduced_composition):
                    has_endpoint = True
            if not has_endpoint:
                return False
        return True

    def has_non_endpoint(self, composition_space, current_population):
        """
        Checks that the current population contains at least one organism not
        at one of the endpoint compositions.

        Returns a boolean indicating whether or not there is an organism not
        at one of the endpoints.

        Args:
            composition_space: the CompositionSpace of the search

            current_population: the list of Organisms comprising the current
                population
        """

        for organism in current_population:
            not_endpoint = True
            for endpoint in composition_space.endpoints:
                if endpoint.almost_equals(
                        organism.composition.reduced_composition):
                    not_endpoint = False
            if not_endpoint:
                return True
        return False
