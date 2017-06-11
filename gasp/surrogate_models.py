# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function


"""
Surrogate Models module:

This module contains classes that correspond to surrogate models, which
estimate the energy per atom and/or relaxed structure of an organism before
its energy is actually computed with an external energy code.

All classes in this module must implement the following three methods:

add_relaxed_organism()
predict_epa()
predict_relaxed_cell()

1. ShreyasModel: predicts energies using Shreyas' machine learning model

"""

from gasp.general import Organism, Cell
# TODO: other imports here


class ShreyasModel(object):
    """
    Shreyas' machine learning surrogate model.
    """

    def __init__(self, model_params):
        """
        Makes a ShreyasModel, and sets default parameter values if necessary.
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
