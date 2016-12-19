# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Utility script for making a plot of the system size.

Usage: python plot_system_size.py /path/to/run_data/file
"""

from gasp.post_processing.plotter import Plotter

import sys

plotter = Plotter(sys.argv[1])
plotter.plot_system_size()
