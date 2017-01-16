# coding: utf-8
# Copywrite (c) Henniggroup.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

"""
Utility script for replacing each tab in a file with four spaces.

Usage: python replace_tabs.py /path/to/file
"""

import sys

with open(sys.argv[1], 'r') as in_file:
    in_file_lines = in_file.readlines()

handled_lines = []
for line in in_file_lines:
    handled_lines.append(line.replace('\t', '    '))

with open(sys.argv[1], 'w') as in_file:
    for line in handled_lines:
        in_file.write(line)
