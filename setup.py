#!/usr/bin/env python

import os

from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='GASP',
        version='0.1',
        description='',
        long_description=open(os.path.join(module_dir, 'README.rst')).read(),
        url='https://github.com/henniggroup/GASP-python',
        author='Benjamin Revard',
        author_email='bcr48@cornell.edu',
        license='MIT',
        packages=find_packages(),
        package_data={},
        zip_safe=False,
        install_requires=['pymatgen>=4.5.2'],
        classifiers=['Programming Language :: Python :: 2.7',
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.5",
                     'Development Status :: 4 - Beta',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        test_suite='nose.collector',
        tests_require=['nose']
    )
