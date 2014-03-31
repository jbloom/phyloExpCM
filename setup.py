"""Setup script for ``phyloExpCM``

This script uses distutils, the standard python mechanism for installing
packages. To build and install the package, use the following
commands::

    python setup.py build
    python setup.py install

If the user does not have permissions to write to the install directory,
the last command may need to be replaced by::

    sudo python setup.py install

to install globally or by::

    sudo python setup.py --user

to install locally.

Written by Jesse Bloom.
"""


import sys
import os
from distutils.core import setup


# main setup command
setup(
    name = 'phyloExpCM', 
    version = '0.2', 
    author = 'Jesse D. Bloom', 
    author_email = 'jbloom@fhcrc.org', 
    url = 'https://github.com/jbloom/phyloExpCM', 
    description = 'Experimentally determined codon models for phylogenetics',
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Mac OS X.',
    packages = ['phyloExpCM'],
    package_dir = {'phyloExpCM':'src'},
    package_data = {'phyloExpCM':['data/*.ibf']},
    scripts = [
            'scripts/phyloExpCM_runcodonPhyML.py',
            'scripts/phyloExpCM_buildHyphyExpCM.py',
            'scripts/phyloExpCM_optimizeHyphyTree.py',
            'scripts/phyloExpCM_ExpModelOptimizeHyphyTree.py',
            'scripts/phyloExpCM_multiHyphyRuns.py',
            'scripts/phyloExpCM_FreqsFromAlignment.py',
            ],
)
