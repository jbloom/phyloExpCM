"""Module for finding the location of the package data.

This module contains a single function called *Path*. It is
used to find the path where the *phyloExpCM* package data
files are installed assuming that the package is installed
with *distutils* using:

    python setup.py build
    python setup.py install

or::

    python setup.py install --user

"""

import os

_here = os.path.dirname(os.path.abspath(__file__))


def Path():
    """Returns absolute path where data files are installed.

    This function returns the absolute path of its own location
    (minus the file name), followed by the suffix ``data/``.
    This is the location where the package data should be installed.
    """
    path = "%s/data/" % _here
    if not os.path.isdir(path):
        raise IOError("Cannot find the expected path of %s -- something is wrong with the way the package or package data was installed.\nDid you use:\npython setup.py install" % path)
    return path

