Installation 
==============================



Dependencies
--------------
* `Python`_ is needed. `phyloExpCM`_ has been tested with `Python`_ 2.6 and 2.7 on Linux and Mac OS X, but will probably work with earlier versions as well.

*  `phyloExpCM`_ uses `mapmuts`_. It has been tested with version 1.0.

* `phyloExpCM`_ interfaces with `codonPhyML`_. It has been tested with version 1.00 201306.18 on Mac OS X and Linux.

* `phyloExpCM`_ interfaces with `HYPHY`_. It has been tested with versions 2.112* for Mac OS X and for Linux.

* Some of the scripts and functions in `phyloExpCM`_ use `numpy`_. It has been tested with `numpy`_ version 1.6.1.

* Some of the scripts and functions in `phyloExpCM`_ use `matplotlib`_. It has been tested with `matplotlib`_ version 1.3.1.

Installation
---------------

To install `phyloExpCM`_, first download the source ZIP repository `on GitHub`_. After unzipping the file, run the following commands::

    cd phyloExpCM
    python setup.py build
    python setup.py install

The last command might require you to use::

    sudo python setup.py install
    
if you do not have privileges to the default installation directory, or::
    
    python setup.py install --user
    
if you want to install locally.
    
These commands install the `Python`_ modules and also install several scripts, which provide the most convenient high-level interface into the package. 


`HYPHY`_ include files
------------------------

Several `HYPHY`_ include batch files (``.ibf`` files) are distributed with this package as they are necessary to run some of the `HYPHY`_ analyses. These are in the ``./src/data/`` subdirectory in the package source. If you install the package using the commands above, these data files should be installed in a location that the scripts are able to find. If you run the package without installing it fully, however, the scripts may not be able to find these files. In that case, you will get errors from `HYPHY`_ indicating that it cannot find specified include files.

.. include:: weblinks.txt
