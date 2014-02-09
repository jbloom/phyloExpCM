.. phyloExpCM documentation master file, created by
   sphinx-quickstart on Fri Aug  9 20:30:31 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=======================================
`phyloExpCM`_ documentation
=======================================

`phyloExpCM`_ is a `Python`_ package for implementing experimentally determined codon models of substitution for phylogenetics.

Written by `Jesse Bloom`_.

Source code is available `on GitHub`_.

`phyloExpCM documentation`_ is available `on GitHub pages`_.


Overview
----------

This package contains scripts for performing codon-level analyses of phylogenetic trees. It is designed for a situation in which you would like to estimate a tree topology (for example with `codonPhyML`_), and then use `HYPHY`_ to compare the extent to which various codon substitution models can describe the evolution of sequences within that tree topology. It is specifically designed to allow you to test experimentally determined codon models working with the output from the `mapmuts`_ package.

This package can be run using a series of scripts which are documented individually. You should be able to run analyses using this package by running those scripts without interfacing directly with the `Python`_ modules.


Contents
------------

.. toctree::
   :maxdepth: 2

   installation
   scripts
   examples
   modules
   acknowledgements



Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`


.. include:: weblinks.txt
