.. _phyloExpCM_plotLogLvsNParams.py:

==========================================
phyloExpCM_plotLogLvsNParams.py
==========================================

This script is designed to help visually display the results of testing several substitution models to infer or optimize phylogenetic trees from the same set of data. The script uses `matplotlib`_ to plot the likelihood versus the number of parameters. You might want to use this script to help display results after several runs of :ref:`phyloExpCM_optimizeHyphyTree.py`.

To run the script, create an `Input file`_  as described below. Then run the script from the command line using the input file as the single argument, as in::

    phyloExpCM_plotLogLvsNParams.py infile.txt

This will create the output plot described in `Output files`_.

Input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 

The input file should contain the following keys:

* *plotfile* should give the name of the created PDF plot. This name should end in the extension ``.pdf``.

* *datafile* should give the name of the existing data file. This is a CSV file which currently does not support entries that themselves contain commas. Essentially, each line is split between commas into entries. Lines in the file that begin with # or are blank are ignored. All other lines must contain at least seven comma-separated entries. Each line should specify a data point. The entries that are important to this script are:

    * Entry 3 should give the log likelihood for that model as a number.

    * Entry 4 should give the number of free parameters for that model as an integer.

    * Entry 7 should be a string that is the "plotting group" to which the data point is assigned. This is the label assigned to the data point in the legend. Each "plotting group" is indicated by a different style point and has its own entry in the legend. So if you make a separate "plotting group" for each data point, then each data point has its own legend entry. On the other hand, if you give several data points the same "plotting group," then they will have the same symbol and legend entry.

  Here is an example of a valid *datafile*::

    #Summary for tree KOSI07.
    #
    #SUBSTITUTION_MODEL, AIC, LOG_LIKELIHOOD, FREE_PARAMETERS, MAXIMUM_LIKELIHOOD_PARAMETERS, EMPIRICAL_PARAMETERS, PLOTTING_GROUP
    experimental, 31000.0, -15000.0, 0, 0, 0, experimental
    KOSI07_F_omega-global-gamma4_rates-one, 32354.2, -16114.1, 63, 3, 60, KOSI07
    KOSI07_F_omega-global-one_rates-gamma4, 32454.8, -16164.4, 63, 3, 60, KOSI07
    KOSI07_F_omega-global-one_rates-one, 33032.5, -16454.3, 62, 2, 60, KOSI07


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_plotLogLvsNParams.py
    plotfile logl_vs_nparams.pdf
    datafile KOSI07_summary.csv


Output files
----------------
The created output file is the PDF plot file specified by *plotfile*. Here is an example:

    **INCLUDE EXAMPLE HERE**

.. include:: weblinks.txt
