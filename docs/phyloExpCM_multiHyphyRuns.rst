.. _phyloExpCM_multiHyphyRuns.py:

==========================================
phyloExpCM_multiHyphyRuns.py 
==========================================

This script implements multiple runs of :ref:`phyloExpCM_optimizeHyphyTree.py` that are all being run on the same tree and sequences but are using different substitution models. Essentially, you can use this to automate multiple runs of the :ref:`phyloExpCM_optimizeHyphyTree.py` to compare different substitution models. For much of the details about the options for this script, see the documentation for :ref:`phyloExpCM_optimizeHyphyTree.py`.

Each run is performed in its own subdirectory, which is created by this script. All of the `HYPHY`_ results will be in this subdirectory for each run, with the exception of the *summaryfile* which is created in the main directory used to run the script.

To run the script, create an input file of the format described below. Then run the script followed by the input file, as in::

    phyloExpCM_runcodonPhyML.py infile.txt

Note that running this script may take quite a while if `HYPHY`_ takes a long time to run. `HYPHY`_ can also start to consume a lot of memory for relatively big input trees. 


Format of the input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 
The entries should appear in the order listed below and shown in the example input file.

The input file should contain the following lines:

* The first seven lines should specify the keys listed below. These keys have exactly the same meaning as for the :ref:`phyloExpCM_optimizeHyphyTree.py` script. These keys are:

    - *hyphypath*

    - *hyphycmdfile* : should be relative, not absolute, file name

    - *hyphyoutfile* : should be relative, not absolute, file name

    - *hyphytreefile* : should be relative, not absolute, file name

    - *hyphydistancesfile* : should be relative, not absolute, file name. Make this *None* if you do not want to create these files. Note that no file will be created for substitution models with branch local parameters.

    - *fastafile*

    - *treefile*

    - *siteslist*

* The key *summaryfile* specifies the name of the output file (CSV format) that summarizes the results from all of the `HYPHY`_ runs. The format is described below under Output files.

* The key *nmulti* specifies the integer number of processes that are run at a single time. Each process is a separate run of :ref:`phyloExpCM_optimizeHyphyTree.py`. If you just want to run the analyses one at a time, set this value to 1. Otherwise if you set it to more than one, that many processes will be run simultaneously. This might make sense to do if you are using a multi-processor machine. However, when choosing the right value, do note that some `HYPHY`_ executables use multiple processes (such as ``HYPHYMP CPU=2``), and also that `HYPHY`_ runs can consume large amounts of RAM which is sometimes limiting.

* The remaining lines should list the different substitution models that are used. Each of these lines should begin with an entry that gives the directory (no spaces) that will be created for the run. It should then be followed by text that provides a valid value for the *model* parameter for input to :ref:`phyloExpCM_optimizeHyphyTree.py`. This script will go through each of these listed substitution models, create the directory, and then run :ref:`phyloExpCM_optimizeHyphyTree.py` within that directory. For example, the entry::

    GY94 GY94_CF3x4_omega-global-one_rates-one

  will lead to the creation of the subdirectory ``./GY94/`` and then the running of :ref:`phyloExpCM_optimizeHyphyTree.py` within that subdirectory using *GY94_CF3x4_omega-global-one_rates-one* as the *model* input parameter. Each of the directory names must be unique. If any of the directories already exist, an exception is raised. If the *experimental* model option is used, the specified ``.ibf`` file giving the experimental `HYPHY`_ substitution model should be in the home directory used to run :ref:`phyloExpCM_multiHyphyRuns.py`. The input file created for each of the runs within its own subdirectory will have the name ``phyloExpCM_optimizeHyphyTree_infile.txt``.


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_multiHyphyRuns.py
    hyphypath HYPHYMP CPU=2
    hyphycmdfile hyphy_cmds.bf
    hyphyoutfile hyphy_output.txt
    hyphytreefile hyphy_tree.newick
    hyphydistancesfile hyphy_distances.txt
    fastafile Human_NPs.fasta
    treefile Human_NPs_codonphyml_tree.newick
    siteslist equilibriumpreferences.txt
    summaryfile hyphy_summary.csv
    nmulti 6
    GY94 GY94_CF3x4_omega-global-one_rates-one
    GY94_gamma-rates GY94_CF3x4_omega-global-one_rates-gamma5
    GY94_gamma-rates_gamma-omega GY94_CF3x4_omega-global-gamma5_rates-gamma5
    GY94_gamma-rates_branchlocal-omega GY94_CF3x4_omega-branchlocal-one_rates-gamma5
    KOSI07 KOSI07_F_omega-global-one_rates-one
    KOSI07_gamma-rates KOSI07_F_omega-global-one_rates-gamma5
    KOSI07_gamma-rates_gamma-omega KOSI07_F_omega-global-gamma5_rates-gamma5
    KOSI07_gamma-rates_branchlocal-omega KOSI07_F_omega-branchlocal-one_rates-gamma5
    experimental_FracTolerated experimental hyphyExpCMs_FracTolerated.ibf
    experimental_HalpernBruno experimental hyphyExpCMs_HalpernBruno.ibf
    experimental_FracTolerated_random1 experimental_randomize1 hyphyExpCMs_FracTolerated.ibf
    experimental_FracTolerated_random2 experimental_randomize2 hyphyExpCMs_FracTolerated.ibf
    experimental_HalpernBruno_random1 experimental_randomize1 hyphyExpCMs_HalpernBruno.ibf
    experimental_HalpernBruno_random2 experimental_randomize2 hyphyExpCMs_HalpernBruno.ibf


Output files
----------------

Some summary output is printed to standard output. In addition, the following output is created:

* For each of the specified models, the specified subdirectory is created with the input file ``phyloExpCM_optimizeHyphyTree_infile.txt``, the log file ``phyloExpCM_optimizeHyphyTree_log.txt``, and the errors file ``phyloExpCM_optimizeHyphyTree_errors.txt`` (this last file is empty if nothing is output to standard error). In addition, the subdirectory contains the output files created by running :ref:`phyloExpCM_optimizeHyphyTree.py`.

* The file *summaryfile* lists the log likelihood, the number of branch lengths, and the number of independently optimized `HYPHY`_ parameters **NOT** including the branch lengths for each model. 

  Note, however, that this parameter count does **NOT** include parameters not explicitly optimized by `HYPHY`_ but that are still estimated from data in an empirical way. For example, for the *GY94* method with *CF3x4*, there are 9 such parameters (three nucleotide frequencies at each of three sites). For the *KOSI07* method with *F*, there are 60 parameters (the frequencies of 60 of the 61 non-stop codons).

  Here is an example of the output in *summaryfile*::

    # model, log_likelihood, nbranchlengths, nparameters, nsharedparameters
    GY94, -4441.19, 45, 2, 2
    KOSI07_gamma-rates, -4389.25, 45, 3, 3
    KOSI07, -4410.96, 45, 2, 2
    GY94_gamma-rates, -4429.59, 45, 3, 3


.. include:: weblinks.txt
