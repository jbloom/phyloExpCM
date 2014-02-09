.. _phyloExpCM_runcodonPhyML.py:

======================================
phyloExpCM_runcodonPhyML.py 
======================================

This script allows you to estimate a tree topology using `codonPhyML`_. It has been tested with `codonPhyML`_ version 1.00 201306.18.

To run the script, create an input file of the format described below. Then run the script followed by the input file, as in::

    phyloExpCM_runcodonPhyML.py infile.txt

Format of the input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). The input file should
contain the following keys:

* *seqfile* should specify the name of a FASTA file with the sequences that we want to analyze. Each must have a unique header of less than 100 characters. These characters should not include any spaces or any of the following: ``:(),"`` -- an exception will be raised if these conditions are not met. In addition, the character ``*`` has special meaning if it is found at the end of a header (it specifies the sequence is an outgroup), so only include this if you know what you are doing. The sequences are assumed to be aligned coding sequences -- an exception will be raised if they are not all of the same length.

* *codonPhyML_path* should specify the path to the `codonPhyML`_ program. For example, this might just be ``codonphyml`` if the program is installed in the search path. Or it could be an absolute path name if the executable is not installed in the search path.

* *seed* specifies the integer random number seed passed to `codonPhyML`_. Setting the same seed for identical runs should guarantee identical output.

* *model* specifies the codon substitution model to use with `codonPhyML`_. Valid values are:

    - *GY94_CF3x4* : The `Goldman and Yang 1994`_ codon model with the equilibrium codon frequencies estimated empirically using the `CF3x4`_ method, with kappa (the transition / transversion ratio) estimated by maximum likelihood, and with omega (the dN/dS ratio) estimated as a single parameter by maximum likelihood.

    - *GY94_CF3x4_omega-gamma4* : The `Goldman and Yang 1994`_ codon model with the equilibrium codon frequencies estimated empirically using the `CF3x4`_ method, with kappa (the transition / transversion ratio) estimated by maximum likelihood, and with omega (the dN/dS ratio) drawn from four gamma-distributed categories with the distribution shape parameter and mean estimated by maximum likelihood.

    - *KOSI07_F_omega-gamma4* : The `Kosiol et al, 2007`_ empirical codon model with the equilibrium codon frequencies estimated empirically using the *F* method, with *kappa(tv)* (the relative decrease in transversions versus transitions) estimated by maximum likelihood, and with *omega* (the elevation in nonsynonymous over synonymous) drawn from four gamma-distributed categories with the distribution shape parameter and mean estimated by maximum likelihood. This is based on the *ECM+F+omega+1kappa(tv)* model described by `Kosiol et al, 2007`_.

* *outprefix* specifies the prefix for the names of files created by this script. Any existing files with these names are overwritten. The names and contents of theese files are specified in the output files section below.


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_runcodonPhyML.py
    seqfile Human_NPs.fasta
    codonPhyML_path codonphyml
    seed 1
    outprefix Human_NPs_codonphyml
    model GY94_CF3x4_omega-gamma4


Output files
----------------

The output is as follows:

* Basic progress of the script is written to standard output.

* Files are created with the prefix *outprefix* and the following suffixes:

    - ``_config.drw`` is the darwin format input file used to run ``codonphyml``.

    - ``_output.txt`` is the output of ``codonphyml``.

    - ``_tree.newick`` is the inferred tree in Newick format. This tree contains branch supports calculated using the `SH-aLRT`_ method.

    - ``_stats.txt`` contains the statistics from the inference, including all parameter values.

    - ``_loglikelihood.txt`` contains a single number representing log likelihood.



.. include:: weblinks.txt
