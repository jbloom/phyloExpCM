.. _phyloExpCM_optimizeHyphyTree.py:

==========================================
phyloExpCM_optimizeHyphyTree.py
==========================================

This script takes a tree of known topology and uses `HYPHY`_ to optimize the branch lengths and any substitution model parameters. It can use several common codon models as well as the experimentally based codon models with known mutation rates. However, for experimentally based evolutionary models it will often be more convenient to instead use the newer script :ref:`phyloExpCM_ExpModelOptimizeHyphyTree.py` -- and use of this latter script is necessary if the mutation rates are unknown.

This script has been tested with `HYPHY`_ versions 2.112*. This script uses several `HYPHY`_ include batch files (``*.ibf`` files) that are distributed with the `phyloExpCM`_ package in the ``./src/data/`` subdirectory, and should be installed in a location that the script can find if you have installed this package as instructed. If `HYPHY`_ fails with messages about being unable to find included files, maybe this didn't work.

To run the script, create an `Input file`_  as described below. Then run the script from the command line using the input file as the single argument, as in::

    phyloExpCM_optimizeHyphyTree.py infile.txt

Note that running this script may take quite a while if `HYPHY`_ takes a long time to run. `HYPHY`_ can also start to consume a lot of memory for relatively big input trees -- so you may need to use a computer with a large amount of RAM.

Note that **if you are running multiple versions of this script then you should place them in separate subdirectories**. This is important because otherwise the scripts may create temporary files with the same names that interfere with / overwrite each other.


Input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 

The input file should contain the following keys:

* *hyphypath* should contain a command that can be used to run `HYPHY`_ in batch mode. If you have installed `HYPHY`_ in the current search path, this might be something like ``HYPHYMP`` or ``HYPHYMP CPU=4`` (the latter command specifying how many CPUs to use for this multi-thread version) or ``HYPHYDEBUG``. If `HYPHY`_ is not installed in the current search path, you will need to specify a full path to the executable.

* *hyphycmdfile* is the name of the created `HYPHY`_ batch mode command file that is used to run `HYPHY`_ by this script.

* *hyphyoutfile* is the name of the created output file that gives the results of the `HYPHY`_ analysis.

* *hyphytreefile* is the name of the created output file containing the `HYPHY`_ optimized tree.

* *hyphydistancesfile* is the name of a created output file file that shows the pairwise distances between all pairs of sequences (those in *fastafile*). Note that these are **not** the distances from the created `HYPHY`_ tree in *hyphytreefile*. Rather, the substitution model has all parameters frozen to the maximum-likelihood values that are determined by `HYPHY`_ during the optimization of *hyphytreefile*. With these optimized substitution model parameters fixed, `HYPHY`_ is then used to estimate the maximum likelihood pairwise distance between each pair of sequences (essentially, the branch length between them in a two-sequence tree). These pairwise distances are then written *hyphydistancesfile* in the format described in the `Output files`_ section below. If you do not want to estimate pairwise distances, then set *hyphydistancesfile* to *None*. 

  Also, note that **pairwise distances can ony be calculated for substitution models in which all parameters are global** to all branches (the reason is that the model parameters cannot all be fixed if they are branch local). So if you are using a substitution model with local branch lengths (for instance *GY94_CF3x4_omega-branchlocal-one_rates-gamma6*), then no pairwise distances can be calculated. In this case, no *hyphydistancesfile* will be created regardless of whether you specify *None* or no file name.

  It is also possible to use *hyphydistancesfile* to calculate pairwise distances to only certain sequences in *fastafile*. This can be particularly advantageous when using one of the *experimental* codon models (specified by the *models* parameter). The reason is that `HYPHY`_ can take a long time to calculate pairwise distances for *experimental* codon models. If you only want to calculate pairwise distances to certain sequences, then include the name of these sequences following they *hyphydistancesfile* key and the file to be created, as in::

    hyphydistancesfile hyphy_distances.txt A/Brevig_Mission/1918_(H1N1) A/Brisbane/10/2007_(H3N2)

  The aforegoing line would cause the script to create the file *hyphy_distances.txt*, but that file would only include pairwise distances to the two indicated sequences from all other sequences in the file (so if there are *N* sequences there would be *2N - 1* distances). Note that the names of the sequences that are specified must exactly match full headers in *fastafile* or an exception is raised.

  Note also that neither the name of the file specified by *hyphydistancesfile* or any of the specified sequences can contain spaces.


* *fastafile* is the name of a FASTA file containing the sequences that comprise the tip nodes. These should be nucleotide coding sequences aligned at the codon level (not at the nucleotide level). Stop codons should be removed, and the aligned sequence length should be a multiple of three. Ambiguous nucleotide identities are not supported. The headers for the sequences should match the names of the tip nodes in *treefile*.

* *treefile* is a Newick-format phylogenetic tree giving the relationship among the sequences in *fastafile*. The names of the tip nodes must match those in *fastafile*. Branch lengths can be specified, but they will be optimized by `HYPHY`_. Branch supports can be present, but they are ignored and removed in the output tree.

* *siteslist* is the name of a file that specifies all of the codon sites that are included in the analysis in 1, 2, ... numbering. The first line of this file should be a header beginning with a *#* character, and is skipped. The remaining lines should begin with an integer >= 1 that gives a site number, followed by whitespace -- anything following this whitespace is ignored. All of the numbers in this first column are taken to be the numbers of the codon sites included. You can use the ``*_equilibriumpreferences.txt`` file generated by the `mapmuts`_ script ``mapmuts_inferpreferences.py`` and used as input to :ref:`phyloExpCM_buildHyphyExpCM.py` as the value for this argument, since it satisfies these formatting requirements.

* *model* specifies the codon substitution model used by `HYPHY`_. Valid values are:

    - *experimental hyphyExpCMs.ibf* specifies that we use an experimentally determined codon model. The first part of this value (the string *experimental*) indicates that we are using an experimentally determined model. The second part of this value should be the name of a file (no spaces) that specifies the codon model in `HYPHY`_ batch format for each site in *sites*. Typically, this would have been generated using the :ref:`phyloExpCM_buildHyphyExpCM.py` script. This file either must exist in the current directory or you must provide a full path name, as the file will be included by `HYPHY`_ in its analysis.

    - *experimental_randomize1 hyphyExpCMs.ibf* specifies that we use a experimentally determined codon model, but that we randomize the assignment of experimentally determined site models to different sites. This is done by randomly shuffling the position indices in the alignment. As a result, each site model will be assigned to some random site -- not necessarily the one that is supposed to match. You would use this is a control. If the experimental models are good, we might expect them to perform well without randomization but poorly with randomization. For this option, the second argument again specifies the codon model in `HYPHY`_ batch file format for each site in *sites* as generated by :ref:`phyloExpCM_buildHyphyExpCM.py`. The number after *randomize* in the first argument is the integer seed used to seed the random number generator -- different seeds will give different results due to different randomization orders.

    - A number of variants of the `Goldman and Yang 1994`_ (GY94) codon model as delineated below. 
    
      The equilibrium codon frequencies are determined empirically from the frequencies of nucleotides at the three nucleotide positions using either the `CF3x4`_ or the F3x4 model. In general, the `CF3x4`_ model should be preferred unless you have a reason to use F3x4, since the former model corrects numerical inconsistencies with the latter. 

      The transition/transversion ratio *kappa* has one global value applied to all sites and branches which is estimated by maximum likelihood.
      
      The non-synonymous / synonymous substitution rate ratio *omega* can be set to either be *global* so that it is shared among all branches, or can be set to be *branchlocal* so that each branch has its own value of *omega*. Taking the branch local approach will obviously greatly increase the number of parameters. In addition, for both *global* and *branchlocal*, it can be set to *one* for one ratio, or *gamma* with an integer indicating the number of rate categories (such as *gamma6*) drawn from a discrete gamma distribution as originally described in `Yang 1994`_. 
      
      In addition, there can be rate variation across sites. Currently, it is only supported to allow rate variation where the rates are drawn with equal probability from a discrete gamma distribution as originally described in `Yang 1994`_. The gamma distribution alpha parameter is estimated by maximum likelihood. This is specified by the word *rates* followed by *one* if there is only one rate (no rate variation). Otherwise *rates* is followed by *gamma* with an integer indicating the number of rate categories, such as *gamma6* (probably a good value).

      Based on these descriptions, here are some valid model specification strings (other combination also are allowed following this format):
      
        - *GY94_CF3x4_omega-global-one_rates-one* is `CF3x4`_ model for equilibrium codon frequencies, one global *omega* (dN/dS) ratio, and one rate for all sites.

        - *GY94_CF3x4_omega-global-one_rates-gamma6* is `CF3x4`_ model of equilibrium codon frequencies, one global *omega* (dN/dS ratio), and rate variation across sites modeled by six categories drawn from a gamma distribution with the shape parameter estimated by maximum likelihood.

        - *GY94_CF3x4_omega-global-gamma6_rates-one* is `CF3x4`_ model of equilibrium codon frequencies, an *omega* (dN/dS) ratio that is drawn from six gamma-distributed categories with both the mean and the shape estimated by maximum likelihood, and one rate for all sites.

        - *GY94_CF3x4_omega-global-gamma6_rates-gamma6* is `CF3x4`_ model of equilibrium codon frequencies, an *omega* (dN/dS) ratio that is drawn from six gamma-distributed categories with both the mean and the shape estimated by maximum likelihood, and rate variation across sites drawn from six gamma-distributed categories with the shape prameter estimated by maximum likelihood.

        - *GY94_CF3x4_omega-branchlocal-one_rates-gamma6* is like *GY94_CF3x4_omega-global-one_rates-gamma6* except now each branch gets its own single dN/dS ratio *omega* that is estimated by maximum likelihood.

        - *GY94_CF3x4_M1a_rates-one* is the *M1a* model described by `Yang et al 2005`_.

        - *GY94_CF3x4_M2a_rates-one* is the *M2a* model described by `Yang et al 2005`_.

        - *GY94_CF3x4_M7_rates-one* is the *M7* model described by `Yang et al 2005`_ with 10 equal probability categories in the beta distribution.

        - *GY94_CF3x4_M8_rates-one* is the *M8* model described by `Yang et al 2005`_ with 10 equal probability categories in the beta distribution plus one additional discrete category.

      Note that it is NOT allowed to combined *branchlocal* with *gamma* for the *omega* as this introduces too many parameters. So the following would NOT be valid: *GY94_CF3x4_omega-branchlocal-gamma6_rates-one*.

    - A number of variants of the KOSI07 empirical codon model described in `Kosiol et al, 2007`_. That reference delineates a large list of model variants. The ones implemented here are all variants of what is described as the *ECM+F+omega+1kappa(tv)* in the original reference. Specifically, the codon exchangeabilities are those given by `Kosiol et al, 2007`_, and allow multi as well as single-nucleotide codon changes. The equilibrium codon frequencies are empirically estimated as the frequencies in the data (the *+F* option). There is a parameter *omega* scaling non-synonymous substitution rates relative to synonymous rates. Codon changes with *ntv* transversions are scaled to a rate of *kappa^ntv* relative to mutations without transversions. Similarly to the *GY94* models, *omega* can be global to all branches with one estimated value, global to all branches with gamma distributed values, or local to branches with one estimated value per branch. The substitution rate can have a single value (no parameters in this case since it just involves branch length scaling) or have gamma distributed values with an estimated shape parameter. *kappa* is always global and estimated from the data. These options are specified with suffixes similar to the *GY94* model, such as:

        - *KOSI07_F_omega-global-gamma6_rates-gamma6*

        - *KOSI07_F_omega-branchlocal-one_rates-one*

* *persitelikelihoods* is an optional argument. If it is *None*, *False*, or not specified, the nothing is done. Otherwise, set it to the name of a file that you would like to create that reports the per-site likelihoods. Specifically, the created file lists the per-site likelihood for each site after the mutation rates and branch lengths have been optimized by maximum likelihood. These parameters and branch lengths are **not** re-optimized for each site -- rather, they are set to their global best values and then the likelihood contribution of each site is computed. The per-site likelihoods are listed in tab-delimited columns with each site on a separate line, and the first column giving the site number and second giving the per-site likelihood. The first line in the file is a header.


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_optimizeHyphyTree.py
    hyphypath HYPHYMP CPU=4
    hyphycmdfile hyphy_cmds.bf
    hyphyoutfile hyphy_output.txt
    hyphytreefile hyphy_tree.newick
    hyphydistancesfile hyphy_distances.txt
    fastafile Human_NPs.fasta
    treefile Human_NPs_codonphyml_tree.newick
    siteslist equilibriumpreferences.txt
    model experimental hyphyExpCMs_FracTolerated.ibf


Output files
----------------

Some summary output is printed to standard output. In addition, the following output files are created. If any of these files already exist, they are overwritten:

* *hyphycmdfile* is a `HYPHY`_ batch file that contains the commands used to run `HYPHY`_ in batch mode.

* *hyphyoutfile* is an output file that contains the results of the `HYPHY`_ analysis. This file may contain many lines dealing with branch lengths and tree topology that you will probably want to ignore in favor of just looking at *hyphytreefile*, but you will want to look at the first few lines of this file, which will have the following format::

    Log likelihood: -671.317
    independent parameters (includes branch lengths): 45
    shared parameters: 0
    number of branch lengths: 45
    number of tip nodes: 24
    number of internal branches: 21


* *hyphytreefile* is an output file that contains the tree with the `HYPHY`_ optimized branch lengths in Newick format. The branch lengths are *t* parameters from `HYPHY`_. For experimentally determined codon models, these are equal to the expected substitutions averaged over all sites. For GY94 models, these are proportional (not equal) to the expected number of substitutions at the site if *omega* is *global* but not if it is *branchlocal*. Any branch supports present in the original tree are removed.

* *hyphydistancesfile* is an output file that shows the pairwise distances between all sequences in *fastafile*. If there are *N* such sequences, then there are :math:`\frac{N (N-1)}{2}` pairwise distances, and so *hyphydistancesfile* contains this many lines. However, if specific headers are specified following *hyphdistancesfile* such that we only estimate distances for a few sequences, then there will be fewer lines -- for instance, if two sequences are specified, then we have *2N - 1* distances. Each line has three entries:
    
    1) The name of the first sequence (the full header in *fastafile*)

    2) The name of the second sequence (the full header in *fastafile*)

    3) The pairwise distance between the sequences.

  These entries are separated by a tab.

  For example, here are a few possible lines::

    Sequence1   Sequence2    0.800143
    Sequence1   Sequence3    0.92947
    Sequence2   Sequence3    2.13376

* If you use the *persitelikelihoods* option, then an output file is created that shows the per-site likelihoods. The format is as follows::

    #SITE   SITE_LOG_LIKELIHOODS
    1   -4.792906363107919
    2   -15.89808216031635
    3   -17.19022748421175

* The script also creates several temporary files that re-map the sequences in *fastafile* and *treefile* to names that are usable by `HYPHY`_ (which only accepts alphanumeric characters and underscores). Specifically, these files begin with the prefix *_codenames_* and then are followed by suffixes which are the names of *fastafile*, *treefile*, and *hyphydistancesfile*. If the script terminates normally, these files are then deleted -- but they will be present when the script is running, and will be preserved if it halts prematurely to aid in debugging.


.. include:: weblinks.txt
