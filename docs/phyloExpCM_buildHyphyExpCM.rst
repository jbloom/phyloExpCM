.. _phyloExpCM_buildHyphyExpCM.py:

========================================
phyloExpCM_buildHyphyExpCM.py 
========================================

This script build experimentally determined codon models for use with `HYPHY`_. The resulting codon models have been tested with `HYPHY`_ versions 2.112*.

This script is designed to take experimentally determined amino-acid preferences as inferred by `mapmuts`_ and then combine those with an experimentally determined nucleotide mutation spectrum to create the codon substitution models. These models are in an include batch file (``.ibf``) in the format to be included in a `HYPHY`_ analysis.

To run the script, create an input file of the format described below. Then run the script followed by the input file, as in::

    phyloExpCM_runcodonPhyML.py infile.txt


Format of the input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 

The input file should contain the following keys:

* *aapreferences* is the name of a file giving the experimentally determined equilibrium amino-acid preferences. These might be the :math:`\pi_{r,a}` values inferred by the `mapmuts_` script ``mapmuts_inferpreferences.py`` into the ``*_equilibriumpreferences.txt`` file. Specifically, :math:`\pi_{r,a}` represents the preference of site *r* for amino acid *a* under the constraint :math:`\sum_a \pi_{r,a} = 1`. This preference is the expected equilibrium amino-acid frequency of amino acid *a* at site *r* for a hypothetical evolving population where every amino-acid is equally likely to mutate to every other amino acid (so all frequencies would be 1/20 = 0.05 in the absence of selection). 

  The format of the file specified by *aapreferences* should match that of the ``*_equilibriumpreferences.txt`` file created by ``mapmuts_inferpreferences.py``. After an initial header line, each line gives the site number (1, 2, ... numbering) for all specified sites, the wildtype amino acid(s) (not used by this script), the site entropy (not used by this script), and then the experimentally determined equilibrium amino-acid frequencies for each amino acid. These equilibrium frequencies can either include or exclude a stop codon (denoted by the * character) as a possible amino acid. However, for the purposes of this script, stop codons are not considered possible amino acid identities. So if there is value for a stop codon, it is set to zero and the remaining :math:`\pi_{r,a}` values are rescaled so that they sum to one (sum over the 20 amino acids *a*. Here are a few example lines::

    #SITE   WT_AA   SITE_ENTROPY    PI_A    PI_C    PI_D    PI_E    PI_F    PI_G    PI_H    PI_I    PI_K    PI_L    PI_M    PI_N    PI_P    PI_Q    PI_R    PI_S    PI_T    PI_V    PI_W    PI_Y    PI_*
    1   M   2.85229 0.0113803   0.0402777   0.0153121   0.0320438   0.013312    0.00936795  0.026916    0.0192017   0.0224047   0.018926    0.554093    0.0445008   0.0138125   0.0212192   0.0131771   0.0152268   0.0237621   0.0102606   0.0369008   0.0367991   0.0211056
    2   A   0.968359    0.878703    0.00384431  0.00390501  0.00815666  0.0048199   0.00244426  0.00333017  0.00258064  0.00391181  0.00094265  0.0248087   0.00238593  0.00064321  0.00288323  0.00610788  0.0012339   0.001455    0.0290568   0.0107545   0.00492073  0.00311186
    3   S   2.93851 0.0295947   0.00304848  0.00227603  0.00668501  0.131168    0.000710502 0.0199653   0.00804841  0.000619715 0.366698    0.0132841   0.0223199   0.0048327   0.0170484   0.00875982  0.090712    0.0171888   0.0102176   0.177257    0.069395    0.000170873
    4   Q   3.93186 0.012565    0.0355955   0.0396076   0.0344002   0.139887    0.0167737   0.0194143   0.0730001   0.0272629   0.131636    0.0424519   0.0840445   0.00063312  0.00533824  0.0554966   0.0647137   0.0153147   0.019934    0.113926    0.0383248   0.0296804

* *mutspectrum* is the name of a file giving the mutational spectrum. These are the mutational rates, and do not include information about selection (that comes from *aapreferences*). Specifically, the entry for *A -> G* represents the probability that a given site mutations from nucleotide *A* to nucleotide *G* in a unit of time given that it was already *A* (in other words, it is :math:`\Pr\left(G \rm{\; at\; time\;} t = 1 | A \rm{\; at\; time\;} t = 0\right)`. Because the mutations are determined by sequencing one nucleic acid strand and it is not known on which strand a mutation originated, mutations and their reverse complements are grouped together as equivalent (i.e. *A -> G* on one strand is the same as *T -> C* on another). 

  There are two options for how you specify *mutspectrum*:


    1) The format below, which gives the rates of all six types of mutations. Note that these rates may frequently not lead to a reversible substitution model, so you may want to consider using the *makereversible* option. The format here is::

        AG, TC, 2.4e-5
        GA, CT, 2.3e-5
        AT, TA, 3.0e-6
        AC, TG, 9.0e-6
        GC, CG, 1.9e-6
        GT, CA, 9.4e-6

       In this file, the first line specifies that the rate of *A -> G* and *T -> C* (which are assumed equal due to the reverse-complement symmetry) is equal to 2.4e-5. There should be six such lines covering all possible mutations. The units for the rate numbers are arbitrary, but if they are small numbers then you may want to adjust *scalefactor* to get them fairly close to one.

    2) The format below, which just gives the rates of five types of mutations::
           
            AC 9.0e-6
            AG 2.4e-5
            AT 3.0e-6
            CA 9.4e-6
            CG 1.9e-6

       In this model, the sixth mutation rate *CT* is calculated from :math:`R_{C \rightarrow T} = \frac{R_{A \rightarrow G} \times R_{C \rightarrow A}}{R_{A \rightarrow C}}` (see :ref:`phyloExpCM_ExpModelOptimizeHyphyTree.py`). Calculating the mutation rate in this way makes the model reversible, so you should be safe setting *makereversible* to *False*. Note that you still may want to make *scalefactor* large if the absolute values of these rates are low.

* *scalefactor* is a number that multiplies all of the substitution matrix entries determined from *mutspectrum* and *aapreferences* to give the values actually placed in the `HYPHY`_ substitution matrices. In a world of perfect numerical accuracy, the value of this number will not affect the **relative** branch lengths (although making it larger will decrease the absolute magnitude of all branch lengths by a constant factor). However, in practice it may be helpful (this is assumed to be the case but has not been carefully verified) to make this number large enough that when it multiplies the typical value in *mutspectrum*, the result is at least fairly close to one. This will prevent the substitution matrices from having lots of very small numbers, which could cause numerical underflow. So for the example values of *mutspectrum* described immediately above, you might want *scalefactor* equal to something like 10000.0.

* *makereversible* specifies that we constrain the mutational spectrum specified by *mutspectrum* to be reversible. Two **important notes**:

    - Note that subsequent optimization with `HYPHY`_ may not function well / properly if this option is *False*, as the *False* option has not been thoroughly tested. The concern is that `HYPHY`_ will not perform well with non-reversible models.

    - Note that this option is really misnamed, and should probably be called *makesymmetric* rather than *makereversible*. This is because it makes the mutation rates symmetric. This is sufficient to make the model reversible -- but there are also less stringent ways to do this that do not require symmetry. See :ref:`phyloExpCM_ExpModelOptimizeHyphyTree.py`, which is a newer and recommended script for using `HYPHY`_ to optimize with experimentally determined models.
    
  If *makereversible* is *False*, then nothing is done. If *True*, then we set the actual mutation values to be the average of the two symmetric exchanges specified. So for instance, let *mutspectrum* be this file::

    AG, TC, 2.4e-5
    GA, CT, 2.3e-5
    AT, TA, 3.0e-6
    AC, TG, 9.0e-6
    GC, CG, 1.9e-6
    GT, CA, 9.4e-6

  This mutation spectrum is not actually quite symmetric, as for example *A -> G* and *G -> A* have slightly different values. With this option, we constrain it to be reversible, giving the following values:

    - *A -> T, T -> A* = 3.0e-6 (no change from specified values)

    - *G -> C, C -> G* = 1.9e-6 (no change from specified values)

    - *A -> G, T -> C, G -> A, C -> T* = (2.4e-5 + 2.3e-5) / 2 = 2.35e-5

    - *A -> C, T -> G, C -> A, G -> T* = (9.0e-6 + 9.4e-6) / 2 = 9.2e-6

  You can see that in this example, the assumption of symmetry seems to be quite good, as the two averaged values are very close. Note that the symmetry enforced by this option is actually stricter than the requirement needed to make things reversible.

  

* *model* specifies how we convert the experimentally determined amino-acid preferences in *aapreferences* into substitution probabilities between amino acids. There are currently two possible values:

    - *FracTolerated* specifies that we compute the substitution probabilities by interpreting the amino-acid preferences in *aapreferences* as the probability (or fraction of genetic backgrounds in which) a substitution is tolerated. Specifically, if *aapreferences* specifies that two amino acids *x* and *y* at site *r* have equilibrium preferences :math:`\pi_{r,x}` and :math:`\pi_{r,y}`, then the substitution probabilities from *x* to *y* is 1 if :math:`\pi_{r,y} > \pi_{r,x}` and is :math:`\pi_{r,y} / \pi_{r,x}` otherwise.

    - *HalpernBruno* specifies that we compute the substitution probabilities provided by `Halpern and Bruno, MBE, 1998`_ (Equation 10 of that paper, although note that their equation contains a typo). Specifically, if *aapreferences* specifies that two amino acids *x* and *y* at site *r* have equilibrium preferences :math:`\pi_{r,x}` and :math:`\pi_{r,y}`, then the substitution probability from *x* to *y* is 1 if :math:`\pi_{r,x} = \pi_{r,y}` and otherwise is :math:`\left[\ln\left(\pi_{r,y} / \pi_{r,x}\right)\right] / \left(1 - \pi_{r,x} / \pi_{r,y}\right)`. This relationship was derived in `Halpern and Bruno, MBE, 1998`_ based on the fixation probabilities of mutations if :math:`\pi_{r,x}` and :math:`\pi_{r,y}` represent the relative fitness values.

* *evolequilfreqs* specifies the name of an output file created by this script which contains the expected equilibrium frequencies of each amino-acid (sum of the frequencies of all of its codons) under the substitution model built by this script. The format of this created file is exactly equivalent to that for *aapreferences*, however the actual frequencies will be different due to the structure of the genetic code, mutational spectrum, and unequal number of codons per amino acid. This file gives our expectation for how often we would observe each amino acid after evolutionary equilibration. Note that this file is in the appropriate format to serve as input for the ``mapmuts_siteprofileplots.py`` script of the `mapmuts`_ package. The wildtype amino acid column is shown as *na* for all sites.

* *hyphyExpCMs* is the name of the created file that contains the `HYPHY`_ substitution models. Its format is described in the `Output files`_ section below.


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_buildHyphyExpCM.py
    aapreferences equilibriumpreferences.txt
    mutspectrum mutspectrum.txt
    scalefactor 10000.0
    makereversible True
    model FracTolerated
    evolequilfreqs FracTolerated_evolutionary_equilibriumfreqs.txt
    hyphyExpCMs hyphyExpCMs_FracTolerated.ibf


Output files
----------------

Some summary output is printed to standard output. In addition, the following output files are created. If any of these files already exist, they are overwritten:

* *evolequilfreqs* is the output file specified by *evolequilfreqs*. It is in the format described in the explanation of that parameter above.

* *hyphyExpCMs* is a created file in `HYPHY`_ include batch format (``.ibf``) for use in `HYPHY`_. This file is designed to be included in `HYPHY`_ batch files using a line such as::

    #include "hyphyExpCMs_FracTolerated.ibf";

  This file creates a substitution model for each site listed in *aapreferences*. For each site *isite* in *aapreferences* (where numbering goes 1, 2, ...) this file defines the following variables for `HYPHY`_ (names are listed below for site *isite = 2* and *isite = 3*, but there would be similar names for all site numbers):

    - *modelmatrix2*, *modelmatrix3*, etc are 61 X 61 matrices giving the exchangeabilities between each non-stop codon, using the standard `HYPHY`_ codon indexing, which lists all of the 61 non-stop codons in alphabetical order. If *makereversible* is *True*, this matrix is symmetric. *modelmatrix2[i, j]* corresponds to the exchangeability from codon *i* to *j* at site *isite = 2* where codons are indexed 0, 1, ...

    - *equilfreqs2*, *equilfreqs3*, etc are 61 X 1 vectors giving the equilibrium frequencies of each of the 61 non-stop codons, again using `HYPHY`_ codon indexing (alphabetical order).

    - *model2*, *model3*, etc are the substitution models, constructed from *modelmatrix2*, *modelmatrix3*, *equilfreqs2*, *equilfreqs3*, etc using the `HYPHY`_ commands::

        Model model2 = (modelmatrix2, equilfreqs2, 1);
        Model model3 = (modelmatrix3, equilfreqs3, 1);

      The third parameter of 1 on each line indicates that the rate matrices are generated as the product of the symmetric exchangeability matrices and the equilibrium frequencies, and so are reversible.. The defined substitution models *model2*, *model3*, etc can then be used in `HYPHY`_ to construct a *LikelihoodFunction* after including the batch file *hyphyExpCMs*. 

  You should apply the models (*model2*, *model3*, etc) to trees. The model entries all multiply a parameter named *t*, which is assumed to be the variable name that you use in `HYPHY`_ to represent branch lengths.



.. include:: weblinks.txt
