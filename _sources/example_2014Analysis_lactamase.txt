===============================================================
2014 analysis of beta-lactamase phylogenetic fit
===============================================================

.. contents::
   :depth: 3

Overview
------------
This is an analysis of the phylogenetic fit of an evolutionary model determined using deep mutational scanning of beta-lactamase. This is the analysis in the paper `An experimentally informed evolutionary model improves phylogenetic fit to divergent lactamase homologs`_.

The results show that an evolutionary model using experimentally determined amino-acid preferences and mutation rates fit from the data outperforms standard phylogenetic evolutionary models.

The data comes from `Firnberg et al 2014`_.

This analysis was performed by `Jesse Bloom`_.

The overall approach is to use the `phyloExpCM`_ package to analyze a beta-lactamase sequence phylogeny (tree topology estimated by `codonPhyML`_ and then fixed) with an experimentally determined evolutionary model based on data from `Firnberg et al 2014`_, and a variety of standard codon substitution models.

Much of the analysis itself is performed using `HYPHY`_, which is run by `phyloExpCM`_ using the scripts detailed in the `phyloExpCM documentation`_. 

The code, data, and results from the analysis are in the ``./examples/2014Analysis_lactamase/`` subdirectory in the main `phyloExpCM`_ package directory. This subdirectory contains all of the necessary input files (see `Input files`_) for running the analysis. It also contains a master script (``run_analysis.py``) that can be used to regenerate the analysis with the command::

    python run_analysis.py

provided that you have installed the software described in `Software used`_.


Sequence numbering
-------------------
Beta-lactamase sequences are often numbered using the `Ambler`_ numbering scheme. This is the numbering scheme used by `Firnberg et al 2014`_. However, the phylogenetic analysis described here converts all numbering to sequential 1, 2, ... numbering (beginning with the N-terminal methionine) of the wildtype lactamase parent used by `Firnberg et al 2014`_. So this is sequential numbering of the protein encoded by the coding sequence in the file ``TEM1_cds.fasta`` described in `Input files`_. 

The file ``sequential_to_Ambler.csv`` converts between the sequential and `Ambler`_ numbering scheme. This file is created from the data provided by `Firnberg et al 2014`_ as described in `Extracting amino-acid preferences`_.

Software used
----------------
The analysis was performed using the following software on the Fred Hutchinson computer cluster using the following Linux distribution::

    Linux version 3.5.0-43-generic (buildd@roseapple) (gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5) ) #66~precise1-Ubuntu SMP Thu Oct 24 14:52:23 UTC 2013

In order to replicate the analysis, you will need to use Linux or Mac OS X, and install the software listed below.

The specific versions of software used are specified. Other versions will likely work as well, but this has not been confirmed.

* `Python`_ (version 2.7.3)

* `phyloExpCM`_ (specifically, `phyloExpCM version 0.3`_)

* `numpy`_ (version 1.6.1)

* `matplotlib`_ (version 1.3.1)

* `Biopython`_ (version 1.58)

* `mapmuts`_ : (`mapmuts version 1.0`_) was used to make the sequence logo plot and to do some basic sequence manipulation in `Python`_.

* `EMBOSS needle`_ (version 6.6.0). The executable with the name ``needle`` is expected to be in the current search path.

* `RAxML`_ (version 7.7.6). The executable with the name ``raxmlHPC-SSE3`` is expected to be in the current search path. This is the SSE3 version, and was built using the command::

    make -f Makefile.SSE3.gcc

  and then copying the executable into the search path.

* `FigTree`_ (version 1.4.0) was used to open tree files for annotation, visualization, and saving of images.

* `Path-O-Gen`_ (version 1.4) was used to manually root the trees for visual display.

* `codonPhyML`_ (version 1.00 201306.18) was used to build phylogenetic trees. The executable with the name ``codonphyml`` is expected to be in the current search path. This executable was compiled with the *omp* options to make it multi-threaded.

* `HYPHY`_ (version HYPHY 2.220131214beta(MP) for Linux on x86_64) was used to optimize the branches and compute the likelihoods under the different models. The multiprocessor (``HYPHYMP``) version of the executable was used. To install this, the following commands were run in the unzipped `HYPHY`_ source code directory after download::

    cmake -DINSTALL_PREFIX=/home/jbloom/.local/ .
    make MP2
    make install

* This analysis was run on the FHCRC scientific computing cluster using `sbatch`_ (version 2.6.2), and requires `sbatch`_ to be installed to run as currently implemented. If you don't want to use `sbatch`_, all you need to do is modify the main ``run_analysis.py`` script to not submit any jobs via `sbatch`_. The use of `sbatch`_ in the current version simply allows more processors to be used simultaneously and so speeds up the analysis.

* `ImageMagick convert`_ 6.6.9-7 2012-08-17 was used to convert PDF files to JPGs.

* `weblogo`_ (version 3.3) was used to make the sequence logo plot.

* `PyPdf`_ (version 1.13) was used to make overlays on the sequence logo plot.

Input files
-----------------------
Here are the input files that are necessary to repeat the analysis. These are the names assigned with the ``./examples/2014Analysis_lactamase/`` subdirectory of the main `phyloExpCM`_ package.

* ``run_analysis.py`` is a `Python`_ script that will perform the analysis assuming that you have all of the appropriate software (see `Software used`_) and input files. To run this script, use::

    python run_analysis.py

  at the command line. 

* ``get_preferences.py`` is a `Python`_ script that is used to extract the amino-acid preferences from ``Firnberg_missense_mutation_fitnesses.csv``.

* ``get_treeseqs.py`` is a `Python`_ script that is used to construct the beta-lactamase sequence alignments.

* ``TEM1_cds.fasta`` is a FASTA file giving the coding DNA sequence for the parent TEM1 beta-lactamase using by `Firnberg et al 2014`_ as reconstructed from ``Data S2`` of their paper.

* ``Firnberg_missense_mutation_fitnesses.csv`` is the Excel spreadsheet ``Data S2`` from `Firnberg et al 2014`_ saved in CSV (comma separated values) format using Excel. This is the table that gives the estimated fitness of each amino-acid mutations, disregarding the possibility that different synonymous codons might have unequal fitnesses.

* Files containing the beta-lactamase sequences used to build the sequence alignments. The beta-lactamase sequences to use were taken from the `Lahey database`_, which contains a curated listing of beta-lactamases. All *TEM* and *SHV* lactamases for which `Genbank`_ nucleotide sequences were available were downloaded. First, the accession numbers for all sequences for which these were provided were parsed from the Excel tables in the `Lahey database`_ on March-9-2014. Then sequences were batch downloaded from `Genbank`_ for all of these accessions for which such sequences were available. The accessions and the batch downloaded sequences are in the following files:

    - ``Lahey_TEM_GenbankAccessions.txt``

    - ``Lahey_TEM_GenbankSequences.gb``

    - ``Lahey_SHV_GenbankAccessions.txt``

    - ``Lahey_SHV_GenbankSequences.gb``

* ``1XPB.pdb`` is PDB file `1XPB`_ of the crystal structure of TEM-1 beta-lactmase, as downloaded from the Protein Data Bank.

* ``1XPB_renumbered.pdb`` is file ``1XPB.pdb`` renumbered using `PDB Goodies`_ so that the numbering is sequential and matches the sequential numbering in ``TEM1_cds.fasta`` in sequential numbering as described in `Sequence numbering`_.

* ``1XPB_renumbered.dssp`` is the result of analyzing ``1XPB_renumbered.pdb`` using the `DSSP webserver`_. This contains the accessible surface area and secondary structure for the crystal structure.

Running the analysis
-----------------------
The analysis in this directory can be repeated simply by running the master script ``run_analysis.py`` from the command line with::

    python run_analysis.py

(or run it using `sbatch`_).

Other than the files delineated under `Input files`_, all data files in this directory can be re-generated with this command. However, the analysis does involve a few manual steps to create nicely formatted and annotated tree images from the raw tree files produced by the script. These manual steps are described in the `Steps in the analysis`_ section below.

There is a variable called *use_existing_output* that is hardcoded into the beginning of this script. This variable is currently set to *True*, which means that if output already exists it is **not** overwritten and regenerated. If you want to regenerate the output (such as after changing input files or altering analysis parameters), then either manually delete the old output or change *use_existing_output* to *False* to force the script to overwrite old output. This script does **not** automatically test for dependencies among files, so if you change a file but keep *use_existing_output* to *True*, then you need to manually remove any old files with dependencies.

Note that ``run_analysis.py`` uses multiple CPUs to speed up some of the analyses. If you don't have as many CPUs on your system, you might want to check that you are not overloading the system -- you can avoid this by using fewer CPUs, which will require you to manually modify the ``run_analysis.py`` script.

Note that ``run_analysis.py`` uses `sbatch`_ by default to run some of its subsidiary programs -- if you don't want to do this, set the *use_sbatch* option to *False* within the program.


Steps in the analysis
---------------------------------------
Here are the steps in the analysis. Except for a few indicated manual steps, all of this is done automatically by the ``run_analysis.py`` script.


Extracting amino-acid preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The approach used here defines the evolutionary model in terms of the amino-acid preferences, where :math:`\pi_{r,a}` denotes the preference of site *r* for amino acid *a*. `Firnberg et al 2014`_ provide their data in terms of the "fitnesses" of different mutations, which are computed from the deep sequencing counts. Some aspects of the procedure used to compute these "fitnesses" from the deep sequencing counts appear to be rather ad hoc. However, here we do not attempt to re-analyze the raw data of `Firnberg et al 2014`_ (which in any case is only partly provided in their publication), but instead simply take their fitness values at face value and convert them to something approximating the amino-acid preferences. This is done by the `Python`_ script ``get_preferences.py``, which can be run using the command::

    python get_preferences.py

This script processes the data from `Firnberg et al 2014`_ in ``Firnberg_missense_mutation_fitnesses.csv`` to generate the file ``amino_acid_preferences.txt``. 

Specifically, the script works as follows:

1) For each residue *r*, we read in the fitness value :math:`w_{r,a}` for amino acid *a* in ``Firnberg_missense_mutation_fitnesses.csv`` if such a value is present. 

2) For the wildtype amino acid at site *r*, we set :math:`w_{r,a} = 1` rather than using the value (usually close to one) provided in ``Firnberg_missense_mutation_fitnesses.csv``. This is done because their method (see supplemental material of `Firnberg et al 2014`_) defines the fitness of the wildtype to be one, and the values in the file that are usually close to but not exactly one are actually for synonymous mutations.

3) Sometimes amino acids don't have a fitness specified. When this is the case, we take the mean of all fitness values at the site (including the wildtype amino acid values set to one in the previous step), and then assign this mean fitness to any amino acids with unknown fitness. This is probably non-ideal, because it seems likely that the missing amino acid fitnesses correspond to very low fitness identities. But because `Firnberg et al 2014`_ fail to provide raw sequencing data or computer code, it is unclear how to do any better.

4) The preference for each amino acid is defined as proportional to its fitness, with the constraint that :math:`\sum_a \pi_{r,a} = 1`. So specifically, 

   .. math::

      \pi_{r,a} = \frac{w_{r,a}}{\sum_{a'} w_{r,a'}}

   where :math:`a'` ranges over all amino acids.

5) The preferences computed using the aforementioned formula are written to the output file ``amino_acid_preferences.txt``.
   Note that the residues in the created ``amino_acid_preferences.txt`` file are numbered sequentially rather than using the `Ambler`_ numbering (see `Sequence numbering`_).

6) The ``get_preferences.py`` script also creates the file ``sequential_to_Ambler.csv``, which converts between the sequential and `Ambler`_ numbering schemes (see `Sequence numbering`_).

Visualizing the amino-acid preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The amino-acid preferences are visualized using a logo plot created by the `mapmuts`_ script ``mapmuts_siteprofileplots.py``, which is detailed in the `mapmuts documentation`_. Briefly, this script uses `weblogo`_ to visualize the preferences by making a plot in which the heights of amino acids are equal to their preferences at each site. Secondary structure and relative solvent accessibility information calculated using the `DSSP webserver`_ from PDB `1XPB`_ is overlaid on the plot. The residues in this plot are numbered using the `Ambler`_ numbering scheme (see `Sequence numbering`_).

The main ``run_analysis.py`` script runs ``mapmuts_siteprofileplots.py`` after creating its input file ``mapmuts_siteprofileplots_infile.txt``. The resulting image is in ``lactamase_site_preferences_logoplot.pdf``. A JPG version of this plot (``lactamase_site_preferences_logoplot.jpg``) is also created using `ImageMagick convert`_. The JPG is shown below:

.. figure:: example_2014Analysis_lactamase_lactamase_site_preferences_logoplot.jpg
   :width: 80%
   :align: center
   :alt: lactamase_site_preferences_logoplot.jpg

   Visual display of amino-acid preferences. Letter heights are proportional to the preference for that amino acid at that site. Secondary structure (SS) and relative solvent accessibility information is overlaid on the plots. Amino-acid letters are colored according to a hydrophobicity scale.


Creation of sequence alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The creation of the sequence alignment ``aligned_lactamases.fasta`` is performed by the `Python`_ script ``get_treeseqs.py`` using the command::

    python get_treeseqs.py

The steps performed by this script are as follows:

1) The CDS sequences are extracted from the Genbank sequences in ``Lahey_TEM_GenbankSequences.gb`` and ``Lahey_SHV_GenbankSequences.gb`` for all records with a single specified CDS that encodes a translatable protein and has not ambiguous nucleotide identities.

2) The CDS sequences are pairwise aligned with the reference sequence in ``TEM1_cds.fasta`` constructing codon alignments from the protein alignments, and any gaps relative to the reference sequence are stripped away.

3) Any sequences that do have less than 60% sequence identity to ``TEM1_cds.fasta`` among alignable sites or have more than 20% gaps are removed. 

4) Aligned sequences are examined for the number of nucleotide differences with other sequences. When there are multiple sequences that have less than four differences from another sequence, only one representative is retained. This reduces the number of highly similar sequences.

5) The remaining set of unique aligned sequences is written to ``aligned_lactamases.fasta``.

Inspection of aligned sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aligned sequences in ``aligned_lactamases.fasta`` are subjected to the following manual steps to check the quality of the data:

1) The ``aligned_lactamases.fasta`` file was manually checked for recombinant sequences using the online interface to `DataMonkey`_ to run the *SBP* and *GARD* programs. Briefly:

        * Using *SBP*, there was no evidence of recombination using the recommended criterion of *cAIC*.

        * Using *GARD*, there was no evidence of recombination.

2) A quick phylogenetic tree was built using `RAxML`_ to make sure that the tree topology seemed reasonable. The `RAxML`_ output was placed in the subdirectory ``./RAxML_output/``. The following command was executed::

        raxmlHPC-SSE3 -w ./RAxML_output/ -n aligned_lactamases -p 1 -m GTRCAT -s aligned_lactamases.fasta

   Note that in the actual command, the full path to ``./RAxML_output/`` is specified. This created the `RAxML`_ tree as ``/RAxML_output/RAxML_bestTree.aligned_lactamases``, as well as several other files in ``./RAxML_output/``. The tree was manually visualized with `FigTree`_ to confirm that it separated into the two expected clades of TEM and SHV sequences. 
    

Construction of phylogenetic trees 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
High-quality phylogenetic trees were then constructed using `codonPhyML`_, which is able to build maximum-likelihood trees with codon substitution models. 

The `phyloExpCM`_ script ``phyloExpCM_runcodonPhyML.py`` (see `phyloExpCM documentation`_ for details) was used to run `codonPhyML`_. Trees were built using two different substitution models, with each tree in its own subdirectory within ``./examples/2014Analysis_lactamase/``. The substitution models used to build the trees:

    * The *GY94* codon model (`Goldman and Yang 1994`_) with:
    
        - The equilibrium codon frequencies estimated empirically using the `CF3x4`_ method.

        - A single transition / transversion ratio (*kappa*) estimated by maximum likelihood.

        - The dN/dS ratio (*omega*) drawn from four discrete gamma-distributed categories (`Yang 1994`_) with the shape parameter and mean estimated by maximum likelihood.

      The tree construction was performed in the ``./GY94_codonPhyML_tree/`` subdirectory. The tree itself is in the file ``./GY94_codonPhyML_tree/codonphyml_tree.newick``.

    * The *KOSI07* empirical codon model (`Kosiol et al 2007`_), using the variant denoted as *ECM+F+omega+1kappa(tv)* in `Kosiol et al 2007`_. Specifically:

        - The equilibrium codon frequencies were estimated empirically using the *F* method (empirical estimation of 60 frequencies for the 61 non-stop codons).

        - The relative decrease in transversions versus transitions, *kappa(tv)*, estimated by maximum likelihood.

        - The relative elevation in nonsynonymous over synonymous mutations (*omega*) drawn from four gamma-distributed categories with shape parameter and mean estimated by maximum likelihood.

      The tree construction was performed in the ``./KOSI07_codonPhyML_tree/`` subdirectory. The tree itself is in the file ``./KOSI07_codonPhyML_tree/codonphyml_tree.newick``.

To visualize these trees, the following **manual** steps were performed for the tree in each subdirectory:

    1) The ``codonphyml_tree.newick`` file was re-rooted using `Path-O-Gen`_ and then opened in `FigTree`_, adjusted for attractive visual display, and saved to the file ``annotated_tree.trees``.

    2) A PDF image of the tree in the ``annotated_tree.trees`` file was saved manually with `FigTree`_, and then converted to a JPG with `ImageMagick convert`_ using the command::

        convert -density 300 annotated_tree.pdf annotated_tree.jpg

       The tree images created by these steps are shown below:

         .. figure:: example_2014Analysis_lactamase_GY94_codonPhyML_tree_annotated_tree.jpg
            :width: 65%
            :align: center
            :alt: GY94_codonPhyML_tree/annotated_tree.jpg

            Image of the tree generated by `codonPhyML`_ for the *GY94* substitution model (file ``./GY94_codonPhyML_tree/annotated_tree.jpg``). The *TEM* sequences are colored red, while the *SHV* sequences are colored blue.

         .. figure:: example_2014Analysis_lactamase_KOSI07_codonPhyML_tree_annotated_tree.jpg
            :width: 65%
            :align: center
            :alt: KOSI07_codonPhyML_tree/annotated_tree.jpg

            Image of the tree generated by `codonPhyML`_ for the *KOSI07* substitution model (file ``./KOSI07_codonPhyML_tree/annotated_tree.jpg``). The *TEM* sequences are colored red, while the *SHV* sequences are colored blue.

The `Robinson-Foulds distance`_ between the trees computed with the *GY94* and *KOSI07* models was computed with `RAxML`_, and the results written to the directory ``./RobinsonFouldsDistances/``. This distance provides a quantitative measure of the similarity between the tree topologies.

Each of these trees was also parsed into the subtrees representing the *TEM* and the *SHV* sequences keeping the tree topology estimated jointly. The subtree parsing is done automatically by the ``run_analysis.py`` script using `BioPython`_. These subtrees have the names:

    * ``./KOSI07_codonPhyML_tree/codonphyml_tree_SHV.newick`` : the *SHV* subtree from ``./KOSI07_codonPhyML_tree/codonphyml_tree.newick``

    * ``./KOSI07_codonPhyML_tree/codonphyml_tree_TEM.newick`` 

    * ``./GY94_codonPhyML_tree/codonphyml_tree_SHV.newick`` 

    * ``./GY94_codonPhyML_tree/codonphyml_tree_TEM.newick`` 

The sequence alignments for these subtrees are in the files

    * ``aligned_GY94_SHV.fasta``

    * ``aligned_GY94_TEM.fasta`` 

    * ``aligned_KOSI07_SHV.fasta``

    * ``aligned_KOSI07_TEM.fasta``


Average amino-acid frequencies in the alignments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The overall average amino-acid frequencies for each alignment (averaged over all sites and sequences) were computed using the `phyloExpCM`_ script ``phyloExpCM_FreqsFromAlignment.py``.

Within the ``./examples/2014Analysis_lactamase/`` these files have the following names:

    * ``avg_aa_freqs_aligned_lactamases.txt`` : the averages over all of the lactamase sequences in ``aligned_lactamases.fasta``

    * ``avg_aa_freqs_aligned_GY94_TEM.txt`` : the averages over just the sequences in ``aligned_GY94_TEM.fasta``

    * ``avg_aa_freqs_aligned_KOSI07_TEM.txt``

    * ``avg_aa_freqs_aligned_GY94_SHV.txt``

    * ``avg_aa_freqs_aligned_KOSI07_SHV.txt``

Although these files list a line of frequencies for each site, the frequencies reported for each site are the **average** over all sites and sequences. These frequencies are used in the models that have the same mutation structure as the experimentally informed models, but replace the preferences with the average amino-acid frequencies. These frequencies are calculated by adding a pseudocount of one to each amino-acid type (a very small correction as there are many sequences and sites).


Optimizing the trees with various substitution models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the heart of the analysis: to use various codon substitution models to optimize the branch lengths and likelihood for the tree topologies estimated with `codonPhyML`_ and the *GY94* or *KOSI07* substitution models as described in `Construction of phylogenetic trees`_. The analyses are also done for the *TEM* and *SHV* subtrees. This analysis compares the experimentally determined codon substitution models with various other substitution models.

This analysis is done with the `phyloExpCM`_ scripts ``phyloExpCM_optimizeHyphyTree.py`` (for standard substitution models) and ``phyloExpCM_ExpModelOptimizeHyphyTree.py`` (for substitution models using the experimentally determined amino-acid preferences), both of which are described in the `phyloExpCM documentation`_. These scripts uses `HYPHY`_ to optimize the branch lengths and any substitution model free parameters while keeping the tree topologies fixed.

The results of the optimizations for the various substitution models are in the subdirectory ``./codonmodel_optimized_trees/``. Within this subdirectory are further subdirectories containing the results of optimizing the branch lengths with `HYPHY`_ using the indicated codon substitution models. For example, ``./codonmodel_optimized_trees/Tree-KOSI07_Model-KOSI07_F_omega-global-gamma4_rates-one/`` contains one such model. Within this subdirectory, the files of interest are:

    * The tree with branch lengths optimized by `HYPHY`_ using the substitution model in question is in a file called either ``hyphy_tree.newick`` or ``optimizedtree.newick``.

    * A summary of the `HYPHY`_ results including the log likelihood is found in a file called either ``hyphy_output.txt`` or ``optimizedtree_results.txt``.

    * ``*_infile.txt`` the input file used to run either ``phyloExpCM_optimizeHyphyTree.py`` or ``phyloExpCM_ExpModelOptimizeHyphyTree.py``.

The models are:

    * *HalpernBruno* : the model described as *HalpernBruno* in the `phyloExpCM documentation`_ for ``phyloExpCM_OptimizeDetectSelection.py``. Briefly, the fixation probabilities :math:`F_{r,xy}` are determined from the amino-acid preferences :math:`\pi_{r,a}` using formula of `Halpern and Bruno 1998`_ as 

        .. math::

           F_{r,xy} = 
           \begin{cases}
           1 & \mbox{if $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
           \frac{\beta \times \ln\left(\pi_{r,\mathcal{A}\left(y\right)} / \pi_{r,\mathcal{A}\left(x\right)}\right)}{1 - \left(\pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}\right)^{\beta}} & \mbox{otherwise.}
           \end{cases}

      The mutation rates involve the fitting of four free parameters as described in the `phyloExpCM documentation`_ for ``phyloExpCM_OptimizeDetectSelection.py``. The stringency parameter :math:`\beta` is either fixed to one or optimized.

    * *FracTolerated* : the model described as *FracTolerated* in the `phyloExpCM documentation`_ for ``phyloExpCM_OptimizeDetectSelection.py``. This model differs from *HalpernBruno* in that the fixation probabilities are now determined as

        .. math::

           F_{r,xy} =
           \begin{cases}
           1 & \mbox{if $\pi_{r,\mathcal{A}\left(y\right)} \ge \pi_{r,\mathcal{A}\left(x\right)}$} \\
           \left(\frac{\pi_{r,\mathcal{A}\left(y\right)}}{\pi_{\mathcal{A}\left(x\right)}}\right)^{\beta} & \mbox{otherwise.}
           \end{cases}

    * *HalpernBrunorandomized* : like *HalpernBruno* but with the experimentally determined amino-acid preferences randomized among sites. This model is not expected to be good because there is no longer any accurate site-specific information.

    * *FracToleratedrandomized* : like *FracTolerated* but with the experimentally determined amino-acid preferences randomized among sites.

    * *avgaafreqsHalpernBruno* : this model uses the mutation model and fixation model of *HalpernBruno*, but does not use the experimentally determined amino-acid preferences. Instead, the amino-acid preferences for all sites are simply set based on the amino-acid frequencies obtained from averaging all sites and sequences in the alignment (see `Average amino-acid frequencies in the alignments`_). So this is a null model that captures the mutation process of the *HalpernBruno* model but isn't site specific. It does contain 60 parameters (the amino-acid frequencies) that are estimated empirically from the seuences, plus the four free mutation parameters used in *HalperBruno*.

    * *avgaafreqsFracTolerated* : like *avgaafreqsHalpernBruno* except that the fixation model is *FracTolerated*.

    * Versions of the six models above where the stringency parameter :math:`\beta` is fixed to one, and versions where :math:`\beta` is treated as a free parameter. These have names like *HalpernBruno* (when :math:`\beta = 1`) and like *fitbetaHalpernBruno* (when :math:`\beta` is a free parameter).

    * *GY94_CF3x4_omega-global-one_rates-one* : the `Goldman and Yang 1994`_ model with the equilibrium frequencies from the `CF3x4`_ method (9 parameters) and a single transition-transversion ratio (*kappa*) and nonsynonymous-synonymous ration (*omega*) estimated by maximum likelihood (2 parameters).

    * *GY94_CF3x4_omega-global-one_rates-gamma4* : like the *GY94_CF3x4_omega-global-one_rates-one* model but with the substitution rate drawn from a 4-category discrete gamma distribution (`Yang 1994`_) with the shape parameter estimated by maximum likelihood (adds 1 parameters).

    * *GY94_CF3x4_omega-global-gamma4_rates-one* : like the *GY94_CF3x4_omega-global-one_rates-one* model but with *omega* drawn from a 4-category discrete gamma distribution (`Yang 1994`_) with the shape parameter estimated by maximum likelihood (adds 1 parameters).

    * *GY94_CF3x4_omega-global-gamma4_rates-gamma4* : like the *GY94_CF3x4_omega-global-one_rates-one* model but with both the rate and *omega* drawn from their own 4-category discrete gamma distributions (`Yang 1994`_) with the shape parameters estimated by maximum likelihood (adds 2 parameters).

    * *GY94_CF3x4_omega-global-branchlocal_rates-gamma4* : like the *GY94_CF3x4_omega-global-one_rates-gamma4* model but with each branch getting its own *omega* value optimized by maximum likelihood (adds a number of parameters equal to the number of branch lengths).

    * Various *KOSI07* models as described in `Kosiol et al 2007`_ with the equilibrium frequencies estimated using the *F* method (60 parameters) and the other options as described immediately above for the *GY94* model.

A helpful way to view the results is to look at the main summary files created by ``run_analysis.py`` in the examples home directory. These files summarize the results for each tree topology.
The results are sorted by the `AIC`_ from best to worst.
Overall, the files contain the following columns:

    1) Description of the substitution model.    

    2) The difference in `AIC`_ from the best model.

    3) The log likelihood

    4) The total number of free parameters in the substitution model. This is the sum of the parameters optimized by maximum likelihood and the parameters estimated empirically from the data.

    5) The total number of parameters in the substitution model optimized by maximum likelihood:

        - For experimental models, this is four or five. There are always four optimized mutation rates, and some models also fit the stringency parameter :math:`\beta` while others constrain it to one.

        - For *GY94* and *KOSI07* models, this is always at least two, as both *omega* and *kappa* are estimated by maximum likelihood. For some of the variants, *omega* is drawn from discrete gamma categories, which adds an additional parameter (there is now a shape and rate parameter). For some of the variants, the rate is also allowed to be drawn from discrete gamma categories, which also adds another parameter. Finally, for some of the variants *omega* is different for each branch, which adds a number of free parameters equal to the number of branches.

    6) The total number of parameters in the substitution model that are empirically estimated from the sequence data:

        - For experimental models, this is zero.

        - For *GY94* models, this is 9 as in the `CF3x4`_ model there are a total of three independent frequencies for each of the three codon positions.

        - For *KOSI07* models, this is 60 as there are a total of 60 independent frequencies for the 61 non-stop codons.

Here are the results for optimization of the combined *TEM* and *SHV* tree built with the *GY94* model in `codonPhyML`_ (this is the file ``GY94_summary.csv``):

    .. include:: example_2014Analysis_lactamase_GY94_summary.csv
       :literal:

Here are the results for optimization of the combined *TEM* and *SHV* tree built with the *KOSI07* model in `codonPhyML`_ (this is the file ``KOSI07_summary.csv``):

    .. include:: example_2014Analysis_lactamase_KOSI07_summary.csv
       :literal:

Here are the results for optimization of the  *TEM*-only tree built with the *GY94* model in `codonPhyML`_ (this is the file ``GY94_TEM_summary.csv``):

    .. include:: example_2014Analysis_lactamase_GY94_TEM_summary.csv
       :literal:

Here are the results for optimization of the  *SHV*-only tree built with the *GY94* model in `codonPhyML`_ (this is the file ``GY94_SHV_summary.csv``):

    .. include:: example_2014Analysis_lactamase_GY94_SHV_summary.csv
       :literal:

Here are the results for optimization of the  *TEM*-only tree built with the *KOSI07* model in `codonPhyML`_ (this is the file ``KOSI07_TEM_summary.csv``):

    .. include:: example_2014Analysis_lactamase_KOSI07_TEM_summary.csv
       :literal:

Here are the results for optimization of the  *SHV*-only tree built with the *KOSI07* model in `codonPhyML`_ (this is the file ``KOSI07_SHV_summary.csv``):

    .. include:: example_2014Analysis_lactamase_KOSI07_SHV_summary.csv
       :literal:


Similar files are created in LaTex formatted tables as ``GY94_summary.tex``, ``KOSI07_summary.tex``, etc.

Comparison of per-site likelihoods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The `phyloExpCM` script ``phyloExpCM_SiteLikelihoodComparison.py`` was used to compare the per-site likelihoods between the best experimentally informed model (the *HalpernBruno* model with a fitted :math:`\beta`) and the best traditional model (the *GY94* model with gamma-distributed rates and :math:`\omega`) values for the tree inferred with *GY94* over the entire lactamase (TEM and SHV) sequences.

This was done by creating the input file ``phyloExpCM_SitelikelihoodComparison_infile.txt`` and then running::

    phyloExpCM_SiteLikelihoodComparison.py phyloExpCM_SiteLikelihoodComparison_infile.txt

This created three output files.

The first, ``sitelikelihoods.txt`` lists the differences in per-site likelihoods by site.

The next two are the following plots, which use Tukey box plots to show the differences in site likelihoods, categorizing sites by the `DSSP webserver`_ computed secondary structures (*SS*) and relative solvent accessibilities (*RSA*) in ``1XPB_renumbered.dssp``. Here are these plots:

.. figure:: example_2014Analysis_lactamase_sitelikelihoodcomparison_bySS.jpg
   :width: 45%
   :align: center
   :alt: sitelikelihoodcomparison_bySS.jpg

   The comparison of site likelihoods with sites categorized by secondary structure (``sitelikelihoodcomparison_bySS.pdf``).

.. figure:: example_2014Analysis_lactamase_sitelikelihoodcomparison_byRSA.jpg
   :width: 45%
   :align: center
   :alt: sitelikelihoodcomparison_byRSA.jpg

   The comparison of site likelihoods with sites categorized by relative solvent accessibility (``sitelikelihoodcomparison_byRSA.pdf``).

These plots show that the experimental model is modestly better for most sites, but there are a small number of sites for which the experimental model is far worse.


Conclusions
--------------
The key conclusions that can be drawn from the tables generated by `Optimizing the trees with various substitution models`_ are as follows:

1) The evolutionary models that use the experimentally determined amino-acid preferences substantially outperform the other models. The results are slightly better if the fixation probabilities are estimated using the *HalpernBruno* method.

2) The superior performance of the models with the experimentally determined preferences is due to the site-specific information. If the preferences are randomized among sites or determined from overall average amino-acid frequencies, the models are far worse than any others.

3) The *GY94* and *KOSI07* models perform similarly -- not as good as the experimentally informed models, but better than the randomized models with the experimental data. Overall *KOSI07* gives better likelihoods, but not enough so to justify the additional empirically determined equilibrium frequencies according to `AIC`_. In general, more parameterized versions of the models give better likelihoods, although the improvements are not enough to support the use of branch-local *omega* values.

4) These results are robust to the exact sequences and model used to the build the tree. They are true for trees built using *GY94* or *KOSI07*, and for trees with the combined *TEM* and *SHV* sequences or for trees with just one sequence subset alone.

5) Incorporation of the stringency parameter :math:`\beta` leads to improvements in the experimentally informed models. The fit values of :math:`\beta` are greater than one, suggesting that real evolution is more stringent than the deep mutational scanning.

**Overall, these results demonstrate that using the experimentally determined amino-acid preferences improves phylogenetic fit over existing approaches.**

.. _`mapmuts`: https://github.com/jbloom/mapmuts
.. _`mapmuts documentation`: http://jbloom.github.io/mapmuts/
.. _`ImageMagick convert` : http://www.imagemagick.org/script/convert.php
.. _`FigTree`: http://tree.bio.ed.ac.uk/software/figtree/
.. _`RAxML`: https://github.com/stamatak/standard-RAxML
.. _`Path-O-Gen`: http://tree.bio.ed.ac.uk/software/pathogen/
.. _`phyloExpCM`: https://github.com/jbloom/phyloExpCM
.. _`phyloExpCM documentation`: http://jbloom.github.io/phyloExpCM
.. _`codonPhyML`: http://sourceforge.net/projects/codonphyml/
.. _`Python`: http://www.python.org/
.. _`sbatch`: https://computing.llnl.gov/linux/slurm/sbatch.html
.. _`CF3x4`: http://www.plosone.org/article/info%3Adoi/10.1371/journal.pone.0011230
.. _`Goldman and Yang 1994`: http://mbe.oxfordjournals.org/content/11/5/725.full.pdf
.. _`Yang 1994`: http://www.ncbi.nlm.nih.gov/pubmed/7932792
.. _`Kosiol et al 2007`: http://www.ncbi.nlm.nih.gov/pubmed/17400572
.. _`HYPHY`: http://www.hyphy.org/w/index.php/Main_Page
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`DSSP webserver`: http://www.cmbi.ru.nl/hsspsoap/
.. _`AIC`: http://en.wikipedia.org/wiki/Akaike_information_criterion
.. _`numpy`: http://www.numpy.org/
.. _`matplotlib`: http://matplotlib.org/
.. _`Firnberg et al 2014`: http://www.ncbi.nlm.nih.gov/pubmed/24567513
.. _`Ambler`: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1151176/
.. _`Lahey database`: http://www.lahey.org/Studies/
.. _`Genbank`: http://www.ncbi.nlm.nih.gov/genbank/
.. _`Biopython`: http://biopython.org/wiki/Main_Page
.. _`DataMonkey`: http://www.datamonkey.org/
.. _`EMBOSS needle`: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
.. _`Halpern and Bruno 1998`: http://www.ncbi.nlm.nih.gov/pubmed/9656490
.. _`1XPB`: http://www.rcsb.org/pdb/explore.do?structureId=1XPB
.. _`PDB Goodies`: http://dicsoft2.physics.iisc.ernet.in/pdbgoodies/inputpage.html
.. _`weblogo`: http://weblogo.threeplusone.com/
.. _`PyPdf`: http://pybrary.net/pyPdf/
.. _`An experimentally informed evolutionary model improves phylogenetic fit to divergent lactamase homologs`: http://dx.doi.org/10.1093/molbev/msu220
.. _`mapmuts version 1.0`: https://github.com/jbloom/mapmuts/tree/v1.0
.. _`Robinson-Foulds distance`: http://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric
.. _`phyloExpCM version 0.3`: https://github.com/jbloom/phyloExpCM/tree/v0.3
