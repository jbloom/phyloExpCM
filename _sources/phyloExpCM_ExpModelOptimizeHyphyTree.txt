.. _phyloExpCM_ExpModelOptimizeHyphyTree.py:

==========================================
phyloExpCM_ExpModelOptimizeHyphyTree.py
==========================================

.. contents::
   :depth: 4

Overview
-------------
This script uses codon substitution models derived from experimental measurements to optimize a tree of known topology. It can use either fully experimentally determined evolutionary models (known amino-acid preferences and known mutation rates), or models for which only the amino-acid preferences are known (and the mutation rates are treated as free parameters).

For most purposes, this script should supersede :ref:`phyloExpCM_optimizeHyphyTree.py` for optimizing trees with experimentally determined amino-acid preferences. For experimentally determined models, it has all of the functionality of that earlier script plus more.

You must provide the script with a phylogenetic tree topology and a sequence alignment, as well as a set of site-specific amino-acid preferences. You can also provide a set of mutation rates, or these mutation rate parameters are optimized as free parameters (in either case, the mutation rates are shared among all sites). The amino-acid preferences and mutation rates define a codon-substitution model, which is then used to optimize the branches of the phylogenetic tree (and also the mutation rates if these are not pre-specified from experimental measurements). 

`HYPHY`_ is used to perform the core analyses, and this script primarily runs `HYPHY`_ and manipulates the input / output.


Substitution models
--------------------

Substitution model
~~~~~~~~~~~~~~~~~~~~~~~
The codon substitution models used here are defined based on experimentally measured parameters. The substitution model is determined by site-specific amino-acid preferences plus mutation rate parameters that are shared across all sites. All selection is assumed to be at the amino-acid level (selection on synonymous sites is not included). So substitution is entirely determined by knowing which mutations arise on the nucleotide level, and then how selection acts on any resulting amino-acid change.

In this model, the rate of substitution :math:`P_{r,xy}` from codon *x* to codon *y* is defined as

.. math::
   :label: Prxy

   P_{r,xy} = 
   Q_{xy} \times F_{r,xy} 

where :math:`Q_{xy}` is the rate of mutation from *x* to *y* and :math:`F_{r,xy}` is the probability that a mutation from *x* to *y* fixes once it arises. Below we describe how these mutation rates and fixation probabilities are determined. Note that this equation assumes that mutation rates are identical for all sites, but that the fixation probabilities are site specific (but site independent).


Mutation rates
~~~~~~~~~~~~~~~~~~~

The codon mutation rates :math:`Q_{xy}` are the result of nucleotide mutations.

The nucleotide mutation rates can either be specified based on experimental measurements, or treated as free parameters that are optimized from the sequence data. The nucleotide mutation rates are assumed to be shared among all sites. Let :math:`R_{m \rightarrow n}` denote the rate of mutation from nucleotide *m* to *n* (more specifically, it is the probability that a nucleotide mutates from *m* to *n* in a small unit of time given that the nucleotide is already *m*), and let :math:`m_c` denote the complement of nucleotide *m* (so for example :math:`T_c = A`). 

We assume that codon mutations happen a single nucleotide at a time, so that the :math:`Q_{xy}` in Equation :eq:`Prxy` are given by

.. math::
   :label: Qxy
   
   Q_{xy} = 
   \begin{cases}
   0 & \mbox{if $x$ and $y$ differ by more than on nucleotide} \\
   R_{m \rightarrow n} & \mbox{if $x$ differs from $y$ by a single-nucleotide change of $m$ to $n$}.
   \end{cases}

Most standard approaches to modeling nucleotide substitution (for example, the *GTR* model) assume some set of "equilibrium" nucleotide frequencies. Here we avoid taking such an approach, because these "equilibrium" frequencies lack a clear biological meaning since nucleotide frequencies are ultimately the result of mutation and selection. We therefore define the mutation rates solely in terms of mutation parameters.

General model
++++++++++++++++
There are four nucleotides, so in principle there are 12 mutation rates. 
However in practice we expect to have the constraint:

.. math::
   :label: Rmn

   R_{m \rightarrow n} = R_{m_c \rightarrow n_c}.
   
This constraint means, for example, that the rate from *A* to *G* is equal to the rate from *T* to *C*. The reason we expect this to be true is biological: a mutation from *m* to *n* on the sequenced strand causes a mutation of :math:`m_c` to :math:`n_c` on the complementary strand, so it is not possible to distinguish between these two types of mutations assuming that both strands are subject to the same mutational processes. This assumption reduces the number of independent mutation parameters to the following six parameters:

    1) :math:`R_{A \rightarrow C}`

    2) :math:`R_{A \rightarrow G}` 

    3) :math:`R_{A \rightarrow T}`

    4) :math:`R_{C \rightarrow A}`

    5) :math:`R_{C \rightarrow G}` 

    6) :math:`R_{C \rightarrow T}` 

Given these mutation, rates, it is useful to determine the expected evolutionary equilibrium codon frequencies in the case where all amino acids and stop codons are selectively neutral, and so evolution is driven entirely by the mutation process defined by Equations :eq:`Rmn` and :eq:`Qxy` (in other words, the expected evolutionary equilibrium codon frequencies in the case where :math:`F_{r,xy} = 1` in Equation :eq:`Prxy`). In this case, the expected evolutionary equilibrium frequencies are determined entirely by the mutation rates as expected biologically. It is straightforward to show that for nonzero mutation rates that satisfy Equation :eq:`Rmn`, the expected nucleotide frequency of a *A* or *T* nucleotide is :math:`\frac{1}{2}\times \frac{R_{C \rightarrow A} + R_{C \rightarrow T}}{R_{A \rightarrow C} + R_{A \rightarrow G} + R_{C \rightarrow A} + R_{C\rightarrow T}}` and the expected nucleotide frequency of a *C* or *G* nucleotide is :math:`\frac{1}{2} \times \frac{R_{A \rightarrow C} + R_{A \rightarrow G}}{R_{A \rightarrow C} + R_{A \rightarrow G} + R_{C \rightarrow A} + R_{C\rightarrow T}}`. The expected evolutionary equilibrium frequency :math:`q_x` of codon *x* in the case where all amino acids and stop codons are selectively neutral are simply the product of these nucleotide frequencies, so 

.. math::
   :label: qx

   q_x = \frac{1}{8} \times \frac{\left(R_{A \rightarrow C} + R_{A \rightarrow G} \right)^{\mathcal{N_{CG}}\left(x\right)}\times \left(R_{C \rightarrow A} + R_{C \rightarrow T} \right)^{\left(3 - \mathcal{N_{CG}}\left(x\right)\right)} }{\left(R_{A \rightarrow C} + R_{A \rightarrow G} + R_{C \rightarrow A} + R_{C\rightarrow T}\right)^3}

where :math:`\mathcal{N_{CG}}\left(x\right)` is the number of *C* or *G* nucleotides in codon *x*.

It is easy to verify that :math:`q_x` as defined by Equation :eq:`qx` specifies the stationary state of a Markov process governed by the transition matrix defined by :math:`Q_{xy}` given in Equation :eq:`Qxy`. For instance, if the diagonal elements of this matrix are defined as

.. math::
   :label: Qxx

   Q_{xx} = 1 - \sum\limits_{y \ne x} Q_{xy}

then the matrix :math:`\left[Q_{xy}\right]` is a stochastic matrix, and so in general will have a unique stationary state with eigenvalue one. It is possible to verify by direct substitution that the :math:`q_x` defined by Equation :eq:`qx` is that eigenvector (normalized to sum to one).

Reversible model
++++++++++++++++++

Unfortunately, the process defined by :math:`Q_{xy}` via Equations :eq:`Rmn` and :eq:`Qxy` is not guaranteed to be `reversible`_. The reversibility condition is

.. math::
   :label: Qxy_reversible

   q_x Q_{xy} = q_y Q_{yx}.

To see that this condition is not satisfied, let :math:`x = AAT` and :math:`y = CAT` so that Equation :eq:`Qxy_reversible` reduces  

.. math::
   :label: Rmn_reversible_algebra

   \left(R_{C \rightarrow A} + R_{C \rightarrow T} \right)^3 \times R_{A \rightarrow C} &=& \left(R_{A \rightarrow C} + R_{A \rightarrow G} \right) \times \left(R_{C \rightarrow A} + R_{C \rightarrow T} \right)^2 \times R_{C \rightarrow A} \implies \\
   \left(R_{C \rightarrow A} + R_{C \rightarrow T} \right) \times R_{A \rightarrow C} &=& \left(R_{A \rightarrow C} + R_{A \rightarrow G} \right) \times R_{C \rightarrow A} \implies \\
   R_{C \rightarrow T} \times R_{A \rightarrow C} &=& R_{A \rightarrow G} \times R_{C \rightarrow A}.

In general this last condition is not guaranteed to be true. Therefore, to make the model `reversible`_, it is necessary to define

.. math::
   :label: Rmn_reversible

   R_{C \rightarrow T} = \frac{R_{A \rightarrow G} \times R_{C \rightarrow A}}{R_{A \rightarrow C}}.

The biological meaning of this constraint is not particularly clear or intuitive, but essentially it means that the path to getting from an *A* to a *T* via an intermediate mutation to *C* is equally likely to the path from getting from an *A* to a *T* via an intermediate mutation to *G*.

With the constraint in Equation :eq:`Rmn_reversible` plus Equation :eq:`Rmn`, there are a total of five unique mutation rates:

    1) :math:`R_{A \rightarrow C}`

    2) :math:`R_{A \rightarrow G}`

    3) :math:`R_{A \rightarrow T}`

    4) :math:`R_{C \rightarrow A}`

    5) :math:`R_{C \rightarrow G}`

The equilibrium codon frequencies with this constraint can be simplified to

.. math::
   :label: qx_reversible

   q_x = \frac{1}{8} \times \frac{\left(R_{A \rightarrow C} + R_{A \rightarrow G} \right)^{\mathcal{N_{CG}}\left(x\right)}\times \left(R_{C \rightarrow A} + \frac{R_{A \rightarrow G} \times R_{C \rightarrow A}}{R_{A \rightarrow C}} \right)^{\left(3 - \mathcal{N_{CG}}\left(x\right)\right)} }{\left(R_{A \rightarrow C} + R_{A \rightarrow G} + R_{C \rightarrow A} + \frac{R_{A \rightarrow G} \times R_{C \rightarrow A}}{R_{A \rightarrow C}}\right)^3}.

Therefore, a reversible codon mutation model can be established in terms of five mutation rate parameters. Now Equation :eq:`Qxy_reversible` is satisfied.

In general, only four of the five mutation rates are really free parameters if the branch lengths are free parameters and the sequences are not date-stamped. In this case, we will use the convention of setting :math:`R_{A \rightarrow C} = 1` and then view the other four mutation rates as relative to this rate of one.

Amino-acid preferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For each codon site *r*, we are assumed to have measured the preference :math:`\pi_{r,a}` for each amino-acid *a*, where :math:`1 = \sum_a \pi_{r,a}`. These must be provided as input to the script based on experimental measurements.

The fixation probabilities are determined from these amino-acid preferences. The general idea is that mutations to higher preference amino acids are more likely to fix than mutations to lower preference amino acids. The exact quantitative relationshps are described below.

Fixation probabilitites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We consider two different models for how the fixation probabilities :math:`F_{r,xy}` depend on the amino-acid preferences :math:`\pi_{r,xy}`. We refer to the first model as the *FracTolerated* model since it conceptually considers the preferences to be related to the fraction of genetic backgrounds that tolerate a given mutation. In this model, 

.. math::
   :label: Frxy_FracTolerated

   F_{r,xy} =
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(y\right)} \ge \pi_{r,\mathcal{A}\left(x\right)}$} \\
   \omega \times \frac{\pi_{r,\mathcal{A}\left(y\right)}}{\pi_{\mathcal{A}\left(x\right)}} & \mbox{otherwise.}
   \end{cases}

where :math:`\mathcal{A}\left(x\right)` denotes the amino acid encoded by codon *x* and the meaning of :math:`\omega` is discussed below.

We refer to the second model as the *HalpernBruno* model since it was originally introduced by `Halpern and Bruno, MBE, 1998`_. This model interprets the amino-acid preferences :math:`\pi_{r,xy}` as being related to selection coefficients, and specifies,

.. math::
   :label: Frxy_HalpernBruno

   F_{r,xy} = 
   \begin{cases}
   1 & \mbox{if $\mathcal{A}\left(x\right) = \mathcal{A}\left(y\right)$} \\
   \omega & \mbox{if $\mathcal{A}\left(x\right) \ne \mathcal{A}\left(y\right)$ and $\pi_{r,\mathcal{A}\left(x\right)} = \pi_{r,\mathcal{A}\left(y\right)}$} \\
   \omega \times \frac{\ln\left(\pi_{r,\mathcal{A}\left(y\right)} / \pi_{r,\mathcal{A}\left(x\right)}\right)}{1 - \pi_{r,\mathcal{A}\left(x\right)} / \pi_{r,\mathcal{A}\left(y\right)}} & \mbox{otherwise.}
   \end{cases}

In both Equations :eq:`Frxy_FracTolerated` and :eq:`Frxy_HalpernBruno`, :math:`\omega` is a parameter that scales the rate of nonsynonymous versus synonymous mutations. A large value of :math:`\omega` means that nonsynonymous mutations are elevated relative to the expectation from the amino-acid preferences, and a small value means that they are depressed. Fixing :math:`\omega = 1` means that the amino-acid preferences naturally describe the relative rates of nonsynonymous and synonymous mutations by capturing all of the selection. This script allows models either with :math:`\omega = 1`, or :math:`\omega` as a free parameter. In general, if the selection in natural evolution is more stringent than that in the experiment used to determine the preferences, then we might expect models with :math:`\omega < 1` might provide better fit by adding extra selection against nonsynoymous mutations.

Note that the definitions in both Equations :eq:`Frxy_FracTolerated` and :eq:`Frxy_HalpernBruno` are `reversible`_ with respect to :math:`\pi_{r,x}`, such that

.. math::
   :label: Frxy_reversible

   \pi_{r,x} F_{r,xy} = \pi_{r,y} F_{r,yx}.


Equilibrium frequencies and reversibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Evaluation of the likelihood of a phylogenetic tree requires specifying expected evolutionary equilibrium frequency for the codons at the root node as well as the transition probabilities :math:`P_{r,xy}` given in Equation :eq:`Prxy`. It is also mathematically convenient if the overall substitution model is `reversible`_.

Here we show that for the choices of :math:`Q_{xy}` and :math:`F_{r,xy}` defined above, the substitution model is `reversible`_ and the equilibrium frequencies can be calculated. 

Specifically, we let :math:`p_{r,x}` denote the evolutionary equilibrium frequency of codon *x* at site *r*. We assume stop codon variants are non-viable (they have :math:`\pi_{r,x} = 0`) and so the evolutionary equilibrium frequency for all stop codons is zero. This reduces us to 61 possible values for *x* at each site, corresponding to the the 61 non-stop codons. 

The evolutionary equilibrium codon frequencies are given by:

.. math::
   :label: prx

   p_{r,x} &=& \frac{\pi_{r,\mathcal{A}\left(x\right)} \times q_x }{\sum\limits_y \pi_{r,\mathcal{A}\left(y\right)} \times q_y} \\


To see that these are the evolutionary equilibrium frequencies, define

.. math::
   :label: Prxx

   P_{r,xx} = -\sum\limits_{y \ne x} P_{r,xy}.

With the definition in Equation :eq:`Prxx`, the matrix :math:`\mathbf{I} + \mathbf{P_r}` is a `stochastic matrix`_, where :math:`\mathbf{I}` is the identity matrix and :math:`\mathbf{P_r} = \left[P_{r,xy}\right]`. For plausible values for the mutation rates and amino-acid preferences, :math:`\mathbf{I} + \mathbf{P_r}` will also be irreducible and acyclic. Therefore, it has a unique (within a scaling constant) principal eigenvector :math:`\mathbf{p_r} = \left[p_{r,x}\right]` with an eigenvalue of one that represents the equilibrium frequencies. In other words, :math:`\mathbf{p_r} = \mathbf{p_r} \left(\mathbf{I} + \mathbf{P_r}\right)`. It can be verified that Equation :eq:`prx` defines such an eigenvector by writing this eigenvector equation in element-wise form. Immediately below, I verify this for the case where :math:`F_{r,xy}` is defined by Equation :eq:`Frxy_FracTolerated`, and :math:`\pi_{r,\mathcal{A}\left(y\right)} > \pi_{r,\mathcal{A}\left(x\right)}`. Other cases can be verified similarly. 

.. math::
   :label: Prx_eigenvector

   p_{r,x} &=& p_{r,x} + \sum\limits_y p_{r,y} P_{r,yx} \\
           &=& p_{r,x} + (\sum\limits_{y\ne x} p_{r,y} P_{r,yx}) + p_{r,x} P_{r,xx} \\
           &=& p_{r,x} + (\sum\limits_{y\ne x} p_{r,y} P_{r,yx}) - p_{r,x} \sum\limits_{y\ne x} P_{r,xy} \\
           &=& p_{r,x} + \sum\limits_{y \ne x} \left(p_{r,y} P_{r,yx} - p_{r,x} P_{r,xy}\right) \\
           &=& p_{r,x} + \sum\limits_{y \ne x} \left(p_{r,y} Q_{yx} F_{r,yx} - p_{r,x} Q_{xy} F_{r,xy}\right) \\
           &=& p_{r,x} + \sum\limits_{y \ne x} \left(\pi_{r,\mathcal{A}\left(y\right)} q_y Q_{yx} F_{r,yx} - \pi_{r,\mathcal{A}\left(x\right)} q_x Q_{xy} F_{r,xy}\right) \\
           &=& p_{r,x} + \sum\limits_{y \ne x} \left(\pi_{r,\mathcal{A}\left(y\right)} q_y Q_{yx} F_{r,yx} - \pi_{r,\mathcal{A}\left(x\right)} q_x Q_{xy} F_{r,xy}\right) \\
           &=& p_{r,x} + \sum\limits_{y \ne x} \left(\pi_{r,\mathcal{A}\left(y\right)} q_y Q_{yx} \frac{\pi_{r,\mathcal{A}\left(x\right)}}{\pi_{r,\mathcal{A}\left(y\right)}} - \pi_{r,\mathcal{A}\left(x\right)} q_x Q_{xy} \right) \\
           &=& p_{r,x} + \pi_{r,\mathcal{A}\left(x\right)} \sum\limits_{y \ne x} \left( q_y Q_{yx} - q_x Q_{xy} \right) \\
           &=& p_{r,x} 

where the last line follows from Equation :eq:`Qxy_reversible`. This verifies that :math:`p_{r,x}` is the stationary state of the Markov process defined by :math:`P_{r,xy}`.

To verify that the substitution model is `reversible`_, it is necessary to show that :math:`0 = p_{r,x} P_{r,xy} - p_{r,y} P_{r,yx}`. This is shown below for the case where :math:`F_{r,xy}` is defined by Equation :eq:`Frxy_FracTolerated`, and :math:`\pi_{r,\mathcal{A}\left(y\right)} > \pi_{r,\mathcal{A}\left(x\right)}`. Other cases can be verified similarly.

.. math::
   :label: Prxy_reversible

   0 &=& p_{r,x} P_{r,xy} - p_{r,y} P_{r,yx} \\
     &=& \pi_{r,\mathcal{A}\left(x\right)} \times q_x  \times Q_{xy} \times F_{r,xy} - \pi_{r,\mathcal{A}\left(y\right)} \times q_y \times Q_{yx} \times F_{r,yx} \\
     &=& \pi_{r,\mathcal{A}\left(x\right)} \times q_x \times Q_{xy} \times F_{r,xy} - \pi_{r,\mathcal{A}\left(y\right)} \times q_y \times Q_{yx} \times F_{r,yx} \\
     &=& \pi_{r,\mathcal{A}\left(x\right)} \times q_x \times Q_{xy} - \pi_{r,\mathcal{A}\left(y\right)} \times q_y \times Q_{yx} \times \frac{\pi_{r, \mathcal{A}\left(x\right)}}{\pi_{r, \mathcal{A}\left(y\right)}} \\
     &=& q_x \times Q_{xy} - q_y \times Q_{yx} \\
     &=& 0

where the last line follows from Equation :eq:`Qxy_reversible`. 
The fact that :math:`P_{r,xy}` defines a `reversible`_ Markov process with stationary state :math:`p_{r,x}` means that it is possible to define a symmetric matrix :math:`\mathbf{S_r}` such that

.. math::
   :label: Sr

   \mathbf{S_r} \rm{diag}\left(\ldots, p_{r,x}, \ldots\right) = \mathbf{P_r}

where :math:`\rm{diag}\left(\ldots, p_{r,x}, \ldots\right)` is the diagonal matrix with :math:`p_{r,x}` along its diagonal. Noting :math:`\mathbf{S_r} = \mathbf{P_r} \rm{diag}\left(\ldots, \frac{1}{p_{r,x}}, \ldots\right)`, we have

.. math::
   :label: Srxy

   S_{r,xy} = 
   \begin{cases}
   \frac{P_{r,xy}}{p_{r,y}} = 0 & \mbox{if $x$ and $y$ differ by more than one nucleotide mutation,} \\
   \frac{P_{r,xy}}{p_{r,y}} = \left( \sum\limits_z \pi_{r,\mathcal{A}\left(z\right)} \times q_z \right) \frac{Q_{xy}}{q_y} \frac{F_{r,xy}}{\pi_{r,\mathcal{A}\left(y\right)}} & \mbox{if $x$ and $y$ differ by one nucleotide,} \\
   \frac{P_{r,xx}}{p_{r,x}} & \mbox{otherwise.} \\
   \end{cases}

This matrix is symmetric since :math:`S_{r,xy} = S_{r,yx}` as can be verified from the fact that :math:`\frac{Q_{xy}}{q_y} = \frac{Q_{yx}}{q_x}` and :math:`\frac{F_{r,xy}}{\pi_{r,\mathcal{A}\left(y\right)}} = \frac{F_{r,yx}}{\pi_{r,\mathcal{A}\left(x\right)}}` as is guaranteed by the reversibility conditions in Equations :eq:`Qxy_reversible` and :eq:`Frxy_reversible`.


Requirements
--------------
* `phyloExpCM`_ : this script is part of the `phyloExpCM`_ package.

* `Python`_ : this script has been tested with version 2.7.3.

* `mapmuts`_ : this script has been tested with version 1.0.

* `HYPHY`_ : this script has been tested with version ``HYPHY 2.220131214beta(MP) for Linux on x86_64``.


Running the script
--------------------

To run the script, create an `Input file`_  as described below. Then run the script from the command line using the input file as the single argument, as in::

    phyloExpCM_ExpModelOptimizeHyphyTree.py infile.txt

Note that running this script may take quite a while if `HYPHY`_ takes a long time to run. `HYPHY`_ can also start to consume a lot of memory for relatively big input trees -- so you may need to use a computer with a large amount of RAM.

**If you are running multiple versions of this script then you should place them in separate subdirectories**. This is important because otherwise the scripts may create temporary or permanent files with the same names that interfere with / overwrite each other.

Input file
----------------------------
The input file specifies `key` / `value` pairs on individual lines in the format::

    key1 value1
    key2 value2
    
Blank lines or lines
that begin with # are ignored (i.e. as comment lines). 

The input file should contain the following keys:

* *hyphypath* should contain a command that can be used to run `HYPHY`_ in batch mode. If you have installed `HYPHY`_ in the current search path, this might be something like ``HYPHYMP`` or ``HYPHYMP CPU=4`` (the latter command specifying how many CPUs to use for this multi-thread version) or ``HYPHYDEBUG``. If `HYPHY`_ is not installed in the current search path, you will need to specify a full path to the executable.

* *fastafile* is the name of a FASTA file containing the sequences that comprise the tip nodes. These should be nucleotide coding sequences aligned at the codon level (not at the nucleotide level). Stop codons should be removed, and the aligned sequence length should be a multiple of three. Ambiguous nucleotide identities are not supported. The headers for the sequences should match the names of the tip nodes in *treefile*.

* *treefile* is a Newick-format phylogenetic tree giving the relationship among the sequences in *fastafile*. The names of the tip nodes must match those in *fastafile*. Branch lengths can be specified, but they will be optimized by `HYPHY`_. Branch supports can be present, but they are ignored and removed in the output tree.

* *aapreferences* is the name of a file giving the experimentally determined equilibrium amino-acid preferences. These might be the :math:`\pi_{r,a}` values inferred by the `mapmuts`_ script ``mapmuts_inferpreferences.py`` into the ``*_equilibriumpreferences.txt`` file. Specifically, :math:`\pi_{r,a}` represents the preference of site *r* for amino acid *a* under the constraint :math:`\sum_a \pi_{r,a} = 1`. This preference is the expected equilibrium amino-acid frequency of amino acid *a* at site *r* for a hypothetical evolving population where every amino-acid is equally likely to mutate to every other amino acid (so all frequencies would be 1/20 = 0.05 in the absence of selection).  

  The format of the file specified by *aapreferences* should match that of the ``*_equilibriumpreferences.txt`` file created by the `mapmuts`_ script ``mapmuts_inferpreferences.py``. After an initial header line, each line gives the site number (1, 2, ... numbering) for all
  specified sites, the wildtype amino acid(s) (not used by this script), the site entropy (not used by this script), and then the experimentally determined equilibrium amino-acid frequencies for each amino acid. These equilibrium frequencies can either include or exclude a stop codon (denoted by the * character) as a possible amino acid.
  However, for the purposes of this script, stop codons are not considered possible
  amino acid identities. So if there is value for a stop codon, it is set to zero and the remaining :math:`\pi_{r,a}` values are rescaled so that they sum to one. Here are a few example lines of the input file format::

      #SITE   WT_AA   SITE_ENTROPY    PI_A    PI_C    PI_D    PI_E    PI_F    PI_G    PI_H    PI_I    PI_K    PI_L    PI_M    PI_N    PI_P    PI_Q    PI_R    PI_S    PI_T    PI_V    PI_W    PI_Y    PI_*    
      1   M   2.85229 0.0113803   0.0402777   0.0153121   0.0320438   0.013312    0.00936795  0.026916    0.0192017   0.0224047   0.018926    0.554093    0.0445008   0.0138125   0.0212192   0.0131771   0.0152268   0.0237621   0.0102606   0.0369008   0.0367991   0.0211056    
      2   A   0.968359    0.878703    0.00384431  0.00390501  0.00815666  0.0048199   0.00244426  0.00333017  0.00258064  0.00391181  0.00094265  0.0248087   0.00238593  0.00064321  0.00288323  0.00610788  0.0012339   0.001455    0.0290568   0.0107545   0.00492073  0.00311186    
      3   S   2.93851 0.0295947   0.00304848  0.00227603  0.00668501  0.131168    0.000710502 0.0199653   0.00804841  0.000619715 0.366698    0.0132841   0.0223199   0.0048327   0.0170484   0.00875982  0.090712    0.0171888   0.0102176   0.177257    0.069395    0.000170873


* *mutationrates* specifies how we determine the mutation rate values. 
  As described in the `Mutation rates`_ section, with the the constraints in Equations :eq:`Rmn` and :eq:`Rmn_reversible`, there are a total of five unique mutation rate parameters:

    1) :math:`R_{A \rightarrow C}`

    2) :math:`R_{A \rightarrow G}`

    3) :math:`R_{A \rightarrow T}`

    4) :math:`R_{C \rightarrow A}`

    5) :math:`R_{C \rightarrow G}`

   There are two options for how these rates can be determined, specified by the choice of a value that you put here for the *mutationrates* key: 

    - The string *freeparameters*. In this case, five independent free parameters are created to describe the mutation rates, and are fitted by maximum likelihood. You should use this option if you do not have experimental data that determines the mutation rates ahead of time.

    - The name of a file specifying experimentally measured mutation rates. Specifically, the entry for *A -> G* represents the probability that a given site mutations from nucleotide *A* to nucleotide *G* in a unit of time given that it was already *A* (in other words, it is :math:`\Pr\left(G \rm{\; at\; time\;} t = 1 | A \rm{\; at\; time\;} t = 0\right)`. The five rates described above should be specified in the file as follows::

          AC 9.0e-6
          AG 2.4e-5
          AT 3.0e-6
          CA 9.4e-6
          CG 1.9e-6

      The rates for other mutations are determined from these values as described by Equations :eq:`Rmn` and :eq:`Rmn_reversible`. 

* *scalefactor* is an option that has meaning only if *mutationrates* specifies the name of a file giving the mutation rates. If you instead set *mutationrates* to *freeparameters*, the *scalefactor* option is ignored if it is provided. The actual substitution model specified by Equation :eq:`Prxy` involves multiplying the mutation rates by the fixation probabilities. So if you provide experimentally determined mutation rates and these are small numbers (as will be the case unless you rescale them yourself), then the entries in the substitution matrix will be very small. This could cause numerical underflow. Therefore, it is suggested that you set *scalefactor* to a number > 0 that is sufficiently large that the product of *scalefactor* and the mutation rates is not dramatically smaller than one. For example, if the mutation rates are in the range of :math:`10^{-5}` to :math:`10^{-6}`, you might set *scalefactor* to 10000.0. Note that if there is perfect numerical accuracy, the value of *scalefactor* does not matter as long as it is > 0.

* *fixationmodel* specifies how we convert the amino-acid preferences into fixation probabilities. It should be one of the following strings:

    - *FracTolerated* specifies that we use the model given by Equation :eq:`Frxy_FracTolerated`.

    - *HalpernBruno* specifies that we use the model given by Equation :eq:`Frxy_HalpernBruno`.

* *siteslist* is the name of a file that specifies all of the codon sites that are included in the analysis in 1, 2, ... numbering. Any lines beginning with a *#* character are skipped. The remaining lines should begin with an integer >= 1 that gives a site number, followed by whitespace -- anything following this whitespace is ignored. All of the numbers in this first column are taken to be the numbers of the codon sites included. You can use the ``*_equilibriumpreferences.txt`` file generated by the `mapmuts`_ script ``mapmuts_inferpreferences.py`` as the value for this argument, since it satisfies these formatting requirements. So potentially, you can set *siteslist* to the same file used for *aapreferences*. All sites listed must be present in the file specified by *aapreferences*.

* *keeptempfiles* specifies whether we create temporary files created by the script. A variety of files are created, primarily for the purpose of running `HYPHY`_ (i.e. many of these are ``*.ibf`` files for `HYPHY`_). These files are detailed in `Temporary output files`_. If *keeptempfiles* is *True*, then these files are retained after conclusion of the script. If *keeptempfiles* is *False*, then these files are deleted if the script terminates normally (they are retained if it ends in error). You might set *keeptempfiles* to *True* to create intermediate files to aid in debugging. You might set it to *False* if you want to avoid creating lots of output files.

* *outfileprefix* gives the name of the prefix pre-pended to the output files. You can also make it *None* if you don't want to pre-pend any prefix. 

* *fitomega* is an optional parameters. If it is *None*, *False*, or is not specified, then the :math:`\omega` parameter defined in Equations :eq:`Frxy_FracTolerated` and :eq:`Frxy_HalpernBruno` is fixed to one. This corresponds to assuming that the amino-acid preferences capture all of the selection. If *fitomega* is instead specified with a value of *freeparameter*, then the value of :math:`\omega` is fit by maximum likelihood. This corresponds to the case where there is some additional selection against (or for) nonsynymous mutations beyond that captured by the amino-acid preferences.

* *randomizepreferences* is an optional argument. If it is not specified or set to a value of *False*, nothing unusual is done. Otherwise it should be set to an integer >= 1 that specifies the random number seed used to randomize the amino-acid preferences among sites. You might want to do this if you are testing that the model relies relies on accurate experimental determination of the amino-acid preferences at each site by scrambling these preferences.


Example input file
---------------------------
Here is an example input file::

    # Example input file for phyloExpCM_OptimizeDetectSelectionHyphy.py
    hyphypath HYPHYMP CPU=4
    fastafile Aligned_NPs.fasta
    treefile codonphyml_tree.newick
    aapreferences equilibriumpreferences.txt
    mutationrates mutations.txt
    scalefactor 10000.0
    fixationmodel FracTolerated
    siteslist siteslist.txt
    keeptempfiles True
    outfileprefix None


Output files
----------------
If you are running multiple instances of this script, **it is very important that you run them in separate subdirectories.** Otherwise they may create output files of the same name that overwrite or interfere with each other. 

Any output files that already exist are overwritten.

Two types of output files are created. Primary output files contain the main results. Temporary output files are created to store intermediate results and operations (primarily related to running `HYPHY`_) -- these files are deleted at the end of the script if *keeptempfiles* is *False*.

Primary output files
~~~~~~~~~~~~~~~~~~~~~~~~~
The following primary output files are created. Each of these files begins with the prefix specified by *outfileprefix*, and then has the file name indicated.

    * ``optimizedtree.newick`` is a Newick file containing the tree with the branch lengths optimized by `HYPHY`_ (using Equation :eq:`Prxy`). The tree topology is unchanged from that specified in *newickfile*.

    * ``optimizedtree_results.txt`` is a text file containing the results of optimizing the tree (this is the tree in ``optimizedtree.newick``). This file gives the log likelihood, the number of branch lengths optimized, the number of free model parameters, and the values of any free parameters. If *mutationrates* specifies known mutation rates then there are no free parameters. If *mutationrates* is set to *freeparameters*, then there are five mutation rate parameters as described above (description of reversible model in `Mutation rates`_ section), but since one of these mutation rates is not independent of the branch lengths, there are really on four free parameters (since we have :math:`R_{A \rightarrow C} = 1` as described above). Here are example contents of this file::

        Log likelihood: -58.6512
        Branch lengths optimized: 27
        Model parameters optimized: 4
        R_AG: 2.86282
        R_AT: 0
        R_CA: 1.62803
        R_CG: 1.8586e-15

      Here *R_AG* is the rate of mutation from *A* to *G* (:math:`R_{A \rightarrow G}`). If the mutation rates are not free parameters, then the last four lines giving these rates will be absent.

      If :math:`\omega` is also being fit via setting *fitomega* to *freeparameter*, then the file also specifies the value of this parameter as in::

        omega 0.03132


Temporary output files
~~~~~~~~~~~~~~~~~~~~~~~~~
The following temporary output files are created, largely for running `HYPHY`_. If *keeptempfiles* is *True*, then these files are kept in the directory. If *keeptempfiles* is *False*, then these files are deleted at the conclusion of the script if it terminates normally (if it ends in error, any created files are retained for debugging purposes). The temporary files are:

    * ``Prxy.ibf`` is a `HYPHY`_ include batch file that specifies the substitution model in Equation :eq:`Prxy`.

    * ``coded_fastafile.fasta`` is a version of the sequence aligment in *fastafile* where the sequence headers have been changed to unique names that conform to `HYPHY`_ naming requirements.

    * ``coded_treefile.newick`` is a version of the tree in *treefile* where the sequence names have been changed to match the renamed headers in ``coded_fastafile.fasta``, so that they conform to `HYPHY`_ naming requirements.

    * ``hyphy_optimizetree_cmds.bf`` is a `HYPHY`_ batch file that performs the operations to optimize the tree branch lengths and any global parameters (such as the mutation rates :math:`Q_{xy}`) over all sites.

    * ``hyphy_optimizetree_output.txt`` is a file that contains the results of running `HYPHY`_ on ``hyphy_optimizetree_cmds.bf``.

    * ``coded_optimizedtree.newick`` is a version of the tree with `HYPHY`_ optimized branch lengths with the names still matching those that have been changed to `HYPHY`_ naming as in ``coded_fastafile.fasta`` and ``coded_treefile.newick``. Except for the names being coded, this matches the primary output file with suffix ``optimizedtree.newick``.

.. include:: weblinks.txt
