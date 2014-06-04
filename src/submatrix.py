"""Module for generating experimentally determined codon substitution models.

This module is designed for constructing codon-level
substitution matrices from the equilibrium amino-acid frequencies
inferred by the ``mapmuts`` package plus additional information about
the nucleotide-level mutational spectrum.

It can generated substitution matrices usable by ``HYPHY``.

Written by Jesse Bloom, 2013.


Dependencies
----------------
This module requires:

* *mapmuts*

* *numpy*


Functions defined in this module
-----------------------------------

* *NonStopCodons* : returns list of all non-stop codons in alphabetical order.

* *CodonTable* : returns codon to amino-acid translation dictionary.

* *ReadMutSpectrum* : reads mutational spectrum.

* *WriteHYPHYMatrices* : writes substitution matrices in ``HYPHY`` format.

* *WriteHYPHYMatrices2* : more versatile version of *WriteHYPHYMatrices*.

* *BuildSubMatrices* : build codon substitution matrices.

* *StationaryStates* : gets stationary states of substitution matrices.

* *DecomposeToExchangeabilityEquilibrium* : decomposes reversible matrix.

"""


import os
import math
import copy
import re
import numpy
import mapmuts.io
import mapmuts.sequtils
import mapmuts.bayesian


def NonStopCodons():
    """Returns a list of all non-stop codons using standard genetic code.

    The codons are returned in a fixed order (alphabetical order)
    as upper-case three-letter strings.

    >>> codons = NonStopCodons()
    >>> len(codons) == 61
    True
    """
    codons = []
    nts = ['A', 'C', 'G', 'T']
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                codon = "%s%s%s" % (nt1, nt2, nt3)
                if CodonTable()[codon] != '*':
                    codons.append(codon)
    return codons



def CodonTable():
    """Returns a dictionary keyed by codons, values amino acids.

    Returns a dictionary of length 64 keyed by all of the codons
    as upper-case strings. The values are the encoded amino acid
    in the one-letter upper-case code, or * for stop codons.

    >>> genetic_code = CodonTable()
    >>> len(genetic_code) == 64
    True
    >>> aas = genetic_code.values()
    >>> aas.count('*') == 3
    True
    >>> len(dict([(aa, True) for aa in aas])) == 21
    True
    >>> aas.count('M') == aas.count('W') == 1
    True
    >>> aas.count('C') == aas.count('D') == aas.count('E') == aas.count('F') == aas.count('H') == aas.count('K') == aas.count('N') == aas.count('Q') == aas.count('Y') == 2
    True
    >>> aas.count('I') == 3
    True
    >>> aas.count('A') == aas.count('G') == aas.count('P') == aas.count('T') == aas.count('V') == 4
    True
    >>> aas.count('L') == aas.count('R') == aas.count('S') == 6
    True
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    return genetic_code




def ReadMutSpectrum(infile):
    """Reads the mutation spectrum.

    On call, *infile* should be the name of a file specifying the rate of
    mutation from each nucleotide to each other nucleotide. For example::

        AG, TC, 2.4e-5
        GA, CT, 2.3e-5
        AT, TA, 3.0e-6
        AC, TG, 9.0e-6
        GC, CG, 1.9e-6
        GT, CA, 9.4e-6

    For example, the first line is interpreted as giving the probability
    that a site that is initially A is mutated to G in a unit time step,
    or the probability that a site that is initially C is mutated to T
    in a unit time step.

    The returned variable is a dictionary keyed by all of the 2-tuples
    found in *mapmuts.sequtils.NTMutTypes()*. The values for each 2-tuple
    is the rate specified by *infile*.
    """
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile of %s" % infile)
    lines = [line.strip() for line in open(infile).readlines() if not line.isspace()]
    d = {}
    for line in lines:
        entries = [x.strip().upper() for x in line.split(',')]
        if len(entries) != 3:
            raise ValueError("line does not contain 3 entries:\n%s" % line)
        tup1 = (entries[0], entries[1])
        tup2 = (entries[1], entries[0])
        r = float(entries[2])
        if tup1 in mapmuts.sequtils.NTMutTypes():
            d[tup1] = r
        elif tup2 in mapmuts.sequtils.NTMutTypes():
            d[tup2] = r
        else:
            raise ValueError("Invalid mutation types in line:\n%s" % line)
    if len(d) != len(mapmuts.sequtils.NTMutTypes()):
        raise ValueError("Did not find right number of mutation types.")
    return d


def WriteHYPHYMatrices(exchangeabilities, stationary_vecs, outfile):
    """Writes substitution matrices in ``HYPHY`` batch file format.

    *exchangeabilities* is a set of symmetrix exchangeability matrices 
    as returned by the function *DecomposeToExchangeabilityEquilibrium*.
    Because this matrix is assumed to be symmetric, only the lower left
    is used.

    *stationary_vecs* is a set of equilibrium frequencies as returned
    by the function *DecomposeToExchangeabilityEquilibrium*

    *outfile* is the file to which we write all of the substitution 
    matrices. This file will have the format of a ``HYPHY`` include
    batch file, so you might typically give it the extension ``.ibf``.
    *outfileprefix* specifies the name of the output files.

    For each site in *exchangeabilities* and *stationary_vecs* (which
    are assumed to contain the same set of sites), this function creates
    a ``HYPHY`` substitution matrix of name *model1*, *model2*, etc
    for sites 1, 2, ... that is a ``HYPHY`` reversible substitution model.
    It also creates the variables *modelmatrix1*, *modelmatrix2*, etc
    to define the exchangeabilities, and *equilfreqs1*, *equilfreqs2*, etc
    to define the equilibrium frequencies. Each matrix entry multiplies
    the variable *t*, which is assumed to be the notation that ``HYPHY``
    uses to denote branch lengths.

    The indexing of the codons for these matrices are indexed according
    to the order of codons defined in *NonStopCodons()*, starting with 0
    and going to 60. This is alphabetical order. This is the same
    standard order that ``HYPHY`` uses to index codons. 
    """
    ncodons = len(NonStopCodons())
    assert ncodons == 61, "Not 61 non-stop codons -- this is likely to cause problems with HYPHY"
    assert len(stationary_vecs) == len(exchangeabilities), "Not the same number of exchangeability and equilibrium vectors"
    f = open(outfile, 'w')
    for isite in stationary_vecs.iterkeys():
        pi = stationary_vecs[isite]
        s = exchangeabilities[isite]
        f.write("modelmatrix%d = {%d, %d};\n" % (isite, ncodons, ncodons))
        for i in range(ncodons):
            for j in range(i + 1, ncodons):
                x = s[i, j]
                f.write("modelmatrix%d[%d][%d] := t * %g;\nmodelmatrix%d[%d][%d] := t * %g;\n" % (isite, i, j, x, isite, j, i, x))
        f.write("\nequilfreqs%d = {%d, 1};\n" % (isite, ncodons))
        for i in range(ncodons):
            f.write("equilfreqs%d[%d] = %g;\n" % (isite, i, pi[i]))
        f.write("\nModel model%d = (modelmatrix%d, equilfreqs%d, 1);\n\n" % (isite, isite, isite))
    f.close()


def WriteHYPHYMatrices2(outfile, sites, aapreferences, fixationmodel, includeselection=True, fitbeta=False):
    """Creates HYPHY include batch file with codon substitution models.

    Created HYPHY include batch file (``*.ibf`` file), include with the command::
        
        #include "myfile.ibf";

    For each site *r*, this batch file defines a HYPHY Model *Prxy* that defines the
    substitution probabilities for site *r* as a 61 X 61 matrix for the rates
    among all non-stop codons, where the codons are indexed from 0 to 60 in the order
    defined by *NonStopCodons()*, which is alphabetical. Element *Prxy[xi, yi]* is
    the substitution rate from codon *NonStopCodons()[xi]* to *NonStopCodons()[yi]*.
    So if there are are two sites 1 and 2, then this file will define the HYPHY Models
    *P1xy* and *P2xy*. Also defines equilibrium frequencies (*prx*) and
    exchangeabilities (*Srxy*) for each site.

    The substitution model is defined by

        .. math::

           P_{r,xy} = 
           \\begin{cases}
           \mu_r \\times Q_{xy} & \mbox{if $\mathcal{A}\left(x\\right) = \mathcal{A}\left(y\\right)$} \\\\
           \mu_r \\times \omega_{r,\mathcal{A}\left(y\\right)} \\times Q_{xy} \\times F_{r,xy} & \mbox{if $\mathcal{A}\left(x\\right) \\ne \mathcal{A}\left(y\\right)$}. \\\\
           \end{cases}

    where :math:`\mathcal{A}\left(x\\right)` is the amino acid encoded by *x*. The
    mutation rate is :math:`Q_{xy}`, the fixation probability is 
    :math:`F_{r,xy}`, selection for amino-acid *a* is :math:`omega_{r,a}`,
    and a scaling factor is represented by :math:`\mu_r`. If you are not trying to
    detect selection, you will want to constrain these last two variables to be one
    in HYPHY, as with::

        mu1 := 1.0;
        omega1A := 1.0;
        omega1C := 1.0;
        omega1D := 1.0;

    If you are detecting *positive selection* for amino-acid change, but not looking 
    for *additional selection* for specific amino acids, constrain :math:`\omega_{r,a}`
    values to have equal values of :math:`\omega_r` with HYPHY, as in::

        omega1A := omega1;
        omega1C := omega1;
        omega1D := omega1;

    Fixation probabilities are defined as numbers in the substition models. 
    These numbers are computed from the amino-acid preferences
    :math:`\pi_{r,a}`, possibly with these preferences taken to the exponent of
    :math:`\\beta` depending on the value of *fitbeta*. There are two different models that can
    be used to calculate these fixation probabilities. The model to be used is specified by 
    *fixationmodel* as one of these:

       * *FracTolerated* specifies that the fixation probabilities are

            .. math::

               F_{r,xy} =
               \\begin{cases}
               1 & \mbox{if $\mathcal{A}\left(x\\right) = \mathcal{A}\left(y\\right)$ or $\pi_{r,\mathcal{A}\left(y\\right)} \ge \pi_{r,\mathcal{A}\left(x\\right)}$} \\\\
               \left(\\frac{\pi_{r, \mathcal{A}\left(y\\right)}}{\pi_{r, \mathcal{A}\left(x\\right)}}\\right)^{\\beta} & \mbox{otherwise.}
               \end{cases}

       * *HalpernBruno* specifies that the fixation probabilities are

            .. math::

               F_{r,xy} = 
               \\begin{cases}
               1 & \mbox{if $\mathcal{A}\left(x\\right) = \mathcal{A}\left(y\\right)$ or $\pi_{r,\mathcal{A}\left(x\\right)} = \pi_{r,\mathcal{A}\left(y\\right)}$} \\\\
               \\frac{\\beta \\times \ln\left(\pi_{r,\mathcal{A}\left(y\\right)} / \pi_{r,\mathcal{A}\left(x\\right)}\\right)}{1 - \left(\pi_{r,\mathcal{A}\left(x\\right)} / \pi_{r,\mathcal{A}\left(y\\right)}\\right)^{\\beta}} & \mbox{otherwise.}
               \end{cases}

    The mutation rates are equal among all sites, and are 

        .. math::

           Q_{xy} = 
           \\begin{cases}
           0 & \mbox{if $x$ and $y$ differ by more than on nucleotide} \\\\
           R_{m \\rightarrow n} & \mbox{if $x$ differs from $y$ by a single-nucleotide change of $m$ to $n$}.
           \end{cases}

    The rates are defined by :math:`R_{m \\rightarrow n}`, which denotes the probability
    that nucleotide *m* mutates to *n* in a unit time given that the identity is already
    *m*. We place the following constraints on these values:

        * :math:`R_{m \\rightarrow n} = R_{m_c \\rightarrow n_c}` where :math:`m_c` is the complement of 
          nucleotide :math:`m` (so :math:`A_c = T`). This constraint assumes that the
          same mutation process operates on the sequenced and complementary strands.

        * :math:`R_{C \\rightarrow T} = \\frac{R_{A \\rightarrow G} \\times R_{C \\rightarrow A}}{R_{A \\rightarrow C}}`. This constraint makes the mutation rates
          reversible. It is defined more for this property than on any
          biological basis.

    With these these constraints, there are four independent mutation rates, denoted
    by the following HYPHY variables::

        RAC
        RAG
        RAT
        RCA
        RCG

    If you want these mutation rates to be free parameters, don't assign them any constraints. If you want them
    to have fixed values, assign these in HYPHY, such as::

        RAC := 1.23e-2;

    CALLING VARIABLES:

    * *outfile* is the created HYPHY include batch file, such as *Prxy.ibf*. This file is
      overwritten if it already exists.

    * *sites* is a list of the integer sites *r* for which we create models. 

    * *aapreferences* is a dictionary that specifies the amino-acid preferences :math:`\pi_{r,a}`.
      For each site *r* in *sites*, there should be a key *r* in *aapreferences*, with
      *aapreferences[r][PI_A]' giving the preference for amino-acid *A* at site *r*, 
      *aapreferences[r][PI_S]* giving the preference for amino-acid *S* at site *r*, etc.
      All preferences are assumed to be > 0.

    * *fixationmodel* is the model used to calculate the fixation probabilities :math:`F_{r,xy}`
      from the amino-acid preferences. It should be one of the two following strings as
      described above:

        * *FracTolerated*

        * *HalpernBruno*

    * *includeselection* is a Boolean switch which is *True* by default. It has 
      the following possible values:
      
        - *True* (default value) means that each site and amino-acid gets
          its own selection value of :math:`\omega_{r,a}` as described
          above.
          
        - *False* means that all :math:`\omega_{r,a}` values are set to
          one. This is what you would want to use if you don't want to
          consider the possibility of any type of selection beyond
          that defined by the amino-acid preferences.
          
         - The string *global_omega* means that there is a single 
           global :math:`\omega` value for all sites and mutations,
           so only one free parameter for all selection. This could be
           useful if there is stronger selection on nonsynymous mutations
           than captured by the amino-acid preferences.

    * *fitbeta* allows the :math:`\\beta` value that scales the preferences as an
      exponent to have values other than one. By default, *fitbeta* is *False*, which
      means that :math:`\\beta` is one. Possible values:

        - *False* (default value) means that :math:`\\beta` is fixed to be one.

        - *global_beta* means that there is a single :math:`\\beta` value for all sites.
          This could be useful if the strength of selection in the experiments is different
          from that in nature.
    """
    assert includeselection in [False, True, 'global_omega'], "Invalid value of includeselection: %s" % includeselection
    if fitbeta == 'global_beta':
        betaexponent = '^beta'
        betamultiplier = 'beta * '
    elif fitbeta == False:
        betaexponent = ''
        betamultiplier = ''
    else:
        raise ValueError("Invalid value of fitbeta: %s" % fitbeta)
    codontable = CodonTable()
    codons = NonStopCodons()
    ncodons = len(codons)
    aminoacids = mapmuts.sequtils.AminoAcids()
    assert ncodons == 61, "Not 61 non-stop codons -- this is likely to cause problems with HYPHY"
    f = open(outfile, 'w')
    f.write("\n//Values for constrained mutation rates defined by others.\n")
    f.write("global RTG := RAC;\nglobal RTC := RAG;\nglobal RTA := RAT;\nglobal RGT := RCA;\nglobal RGC := RCG;\nglobal RCT := RAG * RCA / RAC;\n");
    f.write("\n// Equilibrium codon frequencies for purely mutation-driven evolution, denoted by qx.\n")
    for x in codons:
        ncg = x.count('C') + x.count('G')
        assert 0 <= ncg <= 3
        if ncg == 0:
            f.write('global q%s := RAC^3;\n' % x)
        elif ncg == 3:
            f.write('global q%s := 1.0;\n' % x)
        else:
            f.write('global q%s := RCA^%d;\n' % (x, 3 - ncg))
    for r in sites:
        f.write("\n//Equilibrium frequencies for site %d\np%dx = {%d, 1};\n" % (r, r, ncodons))
        if includeselection == True:
            f.write("global p%dx_denominator := %s;\n" % (r, ' + '.join(["%g%s * q%s * omega%d%s" % (aapreferences[r]['PI_%s' % codontable[x]], betaexponent, x, r, codontable[x]) for x in codons])))
        else:
            f.write("global p%dx_denominator := %s;\n" % (r, ' + '.join(["%g%s * q%s" % (aapreferences[r]['PI_%s' % codontable[x]], betaexponent, x) for x in codons])))
        for xi in range(ncodons):
            x = codons[xi]
            if includeselection == True:
                f.write("p%dx[%d] := (%g%s * q%s * omega%d%s) / p%dx_denominator;\n" % (r, xi, aapreferences[r]['PI_%s' % codontable[x]], betaexponent, x, r, codontable[x], r))
            else:
                f.write("p%dx[%d] := (%g%s * q%s) / p%dx_denominator;\n" % (r, xi, aapreferences[r]['PI_%s' % codontable[x]], betaexponent, x, r))
        f.write("\n// Exchangeabilities for site %d\nS%dxy = {%d, %d};\n" % (r, r, ncodons, ncodons))
        for xi in range(ncodons):
            x = codons[xi]
            for yi in range(xi + 1, ncodons):
                y = codons[yi]
                ntdiffs = [(x[j], y[j]) for j in range(len(x)) if x[j] != y[j]]
                assert ntdiffs > 0, "Found no differences between codons %s and %s" % (x, y)
                if len(ntdiffs) > 1:
                    # zero if amino acids differ by more than one nucleotide
                    f.write("S%dxy[%d][%d] := 0.0;\n" % (r, xi, yi)) 
                    f.write("S%dxy[%d][%d] := 0.0;\n" % (r, yi, xi)) 
                else:
                    ax = codontable[x]
                    ay = codontable[y]
                    pirx = aapreferences[r]['PI_%s' % ax]
                    piry = aapreferences[r]['PI_%s' % ay]
                    assert pirx > 0 and piry > 0, "Preferences must be > 0"
                    if fixationmodel == 'FracTolerated':
                        if piry >= pirx:
                            frxy = '1.0'
                        else:
                            frxy = "%g%s" % (float(piry / pirx), betaexponent)
                    elif fixationmodel == 'HalpernBruno':
                        if piry == pirx:
                            frxy = '1.0'
                        else:
                            frxy = '%s(%g) / (1.0 - %g%s)' % (betamultiplier, math.log(piry / pirx), float(pirx / piry), betaexponent)
                    else:
                        raise ValueError("Invalid fixationmodel of %s" % fixationmodel)
                    qxy = "R%s" % ''.join(ntdiffs[0])
                    if ax != ay: # nonsynonymous
                        if includeselection == 'global_omega':
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s * omega / %g%s;\n" % (r, xi, yi, r, r, qxy, y, frxy, piry, betaexponent))
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s * omega / %g%s;\n" % (r, yi, xi, r, r, qxy, y, frxy, piry, betaexponent))
                        else:
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s;\n" % (r, xi, yi, r, r, qxy, y, frxy, piry, betaexponent))
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s;\n" % (r, yi, xi, r, r, qxy, y, frxy, piry, betaexponent))
                    else: # synonymous
                        if includeselection == True:
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s / omega%d%s;\n" % (r, xi, yi, r, r, qxy, y, frxy, piry, betaexponent, r, ay))
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s / omega%d%s;\n" % (r, yi, xi, r, r, qxy, y, frxy, piry, betaexponent, r, ay))
                        else:
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s;\n" % (r, xi, yi, r, r, qxy, y, frxy, piry, betaexponent))
                            f.write("S%dxy[%d][%d] := t * mu%d * p%dx_denominator * %s / q%s * %s / %g%s;\n" % (r, yi, xi, r, r, qxy, y, frxy, piry, betaexponent))
        f.write("\n//Define substitution model for site %d as exchangeabilities times equilibrium frequencies\nModel P%dxy = (S%dxy, p%dx);\n\n" % (r, r, r, r))
    f.close()


def BuildSubMatrices(mutspectrumfile, equilibriumfreqsfile, model, scalefactor, makereversible):
    """Codon substitution matrix from mutation spectrum and amino-acid frequencies.

    CALLING VARIABLES:

    * *mutspectrumfile* is a file giving the nucleotide mutation spectrum in a
      format that can be read by *ReadMutSpectrum*. Alternatively, it can also
      be a dictionary in which case it should be of the format returned
      by *ReadMutSpectrum*.

    * *equilibriumfreqsfile* is a file giving the equilibrium amino-acid
      frequencies at each site in a format that can be read by
      *mapmuts.io.ReadEntropyAndEquilFreqs*. If stop codon preferences
      are specified, these are set to zero and all other preferences
      rescaled so that they sum to one over the 20 amino acids for
      each site.

    * *model* is a string that specifies the type of relationship between
      the equilibrium amino-acid frequencies and the fixation probability.
      The acceptable values are:

      - *FracTolerated* means that these frequencies are interpreted
        as the fraction of genetic backgrounds in which a mutation is
        tolerated.

      - *HalpernBruno* means that these frequencies are analyzed as
        described by Halpern and Bruno, MBE, 1998.

      - *PreferWildtype_XX* where *XX* is a number > 0 means that the
        wildtype amino acid is preferred over all others by a factor of *XX*.

    * *scalefactor* is a number (> 0) used to multiply all matrix 
      entries. You would want to set this to a number > 1 if you 
      think that the matrix entries would be very small otherwise.
      This is not critical in a world of perfect numerical accuracy,
      but in the real world it might improve things numerically 
      if entries are scaled to not be too small.

    * *makereversible* is a Boolean switch specifying that we make the 
      matrices reversible. All of the allowed *model* settings make the
      amino-acid portion of the selection reversible. So whether the
      matrices are reversible depends on whether the nucleotide mutation
      spectrum is reversible. If it is, then the matrices will be reversible.
      So if *makereversible* is *True*, we the values for *A -> G, T -> C* 
      and *G -> A, C -> T* to both be equal to their average value. We also
      do the same for *A -> C, T -> G* and *C -> A, G -> T*. This is probably
      reasonable to do if these estimated values are just slightly different.
      If they are very different, then you may want to consider the
      plausibility of using reversible models in the first place.

      The alternative is to set *makereversible* to *False*, in which case
      the two pairs of mutation types set above are kept at their original
      values in *mutspectrumfile*, and are not set to be equal. This will
      make the model non-reversible unless these happen to be precisely
      equal to begin with.

    RETURN VARIABLE:

    The returned variable is a dictionary *submatrices*. This dictionary is
    keyed by residue numbers, and *submatrices[isite]* is the substitution 
    matrix for site *isite* (sites numbered as in *equilibriumfreqsfile*, which
    is typically 1, 2, ... numbering). *submatrices[isite]* is a *numpy.array*
    of dimension *(len(NonStopCodons()), len(NonStopCodons()))*. For two different
    non-stop codons *x* and *y*, the substitution probability from *x* to *y* is
    given by *submatrices[isite][xi, yi]* where *xi = NonStopCodons().index(x)*
    and *yi = NonStopCodons().index(y)* for all *x != y*. The diagonal elements
    are given by::

        *submatrices[isite][xi, xi] = -sum([submatrices[isite][xi, yi] for yi in range(len(NonStopCodons())) if yi != xi])

    This means that rows of the matrices sum to zero. This means that the matrices
    plus the identity matrix are right-stochastic matrices.
    """
    if isinstance(mutspectrumfile, str):
        mutspectrum = ReadMutSpectrum(mutspectrumfile)
    elif not isinstance(mutspectrumfile, dict):
        raise ValueError("mutspectrumfile must be a string (file name) or dictionary giving mutation rates.")
    else:
        mutspectrum = copy.deepcopy(mutspectrumfile)
    if makereversible:
        ag = (mutspectrum[('AG', 'TC')] + mutspectrum[('GA', 'CT')]) / 2.0
        mutspectrum[('AG', 'TC')] = ag
        mutspectrum[('GA', 'CT')] = ag
        ac = (mutspectrum[('AC', 'TG')] + mutspectrum[('GT', 'CA')]) / 2.0
        mutspectrum[('AC', 'TG')] = ac
        mutspectrum[('GT', 'CA')] = ac
    equilibriumfreqs = mapmuts.io.ReadEntropyAndEquilFreqs(equilibriumfreqsfile)
    mapmuts.bayesian.PreferencesRemoveStop(equilibriumfreqs)
    codons = NonStopCodons()
    codontable = CodonTable()
    ncodons = len(codons)
    submatrices = {}
    preferwtmatch = re.compile("^PreferWildtype_(?P<preference>\d+(\.\d*){0,1})$")
    m = preferwtmatch.search(model)
    if m:
        preferwt = float(m.group('preference'))
        if preferwt <= 0:
            raise ValueError("Invalid model of %s. If you are using PreferWildtype, the preference must exceed zero." % model)
    else:
        preferwt = False
    for (site, pi_d) in equilibriumfreqs.iteritems():
        wt_aa = pi_d['WT_AA']
        m = numpy.zeros((ncodons, ncodons))
        for xi in range(ncodons):
            x = codons[xi]
            for yi in range(ncodons):
                y = codons[yi]
                diffs = [(x[j], y[j]) for j in range(len(x)) if x[j] != y[j]]
                if len(diffs) == 0:
                    assert x == y
                    continue # don't add diagonal entries yet
                elif len(diffs) > 1:
                    continue # no substitutions between codons differing at > 1 position
                diff = '%s%s' % diffs[0]
                muttype = [tup for tup in mapmuts.sequtils.NTMutTypes() if diff in tup]
                assert len(muttype) == 1, str(muttype)
                muttype = muttype[0] # the mutation type
                qxy = mutspectrum[muttype] * scalefactor
                xaa = codontable[x]
                yaa = codontable[y]
                pi_x = equilibriumfreqs[site]['PI_%s' % xaa]
                pi_y = equilibriumfreqs[site]['PI_%s' % yaa]
                if model == 'HalpernBruno':
                    if xaa == yaa or pi_y == pi_x:
                        m[xi, yi] = qxy
                    else:
                        m[xi, yi] = qxy * math.log(pi_y / pi_x) / (1.0 - pi_x / pi_y)
                elif model == 'FracTolerated':
                    if xaa == yaa or pi_y >= pi_x:
                        m[xi, yi] = qxy
                    else:
                        m[xi, yi] = qxy * pi_y / pi_x
                elif preferwt:
                    if xaa == yaa:
                        m[xi, yi] = qxy # synonymous amino acids are equivalent
                    elif xaa == wt_aa:
                        m[xi, yi] = qxy / preferwt # moves away from WT disfavored by preferwt
                    elif yaa == wt_aa:
                        m[xi, yi] = qxy * preferwt # moves to WT favored by preferWT
                    else:
                        m[xi, yi] = qxy # all non-WT amino acids are equivalent
                else:
                    raise ValueError("Invalid value for model of %s" % model)
        for xi in range(ncodons):
            m[xi, xi] = -sum([m[xi, yi] for yi in range(ncodons) if yi != xi])
        submatrices[site] = m
    return submatrices


def StationaryStates(submatrices, tol=1e-7):
    """Computes stationary states of substitution matrices.

    *submatrices* is a dictionary keyed by site number and with values
    equal to the substitution matrices, as returned by *BuildSubMatrices*.
    Because each of these substitution matrices is a right-stochastic matrix
    minus the identity matrix, the matrices should have one eigenvalue equal
    to zero and all of the rest of the eigenvalues < 0. The left eigenvector
    corresponding to the zero eigenvalue is the stationary state.

    *tol* specifies the amount that the zero eigenvector (stationary state)
    can deviate from zero before an exception is raised.

    The returned variable is a 2-tuple: *(stationary_vecs, pi_d)*:

        * *stationary_vecs* is a dictionary keyed by all site number keys
          in *submatrices*. The values are the left eigenvectors
          of the substitution matrix -- these are *numpy.array* objects
          of length  *len(NonStopCodons())*. These are the eigenvectors
          in codon space. These eigenvectors are scaled so that all entries
          are positive and sum to one.

        * *pi_d* is a dictionary keyed by all site number keys in
          *submatrices*. The values are dictionaries keyed by *'PI_%s' % aa*
          where *aa* is a one-letter amino-acid code for each amino acid.
          The value is the frequency for that amino acid in the stationary
          state, determined as the sum of the frequencies of all of the 
          encoding codons. These are the stationary states in amino-acid
          space.
    """
    pi_d = {}
    stationary_vecs = {}
    codons = NonStopCodons()
    ncodons = len(codons)
    codontable = CodonTable()
    ordered_aas = [codontable[codon] for codon in codons]
    for (site, m) in submatrices.iteritems():
        mt = numpy.transpose(m) # compute eigenvectors of transpose as we want left ones, and eig returns right ones
        (w, v) = numpy.linalg.eig(mt)
        assert len(w) == ncodons
        max_i = 0
        max_w = w[max_i]
        for i in range(1, len(w)):
            if w[i] > max_w:
                max_w = w[i]
                max_i = i
        if abs(max_w) > tol:
            raise ValueError("Maximum eigenvalue is not close to zero: %f" % tol)
        max_v = numpy.transpose(v[:, max_i])
        max_v = max_v / max_v.sum()
        assert numpy.allclose(numpy.zeros(ncodons), abs(numpy.imag(max_v)))
        max_v = max_v.real
        assert numpy.allclose(numpy.zeros(ncodons), numpy.dot(max_v, m)) # should be true since eigenvalue of zero
        stationary_vecs[site] = max_v
        site_d = dict([(aa, 0) for aa in mapmuts.sequtils.AminoAcids()])
        j = 0
        for (x, aa) in zip(max_v, ordered_aas):
            site_d[aa] += x
        pi_d[site] = site_d
    return (stationary_vecs, pi_d)



def DecomposeToExchangeabilityEquilibrium(submatrices, atol=1.0e-8, rtol=1.0e-5, min_p=2.0e-6):
    """Decomposes reversible matrix to exchangeability and equilibrium.

    Let *M* be a substitution matrix of the type returned by *BuildSubMatrices*.
    Assuming that *M* is acyclic and irreducible, then *M* will have one
    eigenvalue of zero and all other eigenvalues < 0 and > -1 (*M + I* is
    a right-stochastic matrix where *I* is the identity matrix). The left
    eigenvector corresponding to the zero eigenvalue is the equilibrium
    (stationary) state. Call this eigenvector *pi*. In other words, we have
    *0 = pi M*.

    If *M* is also reversible, then it is possibly to write *M* as
    *M = S diag(pi)* where *S* is a symmetrix matrix. We call *S* the
    exchangeability matrix.

    This function takes a set of substitition matrices *submatrices*,
    checks that they are reversible, and computes the equilibrium
    vector and the exchangeability matrix for each substitution matrix.
    If it encounters a substitution matrix that is not reversible, it
    raises an exception.

    CALLING VARIABLES:

    * *submatrices* is a dictionary keyed by site number and with values
      equal to the substitution matrices, as returned by *BuildSubMatrices*.
      Because each of these substitution matrices is a right-stochastic matrix
      minus the identity matrix, the matrices should have one eigenvalue equal
      to zero and all of the rest of the eigenvalues < 0. The left eigenvector
      corresponding to the zero eigenvalue is the stationary state.

    * *rtol* and *atol* specify the numerical tolerance for checking detailed
      balance and the closeness of the principal eigenvector of satisfying
      the condition *0 = pi M*. They have the meanings of the same calling
      arguments to *numpy.allclose*.

    * *min_p* specifies the lowest allowed value for an entry in the 
      equilibrium state vector before an error is raised. The reason for this
      argument is that the current method for computing the exchangeability
      matrix *S* will not work if there are entries in *min_p* that get too
      close to zero.

    RETURN VARIABLES:

    If all of the matrices in *submatrices* are reversible within
    the tolerance specified by *atol* and *rtol*, then the 
    returned variable is the 3-tuple: 
    *(stationary_vecs, pi_d, exchangeabilities)*:

        * *stationary_vecs* is a dictionary keyed by all site number keys
          in *submatrices*. The values are the left eigenvectors
          of the substitution matrix -- these are *numpy.array* objects
          of dimension *(1, len(NonStopCodons()))*. These are the eigenvectors
          in codon space. These eigenvectors are scaled so that all entries
          are positive and sum to one.

        * *pi_d* is a dictionary keyed by all site number keys in
          *submatrices*. The values are dictionaries keyed by *'PI_%s' % aa*
          where *aa* is a one-letter amino-acid code for each amino acid.
          The value is the frequency for that amino acid in the stationary
          state, determined as the sum of the frequencies of all of the 
          encoding codons. These are the stationary states in amino-acid
          space.

        * *exchangeabilities* is a dictionary keyed by all site number keys
          in *submatrices*. The values are *numpy.array* objects of the 
          dimension *(len(NonStopCodons()), len(NonStopCodons()))*. These
          are the exchangeability matrices *S*. They are tested to be
          symmetric using the criteria specified by *rtol* and *atol*.

    If we encounter a non-reversible matrix in *submatrices*, an exception
    is raised.
    """
    pi_d = {}
    stationary_vecs = {}
    exchangeabilities = {}
    codons = NonStopCodons()
    ncodons = len(codons)
    codontable = CodonTable()
    ordered_aas = [codontable[codon] for codon in codons]
    for (site, m) in submatrices.iteritems():
        # get equilibrium state
        mt = numpy.transpose(m) # compute eigenvectors of transpose as we want left ones, and eig returns right ones
        (w, v) = numpy.linalg.eig(mt)
        assert len(w) == ncodons
        max_i = 0
        max_w = w[max_i]
        for i in range(1, len(w)):
            if w[i] > max_w:
                max_w = w[i]
                max_i = i
        if not numpy.allclose(max_w, 0.0, rtol=rtol, atol=atol):
            raise ValueError("Maximum eigenvalue is not close to zero for site %d: %f" % (site, max_w))
        max_v = numpy.transpose(v[:, max_i])
        max_v = max_v / max_v.sum()
        assert numpy.allclose(numpy.zeros(ncodons), abs(numpy.imag(max_v)), rtol=rtol, atol=atol), "Substantial imaginary component in principal eigenvalue"
        max_v = max_v.real
        assert numpy.allclose(numpy.zeros(ncodons), numpy.dot(max_v, m), atol=atol, rtol=rtol), "Failure to verify eigenvector of zero by multiplication"
        stationary_vecs[site] = max_v
        site_d = dict([(aa, 0) for aa in mapmuts.sequtils.AminoAcids()])
        j = 0
        for (x, aa) in zip(max_v, ordered_aas):
            site_d[aa] += x
        pi_d[site] = site_d
        # verify reversibility
        pi_rows = numpy.array([max_v for i in range(ncodons)]) # square matrix with each row being stationary vector
        assert pi_rows.shape == (ncodons, ncodons)
        prod = numpy.dot(pi_rows, m) # this product should be symmetrix if reversible
        prod_t = numpy.transpose(prod)
        if not numpy.allclose(prod, prod_t, atol=atol, rtol=rtol):
            raise ValueError("Not reversible for site %d" % site)
        # get symmetrix matrix
        if not numpy.alltrue(max_v > min_p):
            raise ValueError("For site %d, some entries in the stationary state are <= min_p = %g. This will cause problems with the current procedure for computing the symmetric matrix. Vector is: %s" % (site, min_p, str(max_v)))
        pi_inv = numpy.diag(1.0 / max_v)
        s = numpy.dot(m, pi_inv)
        # verify matrix is actually symmetric
        st = numpy.transpose(s)
        if not numpy.allclose(s, st, atol=atol, rtol=rtol):
            raise ValueError("Matrix is not near symmetric for site %d:\n%s" % (site, str(s)))
        exchangeabilities[site] = s
        assert numpy.allclose(m, numpy.dot(s, numpy.diag(max_v)), atol=atol, rtol=rtol), "exchangeability and equilibrium does not recover matrix for site %d" % site
    return (stationary_vecs, pi_d, exchangeabilities)




# Test with doctest
if __name__ == '__main__':
    import doctest
    doctest.testmod()
