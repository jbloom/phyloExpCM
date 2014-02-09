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

* *BuildSubMatrices* : build codon substitution matrices.

* *StationaryStates* : gets stationary states of substitution matrices.

* *DecomposeToExchangeabilityEquilibrium* : decomposes reversible matrix.

"""


import os
import math
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



def BuildSubMatrices(mutspectrumfile, equilibriumfreqsfile, model, scalefactor, makereversible):
    """Codon substitution matrix from mutation spectrum and amino-acid frequencies.

    CALLING VARIABLES:

    * *mutspectrumfile* is a file giving the nucleotide mutation spectrum in a
      format that can be read by *ReadMutSpectrum*.

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
    mutspectrum = ReadMutSpectrum(mutspectrumfile)
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
