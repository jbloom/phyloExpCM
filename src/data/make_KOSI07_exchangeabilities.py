"""Makes ``HYPHY`` include batch file with KOSI07 exchangeabilities.

This script reads in the supplementary information file *ECMunrest.dat* from
Kosiol, Holmes, and Goldman, "An empirical model for protein sequence evolution," 
Mol Biol Evol, 24:1464-1479 (2007)
http://www.ncbi.nlm.nih.gov/pubmed/17400572

It converts the data in that file into a new file *KOSI07_exchangeabilities.ibf*
which is an include batch file for ``HYPHY``.

This file defines a ``HYPHY`` 61x61 matrix called *KOSI07_exchangeabilities* which
contains the exchangeabilities for all non-diagonal entries in the KOSI07 model.
If you include this file in a ``HYPHY`` batch file using::

    #include "KOSI07_exchangeabilities.ibf";

and then then create a 61x1 vector called *codonfreqs* that contains the 
codon equilibrium frequencies and then use the ``HYPHY`` command::

    Model model = (KOSI07_exchangeabilities, codonfreqs, 1);

you will create a ``HYPHY`` substitution model that is reversible.

In *KOSI07_exchangeabilities* (and in the *codonfreqs* variable that you should 
create) the codons are ordered in alphabetical order ("AAA", "AAC", ..., "TTT")
with the exception that the three stop codons ("TAA", "TAG", "TGA")
are excluded. This requires re-ording of the codons in the *ECMunrest.dat*
which orders them using another non-alphabetical ordering scheme.

The exchangeabilities all multiply a branch length denoted *t*, a
rate parameter denoted as *rateparaemter*, and also
multiply a parameter called *omega* if the mutation is non-synonymous
and a parameter called *kappa* a number of times equal to the number
of tranversions. This sets up the model denoted as *ECM+F+w+1k(tv)*
in the original Kosiol et al, 2007 paper::

    KOSI07_exchangeabilities[58][60] := t * rateparameter * 16.0115;
    KOSI07_exchangeabilities[60][58] := t * rateparameter * 16.0115;
    KOSI07_exchangeabilities[57][60] := t * rateparameter * omega * kappa * 2.39582;
    KOSI07_exchangeabilities[60][57] := t * rateparameter * omega * kappa * 2.39582;
"""


import mapmuts.sequtils


def NTransversions(codon1, codon2):
    """Returns the number of transversions that separate two codons."""
    assert len(codon1) == len(codon2) == 3
    ntransversions = 0
    for i in range(3):
        (nt1, nt2) = (codon1[i], codon2[i])
        if nt1 == nt2:
            pass # no mutation
        elif (nt1 == 'A' and nt2 == 'G') or (nt1 == 'G' and nt2 == 'A') or (nt1 == 'C' and nt2 == 'T') or (nt1 == 'T' and nt2 == 'C'):
            pass # transition
        else:
            ntransversions += 1
    return ntransversions


def main():
    """Main body of script."""
    infile = 'ECMunrest.dat'
    outfile = 'KOSI07_exchangeabilities.ibf'
    print "Reading exchangeabilities from %s and creating file %s." % (infile, outfile)
    ncodons = 61
    stopcodons = ['TAA', 'TAG', 'TGA']
    nts = ['A', 'C', 'G', 'T'] # codons in alphabetical order
    hyphy_codons = [] # list of codons in HYPHY (alphabetical) order
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                codon = "%s%s%s" % (nt1, nt2, nt3)
                if codon not in stopcodons:
                    hyphy_codons.append(codon)
    assert len(hyphy_codons) == ncodons
    incodons = 'TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TGT TGC TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG'.split() # list with codons in the order in infile, as taken from that file
    assert len(incodons) == ncodons
    indexmapping = {} # maps infile index to outfile (hyphy) index
    i_incodon = 0
    for incodon in incodons:
        indexmapping[i_incodon] = hyphy_codons.index(incodon)
        i_incodon += 1
    assert len(indexmapping) == ncodons
    assert len(dict([(i, codon) for (codon, i) in indexmapping.iteritems()])) == ncodons
    lines = open(infile).readlines()[ : ncodons - 1]
    exchangeabilities = {} # indexed by (codon_1_index, codon_2_index) hyphy indices
    f = open(outfile, 'w')
    f.write('KOSI07_exchangeabilities = {%d, %d};\n' % (ncodons, ncodons))
    for iline in range(ncodons - 1):
        icodon1 = indexmapping[iline + 1]
        entries = [float(x) for x in lines[iline].split()]
        assert len(entries) == iline + 1
        for ientry in range(iline + 1):
            icodon2 = indexmapping[ientry]
            assert icodon1 != icodon2, "identical codons"
            codon1 = hyphy_codons[icodon1]
            codon2 = hyphy_codons[icodon2]
            aa1 = mapmuts.sequtils.Translate([('head', codon1)])
            aa2 = mapmuts.sequtils.Translate([('head', codon2)])
            if aa1 != aa2:
                omega = "omega * "
            else:
                omega = ""
            kappa = ''.join(["kappa * " for i in range(NTransversions(codon1, codon2))])
            x = float(entries[ientry])
            if x == 0:
                f.write('KOSI07_exchangeabilities[%d][%d] := 0.0;\nKOSI07_exchangeabilities[%d][%d] := 0.0;\n' % (icodon1, icodon2, icodon2, icodon1))
            else:
                f.write('KOSI07_exchangeabilities[%d][%d] := t * rateparameter * %s%s%g;\nKOSI07_exchangeabilities[%d][%d] := t * rateparameter * %s%s%g;\n' % (icodon1, icodon2, omega, kappa, x, icodon2, icodon1, omega, kappa, x))
    f.close()


main() # run the script

