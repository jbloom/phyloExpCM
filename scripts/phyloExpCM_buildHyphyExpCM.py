#!python 

"""Script to build experimentally determined codon models for ``HYPHY``.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import re
import time
import phyloExpCM.io
import phyloExpCM.submatrix
import mapmuts.bayesian


def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_buildHyphyExpCM.py..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    aapreferences = phyloExpCM.io.ParseStringValue(d, 'aapreferences')
    if not os.path.isfile(aapreferences):
        raise IOError("Failed to find aapreferences file of %s" % aapreferences)
    mutspectrum = phyloExpCM.io.ParseStringValue(d, 'mutspectrum')
    if not os.path.isfile(mutspectrum):
        raise IOError("Failed to find mutspectrum file of %s" % mutspectrum)
    scalefactor = phyloExpCM.io.ParseFloatValue(d, 'scalefactor')
    if scalefactor <= 0:
        raise ValueError("scalefactor must be > 0")
    makereversible = phyloExpCM.io.ParseBoolValue(d, 'makereversible')
    model = phyloExpCM.io.ParseStringValue(d, 'model')
    validmodels = ['FracTolerated', 'HalpernBruno']
    if (model not in validmodels) and not re.search("^PreferWildtype_\d+(\.\d*){0,1}$", model):
        raise ValueError("Invalid value of %s for model." % model)
    evolequilfreqs = phyloExpCM.io.ParseStringValue(d, 'evolequilfreqs')
    hyphyExpCMs = phyloExpCM.io.ParseStringValue(d, 'hyphyExpCMs')

    # now build the matrices, exchangeabilities, equilibrium frequencies
    print "\nBuilding substitution matrices for model %s..." % model
    submatrices = phyloExpCM.submatrix.BuildSubMatrices(mutspectrum, aapreferences, model, scalefactor, makereversible)
    print "Computing exchangeabilities and equilibrium frequencies..."
    (stationary_vecs, pi_d, exchangeabilities) = phyloExpCM.submatrix.DecomposeToExchangeabilityEquilibrium(submatrices)

    # write evolequilfreqs
    print "\nNow writing the evolutionary equilibrium frequencies to %s..." % evolequilfreqs
    f = open(evolequilfreqs, 'w')
    f.write('#SITE\tWT_AA\tSITE_ENTROPY')
    for aa in mapmuts.sequtils.AminoAcids():
        f.write('\tPI_%s' % aa)
    f.write('\n')
    sites = pi_d.keys()
    sites.sort()
    for site in sites:
        f.write('%d\tna\t%f' % (site, mapmuts.bayesian.SiteEntropy(pi_d[site])))
        for aa in mapmuts.sequtils.AminoAcids():
            f.write('\t%f' % pi_d[site][aa])
        f.write('\n')
    f.close()

    # write hyphyExpCMs
    print "\nNow writing the experimentally determined substitution models in HYPHY format to %s..." % hyphyExpCMs
    phyloExpCM.submatrix.WriteHYPHYMatrices(exchangeabilities, stationary_vecs, hyphyExpCMs)

    # script done
    print "\nScript complete."



main() # run the script
