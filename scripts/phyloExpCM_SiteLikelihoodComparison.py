#!python 

"""Script to compare per-site likelihoods from two models.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import re
import time
import mapmuts.dssp
import phyloExpCM.io
import phyloExpCM.plot


def GetBin(rsa, nbins):
    """Returns bin assignment string for an *rsa* when there are *nbins* bins."""
    binsize = 1.0 / nbins
    if rsa < binsize:
        return '$< %.2f$' % binsize
    for ibin in range(nbins - 2):
        if binsize * (ibin + 1) <= rsa < binsize * (ibin + 2):
            return '$%.2f - %.2f$' % (binsize * (ibin + 1), binsize * (ibin + 2))
    else:
        return '$> %.2f$' % ((nbins - 1) * binsize)


def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_SiteLikelihoodComparison.py..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    sitelikelihoodfiles = phyloExpCM.io.ParseFileList(d, 'sitelikelihoodfiles')
    assert len(sitelikelihoodfiles) == 2, "Failed to find exactly two sitelikelihoodfiles"
    modelnames = phyloExpCM.io.ParseStringValue(d, 'modelnames')
    modelnames = modelnames.split()
    assert len(modelnames) == len(sitelikelihoodfiles), "Failed to find expected number of modelnames"
    modelnames = [name.strip().replace('_', ' ') for name in modelnames]
    outfileprefix = phyloExpCM.io.ParseStringValue(d, 'outfileprefix')
    if outfileprefix.upper() == 'NONE':
        outfileprefix = ''
    dsspfile = phyloExpCM.io.ParseStringValue(d, 'dsspfile')
    if dsspfile.upper() in ['NONE', 'FALSE']:
        dsspfile = False
    else:
        dsspfile = dsspfile.strip()
        if not os.path.isfile(dsspfile):
            raise IOError("Cannot find dsspfile of %s" % dsspfile)
        dsspchain = phyloExpCM.io.ParseStringValue(d, 'dsspchain')
        if dsspchain.upper() in ['NONE', 'FALSE']:
            dsspchain = None

    # read site likelihoods
    sitelikelihoods = []
    for f in sitelikelihoodfiles:
        print "Reading site likelihoods from %s..." % f
        lines = open(f).readlines()[1 : ]
        print f, len(lines)
        sitelikelihoods.append(dict([(int(line.split()[0]), float(line.split()[1])) for line in lines]))
    sites = sitelikelihoods[0].keys()
    sites.sort()
    for sitedata in sitelikelihoods:
        isites = sitedata.keys()
        isites.sort()
        assert isites == sites, "The specified sites are not the same for all sitelikelihoodfiles: %r and %r" % ([site for site in isites if site not in sites], [site for site in sites if site not in isites])
    print "Read data for %d sites." % len(sites)

    # parse information on sites from dsspfile
    nrsabins = 3 # bin RSA into four bins
    ss_classification_types = ['helix', 'strand', 'loop', 'unknown']
    rsa_classification_types = ['$< 0.33$', '$0.33 - 0.67$', '$> 0.67$', 'unknown']
    if dsspfile:
        dssp = mapmuts.dssp.ReadDSSP(dsspfile, 'Tien2013', chain=dsspchain)
    else:
        dssp = {}
    rsa_classification = {}
    ss_classification = {}
    for site in sites:
        if site in dssp:
            ss_classification[site] = dssp[site]['SS_CLASS']
            rsa_classification[site] = GetBin(dssp[site]['RSA'], nrsabins)
        else:
            rsa_classification[site] = 'unknown'
            ss_classification[site] = 'unknown'

    # make the plot
    for (classification_name, classifications, classificationtypes) in [('RSA', rsa_classification, rsa_classification_types), ('SS', ss_classification, ss_classification_types)]:
        if not dsspfile:
            classificationtypes = ['unknown']
        plotfile = '%ssitelikelihoodcomparison_by%s.pdf' % (outfileprefix, classification_name)
        print "Plotting the data for the %s classification to %s" % (classification_name, plotfile)
        assert len(sitelikelihoods) == 2, "Currently only works for 2 sitelikelihoods"
        xvalues = dict([(site, sitelikelihoods[0][site] - sitelikelihoods[1][site]) for site in sites])
        xlabel = '$\Delta$ (log likelihood):\n%s minus %s' % (modelnames[0], modelnames[1])
        ylabel = {'SS':'secondary structure',
                  'RSA':'relative solvent accessibility'}[classification_name] 
        phyloExpCM.plot.PlotSiteLikelihoods(sites, classifications, classificationtypes, xvalues, xlabel, ylabel, plotfile, symmetrize_axis=True, title='', alpha=0.2)

    # write the text output
    textfile = '%ssitelikelihoods.txt' % outfileprefix
    print "Writing the data in text format to %s" % textfile
    persitedata = []
    for site in sites:
        logl1 = sitelikelihoods[0][site]
        logl2 = sitelikelihoods[1][site]
        diff = logl1 - logl2
        persitedata.append((diff, '%d\t%g\t%g\t%g\t%s\t%s' % (site, diff, logl1, logl2, ss_classification[site], rsa_classification[site])))
    persitedata.sort()
    persitedata = '\n'.join([line for (diff, line) in persitedata])
    f = open(textfile, 'w')
    f.write('#SITE\t%s_likelihood_minus_%s_likelihood\t%s_likelihood\t%s_likelihood\tsecondary_structure\trelative_solvent_accessibility\n%s' % (modelnames[0].replace(' ', '_'), modelnames[1].replace(' ', '_'), modelnames[0].replace(' ', '_'), modelnames[1].replace(' ', '_'), persitedata))
    f.close()
    
    # script done
    print "\nScript complete."


main() # run the script
