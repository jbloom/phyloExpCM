"""Script to check that per-site likelihoods sum to total likelihood.

This is a sanity check. If the site likelihoods are computed correctly, 
they should sum to the total likelihood (perhaps with minor rounding errors).

Written by Jesse Bloom."""


import os
import glob
import phyloExpCM.hyphy


def main():
    """Main body of script."""
    tol = 0.01 # tolerance for difference
    dirs = [os.path.split(dir)[0] for dir in glob.glob('./*/sitelikelihoods.txt')]
    for dir in dirs:
        print "\nAnalyzing %s" % dir
        sitesum = sum([float(line.split()[1]) for line in open('%s/sitelikelihoods.txt' % dir).readlines()[1 : ]])
        if os.path.isfile('%s/hyphy_output.txt' % dir):
            totalfile = '%s/hyphy_output.txt' % dir
        elif os.path.isfile('%s/optimizedtree_results.txt' % dir):
            totalfile = '%s/optimizedtree_results.txt' % dir
        else:
            raise ValueError("Cannot find totalfile for %s" % dir)
        ll = phyloExpCM.hyphy.ExtractLikelihoodAndNParameters(totalfile)[0]
        print "Site sum is %g and total is %g" % (sitesum, ll)
        if abs(sitesum - ll) > tol:
            raise ValueError("Difference exceeds tolerance")


if __name__ == '__main__':
    main() # run the script
