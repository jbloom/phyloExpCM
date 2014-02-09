#!python 

"""Script for running ``codonPhyML``.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import time
import phyloExpCM.io
import phyloExpCM.runphyml


def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_runcodonPhyML..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    seqfile = phyloExpCM.io.ParseStringValue(d, 'seqfile')
    if not os.path.isfile(seqfile):
        raise IOError("Failed to find seqfile of %s" % seqfile)
    path = phyloExpCM.io.ParseStringValue(d, 'codonPhyML_path')
    seed = phyloExpCM.io.ParseIntValue(d, 'seed')
    outprefix = phyloExpCM.io.ParseStringValue(d, 'outprefix')
    model = phyloExpCM.io.ParseStringValue(d, 'model')

    # run codonphyml
    print "\nNow running codonPhyML starting at %s..." % time.asctime()
    start_time = time.time()
    ll = phyloExpCM.runphyml.RunCodonphyml(seqfile, path, seed, outprefix, model)
    print "\nCompleted running codonPhyML at %s, or after %.1f minutes." % (time.asctime(), (time.time() - start_time) / 60.)
    print "\nThe log likelihood is %g." % ll

    print "\nScript complete."



main() # run the script
