#!python 

"""Script that plots log likelihood versus number of parameters using *matplotlib*.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import re
import time
import phyloExpCM.io
import phyloExpCM.plot




def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_plotLogLvsNParams.py..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    plotfile = phyloExpCM.io.ParseStringValue(d, 'plotfile')
    datafile = phyloExpCM.io.ParseStringValue(d, 'datafile')
    if not os.path.isfile(datafile):
        raise IOError("datafile does not exist: %s" % datafile)

    # parse data from datafile into data_d
    print "Now reading data from %s..." % datafile
    data_d = {}
    lines = [line for line in open(datafile).readlines() if not (line.isspace() or line[0] == '#')]
    for line in lines:
        entries = line.split(',')
        if len(entries) < 7:
            raise ValueError("line in %s does not contain at least 7 comma-delimited entries:\n%s" % (datafile, line))
        ll = float(entries[2])
        nparams = int(entries[3])
        plotting_group = entries[6].strip()
        if plotting_group in data_d:
            data_d[plotting_group].append((nparams, ll))
        else:
            data_d[plotting_group] = [(nparams, ll)]

    # make the plot
    print "Now making the plot %s..." % plotfile
    phyloExpCM.plot.PlotLogLvsNParams(plotfile, data_d)

    # script done
    print "\nScript complete."
    sys.stdout.flush()
    time.sleep(0.1)



main() # run the script
