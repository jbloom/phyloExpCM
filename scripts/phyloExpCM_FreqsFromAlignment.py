#!python 

"""Script to get frequencies from sequence alignment

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import re
import time
import phyloExpCM.io
import phyloExpCM.submatrix
import mapmuts.bayesian
import mapmuts.sequtils


def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_FreqsFromAlignment.py..."
    
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    alignmentfile = phyloExpCM.io.ParseStringValue(d, 'alignmentfile')
    if not os.path.isfile(alignmentfile):
        raise IOError("Failed to find alignmentfile of %s" % alignmentfile)
    translateseqs = phyloExpCM.io.ParseBoolValue(d, 'translateseqs')
    includestop = phyloExpCM.io.ParseBoolValue(d, 'includestop')
    pseudocounts = phyloExpCM.io.ParseIntValue(d, 'pseudocounts')
    outputfile = phyloExpCM.io.ParseStringValue(d, 'outputfile')
    requiresubstring = False
    if 'requiresubstring' in d:
        requiresubstring = phyloExpCM.io.ParseStringValue(d, 'requiresubstring')
        if requiresubstring.upper() in ['NONE', 'FALSE']:
            requiresubstring = False

    # read sequences, make sure all of same length
    print "Reading alignment from %s..." % alignmentfile
    seqs = mapmuts.sequtils.ReadFASTA(alignmentfile)
    if not seqs:
        raise ValueError("Failed to find any sequences in alignmentfile")
    if requiresubstring:
        seqs = [(head, seq) for (head, seq) in seqs if requiresubstring in head]
        if not seqs:
            raise ValueError("No sequences in alignmentfile had the requiresubstring substring of %s" % requiresubstring)
    seqlength = len(seqs[0][1])
    for (head, seq) in seqs:
        if len(seq) != seqlength:
            raise ValueError("All sequences are not the same length in alignmentfile. Are you sure they are aligned?")
    if translateseqs:
        if seqlength % 3:
            raise ValueError("Sequences are supposed to be translated, but the lengths are not multiples of 3")
        seqs = mapmuts.sequtils.Translate(seqs, readthrough_stop=True, translate_gaps=True)
        seqlength = len(seqs[0][1])
    count_d = {}
    aminoacids = mapmuts.sequtils.AminoAcids(includestop)
    for r in range(1, seqlength + 1):
        count_d[r] = dict([(aa, pseudocounts) for aa in aminoacids])
    for (head, seq) in seqs:
        for r in range(1, seqlength + 1):
            aa = seq[r - 1]
            if aa == '-':
                continue
            elif aa == '*' and not includestop:
                continue
            else:
                assert aa in aminoacids, "Invalid amino acid %s" % aa
                count_d[r][aa] += 1

    # write outputfile
    print "\nNow writing the frequencies to %s..." % outputfile
    f = open(outputfile, 'w')
    f.write('#SITE\tWT_AA\tSITE_ENTROPY')
    for aa in aminoacids:
        f.write('\tPI_%s' % aa)
    f.write('\n')
    for r in range(1, seqlength + 1):
        pi_d = {}
        n = float(sum(count_d[r].values()))
        if not n:
            raise ValueError("No counts for site %d" % r)
        pi_d = dict([(aa, count_d[r][aa] / n) for aa in aminoacids])
        assert abs(1.0 - sum(pi_d.values())) < 1.0e-7, "Sum of frequencies not close to one for site %d" % r
        f.write('%d\tna\t%f' % (r, mapmuts.bayesian.SiteEntropy(pi_d)))
        for aa in aminoacids:
            f.write('\t%f' % pi_d[aa])
        f.write('\n')
    f.close()
    
    # script done
    print "\nScript complete."


main() # run the script
