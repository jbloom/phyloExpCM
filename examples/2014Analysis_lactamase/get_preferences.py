"""Gets preferences from Firnberg et al amino-acid "fitnesses".

Uses the reported fitnesses except:

* Fitness of wildtype at each site is always set to one.

* If no fitness is provided, uses the average fitness of all amino acids
  with reported fitnesses plus wildtype at that site.

Written by Jesse Bloom, 2014."""


import math
import mapmuts.sequtils
import mapmuts.bayesian


def main():
    """Main body of script."""

    # set up variables
    infile = 'Firnberg_missense_mutation_fitnesses.csv'
    outfile = 'amino_acid_preferences.txt'
    numberingfile = 'sequential_to_Ambler.csv'
    aminoacids = mapmuts.sequtils.AminoAcids()
    residuerange = (1, 286) # first and last residues in sequential numbering

    # read input
    lines = open(infile).readlines()
    plates = [float(c) for c in lines[1].split(',') if c and not c.isspace()]
    fitnesses = {} # keyed by residue r, the amino acid a
    wts = {} # keyed by residue r, value is wildtype
    numbering_f = open(numberingfile, 'w')
    numbering_f.write('#SEQUENTIAL,AMBLER\n')
    for r in range(residuerange[0], residuerange[1] + 1):
        fitnesses[r] = {}
        lineindex = 2 + (r - 1) * len(aminoacids)
        wt = None
        for a in aminoacids:
            line = lines[lineindex]
            lineindex += 1
            entries = line.split(',')
            assert entries[2] == a, "Amino acid mismatch, expected %s on line %s" % (a, line)
            if wt == None:
                numbering_f.write('%d,%d\n' % (r, int(entries[0])))
                wt = entries[1]
                assert wt in aminoacids, "Invalid amino acid %s on line %s" % (wt, line)
            assert wt == entries[1], "Wildtype mismatch, expected %s on line %s" % (wt, line)
            wts[r] = wt
            reportedfitness = entries[4 + len(plates)]
            if reportedfitness:
                reportedfitness = float(reportedfitness)
            else:
                reportedfitness = None
            if a == wt:
                fitnesses[r][a] = 1.0
            else:
                fitnesses[r][a] = reportedfitness
    numbering_f.close()

    # Now compute preferences
    f = open(outfile, 'w')
    f.write("#SITE\tWT_AA\tSITE_ENTROPY\t%s\n" % '\t'.join(['PI_%s' % a for a in aminoacids]))
    for r in range(residuerange[0], residuerange[1] + 1):
        knownfitnesses = [x for x in fitnesses[r].values() if x != None]
        meanknownfitness = sum(knownfitnesses) / float(len(knownfitnesses))
        pi = {}
        for a in aminoacids:
            if fitnesses[r][a] == None:
                fitnesses[r][a] = meanknownfitness # assign unknown fitnesses mean for that residue
            pi[a] = fitnesses[r][a]
        totfitnesses = sum(fitnesses[r].values())
        f.write('%d\t%s\t%.3f\t%s\n' % (r, wts[r], mapmuts.bayesian.SiteEntropy(pi), '\t'.join(['%g' % (fitnesses[r][a] / totfitnesses) for a in aminoacids])))
    f.close()


if __name__ == '__main__':
    main() # run the script
