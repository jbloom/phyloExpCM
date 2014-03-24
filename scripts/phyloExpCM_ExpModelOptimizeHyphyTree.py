#!python 

"""Runs ``HYPHY`` to optimize tree of known topology.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import time
import re
import tempfile
import subprocess
import random
import multiprocessing
import phyloExpCM.io
import phyloExpCM.submatrix
import phyloExpCM.hyphy
import mapmuts.io
import mapmuts.sequtils
import mapmuts.bayesian



def EncodeNames(code_d, infile, outfile, allowmultiple=False):
    """Encodes names in files.

    Designed to translate unique strings to new values.

    *code_d* is a dictionary keyed by strings and with values as string.

    *infile* is the name of an existing file.

    *outfile* is the name of the new file that is created.

    *allowmultiple* specifies that we allow multiple occurrences
    of each key in *code_d* in *infile*. Is *False* by default, 
    meaning that we expect a single occurrence.

    For each key in *code_d*, this function expects to find a single
    unique occurrence of that key in *infile* if *allowmultiple* is *False*,
    or >= one occurrence if *allowmultiple* is *True*. If this is the case, that
    key is replaced with the value corresponding to the key in *code_d*.
    If there is not a single unique occurrence for any key, raises an error.
    """
    text = open(infile).read()
    for (key, value) in code_d.iteritems():
        if allowmultiple:
            if text.count(key) < 1:
                raise ValueError("Failed to find at least one occurrence of %s" % key)
        elif text.count(key) != 1:
            raise ValueError("Failed to find exactly one occurrence of %s in %s: found %d instead" % (key, infile, text.count(key)))
        text = text.replace(key, value)
    open(outfile, 'w').write(text)


def RemoveBranchSupports(treefile):
    """Removes branch supports from *treefile*.

    *treefile* should be a Newick tree. It can optionally include branch
    supports, which are numbers immediately following the closing parentheses
    for a node. The result of this function is that *treefile* is rewritten
    with all of these branch supports removed.

    This function will not work if the tip node names contain parentheses. To
    ensure that this is not the case, it is recommended that you first run
    *EncodeNames* to encode tip names that you can ensure don't have this feature.
    """
    text = open(treefile).read()
    branchsupport = re.compile('\)\d+(\.\d*){0,1}([eE](\-){0,1}\d+){0,1}\:')
    text = branchsupport.sub('):', text)
    open(treefile, 'w').write(text)


def RunHYPHY(hyphypath, cmdfile, outfile):
    """Runs HYPHY.

    *hyphypath* is a path to the HYPHY executable.

    *cmdfile* is the HYPHY command file.

    *outfile* is the created HYPHY output file. An error is raised if
    this file is not created.
    """
    (fd, errorfile) = tempfile.mkstemp()
    os.close(fd)
    print "\nNow running HYPHY (using the command %s) beginning at %s..." % (hyphypath, time.asctime())
    try:
        subprocess.call(hyphypath.split() + [cmdfile], stderr=open(errorfile, 'w')) # run splitting hyphypath in case set to something like HYPHY CPU=2
    except OSError:
        if os.path.isfile(errorfile):
            os.remove(errorfile)
        raise ValueError("Failed to run HYPHY with the following command: %s\nIs HYPHY actually installed with that name in an accessible path?" % hyphypath)
    errors = open(errorfile).read()
    os.remove(errorfile)
    if errors and errors.strip() != 'Check messages.log details of this run.':
        sys.stderr.write("HYPHY generated the following errors:\n%s" % errors)
    if not os.path.isfile(outfile):
        raise ValueError("HYPHY failed to generate the expected output file %s" % outfile)
    print "HYPHY completed execution at %s, and created the output file %s." % (time.asctime(), outfile)



def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_ExpModelOptimizeHyphyTree.py..."
 
    # parse arguments
    args = sys.argv[1 : ]
    if len(args) != 1:
        raise IOError("Script must be called with exactly one argument specifying the name of the input file.")
    infilename = sys.argv[1]
    if not os.path.isfile(infilename):
        raise IOError("Failed to find infile of %s" % infilename)
    d = phyloExpCM.io.ParseInfile(open(infilename))
    print "\nRead the following key / value pairs from infile %s:\n%s" % (infilename, '\n'.join(["%s %s" % tup for tup in d.iteritems()]))
    hyphypath = phyloExpCM.io.ParseStringValue(d, 'hyphypath')
    fastafile = phyloExpCM.io.ParseStringValue(d, 'fastafile')
    if not os.path.isfile(fastafile):
        raise IOError("Cannot find fastafile %s" % fastafile)
    treefile = phyloExpCM.io.ParseStringValue(d, 'treefile')
    if not os.path.isfile(treefile):
        raise IOError("Cannot find treefile %s" % treefile)
    aapreferences = phyloExpCM.io.ParseStringValue(d, 'aapreferences')
    aapreferences = mapmuts.io.ReadEntropyAndEquilFreqs(aapreferences)
    mapmuts.bayesian.PreferencesRemoveStop(aapreferences)
    siteslist = phyloExpCM.io.ParseStringValue(d, 'siteslist')
    if not os.path.isfile(siteslist):
        raise IOError("Cannot file siteslist %s" % siteslist)
    sites = [line.split()[0] for line in open(siteslist).readlines() if not (line.isspace() or line[0] == '#')]
    try:
        sites = [int(site) for site in sites]
    except ValueError:
        raise ValueError("A line in sitesfile does not specify a valid integer")
    if len(sites) != len(dict([(site, True) for site in sites])):
        raise ValueError("Duplicate sites in sitesfile")
    if not sites:
        raise ValueError("No sites specified")
    sites.sort()
    print "\nThe analysis will include the following sites:\n%s" % ', '.join([str(site) for site in sites])
    if sites[0] < 1:
        raise ValueError("The minimum site number specified in sites is < 1")
    for site in sites:
        if site not in aapreferences:
            raise ValueError("site %d is specified in siteslist, but not in aapreferences" % site)
    fixationmodel = phyloExpCM.io.ParseStringValue(d, 'fixationmodel')
    validmodels = ['FracTolerated', 'HalpernBruno']
    if fixationmodel not in validmodels:
        raise ValueError("Invalid fixationmodel of %s\nValid models are:\n%s" % (fixationmodel, '\n'.join(validmodels)))
    mutationrates = phyloExpCM.io.ParseStringValue(d, 'mutationrates')
    if mutationrates == 'freeparameters':
        print "\nMutation rates will be treated as free parameters."
    else:
        if not os.path.isfile(mutationrates):
            raise ValueError("mutationrates must either be freeparameters or specify an existing file.")
        scalefactor = phyloExpCM.io.ParseFloatValue(d, 'scalefactor')
        if scalefactor <= 0:
            raise ValueError("scalefactor must be > 0")
        lines = open(mutationrates).readlines()
        mutationrates = {'AC':None, 'AG':None, 'AT':None, 'CA':None, 'CG':None}
        for line in lines:
            if line.isspace() or line[0] == '#':
                continue
            entries = line.split()
            if len(entries) != 2:
                raise ValueError("Invalid line in mutationrates file: %s" % line)
            mut = entries[0].strip().upper()
            if mut not in mutationrates:
                raise ValueError("Invalid mutation type of %s in mutationrates file\nShould be AC, AG, AT, CA, or CG" % mut)
            if mutationrates[mut] != None:
                raise ValueError("Duplicate entry for mutation rate %s in mutationrates file" % mut)
            try:
                mutationrates[mut] = float(entries[1]) * scalefactor
            except ValueError:
                raise ValueError("mutationrates file does not specify a valid mutation rate for %s" % mut)
    keeptempfiles = phyloExpCM.io.ParseBoolValue(d, 'keeptempfiles')
    outfileprefix = phyloExpCM.io.ParseStringValue(d, 'outfileprefix')
    if outfileprefix.upper() == 'NONE':
        outfileprefix = ''
    optimizedtreefile = '%soptimizedtree.newick' % outfileprefix
    optimizedtreeresults = '%soptimizedtree_results.txt' % outfileprefix
    codeoptimizedtreefile = 'coded_optimizedtree.newick'
    randomizepreferences = False
    if 'randomizepreferences' in d:
        randomizepreferences = phyloExpCM.io.ParseStringValue(d, 'randomizepreferences')
        if randomizepreferences.upper() != 'FALSE':
            randomizepreferences = int(randomizepreferences)
            if randomizepreferences < 1:
                raise ValueError("randomizepreferences must be False or an integer >= 1")
        else:
            randomizepreferences = False
    fitomega = False
    if 'fitomega' in d:
        fitomega = phyloExpCM.io.ParseStringValue(d, 'fitomega')
        if fitomega.upper() in ['NONE', 'FALSE']:
            fitomega = False
        elif fitomega == 'freeparameter':
            fitomega = 'global_omega'
        else:
            raise ValueError("Invalid value of fitomega: %s\nMust be freeparameter, None, or False" % fitomega)

    # create prxy_file which contains the substitution models for HYPHY
    prxy_file = 'Prxy.ibf'
    print "\nCreating the HYPHY include batch file %s which holds the substitution model constructed using the %s fixation model..." % (prxy_file, fixationmodel)
    phyloExpCM.submatrix.WriteHYPHYMatrices2(prxy_file, sites, aapreferences, fixationmodel, includeselection=fitomega)

    # re-map names in treefile and fastafile to HYPHY acceptable formats in codetreefile and codefastafile
    # code_d maps sequence headers to code names
    codetreefile = 'coded_treefile.newick' 
    codefastafile = 'coded_fastafile.fasta' 
    print "\nCreating files containing sequences with recoded names: %s and %s" % (codetreefile, codefastafile)
    seqs = mapmuts.sequtils.ReadFASTA(fastafile)
    code_d = dict([(seqs[i][0], 'TipSeq%d_' % (i + 1)) for i in range(len(seqs))])
    code_d_inv = dict([(y, x) for (x, y) in code_d.iteritems()])
    if len(code_d) != len(code_d_inv):
        raise ValueError("Some the sequence names are not unique")
    EncodeNames(code_d, treefile, codetreefile)
    RemoveBranchSupports(codetreefile)
    EncodeNames(code_d, fastafile, codefastafile)
    if randomizepreferences:
        codeseqs = mapmuts.sequtils.ReadFASTA(codefastafile)
        random.seed(randomizepreferences)
        print "The sequences will be randomized with seed %d" % randomizepreferences
        codeseqs = phyloExpCM.hyphy.RandomizeCodonAlignment(codeseqs)
        mapmuts.sequtils.WriteFASTA(codeseqs, codefastafile)

    # run HYPHY to optimize branch lengths and any global (shared among sites) model parameters, without selection
    aminoacids = mapmuts.sequtils.AminoAcids()
    optimizetree_cmdfile = 'hyphy_optimizetree_cmds.bf'
    optimizetree_outfile = 'hyphy_optimizetree_output.txt'
    print "\nCreating the HYPHY command file %s to optimize branch lengths and any global model parameters shared among all sites." % optimizetree_cmdfile
    constraints = []
    for r in sites:
        constraints.append('global mu%d := 1.0' % r) # constrain to no scaling factors
    if mutationrates != 'freeparameters':
        for (mut, mutrate) in mutationrates.iteritems():
            constraints.append('global R%s := %g' % (mut, mutrate)) # constrain to specified value
    else:
        constraints.append('global RAC := 1.0') # set one mutation rate to one if branch lengths are free parameters and no date stamping
        for mutrate in ['RAG', 'RAT', 'RCA', 'RCG']:
            constraints.append("global %s :> 1.0e-7" % mutrate) # constrain all mutation rates greater than zero
    if fitomega == 'global_omega':
        constraints.append('global omega :> 1.0e-7') # set omega greater than zero
    phyloExpCM.hyphy.CreateHYPHYCommandFile2(optimizetree_cmdfile, optimizetree_outfile, codefastafile, codetreefile, sites, prxy_file, constraints)
    RunHYPHY(hyphypath, optimizetree_cmdfile, optimizetree_outfile)
    # extract tree and write to file
    print "Extracting the optimized tree to files %s and %s." % (optimizedtreefile, codeoptimizedtreefile)
    tree = phyloExpCM.hyphy.ExtractTree(optimizetree_outfile)
    open(codeoptimizedtreefile, 'w').write(tree)
    EncodeNames(code_d_inv, codeoptimizedtreefile, optimizedtreefile)
    # extract output and write to file
    print "Extracting likelihood and results summary to %s." % optimizedtreeresults
    values = {'Log likelihood':None, 'number of branch lengths':None, 'independent parameters (includes branch lengths)':None}
    if mutationrates == 'freeparameters':
        values['RAG'] = None
        values['RAT'] = None
        values['RCA'] = None
        values['RCG'] = None
    if fitomega == 'global_omega':
        values['omega'] = None
    phyloExpCM.hyphy.ExtractValues(optimizetree_outfile, values)
    f = open(optimizedtreeresults, 'w')
    f.write('Log likelihood: %g\n' % values['Log likelihood'])
    f.write('Branch lengths optimized: %d\n' % round(values['number of branch lengths']))
    f.write('Model parameters optimized: %d\n' % round(values['independent parameters (includes branch lengths)'] - values['number of branch lengths']))
    if mutationrates == 'freeparameters':
        for mut in ['AG', 'AT', 'CA', 'CG']:
            f.write('R_%s: %g\n' % (mut, values["R%s" % mut]))
    if fitomega == 'global_omega':
        f.write('omega: %g\n' % values['omega'])
    f.close()
    # extract branch lengths
    nodebranchmatch = re.compile('(?P<node>((TipSeq\d+_)|(Node\d+)))\:(?P<branch>\d+\.*\d*((e|E)\-{0,1}\d+){0,1})[\(\,\)]')
    matches = [m for m in nodebranchmatch.finditer(tree)]
    if len(matches) != values['number of branch lengths']:
        raise ValueError("Failed to find correct number of branch lengths.\nExpected %d, but found %d." % (values['number of branch lengths'], len(matches)))
    branchlengths = [(m.group('node'), float(m.group('branch'))) for m in matches]
    tempfiles = [prxy_file, codetreefile, codefastafile, optimizetree_cmdfile, optimizetree_outfile, codeoptimizedtreefile] # temporary files created so far

    # delete temporary files
    if not keeptempfiles:
        print "Now deleting temporary files."
        for f in tempfiles:
            if os.path.isfile(f):
                print "Deleting %s" % f
                os.remove(f)

    # script done
    print "\nScript complete."
    sys.stdout.flush()
    time.sleep(0.1)



main() # run the script
