#!python 

"""Script that runs ``HYPHY`` to optimize likelihoods for trees of known topology.

This script is part of the ``phyloExpCM`` package.

Written by Jesse Bloom."""


import sys
import os
import time
import re
import tempfile
import subprocess
import random
import phyloExpCM.io
import phyloExpCM.hyphy
import mapmuts.sequtils



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


def ParseModelOptions(line):
    """Parses options from model specification line.

    *line* is a string specifying the model parameters, such as::

        GY94_CF3x4_omega-global-one_rates-one

      or::

        KOSI07_F_omega-global-gamma6_rates-gamma6

    The returned variable is a 5-tuple giving the parsed values as
    *(model, equilfreqs, omega_scope, omega_classes, rate_classes)*.

    This tuple is appropriate for passing to 
    *phyloExpCM.hyphy.CreateHYPHYCommandFile*.
    """
    entries = line.split('_')
    if len(entries) != 4:
        raise ValueError("model must have four underscore separated specifications, this model does not:\n%s" % line)
    modeltype = entries[0]
    if modeltype not in ['GY94', 'KOSI07']:
        raise ValueError("Unrecognized model type of %s" % modeltype)
    equilfreqs = entries[1]
    if not ((modeltype == 'GY94' and equilfreqs in ['CF3x4', 'F3x4']) or (modeltype == 'KOSI07' and equilfreqs == 'F')):
        raise ValueError("Not a valid equilibrium frequencies method of %s for model %s." % (equilfreqs, modeltype))
    omega = entries[2].split('-')
    if omega[0] != 'omega' or len(omega) != 3:
        raise ValueError("Invalid omega specification in model:\n%s" % line)
    omega_scope = omega[1]
    if omega_scope not in ['global', 'branchlocal']:
        raise ValueError("Not a valid omega scope in model:\n%s\nValid values are global or branchlocal" % line)
    omega_classes = omega[2]
    if omega_classes != 'one' and not re.search('^gamma\d+$', omega_classes):
        raise ValueError("Not a valid omega rate class in model:\n%s" % line)
    if omega_scope == 'branchlocal' and omega_classes != 'one':
        raise ValueError("branchlocal cannot be combined with multiple rate classes for omega, so the following is an invalid model\n%s" % line)
    rate_classes = entries[3].split('-')
    if rate_classes[0] != 'rates' or len(rate_classes) != 2:
        raise ValueError("Invalid rates specification in model:\n%s" % line)
    rate_classes = rate_classes[1]
    if rate_classes != 'one' and not re.search('^gamma\d+$', rate_classes):
        raise ValueError("Not a valid rate class in model:\n%s" % line)
    print "\nUsing the %s model with a single global kappa (transition/transversion ratio), equilibrium codon frequencies estimated by %s, omega (dN/dS) %s for branches, %s omega classes, and %s rate classes..." % (modeltype, equilfreqs, omega_scope, omega_classes, rate_classes)
    return (modeltype, equilfreqs, omega_scope, omega_classes, rate_classes)



def main():
    """Main body of script."""
    print "\nRunning phyloExpCM_optimizeHyphyTree.py..."
    
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
    hyphycmdfile = phyloExpCM.io.ParseStringValue(d, 'hyphycmdfile')
    hyphyoutfile = phyloExpCM.io.ParseStringValue(d, 'hyphyoutfile')
    hyphytreefile = phyloExpCM.io.ParseStringValue(d, 'hyphytreefile')
    for (fname, fdescription) in [(hyphyoutfile, 'hyphyoutfile'), (hyphytreefile, 'hyphytreefile')]:
        if os.path.isfile(fname):
            print "Removing existing %s %s" % (fdescription, fname)
            os.remove(fname) 
    fastafile = phyloExpCM.io.ParseStringValue(d, 'fastafile')
    if not os.path.isfile(fastafile):
        raise IOError("Cannot find fastafile %s" % fastafile)
    treefile = phyloExpCM.io.ParseStringValue(d, 'treefile')
    if not os.path.isfile(treefile):
        raise IOError("Cannot find treefile %s" % treefile)
    siteslist = phyloExpCM.io.ParseStringValue(d, 'siteslist')
    if not os.path.isfile(siteslist):
        raise IOError("Cannot file siteslist %s" % siteslist)
    model = phyloExpCM.io.ParseStringValue(d, 'model')
    # parse hyphydistancesfile
    hyphydistancesfile = phyloExpCM.io.ParseStringValue(d, 'hyphydistancesfile')
    if hyphydistancesfile.strip().upper() == 'NONE':
        hyphydistancesfile = None
    elif 'branchlocal' in model:
        print "Cannot compute pairwise distances for substitution models with branch-local parameters. So no hyphydistancesfile of %s will be created for model %s" % (hyphydistancesfile, model)
        hyphydistancesfile = None
    else:
        hyphydistancesfile = hyphydistancesfile.split()
        if len(hyphydistancesfile) == 1:
            hyphydistancesfile = hyphydistancesfile[0]
            print "All pairwise distances will be written to %s" % hyphydistancesfile
            if  os.path.isfile(hyphydistancesfile):
                print "Removing existing hyphydistancesfile %s" % hyphydistancesfile
                os.remove(hyphydistancesfile)
        else:
            distancefilename = hyphydistancesfile[0]
            distanceseqs = [x.strip() for x in hyphydistancesfile[1 : ]]
            hyphydistancesfile = (distancefilename, distanceseqs)
            print "Selected pairwise distances will be written to %s" % distancefilename
            print "Distances will be computed to the following sequences:\n%s" % '\n'.join(distanceseqs)
            if  os.path.isfile(distancefilename):
                print "Removing existing hyphydistancesfile %s" % distancefilename
                os.remove(distancefilename)


    # re-map names in _codenames_* files to HYPHY acceptable format
    # code_d maps sequence headers to code names
    codetreefile = '_codenames_%s' % os.path.basename(treefile)
    codefastafile = '_codenames_%s' % os.path.basename(fastafile)
    seqs = mapmuts.sequtils.ReadFASTA(fastafile)
    code_d = dict([(seqs[i][0], 'TipSeq%d_' % (i + 1)) for i in range(len(seqs))])
    code_d_inv = dict([(y, x) for (x, y) in code_d.iteritems()])
    if len(code_d) != len(code_d_inv):
        raise ValueError("Some the sequence names are not unique")
    EncodeNames(code_d, treefile, codetreefile)
    RemoveBranchSupports(codetreefile)
    EncodeNames(code_d, fastafile, codefastafile)
    if isinstance(hyphydistancesfile, tuple):
        (distancefilename, distanceseqs) = hyphydistancesfile
        encoded_distanceseqs = []
        for x in distanceseqs:
            if x in code_d:
                encoded_distanceseqs.append(code_d[x])
            else:
                raise ValueError("hyphydistancesfile specifies that we compute distances to a sequence that is not a header in fastafile. The invalid sequence name is:\n%s" % x)
        codedistancefilename = '_codenames_%s' % os.path.basename(distancefilename)
        codedistancesfile = (codedistancefilename, encoded_distanceseqs)
    elif hyphydistancesfile:
        assert isinstance(hyphydistancesfile, str)
        codedistancesfile = '_codenames_%s' % os.path.basename(hyphydistancesfile)
    else:
        codedistancesfile = None

    # create HYPHY command file
    sites = [line.split()[0] for line in open(siteslist).readlines()[1 : ] if not line.isspace()]
    try:
        sites = [int(site) for site in sites]
    except ValueError:
        raise ValueError("A line in sitesfile does not specify a valid integer")
    if len(sites) != len(dict([(site, True) for site in sites])):
        raise ValueError("Duplicate sites in sitesfile")
    sites.sort()
    print "\nThe analysis will include the following sites:\n%s" % ', '.join([str(site) for site in sites])
    if sites[0] < 1:
        raise ValueError("The minimum site number specified in sites is < 1")
    if model[ : 12] == 'experimental' and len(model.split()) == 2:
        hyphyExpCM = model.split()[1]
        if not os.path.isfile(hyphyExpCM):
            raise IOError("Cannot find the experimental codon model include batch file %s" % hyphyExpCM)
        print "\nUsing the experimentally determined substitution model defined in %s" % hyphyExpCM
        randomizematch = re.compile('^experimental_randomize(?P<seed>\d+)$')
        m = randomizematch.search(model.split()[0])
        if model.split()[0] == 'experimental':
            pass
        elif m:
            seed = int(m.group('seed'))
            print "\nThe sites in the aligment will be randomized with seed %d." % seed
            random.seed(seed)
            seqs = mapmuts.sequtils.ReadFASTA(codefastafile)
            randseqs = phyloExpCM.hyphy.RandomizeCodonAlignment(seqs)
            mapmuts.sequtils.WriteFASTA(randseqs, codefastafile)
        else:
            raise ValueError("Invalid model specification of %s" % model)
        model = ('experimental', hyphyExpCM)
    elif model[ : 5] == 'GY94_':
        model = ParseModelOptions(model)
    elif model[ : 7] == 'KOSI07_':
        model = ParseModelOptions(model)
    else:
        raise ValueError("Unrecognized model of %s" % model)
    print "\nCreating the HYPHY command file %s..." % hyphycmdfile
    phyloExpCM.hyphy.CreateHYPHYCommandFile(hyphycmdfile, hyphyoutfile, codefastafile, codetreefile, codedistancesfile, sites, model)

    # run HYPHY
    assert not os.path.isfile(hyphyoutfile), "hyphyoutfile should have already been deleted"
    (fd, errorfile) = tempfile.mkstemp()
    os.close(fd)
    print "\nNow running HYPHY (using the command %s) beginning at %s..." % (hyphypath, time.asctime())
    try:
        subprocess.call(hyphypath.split() + [hyphycmdfile], stderr=open(errorfile, 'w')) # run splitting hyphypath in case set to something like HYPHY CPU=2
    except OSError:
        raise
        if os.path.isfile(errorfile):
            os.remove(errorfile)
        raise ValueError("Failed to run HYPHY with the following command: %s\nIs HYPHY actually installed with that name in an accessible path?" % hyphypath)
    errors = open(errorfile).read()
    os.remove(errorfile)
    if errors and errors.strip() != 'Check messages.log details of this run.':
        sys.stderr.write("HYPHY generated the following errors:\n%s" % errors)
    if not os.path.isfile(hyphyoutfile):
        raise ValueError("HYPHY failed to generate the expected output file %s" % hyphyoutfile)
    print "HYPHY completed execution at %s, and created the output file %s." % (time.asctime(), hyphyoutfile)

    # construct the HYPHY tree
    print "\nNow writing the HYPHY optimized tree to %s..." % hyphytreefile
    newicktree = phyloExpCM.hyphy.ExtractTree(hyphyoutfile)
    open(hyphytreefile, 'w').write(newicktree)
    EncodeNames(code_d_inv, hyphytreefile, hyphytreefile)

    # construct hyphydistancesfile from codedistancesfile
    if hyphydistancesfile:
        if isinstance(codedistancesfile, str):
            assert isinstance(hyphydistancesfile, str)
            codedistancefilename = codedistancesfile
            hyphydistancefilename = hyphydistancesfile
        else:
            assert isinstance(codedistancesfile, tuple) and len(codedistancesfile) == 2, "Invalid codedistancesfile: %s" % str(codedistancesfile)
            assert isinstance(hyphydistancesfile, tuple), "Invalid hyphydistancesfile: %s" % str(hyphydistancesfile)
            codedistancefilename = codedistancesfile[0]
            hyphydistancefilename = hyphydistancesfile[0]
        if not os.path.isfile(codedistancefilename):
            raise IOError('HYPHY failed to generate the expected distances file %s' % codedistancefilename)
        EncodeNames(code_d_inv, codedistancefilename, hyphydistancefilename, allowmultiple=True)
        os.remove(codedistancefilename)

    # delete codefastafile and codetreefile
    for f in [codefastafile, codetreefile]:
        if os.path.isfile(f):
            os.remove(f)

    # script done
    print "\nScript complete."
    sys.stdout.flush()
    time.sleep(0.1)



main() # run the script
