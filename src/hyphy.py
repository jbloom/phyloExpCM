"""Module for running ``HYPHY`` analyses.

This module is designed to run the ``HYPHY`` program. It has been 
tested with ``HYPHY`` version 2.112*.

Requires ``mapmuts`` for some basic sequence operations.

Written by Jesse Bloom, 2013.


Functions defined in this module
--------------------------------------

* *RandomizeCodonAlignment* : randomly rearranges sites in an alignment.

* *ExtractTree* : extract tree from ``HYPHY`` output.

* *ExtractLikelihoodAndNParameters* : extracts likelihood and number of parameters from ``HYPHY`` output.

* *CreateHYPHYCommandFile* : creates command file for ``HYPHY`` batch mode.

"""


import re
import os
import random
import mapmuts.sequtils
import phyloExpCM.packagedata



def RandomizeCodonAlignment(seqs):
    """Randomly rearranges sites in a codon alignment.

    *seqs* is a set of aligned coding nucleotide sequences, aligned
    at the codon level (so consecutive triplets are aligned). The
    list consists of *(header, sequence)* 2-tuples. Because the sequences
    are aligned, the *sequence* elements of the 2-tuples should
    all be the same length and should be multiples of three.

    For all sequences, codon sites are randomly interchanged, with the same
    interchange applied to all sequences. The interchange is done
    using the Python *random* module, so for reproducible output seed
    this with *random.seed*.

    The return variable *randomizedseqs* is a copy of *seqs*
    in which the codons have been interchanged randomly.
    """
    length = len(seqs[0][1])
    if length % 3:
        raise ValueError("Alignment lengths are not multiples of three.")
    length = length / 3
    sites = [i for i in range(length)]
    rsites = [i for i in range(length)]
    random.shuffle(rsites)
    rseqs = []
    for (head, seq) in seqs:
        rseq = ['???' for i in range(length)]
        for (i, ri) in zip(sites, rsites):
            rseq[ri] = seq[3 * i : 3 * i + 3]
        rseq = ''.join(rseq)
        assert '?' not in rseq
        rseqs.append((head, rseq))
    return rseqs


def ExtractLikelihoodAndNParameters(hyphyoutfile):
    """Extracts log likelihood and number of parameters from ``HYPHY`` output.

    *hyphyoutfile* is a string giving the name of a ``HYPHY`` output file
    of the type produced by ``phyloExpCM_optimizeHyphyTree.py``.
    It should contain the following lines (along with probably many others)::

        Log likelihood: -671.317
        independent parameters (includes branch lengths): 47
        shared parameters: 0
        number of branch lengths: 45
        number of tip nodes: 24
        number of internal branches: 21

    This function returns the 2-tuple *(loglikelihood, nparameters)*. The entries
    are:

        * *loglikelihood* is the log likelihood, so -671.317 in the example
          output shown above.

        * *nparameters* is the number of free parameters **not** including the
          branch lengths. So in the example above, *nparameters* is
          47 - 45, since the number of branch lengths is subtracted from the
          total number of independent parameters including branch lengths.
          Note that this parameter count will **not** include parameters
          that are empirically estimated from the data (such as the 9
          frequencies in the Goldman-Yang 1994 CF3x4 model).
    """
    llmatch = re.compile ('Log likelihood: (?P<value>\-\d+(\.\d+(e|E\d+){0,1}){0,1})')
    parametersmatch = re.compile('independent parameters \(includes branch lengths\): (?P<value>\d+)')
    branchlengthsmatch = re.compile('number of branch lengths: (?P<value>\d+)')
    d = {'ll':None, 'ntotparameters':None, 'nbranchlengths':None}
    match_d = {'ll':llmatch, 'ntotparameters':parametersmatch, 'nbranchlengths':branchlengthsmatch}
    for line in open(hyphyoutfile):
        for key in d.iterkeys():
            if match_d[key].search(line):
                if d[key] != None:
                    raise ValueError("Found duplicate entry on line:\n%s" % line)
                else:
                    d[key] = match_d[key].search(line).group('value')
    if d['ll']:
        loglikelihood = float(d['ll'])
    else:
        raise ValueError("Failed to parse log likelihood")
    if d['ntotparameters']:
        ntotparameters = int(d['ntotparameters'])
    else:
        raise ValueError("Failed to parse total number of independent parameters")
    if d['nbranchlengths']:
        nbranchlengths = int(d['nbranchlengths'])
    else:
        raise ValueError("Failed to parse total number of branch lengths")
    return (loglikelihood, ntotparameters - nbranchlengths)



def ExtractTree(hyphyoutfile, treename='tree'):
    """Extracts a the Newick tree string from ``HYPHY`` output.

    This function is designed to extract the newick tree string from a
    ``HYPHY`` output file that includes the tree topology and the branch
    lengths. This would be the type of file created by printing the
    optimized *likelihood* object in ``HYPHY`` with 
    *LIKELIHOOD_FUNCTION_OUTPUT = 5*. More specifically, this is the 
    tree that would be present in the *outputfile* created by
    *CreateHYPHYCommandFile*.

    *hyphyoutfile* should specify the name of the ``HYPHY`` output file
    that contains the tree.

    *treename* is a string giving the name used to specify the tree
    in *hyphyoutfile*. By default, this is the string 'tree'.

    This function returns a string giving the Newick string of the tree
    with the correct branch lengths.

    If the same substitution model is used for all sites, it would be
    possible simply to extract the Newick string line from *hyphyoutfile*.
    However, if different sites have different models, then the file will
    contain lots of different tree lines with different branches. To handle
    this situation, this function looks for the newick tree with the name 
    specified by *treename*, a line beginning::

        Tree tree=(((((((((((((((TipSeq13_:0.012291,TipSeq24_:0.0122797)Node14

    and gets that tree topology. However, it does not simply extract that line,
    because the branch lengths are only for the substitution model applied
    to the sites in that tree. It then also looks for entries specifying
    branch lengths, such as::

        tree.Node6.t=5.000000025e-09;

    and uses those to fill in the branch lengths. Once this is done, it
    returns the resulting tree string.
    """
    lines = open(hyphyoutfile).readlines()
    treeline = [line for line in lines if line.split('=')[0] == 'Tree %s' % treename]
    if len(treeline) != 1:
        raise ValueError("Failed to find exactly one tree line")
    treeline = treeline[0]
    blmatch = re.compile('^%s\.(?P<node>\w+)\.t\=(?P<bl>\d+(\.\d+(e(\-){0,1}\d+){0,1}){0,1});\n$' % treename)
    bls = {}
    for line in lines:
        m = blmatch.search(line)
        if m:
            (key, value) = (m.group('node'), float(m.group('bl')))
            if key in bls:
                raise ValueError("Duplicate branch length for node %s" % key)
            bls[key] = value
    if len(bls) != treeline.count(':'):
        raise ValueError("Failed to match a branch length for all colons in treeline: len(bls) = %d, treeline.count(':') = %d" % (len(bls), treeline.count(':')))
    for (key, value) in bls.iteritems():
        nodematch = re.compile("(?P<toreplace>%s\:\d+(\.\d+(e(\-){0,1}\d+){0,1}){0,1})[\)\,\;]" % key)
        ms = nodematch.findall(treeline)
        if len(ms) != 1:
            raise ValueError("Found %d occurrences of %s" % (len(ms), key))
        ms = nodematch.search(treeline)
        treeline = treeline.replace(ms.group('toreplace'), "%s:%g" % (key, value))
    return treeline



def CreateHYPHYCommandFile(cmdfile, outputfile, fastafile, newickfile, distancesfile, sites, model):
    """Creates a ``HYPHY`` command file.

    The created command file *cmdfile* specifies that we run ``HYPHY``
    on a pre-existing tree topology given by *newickfile*, where
    *fastafile* gives the sequences of the tip nodes of the tree.
    The analysis is run using only the sites listed in *sites*, and
    uses the substitution model specified by *model*.
    The tree topology is not changed, but the branch lengths and
    any parameters in the substitution model are optimized.

    This function currently only works for codon sequences and
    codon models.

    CALLING VARIABLES:

    * *cmdfile* is the name of the ``HYPHY`` command file that we create.
      It is overwritten if it already exists. This file can be used
      to run ``HYPHY`` using the command::

        hyphypath cmdfile

      where *hyphypath* gives a path to a ``HYPHY`` executable.

    * *outputfile* is the name of the file that the ``HYPHY`` command
      file *cmdfile* specifies that the final output should be written.

    * *fastafile* is the name of a FASTA file that gives the aligned
      coding sequences (the lengths must be multiples of three,
      and they must be aligned at the codon level). Stop codons
      should be removed or else not included in the listing specified in 
      *sites*. The sequence headers must match the tip node names given
      in *newickfile*. The sequence names must conform to ``HYPHY`` 
      naming conventions, which means only alphanumeric characters
      and underscores.

    * *newickfile* gives the tree in Newick format, with tip names
      matching the sequence names in *fastafile*.

    * *distancesfile* specifies that we create a file that contains the
      pairwise distances between all sequences computed using the model
      parameters obtained after maximizing the likelihood over the whole tree.
      These distances represent the branch length of a two-taxa tree for
      each pair of sequences after fixing the substitution model parameters
      to the maximum-likelihood values for the whole tree. If the substitution
      model contains local parameters (such as an *omega_scope* of "branchlocal"
      for *model*) then *distancesfile* can **not** be used, and must
      be set to *None*. That is because branch local model parameters
      cannot be fixed for the pairwise distance estimations. So if you
      have local parameters and set *distancesfile* to something other than
      *None*, an exception will be raised. 

      The valid values for *distancesfile* are as follows:

        - *None* specifies that no pairwise distances are computed.

        - A string giving the name of the file to create. This file will
          contain the pairwise distances for all possible pairs of
          sequences.

        - A 2-tuple *(filename, distanceseqs)* specifies that we create
          a pairwise-distance file of the name *filename*, but that it only
          contains pairwise distances to all of the sequence listed in
          *distanceseqs*. Then *distanceseqs* should be a non-empty list
          containing one or more sequence names that match the headers in
          *fastafile*. The created *filename* contains the pairwise
          distances from all other sequences to each sequence listed
          in *distanceseqs*.
      
      If a pairwise distances file is create, it contain a line for each
      pair of sequences with three tab-delimited entries: 

        - the first sequence name

        - the second sequence name

        - the pairwise distance between these two sequences


    * *sites* is a listing of all codon sites that we include in the analysis.
      These sites are numbered 1, 2, ... For example, a value
      of *sites = [2, 3, 4]* specifies that we only include
      codon positions 2, 3, and 4. This list must be sorted from
      smallest to largest or an Exception will be raised. Note that
      these are the codon numbers, not the nucleotide numbers.

    * *model* specifies the codon substitution model used. Valid values are:

        - *('experimental', hyphyExpCMs)* : if *model* is set to this
          2-tuple, we use an experimentally determined codon model
          where a different substitution model is applied to each
          site. The first element of the 2-tuple is just the string
          'experimental', while the second element is the name of
          a ``HYPHY`` include batch file specifying the model for
          each site in sites as *model1*, *model2*, etc where the
          numerical suffixes are in 1, 2, ... numbering and
          are present for each site in *sites*. Typically,
          this file would be the file that you create for the argument
          *hyphyExpCMs* using the script ``phyloExpCM_buildHyphyExpCM.py``.

        - *('GY94', equilfreqs, omega_scope, omega_classes, rate_classes)* :
          is a tuple that specifies that we use the Goldman and Yang 1994 
          model. The transition-transversion ratio *kappa* is treated as a
          single global parameter for all branches and sites, and is estimated
          by maximum likelihood. The remaining parameters are as follows:

            - *equilfreqs* specifies how we estimate the equilibrium codon
              frequencies. The supported values are the string *CF3x4* and
              *F3x4*. You should use *CF3x4* unless you have a good reason
              to use *F3x4*. The *CF3x4* is the corrected method for estimating
              the codon frequencies as the products of the individual nucleotides
              at the three positions.

            - *omega_scope* specifies the scope of the dN/dS ratio (*omega*).
              Valid values are the string *global* if the same ratio (or
              set of ratios) are applied to all branches, or *branchlocal*
              if different ratios (or sets of ratios) are applied to each
              branch. The latter approach will greatly increase the number
              of parameters.

            - *omega_classes* specifies how many *omega* values there are
              either for the whole tree (if *omega_scope* is *global*)
              or for each branch (if *omega_scope* is *branchlocal*). This
              variable can be set to the string *one* if there is just one
              *omega* value, or to *gamma6* if the values are drawn from
              six equally probable gamma distributed classes (you can also
              use *gamma2*, *gamma3*, etc). The shape of the gamma distribution
              is estimated by maximum likelihood.

            - *rate_classes* specifies how many rate classes there are for
              the synonymous substitution rate. Valid values are the string
              *one* if there is only one rate class, or *gamma6* if there
              are six rate classes drawn from six equally probably 
              gamma-distributed classes (you can also use *gamma2*, *gamma3*,
              etc). The shape of the gamma distribution is
              estimated by maximum likelihood.

        - *('KOSI07', equilfreqs, omega_scope, omega_classes, rate_classes)* :
          is a tuple that specifies that we use the Kosiol et al, 2007 model.
          This is the variant termed *ECM+F+omega+1kappa(tv)* in that
          original paper. The other variables in the tuple have the same
          meaining as for the *GY94* model except that now the only allowable
          option for *equilfreqs* is *F*.
    """
    seqs = mapmuts.sequtils.ReadFASTA(fastafile)
    if len(seqs) != len(dict([(head, True) for (head, seq) in seqs])):
        raise ValueError("The sequence names in %s are not all identical." % fastafile)
    validnamematch = re.compile('^\w+$')
    for (head, seq) in seqs:
        if not validnamematch.search(head):
            raise ValueError("Invalid sequence name (can only contain letters, numbers, underscore):\n%s" % head)
    if not sites:
        raise ValueError("No sites specified")
    if len(sites) != len(dict([(site, True) for site in sites])):
        raise ValueError("sites contains duplicate entries.")
    if min(sites) < 1:
        raise ValueError("minimum value in sites is less than one")
    sortedsites = [site for site in sites]
    sortedsites.sort()
    skippedsites = max(sites) - len(sites)
    if sortedsites != sites:
        raise ValueError("sites not sorted from smallest to largest")
    ntsites = ','.join(["%d,%d,%d" % (3 * site - 3, 3 * site - 2, 3 * site - 1) for site in sites]) # listing of sites to include for HYPHY DataSetFilter
    includepath = phyloExpCM.packagedata.Path() # path where HYPHY include batch files were installed
    includefiles = ['NTsCodonsAAs.ibf'] # files to always include
    includefiles = ["%s/%s" % (includepath, includefile) for includefile in includefiles] # make paths to include files absolute
    for f in includefiles:
        if not os.path.isfile(f):
            raise IOError("Cannot find HYPHY include file %s. Did you install the package data correctly?" % f)
    # commands to read in codon data
    cmds = [ 
        'INTEGRATION_PRECISION_FACTOR = 5.0e-6;',
        'END_OF_FILE = 0;',
        'LIKELIHOOD_FUNCTION_OUTPUT = 5;',
        'ACCEPT_BRANCH_LENGTHS = 1;',
        '\n'.join(['#include "%s";' % f for f in includefiles]),
        'fprintf(stdout, "Running HYPHY script %s...\\n");' % cmdfile,
        'DataSet data = ReadDataFile("%s");' % fastafile,
        'assert(data.sites % 3 == 0, "Sequence lengths not multiples of 3");',
        'totalcodons = data.sites $ 3;',
        'fprintf(stdout, "Read from %s a set of ", data.species, " seqeunces consisting of ", data.sites, " nucleotides corresponding to ", totalcodons, " codons each.\\n");' % fastafile,
        'fprintf(stdout, "The analysis will include the following %d codon positions (sequential numbering starting with 1):\\n%s\\n");' % (len(sites), ', '.join([str(site) for site in sites])),
        'assert(totalcodons >= %d, "Largest included site exceeds sequence length");' % max(sites),
        'DataSetFilter codonfilter = CreateFilter(data, 3, "%s", "", "TAA,TAG,TGA");' % ntsites,
        'assert(data.species == codonfilter.species, "species number mismatch");',
        'assert(codonfilter.sites == %d, "Codon filtered data does not contain the right number of sites");' % len(sites),
        'fprintf(stdout, "Created a codon filter of ", codonfilter.sites, " sites.\\n");',
        'assert(totalcodons - (totalcodons - %d) - %d == codonfilter.sites, "Codon filtered data is not the expected length. Do sequences contain stop codons?");' % (max(sites), skippedsites),
        'CheckCodonFilter("codonfilter");',
        'fprintf(stdout, "Reading tree string from %s.\\n");' % newickfile,
        'fscanf("%s", String, treestring);' % newickfile,
        ]
    # commands to set up substitution model, apply to tree, define likelihood
    if isinstance(model, tuple) and len(model) == 2 and model[0] == 'experimental':
        hyphyExpCMs = model[1]
        cmds += [
            'fprintf(stdout, "Using the experimentally determined substitution models in %s...\\n");' % hyphyExpCMs,
            '#include "%s";' % hyphyExpCMs,
            'fprintf(stdout, "Now constructing the likelihood function...\\n");',
            ]
        firstsite = True
        for site in sites:
            cmds += [
                'DataSetFilter codonfilter%d = CreateFilter(data, 3, "%d-%d", "", "TAA,TAG,TGA");' % (site, 3 * site - 3, 3 * site - 1),
                'assert(data.species == codonfilter%d.species, "species number mismatch");' % site,
                'assert(1 == codonfilter%d.sites, "codon filter does not contain exactly one site");' % site,
                'CheckCodonFilter("codonfilter%d");' % site,
                'UseModel(model%d);' % site,
                ]
            if firstsite:
                cmds += [
                    'ExecuteCommands("Tree tree = treestring;");',
                    'assert(codonfilter.species == TipCount(tree), "Number of species and number of tips differ");',
                ]
                firstsite = False
            else:
                cmds += [
                    'ExecuteCommands("Tree tree%d = treestring;");' % site,
                    'ReplicateConstraint("this1.?.t := this2.?.t", tree%d, tree);' % site,
                    ]
        cmds.append('LikelihoodFunction likelihood = (codonfilter%d, tree, %s);' % (sites[0], ', '.join(['codonfilter%d, tree%d' % (site, site) for site in sites[1 : ]]))),
    elif isinstance(model, tuple) and len(model) == 5 and model[0] in ['GY94', 'KOSI07']:
        (modeltype, equilfreqs, omega_scope, omega_classes, rate_classes) = model
        if omega_scope == 'branchlocal' and omega_classes != 'one':
            raise ValueError("Cannot combine omega_scope branchlocal and omega_classes != one")
        if omega_scope == 'branchlocal' and distancesfile:
            raise ValueError("Cannot combine omega_scope branchlocal and distancesfile")
        gammamatch = re.compile('^gamma(?P<n>\d+)$')
        if rate_classes == 'one':
            rate_classes = 1
        else:
            rate_classes = gammamatch.search(rate_classes)
            if not rate_classes:
                raise ValueError("Invalid value of rate_classes:\n%s" % str(model))
            rate_classes = int(rate_classes.group('n'))
            if rate_classes < 1:
                raise ValueError("rate_classes must be > 1 for gamma distribution.")
        if omega_classes == 'one':
            omega_classes = 1
        else:
            omega_classes = gammamatch.search(omega_classes)
            if not omega_classes:
                raise ValueError("Invalid value of omega_classes:\n%s" % str(model))
            omega_classes = int(omega_classes.group('n'))
            if omega_classes < 1:
                raise ValueError("omega_classes must be > 1 for gamma distribution.")
        if modeltype == 'GY94':
            cmds += [
                'fprintf(stdout, "Using the Goldman Yang 1994 (GY94) codon model...\\n");',
                '#include "%s/CF3x4.ibf";' % includepath,
                '#include "%s/GY94.ibf";' % includepath,
                'CreateGY94Model("%s", "global", "%s", %d, %d, 1);' % (equilfreqs, omega_scope, omega_classes, rate_classes),
                'UseModel(model);',
                'ExecuteCommands("Tree tree = treestring;")',
                'assert(codonfilter.species == TipCount(tree), "Number of species and number of tips differ");',
                'LikelihoodFunction likelihood = (codonfilter, tree);',
                ]
        elif modeltype == 'KOSI07':
            cmds += [
                'fprintf(stdout, "Using the Kosiol et al 2007 (KOSI07) codon model...\\n");',
                '#include "%s/KOSI07.ibf";' % includepath,
                'CreateKOSI07Model("%s", "ktv", "%s", %d, %d, 1, "%s/KOSI07_exchangeabilities.ibf");' % (equilfreqs, omega_scope, omega_classes, rate_classes, includepath),
                'UseModel(model);',
                'ExecuteCommands("Tree tree = treestring;")',
                'assert(codonfilter.species == TipCount(tree), "Number of species and number of tips differ");',
                'LikelihoodFunction likelihood = (codonfilter, tree);',
                ]
        else:
            raise ValueError("Invalid model type of %s in %s" % (modeltype, str(model)))
    else:
        raise ValueError("Unrecognized value of model: %s" % str(model))
    # now maximize the likelihood
    cmds += [
        'fprintf(stdout, "\\nNow optimizing the likelihood function...\\n");',
        'Optimize(mlestimates, likelihood)',
        'fprintf(stdout, "Completed likelihood optimization. Optimized ", mlestimates[1][1], " indpendent parameters and ", mlestimates[1][2], " shared parameters to obtain a log likelihood of ", mlestimates[1][0], ".\\n");',
        'fprintf(stdout, "Writing the results to %s.\\n");' % outputfile,
        'fprintf("%s", "Log likelihood: ", mlestimates[1][0], "\\nindependent parameters (includes branch lengths): ", mlestimates[1][1], "\\nshared parameters: ", mlestimates[1][2], "\\nnumber of branch lengths: ", TipCount(tree) + BranchCount(tree), "\\nnumber of tip nodes: ", TipCount(tree), "\\nnumber of internal branches: ", BranchCount(tree), "\\n",likelihood);' % outputfile,
        ]
    # potentially add command to compute pairwise distances
    if distancesfile:
        if isinstance(distancesfile, str):
            distanceseqs = 'ALL'
        elif isinstance(distancesfile, tuple) and len(distancesfile) == 2:
            (distancesfile, distanceseqs) = distancesfile
            if not (isinstance(distancesfile, str) and distancesfile):
                raise ValueError("distancesfile 2-tuple must specify non-empty string as first entry")
            if not (isinstance(distanceseqs, list) and len(distanceseqs) >= 1):
                raise ValueError("distancesfile 2-tuple must specify non-empty list as second entry")
            heads = dict(seqs)
            for x in distanceseqs:
                if x not in heads:
                    raise ValueError("distancesfile 2-tuple's second entry specifies a sequence header of %s, which is not in fastafile" % x)
        else:
            raise ValueError("distancesfile must either be a 2-tuple or a string, but got a value of: %s" % str(distancesfile))
        cmds += [
            'fprintf(stdout, "\\nNow computing pairwise distances.\\n");',
            'fprintf(stdout, "\\nFirst fixing all global variables to the maximum-likelihood values estimated on the entire tree.\\n");',
            'GetString(associativearray, likelihood, -1);',
            'globalindependentvariables = associativearray["Global Independent"];',
            'for (ivariable=0; ivariable<Columns(globalindependentvariables); ivariable=ivariable+1) {',
            '  variable = globalindependentvariables[ivariable];',
            '  cmdstring = variable + " := " + Format(variable, 0, 30) + ";";',
            '  fprintf(stdout, "\\nFixing variable as follows: ", cmdstring, "\\n");',
            '  ExecuteCommands(cmdstring);',
            '}',
            'pairwisedistancesfile = "%s";' % distancesfile,
            ]
        if distanceseqs == 'ALL':
            cmds += [
                'fprintf(stdout, "\\nNow computing all pairwise distances and writing to ", pairwisedistancesfile, "...\\n");',
                'for (seq1=0; seq1<codonfilter.species; seq1=seq1+1) {',
                '  for (seq2=seq1+1; seq2<codonfilter.species; seq2=seq2+1) {',
                '    pairstring = Format(seq1, 0, 0) + "," + Format(seq2, 0, 0);',
                '    DataSetFilter pairfilter = CreateFilter(codonfilter, 3, "", pairstring, "TAA,TAG,TGA");',
                '    GetString(seq1name, pairfilter, 0);',
                '    GetString(seq2name, pairfilter, 1);',
                '    pairtreestring = "(" + seq1name + "," + seq2name + ")";',
                ]
        else:
            assert isinstance(distanceseqs, list) and len(distanceseqs) >= 1
            cmds += [
                'fprintf(stdout, "\\nNow computing pairwise distances from all sequences to the %d sequences %s and writing to ", pairwisedistancesfile, "...\\n");' % (len(distanceseqs), ', '.join(distanceseqs)),
                'for (seq1=0; seq1<codonfilter.species; seq1=seq1+1) {',
                '  for (seq2=0; seq2<codonfilter.species; seq2=seq2+1) {',
                '    if (seq1 == seq2) {',
                '      continue;',
                '    }',
                '    pairstring = Format(seq1, 0, 0) + "," + Format(seq2, 0, 0);',
                '    DataSetFilter pairfilter = CreateFilter(codonfilter, 3, "", pairstring, "TAA,TAG,TGA");',
                '    GetString(seq1name, pairfilter, 0);',
                '    if (%s) {' % (' && '.join(['seq1name != "%s"' % x for x in distanceseqs])),
                '      continue;',
                '    }',
                '    GetString(seq2name, pairfilter, 1);',
                '    pairtreestring = "(" + seq1name + "," + seq2name + ")";',
                ]
        if isinstance(model, tuple) and len(model) == 2 and model[0] == 'experimental':
            firstsite = True
            for site in sites:
                cmds += [
                    '    DataSetFilter pairfilter%d = CreateFilter(data, 3, "%d-%d", pairstring, "TAA,TAG,TGA");' % (site, 3 * site - 3, 3 * site - 1),
                    '    assert(2 == pairfilter%d.species, "pairfilter does not have exactly two species");' % site,
                    '    assert(1 == pairfilter%d.sites, "pairfilter does not contain exactly one site");' % site,
                    '    CheckCodonFilter("pairfilter%d");' % site,
                    '    UseModel(model%d);' % site,
                    ]
                if firstsite:
                    cmds += [
                        '    Tree pairtree = pairtreestring;',
                ]
                    firstsite = False
                else:
                    cmds += [
                        '    Tree pairtree%d = pairtreestring;' % site,
                        '    ReplicateConstraint("this1.?.t := this2.?.t", pairtree%d, pairtree);' % site,
                        ]
            cmds.append('    LikelihoodFunction pairlikelihood = (pairfilter%d, pairtree, %s);' % (sites[0], ', '.join(['pairfilter%d, pairtree%d' % (site, site) for site in sites[1 : ]]))),
        elif isinstance(model, tuple) and len(model) == 5 and model[0] in ['GY94', 'KOSI07']:
            cmds += [
                '    assert(pairfilter.sites == codonfilter.sites, "pairfilter and codonfilter have different numbers of sites");',
                '    assert(pairfilter.species == 2, "pairfilter does not have exactly two species");',
                '    UseModel(model);',
                '    Tree pairtree = pairtreestring;',
                '    LikelihoodFunction pairlikelihood = (pairfilter, pairtree);',
                ]
        else:
            raise ValueError("Unrecognized value of model: %s" % str(model))
        cmds += [
            '    fprintf(stdout, "Computing distance between sequences ", seq1 + 1, " and ", seq2 + 1, "... ");',
            '    Optimize(pairestimates, pairlikelihood);',
            '    assert(pairestimates[1][1] == 1, "Found more than one independent variable that was optimized. You are probably trying to compute pairwise distances using a model that has local parameters (such as branchlocal omega values). This is NOT allowed. You can only compute pairwise distances on models with all global parameters.");',
            '    assert(pairestimates[1][2] == 0, "At least one global parameter was optimized when computing pairwise distances. Pairwise distances should only be computed after fixing the global model parameters.");',
            '    pairdistance = pairestimates[0][0];',
            '    assert(pairdistance >= 0, "Pairwise distance between sequences is negative. Something is wrong!");',
            '    fprintf(stdout, "computed distance of ", pairdistance, ".\n");',
            '    fprintf(pairwisedistancesfile, seq1name, "\t", seq2name, "\t", pairdistance, "\n");',
            '  }',
            '}',
            ]
    # last command
    cmds += [
        'fprintf(stdout, "Completed HYPHY script %s.\\n");' % cmdfile,
        ]
    # write the file
    open(cmdfile, 'w').write('\n'.join(cmds))


if __name__ == '__main__':
    import doctest
    doctest.testmod()

