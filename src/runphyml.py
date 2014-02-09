"""Python module for running ``phyml`` and ``codonphyml``.

Requires the ``mapmuts`` package.

Written by Jesse Bloom, 2013.


Functions defined in this module
-----------------------------------

* *WritePhylipSequenceFile* : writes sequences in PHYLIP format

* *RunCodonphyml* : runs ``codonphyml``

"""


import os
import re
import tempfile
import subprocess
import cStringIO
import mapmuts.sequtils


def WritePhylipSequenceFile(sequences, f, add_mm_line=False):
    """Writes a sequence file in sequential relaxed Phylip format.

    *sequences* is a list of 2-tuples *(head, seq)* where *head* is the
        sequence header and *seq* is the sequence (both strings). The sequences
        should be aligned, and so should all be the same length.
        
    *f* is a writeable file-like object to which the Phylip formatted
        sequences are written.

    *add_mm_line* specifies that we add the final line necessary for using
    ``codonPhyML`` with multi-models. If this option has its default value
    of *None*, then no such line is added. Otherwise, *add_mm_line* should
    be a string of length equal to the length of the sequences in *sequences*.
    The characters in this line are typically taken to specify the model 
    assigned to each site, and are written with the prefix ``#=GR mods``, 
    as in::

        #=GR mods AAAAABBAAAAAAAAAA

    Phylip file format is defined here: http://www.phylo.org/tools/phylip.html

    This function modifies that Phylip format as follows:

        * the sequence name (*head*) can be up to 100 characters long

        * there is always a space between the sequence name and the sequence

        * the following symbols are not allowed, blanks and any of :(),"

    If any of the headers violate these restrictions, or are not unique,
    an exception is raised.

    >>> sequences = [('head1', 'ATGACT'), ('head2', 'ATGCTA')]
    >>> f = cStringIO.StringIO()
    >>> WritePhylipSequenceFile(sequences, f)
    >>> f.seek(0)
    >>> print f.read()
    2 6
    head1 ATGACT
    head2 ATGCTA
    <BLANKLINE>

    >>> sequences = [('head1', 'ATGACT'), ('head1', 'ATGCTA')]
    >>> f = cStringIO.StringIO()
    >>> WritePhylipSequenceFile(sequences, f)
    Traceback (most recent call last):
        ...
    ValueError: Duplicate name of head1

    """
    maxnamelength = 100
    notallowed = ' :(),"'
    names_used = {}
    length = len(sequences[0][1])
    f.write("%d %d\n" % (len(sequences), length))
    for (name, seq) in sequences:
        if length != len(seq):
            raise ValueError("Sequences not all of the same length")
        if len(name) >= maxnamelength:
            raise ValueError("Name is too long:\n%s" % name)
        for x in notallowed:
            if x in name:
                raise ValueError("Name contains dis-allowed character %s" % x)
        if name in names_used:
            raise ValueError("Duplicate name of %s" % name)
        names_used[name] = True
        f.write("%s %s\n" % (name, seq))
    if add_mm_line:
        if not (isinstance(add_mm_line, str) and len(add_mm_line) == length):
            raise ValueError("add_mm_line is not a string of the same length as the sequences: %d and %d" % (len(add_mm_line), length))
        f.write("#=GR mods %s\n" % add_mm_line)


def RunCodonphyml(seqfile, path, seed, outprefix, model, add_mm_line=None):
    """Runs ``codonphyml``, returns log likelihood of tree.

    CALLING VARIABLES:

    * *seqfile* is a FASTA file with the sequences that we want to
      analyze. Each must have a unique header of less than 100 characters
      that does not include a blank or the characters :(),"

    * *path* is the path to the ``codonphyml`` executable. 

    * *seed* is the integer random number seed.

    * *outprefix* is the prefix for the files created by this function.
      These files have that prefix with the following suffixes. 
      Any existing files with these names are overwritten:

        - ``_config.drw`` is the darwin format input file used to run
          ``codonphyml``.

        - ``_output.txt`` is the output of ``codonphyml``.

        - ``_tree.newick`` is the inferred tree in Newick format. Branch
          supports are shown, calculated using the SH-aLRT approach
          (http://sysbio.oxfordjournals.org/content/60/5/685.long)

        - ``_stats.txt`` contains the statistics from the inference,
          including all parameter values.

        - ``_loglikelihood.txt`` contains a single number representing
          the log likelihood.

    * *model* specifies the codon substitution model used. Valid values are:
   
        - *GY94_CF3x4* : The Goldman-Yang codon model with 
          the equilibrium codon frequencies estimated empirically using 
          the CF3x4 method, with kappa (the transition / transversion ratio) 
          estimated by maximum likelihood, and with omega (the dN/dS ratio) 
          estimated as a single parameter by maximum likelihood

        - *GY94_CF3x4_omega-gamma4* : The Goldman-Yang codon model with 
          the equilibrium codon frequencies estimated empirically using 
          the CF3x4 method, with kappa (the transition / transversion ratio) 
          estimated by maximum likelihood, and with omega (the dN/dS ratio) 
          drawn from four gamma-distributed categories with the distribution 
          shape parameter (alpha) estimated by maximum likelihood.

        - *KOSI07_F_omega-gamma4 : The Kosiol et al, 2007 empirical codon model 
          with the equilibrium codon frequencies estimated empirically using 
          the *F* method, with *kappa(tv)* (the relative decrease in transversions 
          versus transitions) estimated by maximum likelihood, and with *omega*
          (the elevation in nonsynonymous over synonymous) drawn from four 
          gamma-distributed categories with the distribution shape parameter 
          and mean estimated by maximum likelihood. This is based on 
          the *ECM+F+omega+1kappa(tv)* model described by Kosiol et al, 2007.

    * *add_mm_line* is *None* by default. If it is set to another
      value, it should be set to a string specifying a multi-model
      line in the format described in the documentation for the 
      parameter of the same name in *WritePhylipSequenceFile*.

    RETURN VARIABLE:

    This function returns a number giving the calculated log likelihood.

    It also creates the files specified by *outprefix*.

    """
    # define file names
    seqs = mapmuts.sequtils.ReadFASTA(seqfile)
    phylipfile = tempfile.mkstemp(suffix='phylip')[1]
    WritePhylipSequenceFile(seqs, open(phylipfile, 'w'), add_mm_line=add_mm_line)
    statsfile = "%s_codonphyml_stats.txt" % phylipfile
    treefile = "%s_codonphyml_tree.txt" % phylipfile
    outfile = "%s_output.txt" % outprefix
    outtreefile = "%s_tree.newick" % outprefix
    outstatsfile = "%s_stats.txt" % outprefix
    llfile = "%s_loglikelihood.txt" % outprefix
    configfile = "%s_config.drw" % outprefix
    # create darwin format config file
    configlines = [
            "cpconfig := table();", # options in table called cpconfig
            "cpconfig['inputfile'] := '%s';" % phylipfile, # sequences
            "cpconfig['oformat'] := 'txt';", # output if txt format
            "cpconfig['sequential'] := true;", # sequences sequential phylip format
            "cpconfig['multiple'] := 1;", # number of datasets to analyze
            "cpconfig['bootstrap'] := -4;", # SH-aLRT supports
            "cpconfig['search'] := 'NNI';", # tree topology search
            "cpconfig['optimize'] := 'tlr';", # optimize topology (t), branch length (l), and rates (r)
            "cpconfig['quiet'] := true;", # no interactive questions asked
            "cpconfig['modrates'] := true;", # not sure what this option does
            "cpconfig['r_seed'] := %d;" % seed, # random number seed
            "cpconfig['datatype'] := 'codon';", # codon data
            ]
    if model == 'GY94_CF3x4_omega-gamma4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'CF3X4';", # CF3X4 codon frequency model
            "cpconfig['kappa'] := 'e';", # ML estimation of transition/transversion ratio
            "cpconfig['omega'] := 'DGAMMA';", # discrete gamma model for dN/dS
            "cpconfig['wclasses'] := 4;", # number of omega rate categories
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    elif model == 'GY94_CF3x4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'CF3X4';", # CF3X4 codon frequency model
            "cpconfig['kappa'] := 'e';", # ML estimation of transition/transversion ratio
            "cpconfig['omega'] := 'DM0';", # one gamma estimated by ML
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    elif model == 'KOSI07_F_omega-gamma4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model for codons
            "cpconfig['qrates'] := 'KOSI07';", # Kosiol 2007 empirical codon model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'F1XCODONS';", # F codon frequency model
            "cpconfig['kappa'] := 'KAP3';", # ML estimation of kappa(tv)
            "cpconfig['omega'] := 'DGAMMA';", # discrete gamma model for dN/dS
            "cpconfig['wclasses'] := 4;", # number of omega rate categories
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    else:
        raise ValueError("Unrecognized model specification of %s" % model)
    open(configfile, 'w').write('\n'.join(configlines))
    # now run codonphyml
#    cmds = [path, '-i', phylipfile, '-q', '-d', 'codon', '--r_seed', str(seed), '-b', '0', '--quiet']
    cmds = [path, '--darwinconfig', configfile]
    try:
        try:
            p = subprocess.Popen(cmds, stdout=open(outfile, 'w'), stderr=subprocess.PIPE, stdin=subprocess.PIPE)
            (stdoutdata, stderrdata) = p.communicate('Y\n')
        except OSError:
            raise ValueError("codonphyml failed to run, probably the executable is not callable using the following path: %s\n%s" % (path, str(stderrdata)))
        # capture output
        if stderrdata:
            raise ValueError("codonphyml encountered errors:\n%s" % stderrdata)
        if not os.path.isfile(statsfile):
            raise IOError("Failed to find codonphyml stats file %s" % (statsfile))
        if not os.path.isfile(treefile):
            raise IOError("Failed to find codonphyml tree file %s" % (treefile))
        open(outtreefile, 'w').write(open(treefile).read())
        stats = open(statsfile).read()
        open(outstatsfile, 'w').write(stats)
        llmatch = re.compile('\. Log-likelihood: (?P<ll>\-\d+(\.\d+){0,1})')
        m = llmatch.search(stats)
        if not m:
            raise ValueError("Failed to find Log-likelihood in:\n%s" % stats)
        ll = m.group('ll')
        open(llfile, 'w').write(ll)
        ll = float(ll)
    finally:
        for f in [phylipfile, statsfile, treefile]:
            if os.path.isfile(f):
                os.remove(f)
    return ll



if __name__ == '__main__':
    import doctest
    doctest.testmod()
